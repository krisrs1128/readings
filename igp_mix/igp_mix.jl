#!/usr/bin/env julia

# File description -------------------------------------------------------------
#
# Functions for using HMC to optimize hyperparameters in a GP, following
# http://mlg.eng.cam.ac.uk/pub/pdf/WilRas96.pdf
#
# This is useful to think about in the context of mixtures of GPs, because it's
# step 3 in http://mlg.eng.cam.ac.uk/zoubin/papers/iMGPE.pdf
#
# author: sankaran.kris@gmail.com
# date: 08/29/2017

using Distributions
using Mamba

###############################################################################
#                        Define types for GP inference                        #
###############################################################################

type KernelParam
  l::Float64
  v0::Float64
  v1::Float64
end


type GPHyper
  log_l2
  log_v0
  log_v1
end


type GPModel
  theta::KernelParam
  x_train::Matrix
  y_train::Vector
end

function param_from_theta(theta::Vector)
  KernelParam(theta[1], theta[2], theta[3])
end

###############################################################################
#                                Define kernel                                #
###############################################################################
"""Covariance Kernel with Noise

# Arguments
x::Matrix: Each row gives a value for rows x_{i} in the matrix K(x_{i}, y_{j})
y::Matrix: Each row gives a value for columns y_{j} in the matrix K(x_{i}, y_{j})
theta::KernelParam: Parameters for the kernel (bandwidth l, amplitude v0, and
  noise v1)

# Examples
```julia-repl
theta = KernelParam(1, .2, 1)
kernel(rand(10, 2), rand(10, 2), theta)
```
"""
function kernel(x::Matrix, y::Matrix, theta::KernelParam)
  n1 = size(x, 1)
  n2 = size(y, 1)
  K = zeros(n1, n2)

  for i = 1:n1
    for j = 1:n2
      K[i, j] = theta.v0 * exp(- 1 / (2 * theta.l ^ 2) * norm(x[i, :] - y[j, :]) ^ 2)
      if i == j
        K[i, j] = theta.v1 + K[i, j]
      end
    end
  end

  K
end


"""Simulate toy data

# Examples
```julia-repl
julia> n = 100
julia> theta = KernelParam()
julia> simulate(n, theta)
```
"""
function simulate(n::Int64, theta::KernelParam)
  x = rand(n, 1)
  K = kernel(x, x, theta)
  y = rand(MultivariateNormal(zeros(n), K))
  x, y
end


function gp_posterior(x_new::Matrix,
                      gp::GPModel,
                      epsilon::Float64 = 1e-8)
  theta_noiseless = deepcopy(gp.theta)
  theta_noiseless.v1 = 0

  k_x_xnew = kernel(gp.x_train, x_new, theta_noiseless)
  k_xnew_xnew = kernel(x_new, x_new, theta_noiseless)
  k_train = kernel(gp.x_train, gp.x_train, gp.theta)
  inv_k_train = inv(k_train)

  mu = k_x_xnew' * inv_k_train * gp.y_train
  Sigma = k_xnew_xnew - k_x_xnew' * inv_k_train * k_x_xnew + epsilon * eye(size(x_new, 1))
  Distributions.MvNormal(mu, Hermitian(Sigma))
end


###############################################################################
#                   Define kernel derivatives and GP sampelr                  #
###############################################################################
"""Gradient of Marginal Likelihood

Computes the gradient of the marginal likelihood, given dK / dtheta

# Arguments
k_theta::Matrix: The kernel matrix evaluated at some coordinates x
y::Vector: The response at the coordinates x
k_deriv: The matrix dK / dtheta
"""
function marginal_grad(k_theta::Matrix, y::Vector, k_deriv::Matrix)
  alpha = inv(k_theta) * y
  return trace((alpha * alpha' - inv(k_theta)) * k_deriv)
end


"""Derivative of Kernel w.r.t log-amplitude log_v0

# Arguments
x::Matrix: Coordinates at which the kernel matrix was evaluated
theta::KernelParam: Parameters for the kernel (bandwidth l, amplitude v0, and
  noise v1)
"""
function k_deriv_logv0(x::Matrix, theta::KernelParam)
  theta.v1 = 0
  kernel(x, x, theta)
end


"""Derivative of Kernel w.r.t. log-bandwidth log_l2

# Arguments
x::Matrix: Coordinates at which the kernel matrix was evaluated
theta::KernelParam: Parameters for the kernel (bandwidth l, amplitude v0, and
  noise v1)
"""
function k_deriv_logl2(x::Matrix, theta::KernelParam)
  n = size(x, 1)
  k_deriv = zeros(n, n)
  for i = 1:n
    for j = 1:n
      norm_sq = norm(x[i, :] - x[j, :]) ^ 2
      k_deriv[i, j] = theta.v0 * exp(- 1 / (2 * theta.l ^ 2) * norm_sq) *
        norm_sq / (2 * theta.l ^ 2)
    end
  end

  k_deriv
end


"""Gradient Functions at given Coordinates

# Arguments
y::Vector: The response at the coordinates x
x::Matrix: Coordinates at which the kernel matrix was evaluated
a::GPHyper: Priors for parameters in kernel function
"""
function gradient_generator(y::Vector, x::Matrix, a::GPHyper)

  function log_posterior_grad(log_theta::Vector)
    log_l2 = log_theta[1]
    log_v0 = log_theta[2]
    log_v1 = log_theta[3]
    theta = KernelParam(sqrt(exp(log_l2)), exp(log_v0), exp(log_v1))

    k_theta = kernel(x, x, theta)
    n = length(y)

    ## compute (unnormalized) posterior
    unnorm_posterior = logpdf(MvNormal(zeros(n), k_theta), y) +
      logpdf(a.log_l2, log_v0) +
      logpdf(a.log_v0, log_v1) +
      logpdf(a.log_v1, log_l2)

    ## compute gradient
    grad_logl2 = marginal_grad(k_theta, y, k_deriv_logl2(x, theta)) +
      gradlogpdf(a.log_l2, log_l2)
    grad_logv0 = marginal_grad(k_theta, y, k_deriv_logv0(x, theta)) +
      gradlogpdf(a.log_v0, log_v0)
    grad_logv1 = marginal_grad(k_theta, y, exp(log_v1) * eye(n)) +
      gradlogpdf(a.log_v1, log_v1)

    unnorm_posterior, [grad_logl2, grad_logv0, grad_logv1]
  end
end

function GPSampler(x::Matrix,
                   y::Vector,
                   a::GPHyper,
                   n_samples::Int64 = 1000,
                   L::Int64 = 5,
                   epsilon::Float64 = 0.005,
                   theta0::Vector = zeros(3))
  logf_grad = gradient_generator(y, x, a)

  samples = zeros(n_samples, 3)
  sampler = HMCVariate(theta0, epsilon, L, logf_grad)
  for i = 1:n_samples
    samples[i, :] = exp.(sampler)
    println(i, "\t", samples[i, :])
    sample!(sampler)
  end

  samples
end

