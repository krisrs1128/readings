#!/usr/bin/env julia

# File description -------------------------------------------------------------
#
# Experiment using HMC to optimize hyperparameters in a GP, following
# http://mlg.eng.cam.ac.uk/pub/pdf/WilRas96.pdf
#
# This is useful to think about in the context of mixtures of GPs, because it's
# step 3 in http://mlg.eng.cam.ac.uk/zoubin/papers/iMGPE.pdf
#
# author: sankaran.kris@gmail.com
# date: 08/24/2017

using Distributions
using Gadfly
using Mamba


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


type GP
  theta::KernelParam
  a::GPHyper
end


# Gaussian Kernel matrix (with noise)
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

  return K
end


function marginal_grad(k_theta::Matrix, y::Vector, k_deriv::Matrix)
  alpha = inv(k_theta) * y
  return trace((alpha * alpha' - inv(k_theta)) * k_deriv)
end


function k_deriv_logv0(x::Matrix, theta::KernelParam)
  theta.v1 = 0
  return kernel(x, x, theta)
end


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

  return k_deriv
end


function gradient_generator(y::Vector, x::Matrix, a::GPHyper)

  return function log_posterior_grad(log_theta::Vector)
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

    return unnorm_posterior, [grad_logl2, grad_logv0, grad_logv1]
  end
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

  return x, y
end

theta = KernelParam(sqrt(0.05), 10.0, 0.5)
x, y = simulate(50, theta)
plot(x = x[:, 1], y = y)

n_samples = 1000
samples = zeros(n_samples, 3)
a = GPHyper(
  Distributions.TDist(1),
  Distributions.TDist(1),
  Distributions.TDist(1)
)
logf_grad = gradient_generator(y, x, a)
L = 5
epsilon = 0.005
sampler = HMCVariate(samples[1, :], epsilon, L, logf_grad)
for i = 1:n_samples
  sample!(sampler)
  samples[i, :] = exp.(sampler)
  print(i)
  println(samples[i, :])
end

plot(x = samples[:, 1], y = samples[:, 2])
plot(x = 1:n_samples, y = samples[:, 1])
plot(x = 1:n_samples, y = samples[:, 2])
plot(x = 1:n_samples, y = samples[:, 3])
