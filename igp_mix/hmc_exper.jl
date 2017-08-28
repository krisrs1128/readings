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

function marginal_grad(k_theta::Matrix, y::Vector, k_deriv::Matrix)
  alpha = inv(k_theta) * y
  return trace((alpha * alpha' - inv(k_theta)) * k_deriv)
end

# Gaussian Kernel matrix (with noise)
function kernel(x::Matrix, y::Matrix, l::Float64, v0::Float64, v1::Float64)
  n1 = size(x, 1)
  n2 = size(y, 1)
  K = zeros(n1, n2)

  for i = 1:n1
    for j = 1:n2
      K[i, j] = v0 * exp(- 1 / (2 * l ^ 2) * norm(x[i, :] - y[j, :]) ^ 2)
      if i == j
        K[i, j] = v1 + K[i, j]
      end
    end
  end

  return K
end

function k_deriv_logv0(x::Matrix, log_v0::Float64, log_l2::Float64)
  return kernel(x, x, sqrt(exp(log_l2)), exp(log_v0), 0.0)
end

function k_deriv_logl2(x::Matrix, log_v0::Float64, log_l2::Float64)
  n = size(x, 1)
  k_deriv = zeros(n, n)
  for i = 1:n
    for j = 1:n
      norm_sq = norm(x[i, :] - x[j, :]) ^ 2
      k_deriv[i, j] = exp(log_v0 - (1 / 2) * exp(-log_l2) * norm_sq) *
        (exp(-log_l2) * norm_sq / 2)
    end
  end

  return k_deriv
end


function gradient_generator(y::Vector, x::Matrix, a::Float64 = 20.0)

  return function log_posterior_grad(theta::Vector)
    log_l2 = theta[1]
    log_v0 = theta[2]
    log_v1 = theta[3]

    n = size(y, 1)
    k_theta = kernel(x, x, sqrt(exp(log_l2)), exp(log_v0), exp(log_v1))

    ## compute (unnormalized) posterior
    n0a = Normal(0, a)
    unnorm_posterior = logpdf(MvNormal(zeros(n), k_theta), y) +
      logpdf(n0a, log_v0) +
      logpdf(n0a, log_v1) +
      logpdf(n0a, log_l2)

    ## compute gradient
    grad_logl2 = marginal_grad(k_theta, y, k_deriv_logl2(x, log_v0, log_l2)) +
      gradlogpdf(n0a, log_l2)
    grad_logv0 = marginal_grad(k_theta, y, k_deriv_logv0(x, log_v0, log_l2)) +
      gradlogpdf(n0a, log_v0)
    grad_logv1 = marginal_grad(k_theta, y, exp(log_v1) * eye(n)) +
      gradlogpdf(n0a, log_v1)

    return unnorm_posterior, [grad_logl2, grad_logv0, grad_logv1]
  end
end

"""Simulate toy data

# Examples
```julia-repl
julia> n = 100
julia> l = 0.1
julia> v0 = 1.0
julia> v1 = 0.02
julia> simulate(n, l, v0, v1)
```
"""
function simulate(n::Int64, l::Float64, v0::Float64, v1::Float64)
  x = rand(n, 1)
  K = kernel(x, x, l, v0, v1)
  y = rand(MultivariateNormal(zeros(n), K))

  return x, y
end

x, y = simulate(50, sqrt(0.05), 10.0, 0.5)
plot(x = x[:, 1], y = y)

n_samples = 4000
samples = zeros(n_samples, 3)
logf_grad = gradient_generator(y, x)
L = 10
epsilon = 0.005
theta = HMCVariate(samples[1, :], epsilon, L, logf_grad)
for i = 1:n_samples
  sample!(theta)
  samples[i, :] = exp.(theta)
  print(i)
  println(samples[i, :])
end

plot(x = samples[:, 1], y = samples[:, 2])
plot(x = 1:n_samples, y = samples[:, 1])
plot(x = 1:n_samples, y = samples[:, 2])
plot(x = 1:n_samples, y = samples[:, 3])
