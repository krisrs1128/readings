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
      K[i, j] = v0 * exp(- 1 / (2 * l^2) * norm(x[i, :] - y[j, :]) ^ 2)
      if i == j
        K[i, j] += v1
      end
    end
  end

  return K
end


function k_deriv_v0(x::Matrix, l::Float64)
  return kernel(x, x, l, 1.0, 0.0)
end


function k_deriv_l(x::Matrix, v0::Float64, l::Float64)
  n = size(x, 1)
  k_deriv = zeros(n, n)
  for i = 1:n
    for j = 1:n
      norm_sq = norm(x[i, :] - x[j, :]) ^ 2
      k_deriv[i, j] = v0 * exp(- 1 / (2 * l ^ 2) * norm_sq) * 1 / (2 * l^3) * norm_sq
    end
  end

  return k_deriv
end

function gradient_generator(y::Vector, x::Matrix)

  return function log_posterior_grad(theta::Vector)
    v0 = theta[1]
    v1 = theta[2]
    l = theta[3]

    k_theta = kernel(x, x, l, v0, v1)

    ## compute (unnormalized) posterior
    n = size(y, 1)
    unnorm_posterior = logpdf(MvNormal(zeros(n), k_theta), y)

    ## compute gradient
    grad_v0 = marginal_grad(k_theta, y, k_deriv_v0(x, l))
    grad_v1 = marginal_grad(k_theta, y, eye(n))
    grad_l = marginal_grad(k_theta, y, k_deriv_l(x, v0, l))

    return unnorm_posterior, [grad_v0, grad_v1, grad_l]

  end
end

function marginal_grad_wrapper(v0, v1, l)
  marginal_grad(k_theta, y, eye(n))
  marginal_grad(k_theta, y, k_deriv_v0(x, l))
  marginal_grad(k_theta, y, k_deriv_l(x, v0, l))
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
function simulate(n::Int, l::Float64, v0::Float64, v1::Float64)
  x = rand(n, 1)
  K = kernel(x, x, l, v0, v1)
  y = rand(MultivariateNormal(zeros(n), K))

  return x, y
end

epsilon = 0.1
L = 50
Sigma = eye(3)
theta1 = HMCVariate([0.0, 0.0, 0.0], epsilon, L, logfgrad)
