#!/usr/bin/env julia

# File description -------------------------------------------------------------
#
# Experiment using svgd for normal mixture problem.
#
# author: sankaran.kris@gmail.com
# date: 12/04/2017

using PyPlot
using Distributions
include("svgd.jl")

"""Gradient of Log Likelihood / Log Posterior

This gradient can be computed tractably, and is required for inference with
SVGD.

# Arguments
x::Array{Float64}:
mu2::Float64:
sigma2::Float64:
"""
function grad_logp(x::Array{Float64}, mu2::Float64 = 10.0, sigma2::Float64 = 10.0)
  p1 = Normal()
  p2 = Normal(mu2, sigma2)

  (-x .* pdf(p1, x) - (x - mu2) / sigma2 ^ 2 .* pdf(p2, x)) ./
    (pdf(p1, x) + pdf(p2, x))
end

"""Gradient of Log Likelihood / Log Posterior
"""
function grad_logp(x::Float64, mu2::Float64 = 10.0, sigma2::Float64 = 1.0)
  x_vec = [x]
  grad_p([x], mu2, sigma2)
end


"""Kernel Expansion around Observed Data
"""
function kernel(x::Array{Float64}, z::Float64, h::Float64 = 1.0)
  n = length(x)
  pdfs = zeros(n)
  for i = 1:n
    pdfs[i] = pdf(Normal(x[i], h), z)
  end

  pdfs
end

"""Gradient of RBF Kernel
"""
function grad_kernel(x::Array{Float64}, z::Float64, h::Float64 = 1.0)
  - (x - z) / h .* pdf(Normal(z, h), x)
end

x0 = rand(Uniform(-10, 20), 400)
eps = collect(4:0.1:100) .^ (-1.5)
p_fit = svgd(x0, eps, grad_logp, kernel, grad_kernel)

plt[:hist](p_fit, 100)
