#!/usr/bin/env julia

# File description -------------------------------------------------------------
#
# Some functions used in stochastic variational gradient descent. Just some
# experiments trying to better understand the method.
#
# Basically just me rehashing
# https://github.com/DartML/Stein-Variational-Gradient-Descent
#
# author: sankaran.kris@gmail.com
# date: 12/04/2017

function svgd(x::Array{Float64},
              eps::Array{Float64},
              grad_logp::Function,
              kernel::Function,
              grad_kernel::Function)
  for iter = 1:length(eps)
    if iter % 10 == 0
      println("Iteration ", iter)
    end

    phi_hat = estimate_phi(x, grad_logp, kernel, grad_kernel)
    for i = 1:length(x)
      x[i] = x[i] + eps[iter] * phi_hat(x[i])
    end
  end

  x
end

function estimate_phi(x::Array{Float64},
                      grad_logp::Function,
                      kernel::Function,
                      grad_kernel::Function)
  return function(z::Float64)
    phi_vec = dot(kernel(x, z), grad_logp(x)) +
      grad_kernel(x, z)
    mean(phi_vec)
  end
end
