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

function marginal_grad(k_theta::Matrix, y::Vector, k_deriv::Matrix)
  alpha = inv(k_theta) * y
  return trace((alpha * alpha' - inv(k_theta)) * k_deriv)
end

# Gaussian Kernel matrix (with noise)
function kernel(X::Matrix, Y::Matrix, l::float, v0::float, v1::float)
  n1 = size(X, 1)
  n2 = size(Y, 1)
  K = zeros(n1, n2)

  for i = 1:n1
    for j = 1:n2
      K[i, j] = v0 * exp(- 1 / (2 * l^2) * norm(X[i, :] - Y[j, :]) ^ 2)
      if i == j
        K[i, j] += v1
      end
    end
  end

  return K
end


function k_deriv_v0(X::Matrix, l::float)
  return kernel(X, X, l, 1, 0)
end


function k_deriv_l(X::Matrix, v0::float, l;;float)
  n = size(X, 1)
  k_deriv = zeros(n, n)
  for i = 1:n1
    for j = 1:n2
      norm_sq = norm(X[i, :] - X[j, :]) ^ 2
      k_deriv[i, j] = v0 * exp(- 1 / (2 * l ^ 2) * norm_sq) * 1 / (2 * l^3) * norm_sq
    end
  end

  return k_deriv
end
