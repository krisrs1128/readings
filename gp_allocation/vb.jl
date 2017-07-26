#!/usr/bin/env julia

# File description -------------------------------------------------------------
# Code for allocation of GP experts, based on https://github.com/trungngv/fgp/
#
# author: sankaran.kris@gmail.com
# date: 07/25/2017

using Distributions

# generated in the simulate.R file
y = readcsv("data/y.csv", header = true)[1];

function kernel(x::Array{Float64, 2},
                y::Array{Float64, 2},
                sigma::Float64,
                l::Float64)
  print(x)
  n_x = size(x, 1)
  n_y = size(y, 1)
  kernel_matrix = zeros(n_x, n_y)

  for i = 1:n_x
    for j = 1:n_y
      kernel_matrix[i, j] = sigma ^ 2 * exp(-1 / (2 * l ^ 2) * (x[i] - y[j]) ^ 2)
    end
  end

  return kernel_matrix
end

function update_log_rho(x::Array{Float64, 2},
                        y::Vector{Float64},
                        r::Array{Float64, 2},
                        m::Array{Float64, 2},
                        V::Array{Float64, 2},
                        u::Array{Array{Float64, 2}, 1},
                        sigma::Array{Float64, 1},
                        l::Float64,
                        a::Float64)
  n, d = size(x)
  _, K = size(m)
  log_rho = zeros(n, K)

  kappa = function(x, y)
    kernel(x, y, l, a)
  end

  for i = 1:n
    for k = 1:K
      x_distn = MvNormal(m[:, k], V)

      lambda_k = diagm(diag(kappa(x, x) - kappa(x, u[k]) * inv(kappa(u[k], u[k])) * kappa(x, u[k])'))
      gamma_k = diagm(r[:, k]) .* inv(lambda_k + (sigma[k] ^ 2) * eye(n))
      psi_k = kappa(u[k], u[k]) + kappa(u[k], x) * gamma_k * kappa(x, u[k])
      mu_gk = kappa(u[k], u[k]) *  inv(psi_k) * kappa(u[k], x) * gamma_k * y
      mu_y = kappa(x[i:i, :], u[k]) * inv(kappa(u[k], u[k])) * mu_gk
      sigma_y = lambda_k[i, i] + sigma[k] ^ 2

      y_distn = Normal(mu_y[1, 1], sigma_y)
      log_rho[i, k] = normalize_log_space(logpdf(x_distn) + logpdf(y_distn))
    end
  end

end

# x = randn(10, 2)
# y = randn(10, 1)
# r = rand(Uniform(), 10)
# r = [r 1 - r]
# m = rand(2, 3)
# V = [1 0; 0 1]
# u = [rand(2, 2), rand(2, 2), rand(2, 2)]
# sigma = [1, 1, 1]
