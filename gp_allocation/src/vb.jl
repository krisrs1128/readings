#!/usr/bin/env julia

# File description -------------------------------------------------------------
# Code for allocation of GP experts, based on https://github.com/trungngv/fgp/
#
# author: sankaran.kris@gmail.com
# date: 07/25/2017

using Distributions

###############################################################################
## helper functions used throughout the algorithm
###############################################################################

"""
    kernel(x::Array{Float64,2},
           y::Array{Float64,2},
           sigma::Float64 = 1.0,
           l::Float64 = 1.0)

Compute the kernel matrix between two matrices x and y.

# Arguments
- `x::Array{Float64, 2}`: Each row is an ``x_{i}`` when computing ``k\left(x_{i}, y_{j}\right)``.
- `y::Array{Float64, 2}`: Each row is a ``y_{j}`` when computing ``k\left(x_{i}, y_{j}\right)``.
- `sigma::Float64 = 1.0`: The bandwidth of the kernel.
- `l::Float64 = 1.0`: The scale of the kernel.

```julia-repl
kernel(randn(10, 2), randn(5, 2))
```
"""
function kernel(x::Array{Float64, 2},
                y::Array{Float64, 2},
                sigma::Float64 = 1.0,
                l::Float64 = 1.0)
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

"""
    lse(x::Array{Float64, 1})

Compute the log-sum-exp function.

```jldoctest
a = lse([log(0.1), log(0.4)]);
b = log(0.5);

a == b
true
```
"""
function lse(x::Array{Float64, 1})
  return log(sum(exp(x)))
end

function normalize_log_space(log_x::Array{Float64, 1})
  log_x - lse(log_x)
end

###############################################################################
## E-step for variational bayes updates
###############################################################################

function lambda(x::Array{Float64, 2},
                uk::Array{Float64, 2},
                l::Float64,
                a::Float64)

"""
  Compute the log marginal likelihood

This is the log marginal likelihood, as defined in equation 22 of "Fast
Allocation of Gaussian Process Experts"

"""
  kappa = function(x, y)
    kernel(x, y, l, a)
  end

  return diagm(
    diag(kappa(x, x) - kappa(x, uk) * inv(kappa(uk, uk)) * kappa(x, uk)')
  )
end

"""
    update_log_rho(x::Array{Float64, 2},
                   y::Vector{Float64},
                   r::Array{Float64, 2},
                   m::Array{Float64, 2},
                   V::Array{Float64, 2},
                   u::Array{Array{Float64, 2}, 1},
                   sigma::Array{Float64, 1},
                   l::Float64 = 1.0,
                   a::Float64 = 1.0)

Compute the rho's necessary for the E-step

# Arguments
- `x::Array{Float64, 2}`: The array of observed covariates.
- `y::Vector{Float64}`: The vector of observed response values. The i^{th}
  row of x corresponds to the i^{th} entry of y.
- `r::Array{FLoat64, 2}`: An n x K array of responsibilities computed in the
  previous update. The i^{th} row corresponds to the i^{th} sample, and the
  k^{th} column corresponds to the k^{th} cluster component
- `m::Array{Float64, 2}`: An P x K array whose k^{th} column is the mean for the
  k^{th} variational approximation distribution.
- `V::Array{Float64, 2}`: A P x P covariance for the variational approximation
  distributions. This is shared across all K clusters.
- `u::Array{Array{Float64, 2}, 2}`: The positions of the inducing points. The
  k^{th} element of the outer array corresponds to the k^{th} cluster. The Nk x
  P array within this element gives the actual inducing points for this cluster.
- `sigma::Array{Float64, 1}`: The noise of observed data y around the GP mean f.
- `l::Float64 = 1.0`: The scale of the kernel.
- `a::Float64 = 1.0`: The bandwidth of the kernel.

```julia-repl
# generated in the simulate.R file
x = randn(10, 2)
y = randn(10)
r = rand(Uniform(), 10)
r = [r/2 r/2 1 - r]
m = rand(3, 2)
V = [1.0 0; 0 1.0]
u = [rand(2, 2), rand(2, 2), rand(2, 2)]
sigma = [1.0, 1.0, 1.0]
update_log_rho(x, y, r, m, V, u, sigma, 1.0, 1.0)
```
"""
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
  K, _ = size(m)
  log_rho = zeros(n, K)

  kappa = function(x, y)
    kernel(x, y, l, a)
  end

  for i = 1:n
    for k = 1:K
      x_distn = MvNormal(m[k, :], V)
      lambda_k = lambda(x, u[k], l, a)
      gamma_k = diagm(r[:, k]) .* inv(lambda_k + (sigma[k] ^ 2) * eye(n))
      psi_k = kappa(u[k], u[k]) + kappa(u[k], x) * gamma_k * kappa(x, u[k])
      mu_gk = kappa(u[k], u[k]) *  inv(psi_k) * kappa(u[k], x) * gamma_k * y
      mu_y = kappa(x[i:i, :], u[k]) * inv(kappa(u[k], u[k])) * mu_gk
      sigma_y = lambda_k[i, i] + sigma[k] ^ 2

      y_distn = Normal(mu_y[1, 1], sigma_y)
      log_rho[i, k] = logpdf(x_distn, x[i, :]) + logpdf(y_distn, y[i])
    end
    log_rho[i, :] = normalize_log_space(log_rho[i, :])
  end

  return (log_rho)
end

function Q(xk::Array{Float64, 2}, uk::Array{Float64, 2})
  kappa = function(x, y)
    kernel(x, y, l, a)
  end

  return kappa(xk, uk) * inv(kappa(uk, uk)) * kappa(uk, xk)
end


"""
    log_marginal(y::Vector{Float64},
                 x::Array{Float64, 2},
                 z::Vector{Int64},
                 u::Array{Array{Float64, 2}, 1},
                 sigma::Array{Float64, 1},
                 l::Float64,
                 a::Float64)

Compute the log marginal probability of y, as defined in equation (22) of "Fast
Allocation of Gaussian Processes"

```julia-repl
y = randn(10)
x = randn(10, 2)
z = [1, 1, 1, 1, 2, 2, 2, 3, 3, 1]
u = [rand(2, 2), rand(2, 2), rand(2, 2)]
sigma = [1.0, 1.0, 1.0]
log_marginal(y, x, z, u, sigma, 1.0, 1.0)
```
"""
function log_marginal(y::Vector{Float64},
                      x::Array{Float64, 2},
                      z::Vector{Int64},
                      u::Array{Array{Float64, 2}, 1},
                      sigma::Array{Float64, 1},
                      l::Float64,
                      a::Float64)
  n, p = size(x)
  K = length(u)
  log_liks = zeros(K)

  for k = 1:K
    xk = x[z .== k, :]
    yk = y[z .== k]
    Q_k = Q(xk, u[k])
    lambda_k = lambda(xk, u[k], l, a)
    marginal_cov = Q_k + lambda_k + diagm(sigma[k] * ones(size(xk, 1)))
    log_liks[k] = -0.5 * (logdet(marginal_cov) + yk' * inv(marginal_cov) * yk)[1]
  end

  return log_liks
end
