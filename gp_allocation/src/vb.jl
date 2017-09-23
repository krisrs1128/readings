#!/usr/bin/env julia

# File description -------------------------------------------------------------
# Code for allocation of GP experts, based on https://github.com/trungngv/fgp/
#
# author: sankaran.kris@gmail.com
# date: 07/25/2017

using Distributions
using Clustering
using Optim

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
- `x::Matrix{Float64}`: Each row is an ``x_{i}`` when computing ``k\left(x_{i}, y_{j}\right)``.
- `y::Matrix{Float64}`: Each row is a ``y_{j}`` when computing ``k\left(x_{i}, y_{j}\right)``.
- `sigma::Float64 = 1.0`: The bandwidth of the kernel.
- `l::Float64 = 1.0`: The scale of the kernel.

```julia-repl
kernel(randn(10, 2), randn(5, 2))
```
"""
function kernel(x::Matrix{Float64},
                y::Matrix{Float64},
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
    lse(x::Vector{Float64})

Compute the log-sum-exp function.

```jldoctest
a = lse([log(0.1), log(0.4)]);
b = log(0.5);

a == b
true
```
"""
function lse(x::Vector{Float64})
  return log(sum(exp(x)))
end

function normalize_log_space(log_x::Vector{Float64})
  log_x - lse(log_x)
end

###############################################################################
## E-step for variational bayes updates
###############################################################################

"""
    lambda(x::Matrix{Float64},
           uk::Matrix{Float64},
           l::Float64,
           a::Float64)

Compute the covariance for the variational approximation

This is defined in equation (8) of "Fast Allocation of Gaussian Process Experts"
"""
function lambda(x::Matrix{Float64},
                uk::Matrix{Float64},
                l::Float64,
                a::Float64)
  kappa = function(x, y)
    kernel(x, y, l, a)
  end

  return diagm(
    diag(kappa(x, x) - kappa(x, uk) * inv(kappa(uk, uk)) * kappa(x, uk)')
  )
end

function update_m(u::Vector{Matrix{Float64}})
  K = length(u)
  _, p = size(u[1])
  m = zeros(p, K)

  for k = 1:K
    m[:, k] = mean(u[k], 1)
  end

  return m
end

function update_v(u::Vector{Matrix{Float64}})
  K = length(u)
  _, p = size(u[1])
  vk = zeros(p, K)

  for j = 1:p
    for k = 1:K
      vk[j, k] = var(u[k][:, j])
    end
  end

  return Matrix(Diagonal(mean(vk, 2)[:, 1]))
end

"""
    update_log_rho(x::Matrix{Float64},
                   y::Vector{Float64},
                   r::Matrix{Float64},
                   m::Matrix{Float64},
                   V::Matrix{Float64},
                   u::Array{Matrix{Float64}, 1},
                   sigma::Vector{Float64},
                   l::Float64 = 1.0,
                   a::Float64 = 1.0)

Compute the rho's necessary for the E-step

# Arguments
- `x::Matrix{Float64}`: The array of observed covariates.
- `y::Vector{Float64}`: The vector of observed response values. The i^{th}
  row of x corresponds to the i^{th} entry of y.
- `r::Array{FLoat64, 2}`: An n x K array of responsibilities computed in the
  previous update. The i^{th} row corresponds to the i^{th} sample, and the
  k^{th} column corresponds to the k^{th} cluster component
- `m::Matrix{Float64}`: An p x K array whose k^{th} column is the mean for the
  k^{th} variational approximation distribution.
- `V::Matrix{Float64}`: A p x p covariance for the variational approximation
  distributions. This is shared across all K clusters.
- `u::Array{Matrix{Float64}, 2}`: The positions of the inducing points. The
  k^{th} element of the outer array corresponds to the k^{th} cluster. The Nk x
  p array within this element gives the actual inducing points for this cluster.
- `sigma::Vector{Float64}`: The noise of observed data y around the GP mean f.
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
function update_log_rho(x::Matrix{Float64},
                        y::Vector{Float64},
                        r::Matrix{Float64},
                        m::Matrix{Float64},
                        V::Matrix{Float64},
                        u::Vector{Matrix{Float64}},
                        sigma::Vector{Float64},
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
      x_distn = MvNormal(m[:, k], V)
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

  return log_rho
end

"""
    Q(xk::Matrix{Float64}, uk::Matrix{Float64})

Compute the marginal covariance for the k^{th} cluster, as defined just after
equation (22) of "Fast Allocation of Gaussian Process Experts"
"""
function Q(xk::Matrix{Float64}, uk::Matrix{Float64})
  kappa = function(x, y)
    kernel(x, y, l, a)
  end

  return kappa(xk, uk) * inv(kappa(uk, uk)) * kappa(uk, xk)
end


"""
    log_marginal(y::Vector{Float64},
                 x::Matrix{Float64},
                 z::Vector{Int64},
                 u::Array{Matrix{Float64}, 1},
                 sigma::Vector{Float64},
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
                      x::Matrix{Float64},
                      z::Vector{Int64},
                      u::Array{Matrix{Float64}, 1},
                      sigmas::Vector{Float64},
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
    marginal_cov = Q_k + lambda_k + diagm(sigmas[k] * ones(size(xk, 1)))
    println(eig(marginal_cov)[1])
    log_liks[k] = -0.5 * (logdet(marginal_cov) + yk' * inv(marginal_cov) * yk)[1]
  end

  return log_liks
end

###############################################################################
## Full VB iteration
###############################################################################
function z_assignment(rho::Matrix{Float64})
  n, _ = size(rho)
  z = Vector{Int64}(n)
  for i = 1:n
    _, z[i] = findmax(rho[i, :])
  end

  return z
end

function vb(y::Vector{Float64},
            x::Matrix{Float64},
            K::Int64,
            n_inducing::Int64,
            n_iter::Int64)

  ## Initialization
  n = size(x, 1)
  inducing_init = kmeans(x', n_inducing).centers'
  u = Vector{Matrix{Float64}}(K)
  ix = round(linspace(1, K, n_inducing))
  for k = 1:K
    u[k] = inducing_init[ix .== k, :]
  end

  log_rho = log(rand(Dirichlet(ones(K)), n))'
  log_sigmas = log(var(y) * ones(K))
  l = 1.0
  a = 1.0

  ## get size is ns
  ns = zeros(Int64, K)
  for k = 1:K
    ns[k] = size(u[k], 1)
  end

  for iter = 1:n_iter
    ## E-step
    println(iter)
    m = update_m(u)
    V = update_v(u)
    log_rho = update_log_rho(x, y, exp(log_rho), m, V, u, exp(log_sigmas), l, a)
    z = z_assignment(exp(log_rho))

    ## M-step
    p = size(u[1], 2)
    f = function(hyper::Vector{Float64})
      u2, log_sigmas2, l2, a2 = unstack_params(hyper, ns, p)
      val = -sum(log_marginal(y, x, z, u2, exp(log_sigmas2), l2, a2))
     return val
    end

    hyper = stack_params(u, log_sigmas, l, a)
    opt = optimize(f, hyper, NelderMead(), Optim.Options(iterations = 10))
    u, log_sigmas, l, a = unstack_params(opt.minimizer, ns, p)
  end

end

function stack_params(u::Vector{Matrix{Float64}},
                      sigmas::Vector{Float64},
                      l::Float64,
                      a::Float64)
  K = length(u)
  p = size(u[1], 2)
  hyper = [sigmas; [l, a]]

  for k in 1:K
    for j in 1:p
      hyper = [hyper; u[k][:, j]]
    end
  end

  return hyper
end

function unstack_params(hyper::Vector{Float64}, ns::Vector{Int64}, p::Int64)
  chyper = copy(hyper)
  K = length(ns)
  log_sigmas = chyper[1:K]
  l = chyper[K + 1]
  a = chyper[K + 2]

  deleteat!(chyper, 1:(K + 2))
  u = Vector{Matrix{Float64}}(K)
  for k = 1:K
    cur_u = zeros(ns[k], p)
    for j = 1:p
      cur_u[:, j] = chyper[1:ns[k]]
      deleteat!(chyper, 1:ns[k])
    end
    u[k] = cur_u
  end

  return u, log_sigmas, l, a
end

y = randn(200)
x = randn(200, 2)
n_iter = 100
n_inducing = 50
K = 3
