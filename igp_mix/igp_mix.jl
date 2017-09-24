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
# date: 09/06/2017

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

###############################################################################
#                           generic helper functions                          #
###############################################################################
function param_from_theta(theta::Vector)
  KernelParam(theta[1], theta[2], theta[3])
end

function lse(log_x::Vector)
  m = maximum(log_x)
  m + log(sum(exp.(log_x - m)))
end

function normalize_log_space(log_x::Vector{Float64})
  log_x - ones(length(log_x)) * lse(log_x)
end

function normalize_log_space(log_x::Dict{Int64, Float64})
  dict_keys = collect(keys(log_x))
  new_vals = normalize_log_space(collect(values(log_x)))

  normalized_log_x = Dict{Int64, Float64}()
  for k in 1:length(dict_keys)
    normalized_log_x[dict_keys[k]] = new_vals[k]
  end

  return normalized_log_x
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

function rand_kernel(a::GPHyper)
  KernelParam(
    sqrt(exp(rand(a.log_l2))),
    exp(rand(a.log_v0)),
    exp(rand(a.log_v1))
  )
end

"""Simulate a Mixture of GPs

# Arguments

# Examples
```julia-repl
a = GPHyper(
  Distributions.Normal(),
  Distributions.Normal(),
  Distributions.Normal(),
)
thetas = Dict(1 => rand_kernel(a), 2 => rand_kernel(a), 3 => rand_kernel(a))
c, x, y = simulate_mix(200, thetas)

using Gadfly
plot(x = x[:, 1], y = y, color = c)
```
"""
function simulate_mix(n::Int64, thetas::Dict{Int64, KernelParam})
  x = zeros(n, 1)
  y = zeros(n)
  K = length(thetas)
  c = rand(1:K, n)

  for k = 1:K
    x[c .== k, :], y[c .== k] = simulate(sum(c .== k), thetas[k])
    x[c .== k, :] += k / 2
  end

  c, x, y
end

function gp_posterior(x_new::Matrix,
                      gp::GPModel,
                      epsilon::Float64 = 1e-5)
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

"""Log PDF for a GP Model

During sampling, it is useful at various points to evaluate the log probability
of a single GP model.
"""
function gp_logpdf(gp::GPModel)
  k_theta = kernel(gp.x_train, gp.x_train, gp.theta)
  gsn = MvNormal(k_theta)
  logpdf(gsn, gp.y_train)
end

"""Log PDF for Independent GPs

This just adds together the log PDFs for several GP models.
"""
function mix_gp_logpdf(gps::Vector{GPModel})
  K = length(gps)
  logpdf = 0
  for k = 1:K
    logpdf += gp_logpdf(gps[k])
  end
end

"""Log EPPF for Dirichlet Process

See page 23 in this tutorial, for example
https://www.stats.ox.ac.uk/~teh/teaching/npbayes2012/dp.pdf

```julia-repl
nk = [12, 5, 24]
alpha = 2.0
dp_log_eppf(nk, alpha)
```
"""
function dp_log_eppf(nk::Vector, alpha::Float64)
  K = length(nk)
  K * log(alpha) + lgamma(alpha) - lgamma(alpha + sum(nk)) + sum(lgamma.(nk))
end

"""Table Probabilities for Next Customer in CRP

```julia-repl
c = [1, 2, 2, 1, 1, 3]
alpha = 2.0
crp_log_prob(c, alpha)
# Case with no member in class 3
c = [1, 2, 2, 1, 1, 4]
crp_log_prob(c, alpha)
```
"""
function crp_log_prob(c::Vector{Int64}, alpha::Float64, K::Int64)
  n = length(c) + 1

  log_probs = Dict{Int64, Float64}()
  log_probs[maximum(c) + 1] = log(alpha) - log(n + alpha - 1)
  for k in unique(c)
    log_probs[k] = log(sum(c .== k)) - log(n + alpha - 1)
  end

  normalize_log_space(log_probs)
end

###############################################################################
#        Conditional probabilities for sampling class memberships c[i]         #
###############################################################################
function substitute_probs(update_ix::Int64,
                          c0::Vector{Int64},
                          x::Matrix,
                          y::Vector,
                          thetas::Dict{Int64, KernelParam})
  K = length(thetas)
  c = deepcopy(c0)

  ## add that sample into different sets
  sub_probs = Dict{Int64, Float64}()
  for k in collect(keys(thetas))
    c[update_ix] = k
    sub_probs[k] = gp_logpdf(
      GPModel(thetas[k], x[c .== k, :], y[c .== k])
    )
  end

  sub_probs
end

function gp_logpdf_wrapper(c::Vector{Int64},
                           x::Matrix,
                           y::Vector,
                           thetas::Dict{Int64, KernelParam})
  K = length(thetas)
  ref_probs = Dict{Int64, Float64}()

  ## Get probs for each component
  for k in collect(keys(thetas))
    ref_probs[k] = gp_logpdf(
      GPModel(thetas[k], x[c .== k, :], y[c .== k])
    )
  end

  ref_probs
end

function class_conditional(update_ix::Int64,
                           c::Vector{Int64},
                           x::Matrix,
                           y::Vector,
                           thetas::Dict{Int64, KernelParam},
                           alpha::Float64,
                           a::GPHyper)
  keep_ix = 1:size(x, 1) .!= update_ix
  prior = crp_log_prob(c[keep_ix], alpha, maximum(keys(thetas)))
  if !(c[update_ix] in keys(prior))
    prior[c[update_ix]] = 0
  end

  sub_probs = substitute_probs(update_ix, c, x, y, thetas)
  ref_probs = gp_logpdf_wrapper(c[keep_ix], x[keep_ix, :], y[keep_ix], thetas)

  liks = Dict{Int64, Float64}()
  ref_sum = sum(collect(values(ref_probs)))
  for k in unique(c)
    liks[k] = ref_sum + sub_probs[k] - ref_probs[k]
  end

  theta_new = rand_kernel(a)
  liks[maximum(keys(prior))] = ref_sum +
    gp_logpdf(GPModel(theta_new, x[[update_ix], :], y[[update_ix]]))

  post = merge(+, b, e)
  normalize_log_space(post), theta_new
end

function class_counts(c::Vector{Int64})
  K = maximum(c)
  counts = zeros(K)
  for i = 1:length(c)
    counts[c[i]] += 1
  end
  counts
end

function joint_log_prob(c::Vector{Int64},
                        x::Matrix,
                        y::Vector,
                        thetas::Dict{Int64, KernelParam},
                        alpha::Float64)
  liks = gp_logpdf_wrapper(c, x, y, thetas)
  nk = class_counts(c)
  eppf = dp_log_eppf(nk, alpha)
  sum(liks) + eppf
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
    unnorm_posterior = logpdf(MvNormal(k_theta), y) +
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
    sample!(sampler)
  end

  samples
end

###############################################################################
#                             Mix-GP Gibbs sampler                            #
###############################################################################

type MixGPState
  c::Vector{Int64}
  thetas::Dict{Int64, KernelParam}
end

function write_state(iter::Int64, state::MixGPState, out_path::String)
  c_path = string(out_path, "c.csv")
  thetas_path = string(out_path, "thetas.csv")

  n = length(state.c)
  open(c_path, "a") do x
    writecsv(x, [iter * ones(n) collect(1:n) state.c])
  end

  open(thetas_path, "a") do x
    theta_array = Matrix(state.thetas)
    writecsv(x, [iter * ones(size(theta_array, 1)) theta_array])
  end
end

function Vector(theta::KernelParam)
  [theta.l theta.v0 theta.v1]
end

function Matrix(thetas::Dict{Int64, KernelParam})
  theta_array = zeros(length(thetas), 4)
  dict_keys = collect(keys(thetas))
  for k = 1:length(dict_keys)
    theta_vec = Vector(thetas[dict_keys[k]])
    theta_array[k, 1] = dict_keys[k]
    theta_array[k, 2:4] = theta_vec
  end

  theta_array
end

function sweep_indicators!(state::MixGPState,
                           x::Matrix,
                           y::Vector,
                           alpha::Float64,
                           a::GPHyper)
  n = length(state.c)

  for i = 1:n
    c_logprob, new_kernel = class_conditional(
      i, state.c, x, y, state.thetas, alpha, a
    )

    log_probs = collect(values(c_logprob))
    state_index = rand(Categorical(exp.(log_probs)))
    state.c[i] = collect(keys(c_logprob))[state_index]
    if !(state.c[i] in collect(keys(state.thetas)))
      state.thetas[state.c[i]] = new_kernel
    end
  end

  state
end

function add_kernels!(state::MixGPState,
                      new_cs::Vector{Int64},
                      a::GPHyper)
  for c in new_cs
    state.thetas[string(c)] = rand_kernel(a)
  end
end

function MixGPSampler(x::Matrix,
                      y::Vector,
                      alpha::Float64,
                      a::GPHyper,
                      out_path::String,
                      n_iter::Int64 = 20,
                      thin::Int64 = 20,
                      n_hmc::Int64 = 10,
                      L::Int64 = 5,
                      epsilon::Float64 = 0.005)
  if ispath(out_path)
    rm(out_path, recursive = true)
  end
  mkpath(out_path)

  ## initialize the sampling state
  n = length(y)
  state = MixGPState(
    ones(n),
    Dict(1 => rand_kernel(a))
  )

  for iter = 1:n_iter
    println("iter ", iter)
    if (iter - 1) % thin == 0
      write_state(iter, state, out_path)
    end

    ## resample all the cs
    sweep_indicators!(state, x, y, alpha, a)

    ## resample the kernel hyperparameters
    for k in collect(keys(state.thetas))
      if !any(state.c .== k)
        println("deleting state")
        delete!(state.thetas, k)
        print(k)
        print(keys(state.thetas))
        continue
      end

      theta0 = log.([state.thetas[k].l ^ 2, state.thetas[k].v0, state.thetas[k].v1])
      theta_samples = GPSampler(x[state.c .== k, :], y[state.c .== k], a, n_hmc, L, epsilon, theta0)
      state.thetas[k] = param_from_theta(
        [[sqrt(theta_samples[end, 1])];
         theta_samples[end, 2:3]]
      )
    end
  end

  state
end

function mix_posteriors(x_new::Matrix, state::MixGPState)
  post = Dict{Int64, Distributions.MvNormal}()
  for k in keys(state.thetas)
    gp = GPModel(state.thetas[k], x, y)
    post[k] = gp_posterior(x_new, gp)
  end

  post
end

function mix_posteriors(x_new::Matrix, states::Dict{Int64, MixGPState})
  posteriors = Dict{Int64, Dict{Int64, Distributions.MvNormal}}()
  for k in keys(states)
    println("Processing ", k)
    posteriors[k] = mix_posteriors(x_new, states[k])
  end

  posteriors
end

function write_posteriors(output_path::String,
                          x_new::Matrix,
                          post::Dict{Int64, Distributions.MvNormal})
  if isfile(output_path)
    rm(output_path)
  end

  open(output_path, "a") do x
    for k in keys(post)
      append_component(x, k, x_new, post[k])
    end
  end
end

function write_posteriors(output_path::String,
                          x_new::Vector,
                          posteriors::Dict{Int64, Dict{Int64, Distributions.MvNormal}})
  if isfile(output_path)
    rm(output_path)
  end

  open(output_path, "a") do x
    for i in keys(posteriors)
      for k in keys(posteriors[i])
        append_component(i, x, k, x_new, posteriors[i][k])
      end
    end
  end
end

function append_component(i::Int64,
                          x::IOStream,
                          k::Int64,
                          x_new::Vector,
                          distn::Distributions.MvNormal)
  ones_d = ones(length(distn))
  writecsv(x, [i * ones_d k * ones_d x_new mean(distn)])
end

"""Read States from File

This reads the thetas and cs written to file by write_state().
"""
function read_states(thetas_path::String, c_path::String)
  thetas_array = readcsv(thetas_path)
  c_array = readcsv(c_path)
  iters = unique(thetas_array[:, 1])
  states = Dict{Int64, MixGPState}()

  for iter in iters
    cur_thetas = thetas_array[thetas_array[:, 1] .== iter, :]
    cur_c = c_array[c_array[:, 1] .== iter, 3]

    ## construct kernel hyperparameter dictionary
    cur_param = Dict{Int64, KernelParam}()
    for j in 1:size(cur_thetas, 1)
      cur_param[j] = param_from_theta(cur_thetas[j, 3:5])
    end

    states[iter] = MixGPState(
      [Int64(ci) for ci in cur_c],
      cur_param
    )
  end

  states
end
