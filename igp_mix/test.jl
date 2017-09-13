#!/usr/bin/env julia

# File description -------------------------------------------------------------
# Test functions in hmc_exper.jl. In particular, check that the analytical
# derivatives required for Hamiltonian Monte Carlo agree with the finite
# difference approximation.
#
# author: sankaran.kris@gmail.com
# date: 08/29/2018

using Base.Test
using Calculus


"""Coordinatewise Projections of Kernel

This returns a dictionary of functions that evaluate the kernel varying one
argument at a time.

# Arguments

x::Matrix: The data on which the kernel is evaluated
theta::KernelParam: The bandwidth / noise parameters of the kernel
i::Int64 = 1: The row-coordinate of the kernel matrix to return
j::Int64 = 2: The column-coordinate of the kernel matrix to return
"""
function k_coords(x::Matrix,
                  theta::KernelParam,
                  i::Int64 = 1,
                  j::Int64 = 2)
  f_list = Dict{String, Function}()

  f_list["l"] = function(log_l2::Float64)
    theta.l = sqrt(exp(log_l2))
    kernel(x, x, theta)[i, j]
  end

  f_list["v0"] = function(log_v0::Float64)
    theta.v0 = exp(log_v0)
    kernel(x, x, theta)[i, j]
  end

  f_list["v1"] = function(log_v1::Float64)
    theta.v1 = exp(log_v1)
    kernel(x, x, theta)[i, j]
  end

  f_list
end


"""Analytical Derivatives of Kernel

Return a dictionary of the derivatives of the kernel function, one parameter at
a time.

# Arguments
x::Matrix: The data on which the kernel derivative is evaluated
theta::KernelParam: The bandwidth / noise parameters of the kernel
i::Int64 = 1: The row-coordinate of the kernel derivative matrix to return
j::Int64 = 2: The column-coordinate of the kernel derivative matrix to return
"""
function k_deriv_coords(x::Matrix,
                        theta::KernelParam,
                        i::Int64 = 1,
                        j::Int64 = 2)
  f_list = Dict{String, Function}()

  f_list["log_l2"] = function(log_l2::Float64)
    theta.l = sqrt(exp(log_l2))
    k_theta_deriv = k_deriv_logl2(x, theta)
    k_theta_deriv[i, j]
  end

  f_list["log_v0"] = function(log_v0::Float64)
    theta.v0 = exp(log_v0)
    k_theta_deriv = k_deriv_logv0(x, theta)
    k_theta_deriv[i, j]
  end

  f_list["log_v1"] = function(log_v1::Float64)
    k_theta_deriv = exp(log_v1) * eye(size(x, 1))
    k_theta_deriv[i, j]
  end

  f_list
end

a = GPHyper(
  Distributions.Logistic(-1, 4),
  Distributions.Logistic(-1, 4),
  Distributions.Logistic(-1, 4)
)

theta = rand_kernel(a)
x, y = simulate(100, theta)
coord_f = k_coords(x, theta)
deriv_f = k_deriv_coords(x, theta)

@testset "Derivative tests" begin
  for u in rand(Uniform(-10, 10), 10)
    @test check_derivative(coord_f["l"], deriv_f["log_l2"], u) ≈ 0 atol = 1e-5
    @test check_derivative(coord_f["v0"], deriv_f["log_v0"], u) ≈ 0 atol = 1e-5
    @test check_derivative(coord_f["v1"], deriv_f["log_v1"], u) ≈ 0 atol = 1e-5
  end
end

@testset "Gibbs conditionals" begin
  for i = 1:10

    ## simulate toy data
    n = 100
    K = 3
    c = rand(1:K, n)
    update_ix = rand(1:n)
    alpha = 2.0
    thetas = Dict{Int64, KernelParam}()

    for k = 1:K
      thetas[k] = rand_kernel(a)
    end
    c, x, y = simulate_mix(n, thetas)

    ## swap one element
    c_prime = deepcopy(c)
    c_prime[update_ix] = rand(1:K)
    joint_diff = joint_log_prob(c, x, y, thetas, alpha) -
      joint_log_prob(c_prime, x, y, thetas, alpha)
    conditional_probs = class_conditional(update_ix, c, x, y, thetas, alpha, a)[1]
    conditional_diff = conditional_probs[c[update_ix]] -
      conditional_probs[c_prime[update_ix]]
    @test joint_diff - conditional_diff ≈ 0 atol = 1e-10
  end
end
