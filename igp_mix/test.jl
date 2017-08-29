#!/usr/bin/env julia

# File description -------------------------------------------------------------
# Test functions in hmc_exper.jl. In particular, check that the analytical
# derivatives required for Hamiltonian Monte Carlo agree with the finite
# difference approximation.
#
# author: sankaran.kris@gmail.com
# date: 08/29/2018

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

coord_f = k_coords(x, theta)
deriv_f = k_deriv_coords(x, theta)

for u in rand(Uniform(-10, 10), 10)
  println(check_derivative(coord_f["l"], deriv_f["log_l2"], u))
  println(check_derivative(coord_f["v0"], deriv_f["log_v0"], u))
  println(check_derivative(coord_f["v1"], deriv_f["log_v1"], u))
end
