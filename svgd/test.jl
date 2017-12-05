#!/usr/bin/env julia

# File description -------------------------------------------------------------
#
# Some tests.
#
# author: sankaran.kris@gmail.com
# date:

include("normal_mixture.jl")
using Calculus

function log_pdf(x::Float64)
  p1 = Normal()
  p2 = Normal(10.0, 3.0)
  log(pdf(p1, x) + pdf(p2, x))
end

check_derivative(log_pdf, grad_logp, -10.0)
