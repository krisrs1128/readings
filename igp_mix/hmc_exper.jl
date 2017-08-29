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

using Distributions
include("igp_mix.jl")

theta = KernelParam(sqrt(0.05), 10.0, 0.5)
x, y = simulate(70, theta)

a = GPHyper(
  Distributions.TDist(1),
  Distributions.TDist(1),
  Distributions.TDist(1)
)
samples = GPSampler(x, y, a, 20, 5, 0.001, [-1, 2, -1])

x_new = collect(0:0.01:1)[:, :]
post_data = zeros(0, 4)
for i in 1:size(samples, 1)
  if i % 50 != 0
    next
  end

  theta_fit = param_from_theta(samples[i, :])
  post = gp_posterior(x_new, GPModel(theta_fit, x, y))
  post_data = [
    post_data;
    i * ones(size(x_new, 1)) x_new mean(post) var(post)
  ]
end

writedlm("train_data.csv", [x y])
writedlm("samples.csv", samples)
writedlm("post_data.csv", post_data)
