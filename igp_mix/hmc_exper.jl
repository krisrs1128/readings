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
  Distributions.Logistic(-1, 4),
  Distributions.Logistic(-1, 4),
  Distributions.Logistic(-1, 4)
)
samples = GPSampler(x, y, a, 5000, 5, 0.005, zeros(3))

x_new = collect(0:0.01:1)[:, :]
post_data = zeros(0, 4)
for i in 1:size(samples, 1)
  if (i - 1) % 50 != 0
    continue
  end

  println(i)
  theta_fit = param_from_theta(samples[i, :])
  post = gp_posterior(x_new, GPModel(theta_fit, x, y))
  post_data = [
    post_data;
    i * ones(size(x_new, 1)) x_new mean(post) var(post)
  ]
end

writedlm("data/train_data.tsv", [x y])
writedlm("data/samples.tsv", samples)
writedlm("data/post_data.tsv", post_data)
