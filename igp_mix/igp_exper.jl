#!/usr/bin/env julia

## File description -------------------------------------------------------------
##
## Experiment using the mixture of GPs sampler on data drawn from the assumed
## model.
##
## author: sankaran.kris@gmail.com
## date: 09/13/2017


## simulate toy data
n = 50
K = 3
c = rand(1:K, n)
update_ix = rand(1:n)
alpha = 2.0
thetas = Dict{Int64, KernelParam}()

a = GPHyper(
  Distributions.Logistic(-1, 4),
  Distributions.Logistic(-1, 4),
  Distributions.Logistic(-1, 4)
)

for k = 1:K
  thetas[k] = rand_kernel(a)
end
c, x, y = simulate_mix(n, thetas)

## run the sampler, and keep last iteration
samples = MixGPSampler(x, y, alpha, a, 10)

## get posterior estimates for each component
x_new = collect(minimum(x):0.1:maximum(x))[:, :]
post = Dict{Int64, Distributions.MvNormal}()
for k = 1:maximum(samples.c)
  if sum(samples.c .== k) == 0
    continue
  end

  gp = GPModel(samples.thetas[k], x, y)
  post[k] = gp_posterior(x_new, gp)
end

## write posteriors to file
open("data/post_mix.csv", "a") do x
  for k in keys(post)
    writedlm(
      x,
      [k * ones(length(post[k])) mean(post[k])]
    )
  end
end
