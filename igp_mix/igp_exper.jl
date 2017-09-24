#!/usr/bin/env julia

## File description -------------------------------------------------------------
##
## Experiment using the mixture of GPs sampler on data drawn from the assumed
## model.
##
## author: sankaran.kris@gmail.com
## date: 09/15/2017

## simulate toy data
include("igp_mix.jl")
srand(09142017)
n = 60
n_iter = 10000
K = 3
c = rand(1:K, n)
update_ix = rand(1:n)
alpha = 1.0
thetas = Dict{Int64, KernelParam}()

a = GPHyper(
  Distributions.Logistic(-1, 1),
  Distributions.Logistic(-1, 1),
  Distributions.Logistic(-1, 1)
)

for k = 1:K
  thetas[k] = rand_kernel(a)
end
c, x, y = simulate_mix(n, thetas)

## run the sampler, and keep last iteration
state = MixGPSampler(x, y, alpha, a, "data/samples/sim0914/", n_iter)

## get posterior estimates for each component
x_new = collect(minimum(x):0.01:maximum(x))[:, :]
post = mix_posteriors(x_new, state)
write_posteriors("data/mix_post.csv", post)
writecsv("data/mix_data.csv", [c x y])

## consider instead a zeros + departures dataset
x = rand(n, 1)
function f(x)
  if x < 0.2
    return 1, 0
  elseif x < 0.4
    return 2, -(x - 0.4) * (x - 0.2)
  elseif x < 0.8
    return 3, 0
  elseif x <= 1
    return 4, -(x - 0.8) * (x - 1.2)
  end
end

## prepare data
c = [f(z)[1] for z in x[:, 1]]
y = [f(z)[2] for z in x[:, 1]]
y += 0.01 * rand(length(y))
y -= mean(y)

## sample and write to file
alpha = 3.0
state = MixGPSampler(x, y, alpha, a, "data/samples/bump0914/", n_iter)

x_new = collect(minimum(x):0.01:maximum(x))[:, :]
post = mix_posteriors(x_new, state)
write_posteriors("data/bump_posteriors.csv", post)
writecsv("data/bump_data.csv", [c x y])

## consider real microbiome series
alpha = 1.0
y = readcsv("data/unc063x1_data.csv")[:]
x = collect(linspace(0, 1, length(y)))[:, :]
MixGPSampler(x, y, alpha, a, "data/samples/unc063x1/", n_iter)

states = read_states(
  "data/samples/unc063x1/thetas.csv",
  "data/samples/unc063x1/c.csv"
)

x_new = collect(minimum(x):0.01:maximum(x))[:, :]
posteriors = mix_posteriors(x_new, states)
write_posteriors("data/unc063x1_posteriors.csv", x_new[:], posteriors)
writecsv("data/unc063x1_data.csv", [zeros(length(y)) x y])
