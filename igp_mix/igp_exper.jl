#!/usr/bin/env julia

## File description -------------------------------------------------------------
##
## Experiment using the mixture of GPs sampler on data drawn from the assumed
## model.
##
## author: sankaran.kris@gmail.com
## date: 10/08/2017

## simulate toy data
include("igp_mix.jl")
srand(09142017)
n = 60
n_iter = 15
K = 3
c = rand(1:K, n)
update_ix = rand(1:n)
alpha = 1.0
thetas = Dict{Int64, KernelParam}()

a = GPHyper(
  Distributions.Logistic(-3, 2),
  Distributions.Logistic(0, 0.1),
  Distributions.Logistic(-1, 1)
)

for k = 1:K
  thetas[k] = rand_kernel(a)
end
c, x, y = simulate_mix(n, thetas)

## run the sampler, and keep last iteration
MixGPSampler(x, y, alpha, a, "data/sim0914/samples/", n_iter)

## get posterior estimates for each component
states = read_states(
  "data/sim0914/samples/thetas.csv",
  "data/sim0914/samples/c.csv"
)
x_new = collect(minimum(x):0.01:maximum(x))[:, :]
post = mix_posteriors(x_new, states)
write_posteriors("data/sim0914/post.csv", x_new, post)
writecsv("data/sim0914/data.csv", [c x y])

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
MixGPSampler(x, y, alpha, a, "data/bump0914/samples/", n_iter)

states = read_states(
  "data/bump0914/samples/thetas.csv",
  "data/bump0914/samples/c.csv"
)

x_new = collect(minimum(x):0.01:maximum(x))[:, :]
posteriors = mix_posteriors(x_new, states)
write_posteriors("data/bump0914/posteriors.csv", x_new, posteriors)
writecsv("data/bump0914/data.csv", [c x y])

## consider real microbiome series
alpha = 0.1
y = readcsv("data/unc063x1/raw_data.csv")[:]
y += 0.01 * rand(length(y))
x = collect(linspace(0, 1, length(y)))[:, :]
MixGPSampler(x, y, alpha, a, "data/unc063x1/samples/", n_iter, 2)

states = read_states(
  "data/unc063x1/samples/thetas.csv",
  "data/unc063x1/samples/c.csv"
)

x_new = collect(minimum(x):0.005:maximum(x))[:, :]
posteriors = mix_posteriors(x_new, states)
write_posteriors("data/unc063x1/posteriors.csv", x_new, posteriors)
writecsv("data/unc063x1/data.csv", [zeros(length(y)) x y])
