#!/usr/bin/env julia

## File description -------------------------------------------------------------
##
## Experiment using the mixture of GPs sampler on data drawn from the assumed
## model.
##
## author: sankaran.kris@gmail.com
## date: 09/13/2017


## simulate toy data
srand(09142017)
n = 100
K = 4
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
state = MixGPSampler(x, y, alpha, a, "data/samples/sim0914/", 500)

## get posterior estimates for each component
x_new = collect(minimum(x):0.01:maximum(x))[:, :]
post = mix_posteriors(x_new, state)
write_posteriors("data/mix_post.csv", post)
writecsv("data/mix_data.csv", [c x y])

## consider instead a zeros + departures dataset
x = rand(150, 1)
function f(x)
  if x < 0.2
    return 0
  elseif x < 0.4
    return -(x - 0.4) * (x - 0.2)
  elseif x < 0.8
    return 0
  elseif x <= 1
    return -(x - 0.8) * (x - 1.2)
  end
end

y = [f(z) for z in x[:, 1]]
y += 0.01 * rand(length(y))
y -= mean(y)

state = MixGPSampler(x, y, alpha, a, 100)

x_new = collect(minimum(x):0.01:maximum(x))[:, :]
post = mix_posteriors(x_new, state)
write_posteriors("data/bump_posteriors.csv", post)
writecsv("data/bump_data.csv", [x y])
