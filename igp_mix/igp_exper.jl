#!/usr/bin/env julia

## File description -------------------------------------------------------------
##
## Experiment using the mixture of GPs sampler on data drawn from the assumed
## model.
##
## author: sankaran.kris@gmail.com
## date: 09/13/2017


## simulate toy data
n = 150
K = 3
c = rand(1:K, n)
update_ix = rand(1:n)
alpha = 2.0
thetas = Dict{Int64, KernelParam}()

a = GPHyper(
  Distributions.Normal(),
  Distributions.Normal(),
  Distributions.Normal()
)

for k = 1:K
  thetas[k] = rand_kernel(a)
end
c, x, y = simulate_mix(n, thetas)

## run the sampler, and keep last iteration
state = MixGPSampler(x, y, alpha, a, 30)

## get posterior estimates for each component
x_new = collect(minimum(x):0.01:maximum(x))[:, :]
post = mix_posteriors(x_new, state)
write_post("data/mix_post.csv", post)
writedlm("data/mix_data.csv", [x y])

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
writedlm("data/bump_data.csv", [x y])
