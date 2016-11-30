#! /usr/bin/env Rscript

# File description -------------------------------------------------------------

## ---- setup ----
source("./dtm.R")
library("rstan")

## ---- generate_data ----
N <- 100
V <- 15
sim_data <- dtm_data(N, V, rep(1000, N), sigma = 0.1, alpha = rep(.01, 3))

## ---- fit_model ----
stan_data <- list(
  N = N,
  V = V,
  T = N,
  K = 5,
  sigma = 0.1,
  delta = 0.1,
  times = 1:N,
  times_mapping = 1:N,
  X = sim_data$X
)

#m <- stan_model("./dtm.stan")
stan_fit <- vb(m, data = stan_data)
samples <- extract(stan_fit)

alpha_hat <- apply(samples$alpha, c(2, 3), mean)
pairs(cbind(sim_data$theta, t(apply(alpha_hat, 1, softmax))))

mu_hat <- apply(samples$mu, c(2, 3, 4), mean)
beta_hat <- apply(samples$beta, c(2, 3, 4), mean)
