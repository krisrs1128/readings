#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Simulate and estimate simple PPCA model.
##
## author: sankaran.kris@gmail.com
## date: 2/18/2018

## libraries
library("rstan")

## simulate data
N <- 500
p <- 10
k <- 2
sigma <- 0.2

W <- matrix(rnorm(p * k), p, k)
z <- matrix(rnorm(N * k), N, k)

x <- matrix(nrow = N, ncol = p)
for (i in seq_len(N)) {
  x[i, ] <- W %*% z[i, ] + rnorm(p, 0, sigma)
}

## prepare stan model
stan_data <- list("N" = N, "p" = p, "k" = k, "x" = x)
fit <- vb(stan_model("ppca.stan"), stan_data)

fit_data <- extract(fit)
dim(fit_data$z)

plot(-fit_data$z[, 2:1, 1])
points(z[1:2, ], col = "red")

w_means <- apply(fit_data$w, c(2,3 ), mean)
w_means
W
plot(w_means[, 1], W[, 2])
plot(w_means[, 1], W[, 2])
pairs(cbind(w_means, W))
