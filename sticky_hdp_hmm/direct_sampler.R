#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Code for the direct assignment collapsed gibbs sampler for the sticky
## HDP-HMM, as described in Algorithm 8 of Emily Fox's thesis
##
## author: Kris Sankaran
## date: 06/19/2017

## ---- libraries ----
library("mvtnorm")

## ---- utils ----
#' @examples
#' z <- c(1, 1, 2, 1, 1, 3, 3, 3, 1)
#' transition_counts(z)
transition_counts <- function(z) {
  K <- max(z)
  n <- matrix(0, nrow = K, ncol = K)
  time_len <- length(z)

  for (i in seq_len(time_len - 1)) {
    n[z[i], z[i + 1]] <- n[z[i], z[i + 1]] + 1
  }
  n
}

decrement_count <- function(n, z_prev, z_cur, z_next) {
  n[z_prev, z_cur] <- n[z_prev, z_cur] - 1
  n[z_cur, z_next] <- n[z_cur, z_next] - 1
}

increment_count <- function(n, z_prev, z_cur, z_next) {
  n[z_prev, z_cur] <- n[z_prev, z_cur] + 1
  n[z_cur, z_next] <- n[z_cur, z_next] + 1
}

f_vector <- function(yt, z_prev, z_next, n, alpha, beta, kappa, theta) {
  K <- nrow(n)
  f <- vector(length = K)

  for (k in seq_len(K)) {
    f[k] <- (alpha * beta[k] + n[k, z_next] + kappa * (z_next == k)) / (alpha + sum(n[k, ]) + kappa)
    f[k] <- f[k] * dmvt(yt, theta[k]$mu, theta[k]$sigma_sq, theta[k]$nu)
  }
  f
}
