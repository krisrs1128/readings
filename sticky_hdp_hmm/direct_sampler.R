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

prior_predictive <- function(z_prev, z_next, n, alpha, beta, kappa) {
  K <- nrow(n)
  pzt <- vector(length = K + 1)

  ## equation (A10, page 216)
  for (k in seq_len(K)) {
    pzt[k] <- (alpha * beta[k] + n[z_prev, k] + kappa * (z_prev == k)) *
      (alpha * beta[z_next] + n[k, z_next] + kappa * (z_next == k) + (z_prev == k) * (z_next == k)) *
      (alpha + sum(n[k, ]) + kappa + (z_prev == k)) ^ (-1)
  }
  pzt[K + 1] <- (alpha ^ 2 * beta[K + 1] * beta[z_next]) / (alpha + kappa)

  pzt
}

update_f <- function(yt, pzt, emission) {
  K <- length(pzt) - 1
  for (k in seq_len(K + 1)) {
    f[k] <- f[k] * likelihood_reweight(yt, emission[k])
  }
  f
}

grow_beta <- function(beta, gamma) {
  K <- length(beta) - 1
  beta_new <- rbeta(1, gamma) * beta[K + 1]
  c(beta[-(K + 1)], beta_new, 1 - sum(beta[-(K + 1)]) - beta_new)
}

#' single term emission updates for the normal inverse wishart prior
#' see page 217 of Emily Fox's thesis
decrement_emission <- function(emission_k, yt) {
  emission_k$zeta <- emission_k$zeta - 1
  emission_k$nu <- emission_k$nu - 1
  old_zeta_theta <- emission_k$zeta_theta
  emission_k$zeta_theta <- emission_k$zeta_theta - yt
  emission_k$nu_delta <- emission_k$nu_delta - yt %*% t(yt) + old_zeta_theta -
    emission_k$zeta_theta %*% t(emission_k$zeta_theta) / emission_k$zeta

  emission_k
}

increment_emission <- function(emission_k, yt) {
  emission_k$zeta <- emission_k$zeta + 1
  emission_k$nu <- emission_k$nu + 1
  old_zeta_theta <- emission_k$zeta_theta
  emission_k$zeta_theta <- emission_k$zeta_theta + yt
  emission_k$nu_delta <- emission_k$nu_delta + yt %*% t(yt) + old_zeta_theta -
    emission_k$zeta_theta %*% t(emission_k$zeta_theta) / emission_k$zeta

  emission_k
}

likelihood_reweight <- function(y, emission_k) {
  theta <- (1 / emission_k$zeta) * emission_k$zeta_theta
  d <- nrow(emission_k$theta)
  nu_delta_coef <- (emission_k$zeta + 1) / (emission_k$zeta * (emission_k$nu - d - 1))
  dmvt(y, theta, nu_delta_coef * emission_k$nu_delta)
}
