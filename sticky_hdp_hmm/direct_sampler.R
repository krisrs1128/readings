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

K <- 20
gamma <- 2
alpha <- 2
kappa <- 4

z <- markov_chain(trans_mat(gamma, gamma, K, kappa))
lambda <- list(zeta = 2, theta = c(0, 0), nu = 10, Delta = diag(c(0.01, 0.01)))
theta <- emission_parameters(K, lambda)
y <- emissions(z, theta)

z <- rep(1, length(y))
emission <- lapply(1:2, function(i) emission_prior(lambda))
beta <- c(0.1, 0.9)
#sample_z(y, rep(1, length(y)), emission, alpha, beta, kappa, gamma, lambda)

sample_z <- function(y, z, emission, alpha, beta, kappa, gamma, lambda) {
  n <- transition_counts(z)
  time_len <- length(z)
  for (i in seq(2, time_len - 1)) {
    K <- max(z)

    ## decrements
    n <- update_count(n, z[i - 1], z[i], z[i + 1], -1)
    emission[[z[i]]] <- update_emission(emission[[z[i]]], y[i, ], -1)

    ## update according to prior predictive * likelihood
    pzt <- prior_predictive(z[i - 1], z[i + 1], n, alpha, beta, kappa)
    f <- update_f(y[i, ], pzt, emission)
    z[i] <- sample(seq_along(f), 1, prob = f)

    ## grow K, if sample new state
    if (max(z) > K) {
      beta <- grow_beta(beta, gamma)
      emission <- c(emission, list(emission_prior(lambda)))
    }

    ## increment
    n <- update_count(n, z[i - 1], z[i], z[i + 1], 1)
    emission[[z[i]]] <- update_emission(emission[[z[i]]], y[i, ], 1)
  }

  list("z" = z, "beta" = beta, "emission" = emission)
}

#' @examples
#' z <- c(1, 1, 2, 1, 1, 3, 3, 3, 1)
#' transition_counts(z)
transition_counts <- function(z) {
  K <- max(z)
  n <- matrix(0, nrow = K + 1, ncol = K + 1)
  time_len <- length(z)

  for (i in seq_len(time_len - 1)) {
    n[z[i], z[i + 1]] <- n[z[i], z[i + 1]] + 1
  }
  n
}

update_count <- function(n, z_prev, z_cur, z_next, s) {
  n[z_prev, z_cur] <- n[z_prev, z_cur] + s
  n[z_cur, z_next] <- n[z_cur, z_next] + s
  n
}

prior_predictive <- function(z_prev, z_next, n, alpha, beta, kappa) {
  K <- nrow(n) - 1
  pzt <- vector(length = K + 1)

  ## equation (A10, page 216)
  for (k in seq_len(K)) {
    pzt[k] <- (alpha * beta[k] + n[z_prev, k] + kappa * (z_prev == k)) *
      (alpha * beta[z_next] + n[k, z_next] + kappa * (z_next == k) + (z_prev == k) * (z_next == k)) *
      (alpha + sum(n[k, ]) + kappa + (z_prev == k)) ^ (-1)
  }
  pzt[K + 1] <- (alpha ^ 2 * beta[K + 1] * beta[z_next]) / (alpha + kappa)

  pzt / sum(pzt)
}

update_f <- function(yt, pzt, emission) {
  K <- length(pzt) - 1
  f <- vector(length = K + 1)
  for (k in seq_len(K + 1)) {
    f[k] <- pzt[k] * likelihood_reweight(yt, emission[[k]])
  }
  f / sum(f)
}

grow_beta <- function(beta, gamma) {
  K <- length(beta) - 1
  beta_new <- rbeta(1, 1, gamma) * beta[K + 1]
  c(beta[-(K + 1)], beta_new, 1 - sum(beta[-(K + 1)]) - beta_new)
}

#' single term emission updates for the normal inverse wishart prior
#' see page 217 of Emily Fox's thesis
update_emission <- function(emission_k, yt, s) {
  emission_k$zeta <- emission_k$zeta + s
  emission_k$nu <- emission_k$nu + s
  old_zeta_theta <- emission_k$zeta_theta
  emission_k$zeta_theta <- emission_k$zeta_theta + s * yt
  emission_k$nu_delta <- emission_k$nu_delta + s * yt %*% t(yt) + old_zeta_theta -
    emission_k$zeta_theta %*% t(emission_k$zeta_theta) / emission_k$zeta

  emission_k
}

emission_prior <- function(lambda) {
  list(
    "zeta" = lambda$zeta,
    "nu" = lambda$nu,
    "zeta_theta" = lambda$zeta * lambda$theta,
    "nu_delta" = lambda$nu * lambda$Delta
  )
}

likelihood_reweight <- function(y, emission_k) {
  theta <- (1 / emission_k$zeta) * emission_k$zeta_theta
  d <- length(y)
  nu_delta_coef <- (emission_k$zeta + 1) / (emission_k$zeta * (emission_k$nu - d - 1))
  dmvt(y, theta, nu_delta_coef * emission_k$nu_delta, log = FALSE)
}
