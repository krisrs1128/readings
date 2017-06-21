#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Code for the direct assignment collapsed gibbs sampler for the sticky
## HDP-HMM, as described in Algorithm 8 of Emily Fox's thesis
##
## author: Kris Sankaran
## date: 06/19/2017

## ---- libraries ----
library("mvtnorm")
library("jsonlite")

## ---- utils ----
#' @examples
#' K <- 20
#' gamma <- 2
#' alpha <- 2
#' kappa <- 3
#'
#' z <- markov_chain(trans_mat(gamma, gamma, K, kappa))
#' lambda <- list(zeta = .1, theta = c(0, 0), nu = 10, Delta = diag(c(0.1, 0.1)))
#' theta <- emission_parameters(K, lambda)
#' y <- emissions(z, theta)
#' plot(y[, 1], col = z)
#'
#' direct_sampler(y, alpha, kappa, gamma, lambda)
#'
#' emission <- lapply(1:2, function(i) emission_prior(lambda))
#' emission[[1]]$zeta <- emission[[2]]$zeta + nrow(y)
#' emission[[1]]$nu <- emission[[2]]$nu + nrow(y)
#' beta <- c(0.1, 0.9)
#' sample_z(y, rep(1, nrow(y)), emission, alpha, beta, kappa, gamma, lambda)
direct_sampler <- function(y, alpha, kappa, gamma, lambda, n_iter = 1000) {
  ## initialize markov chain state space
  z <- rep(1, nrow(y))
  emission <- lapply(1:2, function(i) emission_prior(lambda))
  emission[[1]]$zeta <- emission[[2]]$zeta + nrow(y)
  emission[[1]]$nu <- emission[[2]]$nu + nrow(y)
  beta <- c(0.1, 0.1)

  for (i in seq_len(n_iter)) {
    ## if (i %% 50 == 0) {
      cat(sprintf("iteration %s\n", i))
    ## }

    ## step 1: sample the state sequence
    z_updates <- sample_z(y, z, emission, alpha, beta, kappa, gamma, lambda)
    z <- z_updates$z
    beta <- z_updates$beta
    emission <- z_updates$emission

    ## step 2: delete unused modes
    deleted_modes <- delete_unused_modes(z, beta, emission)
    z <- deleted_modes$z
    emission <- deleted_modes$emission
    beta <- deleted_modes$beta

    ## step 3: sample auxiliary varaibles
    m <- sample_m(z, alpha, beta, kappa)
    w <- sample_override(diag(m), kappa / (kappa + alpha), beta)

    ## step 4: sample global transition distribution
    beta <- sample_beta(m, w, gamma)

    state <- list("z" = z, "beta" = beta, "m" = m, "w" = w, "emission" = emission)
    write_state("sampler_data", state, i)
  }

  state
}

sample_z <- function(y, z, emission, alpha, beta, kappa, gamma, lambda) {
  time_len <- length(z)
  for (i in seq(2, time_len - 1)) {
    K <- length(emission) - 1
    n <- transition_counts(z)

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
  modes <- c(unique(z), "new")
  K <- length(modes) - 1
  n <- matrix(0, nrow = K + 1, ncol = K + 1, dimnames = list(modes, modes))
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
  modes <- as.numeric(names(beta[1:K]))
  beta_new <- rbeta(1, 1, gamma) * beta["new"]

  beta_update <- c(beta[1:K], beta_new, 1 - sum(beta[1:K]) - beta_new)
  names(beta_update) <- c(modes, max(modes) + 1, "new")
  beta_update
}

#' single term emission updates for the normal inverse wishart prior
#' see page 217 of Emily Fox's thesis
update_emission <- function(emission_k, yt, s) {
  old_zeta_theta <- emission_k$zeta_theta
  old_zeta <- emission_k$zeta

  emission_k$zeta <- emission_k$zeta + s
  emission_k$nu <- emission_k$nu + s
  emission_k$zeta_theta <- emission_k$zeta_theta + s * yt
  emission_k$nu_delta <- emission_k$nu_delta + s * yt %*% t(yt) + old_zeta_theta %*% t(old_zeta_theta) / old_zeta-
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
  t_dens <- dmvt(y, theta, nu_delta_coef * emission_k$nu_delta, log = FALSE)
  ifelse(is.na(t_dens), 0, t_dens)
}

#' @examples
#' z <- c(1, 2, 3, 1, 5)
#' beta <- c(.1, .1, .2, .1, .1, .4)
#' emission <- rep(list("emission_params"), 6)
#' names(emission) <- as.character(1:6)
#' delete_unused_modes(z, beta, emission)
delete_unused_modes <- function(z, beta, emission) {
  z_fact <- factor(z, levels = names(emission))
  z_counts <- table(z_fact)
  unused_modes <- names(z_counts)[z_counts == 0]

  list(
    "beta" = beta[setdiff(names(beta), unused_modes)]
    "emission" = emission[setdiff(names(emission), unused_modes)]
  )
}

#' From MCMCpack library
rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x / as.vector(sm)
}

sample_beta <- function(m, w, gamma) {
  m_bar <- m
  diag(m_bar) <- diag(m_bar) - w
  K <- length(w)
  rdirichlet(1, c(colSums(m_bar), gamma))[1, ]
}
