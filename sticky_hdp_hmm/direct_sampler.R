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
K <- 20
gamma <- 2
alpha <- 2
kappa <- 3

z <- markov_chain(trans_mat(gamma, gamma, K, kappa))
lambda <- list(zeta = .1, theta = c(0, 0), nu = 10, Delta = diag(c(0.1, 0.1)))
theta <- emission_parameters(K, lambda)
y <- emissions(z, theta)
plot(y[, 1], col = z)

direct_sampler(y, alpha, kappa, gamma, lambda)

## emission <- lapply(1:2, function(i) emission_prior(lambda))
## emission[[1]]$zeta <- emission[[2]]$zeta + nrow(y)
## emission[[1]]$nu <- emission[[2]]$nu + nrow(y)
## beta <- c(0.1, 0.9)
## sample_z(y, rep(1, nrow(y)), emission, alpha, beta, kappa, gamma, lambda)

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

    if (any(is.na(f))) {
      f <- rep(1 / (K + 1), K + 1) ## somehow this is broken...
    }
    print(f)
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
#' delete_unused_modes(z, beta, emission)
delete_unused_modes <- function(z, beta, emission) {
  K <- length(emission) - 1
  mapping <- cbind(sort(c(unique(z), K + 1)), seq_along(c(K + 1, unique(z))))
  z_new <- vector(length = length(z))
  for (k in seq_len(K + 1)) {
    if (k != K + 1 && !(k %in% mapping[, 1])) {
      emission[[k]] <- "delete_flag"
      beta[K + 1] <- beta[K + 1] + beta[k]
      beta[k] <- NA
    }

    z_new[z == k] <- mapping[mapping[, 1] == k, 2]
  }

  list(
    "z" = z,
    "beta" = beta[!is.na(beta)],
    "emission" = emission[emission != "delete_flag"]
  )
}

sample_m_coord <- function(n_jk, alpha, beta_k, kappa, k, j) {
  m_jk <- 0
  for (n in seq_len(n_jk)) {
    p <- (alpha * beta_k + kappa * (k == j)) /
      (n + alpha * beta_k + kappa * (k == j))
    if (runif(1) < p) {
      m_jk <- m_jk + 1
    }
  }
  m_jk
}

sample_m <- function(z, alpha, beta, kappa) {
  K <- length(beta) - 1
  m <- matrix(0, K, K)
  n <- transition_counts(z)

  for (j in seq_len(K)) {
    for (k in seq_len(K)) {
      m[j, k] <- sample_m_coord(
        n[j, k], alpha, beta[k], kappa, k, j
      )
    }
  }

  m
}

sample_override <- function(m_diag, rho, beta) {
  K <- length(m_diag)
  w <- vector(length = K)
  for (k in seq_len(K)) {
    w[k] <- rbinom(1, m_diag[k], rho / (rho + beta[k] * (1 - rho)))
  }

  w
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

append_to_file <- function(name, x) {
  write.table(
    x, name,
    append = TRUE,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(name)
  )
}

write_state <- function(out_dir, state, iter) {
  dir.create(out_dir, recursive = TRUE)
  z_mat <- cbind(iter, t(as.matrix(state$z)))
  ## colnames(z_mat) <- c("iter", paste0("time_", seq_along(state$z)))
  append_to_file(file.path(out_dir, "z.csv"), z_mat)

  beta_mat <- cbind(iter, t(as.matrix(state$beta)))
  ## colnames(z_mat) <- c("iter", paste0("k_", seq_along(state$beta)))
  append_to_file(file.path(out_dir, "beta.csv"), beta_mat)

  w_mat <- cbind(iter, t(as.matrix(state$w)))
  ## colnames(w_mat) <- c("iter", paste0("k_", seq_along(state$w)))
  append_to_file(file.path(out_dir, "w.csv"), w_mat)

  em <- state$emission
  em$iter <- iter
  cat(toJSON(em), file = file.path(out_dir, "emissions.json"), append = TRUE)

  m_mat <- cbind(iter, t(as.matrix(state$m)))
  ## colnames(m_mat) <- c("iter", paste0("k_", seq_len(ncol(m_mat))))
  append_to_file(file.path(out_dir, "m.csv"), m_mat)
}
