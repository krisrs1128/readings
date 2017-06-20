#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulation to use when checking sticky hdp hmm sampler.
##
## author: sankaran.kris@gmail.com
## date: June 19, 2017

## ---- libraries ----
library("tidyverse")


## ---- utils ----
tidy_melt <- function(x, row_var = "row", col_var = "col", value_name = "value") {
  x <- x %>%
    as_data_frame() %>%
    rownames_to_column("row") %>%
    gather(col, value, -row) %>%
    mutate(
      row = as.numeric(gsub("V", "", row)),
      col = as.numeric(gsub("V", "", col))
    )
  colnames(x) <- c(row_var, col_var, value_name)
  x
}

#' @examples
#' pi_df <- t(sapply(1:1000, function(i) sbp(10, 50))) %>%
#'   tidy_melt("i", "k")
#' ggplot(pi_df) +
#'   geom_point(aes(x = k, y = value), alpha = 0.05, size = 0.3)
sbp <- function(alpha, K) {
  stick_remaining <- 1
  pi <- vector(length = K + 1)
  for (k in seq_len(K)) {
    pi[k] <- rbeta(1, 1, alpha) * stick_remaining
    stick_remaining <- stick_remaining - pi[k]
  }
  pi[K + 1] <- stick_remaining

  pi
}

#' @examples
#' beta0 <- sbp(1, 10) ## weights associated with draw from G0
#' pi_j <- sticky_hdp_masses(5, 0.5, 3, beta0)
sticky_hdp_given_G0 <- function(alpha, kappa, j, beta0) {
  K <- length(beta0) - 1
  if (j > K) {
    stop("sticky transition index must be smaller than K")
  }
  beta1 <- sbp(alpha + kappa, K)
  pi_j <- vector(length = K + 1)

  for (k in seq_len(K)) {
    sticky_indic <- runif(1) < kappa / (alpha + kappa)
    if (sticky_indic) {
      pi_j[j] <- pi_j[j] + beta1[k]
    } else {
      cur_atom <- sample(seq_len(K + 1), 1, prob = beta0) ## a draw from G_{0}
      pi_j[cur_atom] <- pi_j[cur_atom] + beta1[k]
    }
  }

  pi_j
}

#' @examples
#' P <- trans_mat(3, 3, 15, 1.5)
trans_mat <- function(gamma, alpha, K, kappa) {
  P <- matrix(0, K + 1, K + 1)
  beta0 <- sbp(gamma, K) ## weights from G0
  for (k in seq_len(K)) {
    P[k, ] <- sticky_hdp_given_G0(alpha, kappa, k, beta0)
  }
  P
}

#' @examples
#' P <- trans_mat(3, 3, 15, 1.5)
#' z <- markov_chain(P)
markov_chain <- function(P, time_len = 100) {
  ## start in a random state
  K <- nrow(P) - 1
  z <- c(sample(seq_len(K), 1), vector(length = time_len - 1))
  for (i in seq_len(time_len - 1)) {
    if (z[i] <= K) {
      z[i + 1] <- sample(seq_len(K + 1), size = 1, prob = P[z[i], ])
    } else {
      z[i + 1] <- z[1] ## restart chain when escaped truncation level
      warning("escaped state space 1 <= k <= K, restarting from z[0]")
    }
  }
  z
}

#' See definition on page 54 of Emily Fox's thesis
#' @examples
#' rnorm_inv_wishart(2, c(0, 0), 3, diag(c(.1, .1)))
rnorm_inv_wishart <- function(zeta, theta, nu, Delta) {
  Sigma <- rWishart(1, nu, Delta)[,, 1]
  mu <- rmvnorm(1, theta, (1 / zeta) * Sigma)
  list("mu" = mu, "Sigma" = Sigma)
}

emission_parameters <- function(K, lambda) {
  emissions <- list()
  for (k in seq_len(K)) {
    emissions[[k]] <- rnorm_inv_wishart(lambda$zeta, lambda$theta, lambda$nu, lambda$Delta)
  }
  emissions
}

#' @examples
#' K <- 50
#' z <- markov_chain(trans_mat(3, 3, K, 2.5))
#' lambda <- list(zeta = 2, theta = c(0, 0), nu = 3, Delta = diag(c(0.01, 0.01)))
#' theta <- emission_parameters(K, lambda)
#' y <- emissions(z, theta)
#' plot(y, col = z, asp = 1)
#' lines(y, col = "grey")
#' points(y, col = z)
#' plot(y[, 1], col = z)
emissions <- function(z, theta) {
  time_len <- length(z)
  K <- length(theta)
  y <- matrix(nrow = time_len, ncol = length(theta[[1]]$mu))
  for (k in seq_len(K)) {
    n_in_state <- sum(z == k)
    if (n_in_state > 0) {
      y[z == k, ] <- rmvnorm(n_in_state, theta[[k]]$mu, theta[[k]]$Sigma)
    }
  }
  y
}
