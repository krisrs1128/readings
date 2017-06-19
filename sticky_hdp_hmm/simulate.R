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
#' P <- trans_mat(10, 50)
trans_mat <- function(alpha, K) {
  P <- matrix(0, K, K)
  for (k in seq_len(K)) {
    P[k, ] <- sbp(alpha, K)
  }
  P
}

P <- trans_mat(1, 10)
markov_chain <- function(P, time_len = 100) {
  ## start in a random state
  K <- nrow(P)
  z <- c(sample(seq_len(K), 1), vector(length = time_len - 1))
  for (i in seq_len(time_len - 1)) {
    z[i + 1] <- sample(seq_len(K), size = 1, prob = P[z[i], ])
  }
  z
}
