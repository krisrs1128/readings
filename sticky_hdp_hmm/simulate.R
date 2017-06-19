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
  pi <- vector(length = K)
  for (k in seq_len(K - 1)) {
    pi[k] <- rbeta(1, 1, alpha) * stick_remaining
    stick_remaining <- stick_remaining - pi[k]
  }
  pi[K] <- stick_remaining

  pi
}

#' @examples
#' P <- trans_mat(10, 50)
trans_mat <- function(alpha, K, kappa = 0) {
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
