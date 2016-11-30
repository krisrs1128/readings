#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# Functions to simulate data according to a dynamic topic model

source("../lda/lda_counts.R")
source("unigram.R")

#' @examples
#' beta <- dtm_beta(15, 4, 50)
dtm_beta <- function(V, K, n_times, sigma = 1, beta0 = NULL) {
  beta <- list()
  for (k in seq_len(K)) {
    beta[[k]] <- topic_params(V, n_times, sigma, beta0)
  }
  beta <- do.call(abind::abind, c(beta, list(along = 0)))
  aperm(beta, c(2, 3, 1))
}
