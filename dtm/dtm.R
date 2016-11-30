#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# Functions to simulate data according to a dynamic topic model

source("../lda/lda_counts.R")
source("unigram.R")

#' @examples
#' beta <- dtm_probs(15, 4, 50)
dtm_probs <- function(V, K, n_times, sigma = 0.1, beta0 = NULL) {
  beta <- list()
  for (k in seq_len(K)) {
    beta[[k]] <- topic_params(V, n_times, sigma, beta0)
  }
  beta <- do.call(abind::abind, c(beta, list(along = 0)))
  aperm(beta, c(2, 3, 1))
}

#' @examples
#' N <- 100
#' document_lengths <- rep(1000, N)
#' V <- 15
#' X <- dtm_data(N, V, document_lengths, sigma = 0.1)
dtm_data <- function(N, V, document_lengths, alpha = NULL, sigma = 0.1,
                     beta0 = NULL) {
  if (is.null(alpha)) {
    alpha <- rep(.01, 5) # topics prior
  }

  K <- length(alpha)
  theta <- topic_params(K, N, sigma) %>%
    apply(1, softmax) %>%
    t()

  beta <- dtm_beta(V, K, N, sigma, beta0)
  X <- matrix(nrow = N, ncol = V)
  for (i in seq_len(N)) {
    probs <- softmax(beta[i,, ] %*% theta[i, ])
    X[i, ] <- rmultinom(1, document_lengths[i], probs)
  }

  list(X = X, beta = beta, theta = theta)
}
