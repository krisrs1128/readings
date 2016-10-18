
###############################################################################
# Dynamic unigram model
#
# This is a simple model used to describe the essential ideas in Blei and
# Lafferty's Dynamic Topic Models.
###############################################################################

#' @title Generate topic parameters before putting in logit
#' @examples
#' beta <- topic_params(15, 1000)
#' plot(beta[, 2])
topic_params <- function(V, n_times, sigma = 1, beta0 = NULL) {
  beta <- matrix(nrow = n_times, ncol = V)
  if (is.null(beta0)) {
    beta0 <- rep(0, V)
  }

  beta[1, ] <- beta0
  for (i in seq_len(n_times - 1)) {
    beta[i + 1, ] <- beta[i, ] + rnorm(V, 0, sigma)
  }

  beta
}

softmax <- function(mu) {
  exp(mu) / sum(exp(mu))
}

#' @examples
#' beta <- topic_params(10, 200, sigma = 0.1)
#' counts <- word_counts(beta)
#' library("reshape2")
#' library("ggplot2")
#' ggplot(melt(counts)) +
#'   geom_line(aes(x = Var1, y = value, group = Var2))
word_counts <- function(beta, document_lengths = NULL) {
  if (is.null(document_lengths)) {
    document_lengths <- rep(1000, nrow(beta))
  }

  n_times <- nrow(beta)
  V <- ncol(beta)
  counts <- matrix(nrow = n_times, ncol = V)

  for (i in seq_len(n_times)) {
    counts[i, ] <- rmultinom(1, document_lengths[i], softmax(beta[i, ]))
  }

  counts
}
