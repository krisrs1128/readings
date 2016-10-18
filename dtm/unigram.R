
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
