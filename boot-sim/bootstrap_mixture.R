#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# Helper functions to simulate nonparametric mixtures using the bootstrap. This
# is interesting for generating microbiome-like data (zero-inflated, heavy
# tailed) with known underlying structure (like clusters or gradients).

#' Nonparametric bootstrap samples mixing across distributions
#'
#' @param Xs [list of K vectors] Each vector in the list is a sample
#'   from some true F_{k}. We will return samples From \sum_{k}
#'   \theta_{k}\hat{F}_{k}^{ast}
#' @param theta [B x K matrix] A matrix giving mixing proportions
#'   between the \hat{F}_{k}, for each of the B bootstrap samples.
#' @example
#' Xs <- list(
#'   rnorm(1000, 0, 1),
#'   rnorm(820, 3, 1)
#' )
#' theta <- matrix(rgamma(1000 * 2, 1), 1000, 2)
#' theta <- (1 / rowSums(theta)) * theta
#' x_star <- bootstrap_mixture(Xs, theta)
#' library("ggplot2")
#' ggplot(data.frame(theta = theta, x = x_star)) +
#'   geom_point(aes(x = x, y = theta.1))
bootstrap_mixture <- function(Xs, theta) {
  stopifnot(ncol(theta) == length(Xs))

  B <- nrow(theta)
  x_star <- vector(length = B)
  for (k in seq_along(Xs)) {
    x_star <- x_star + theta[, k] * sample(Xs[[k]], size = B, replace = TRUE)
  }
  x_star
}
