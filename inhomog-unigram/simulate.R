#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulating data from unigram models that are nonstationary (e.g., evolve
## faster at some points than others).
## author: kriss1@stanford.edu

#' Multidimensional Gaussian RW
#'
#' @param mu0 [length D numeric vector] The initial value of the RW.
#' @param ts [length T numeric vector] The timepoints at which to sample points.
#'   Must be sorted.
#' @return mus [T x D numeric matrix] The sampled values of the random walk.
#' @examples
#' mu0 <- rnorm(10)
#' ts <- sort(10 * rbeta(50, 1, 4))
#' mus <- gaussian_rw(mu0, ts)
gaussian_rw <- function(mu0, ts) {
  n_times <- length(ts)
  V <- length(mu0)
  mus <- matrix(0, n_times, V)

  for (i in seq_along(ts)[-1]) {
    mus[i, ] <- sqrt(ts[i] - ts[i - 1]) * rnorm(V, mean = mus[i - 1, ])
  }
  mus
}

