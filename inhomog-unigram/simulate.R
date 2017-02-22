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
#' @param sigmas [length T - 1 numeric vector] Factors that induce larger or
#'   smaller variances at specified timepoints. Specifically, the variance of
#'   the change in values from time t[i] to t[i + 1] is given by sigma[i] *
#'   sqrt(t[i + 1] - t[i]). If not provided, defaults to vector of all 1s.
#' @return mus [T x D numeric matrix] The sampled values of the random walk.
#' @examples
#' ## standard RW
#' mu0 <- rnorm(10)
#' ts <- seq(0, 1, by = 0.02)
#' mus <- gaussian_rw(mu0, ts)
#'
#' ## changepoint in variances
#' sigmas <- c(rep(1, 20), rep(4, 5), rep(1, 25))
#' mus <- gaussian_rw(mu0, ts, sigmas)
#' plot(ts, mus[, 1])
gaussian_rw <- function(mu0, ts, sigmas = NULL) {
  n_times <- length(ts)
  V <- length(mu0)
  mus <- matrix(0, n_times, V)

  if (is.null(sigmas)) {
    sigmas <- rep(1, n_times - 1)
  }
  print(sigmas)

  for (i in seq_along(ts)[-1]) {
    mus[i, ] <- sigmas[i - 1] * sqrt(ts[i] - ts[i - 1]) * rnorm(V, mean = mus[i - 1, ])
  }
  mus
}
