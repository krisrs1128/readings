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
#' sigmas <- c(rep(1, 20), rep(6, 10), rep(1, 19))
#' mus <- gaussian_rw(mu0, ts, sigmas)
#' plot(ts, mus[, 1])
#'
#' ## random walk on sigmas, in log space
#' log_sigmas <- gaussian_rw(0, ts)
#' plot(5 * exp(log_sigmas))
#' mus <- gaussian_rw(mu0, ts, 5 * exp(log_sigmas))
#' plot(ts, mus[, 1])
#'
#' ## gaussian process on sigmas
#' library("simData") # github.com/krisrs1128/simData
#' sigmas <- abs(gp_data(ts, bandwidth = 5e-3))
#' mus <- gaussian_rw(mu0, ts, sigmas)
#' plot(ts, mus[, 1])
#'
#' ## iid inverse gammas for sigma
#' sigmas <- 1 / rgamma(length(ts), 5, 5)
#' plot(sigmas)
#' mus <- gaussian_rw(mu0, ts, sqrt(sigmas))
#' plot(ts, mus[, 1])
gaussian_rw <- function(mu0, ts, sigmas = NULL) {
  n_times <- length(ts)
  p <- length(mu0)
  deltas <- matrix(
    rnorm((n_times - 1) * p),
    n_times - 1,
    p
  )

  if (is.null(sigmas)) {
    sigmas <- rep(1, n_times - 1)
  }

  for (i in seq_len(nrow(deltas))) {
    deltas[i, ] <- sqrt(ts[i + 1] - ts[i]) * sigmas[i] * deltas[i, ]
  }

  rw(mu0, deltas)
}

#' Standard RW
#'
#' @param x0 [length p numeric vector] Initial value of random walk
#' @param deltas [dimension T x p matrix] Value of changes for each (of T) steps
#'   in the RW
#' @return x [dimension T + 1 x p matrix] The summed deltas at each time point
rw <- function(x0, deltas) {
  n_times <- 1 + nrow(deltas)
  p <- ncol(deltas)
  x <- matrix(0, n_times, p)

  for (i in seq_len(n_times)[-1]) {
    x[i, ]  <- x[i - 1] + deltas[i - 1, ]
  }

  x
}

#' Order 1 ARCH
#'
#' @param ts [length n_times vector] A sorted vector of times
#' @param a0 [scalar] The initial value of the arch process
#' @param alpha_0 [scalar] The baseline noise in the ARCH
#' @param alpha_1 [scalar] How strongly the variances persist across timepoints
#'
#' @examples
#' mus <- arch(ts, alpha_0 = 1, alpha_1 = 0.95)
#' plot(ts, mus)
arch <- function(ts, a0 = 0.1, alpha_0 = 1, alpha_1 = 0.7) {
  n_times <- length(ts)
  eps_sq <- rnorm(n_times - 1) ^ 2
  a_sq <- c(a0 ^ 2, rep(NA, n_times - 1))

  for (i in seq_along(a_sq)[-1]) {
    a_sq[i] <- (alpha_0 + alpha_1 * a_sq[i - 1]) * (ts[i] - ts[i - 1]) * eps_sq[i - 1]
  }

  sqrt(a_sq)
}

