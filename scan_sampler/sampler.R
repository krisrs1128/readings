#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Studying the scan sampling algorithm, as described in
## De Jong, Piet. "The scan sampler for time series models." Biometrika 84.4
## (1997): 929-937.
## De Jong, Piet. "Smoothing and interpolation with the state-space model."
## Journal of the American Statistical Association 84.408 (1989): 1085-1088.
## APA
##
## author: sankaran.kris@gmail.com
## date: 07/17/2017

###############################################################################
## Necessary libraries
###############################################################################
library("expm")

###############################################################################
## simulate data from the kalman filter
###############################################################################


#' @examples
#' ## perturbation like example
#' A <- diag(0.9, nrow = 1)
#' C <- diag(1, nrow = 1)
#' Q <- diag(4, nrow = 1)
#' R <- diag(1, nrow = 1)
#' res <- simulate(A, C, 0, Q, R)
#' plot(res$y)
#' lines(res$x, col = "blue")
simulate <- function(A ,C, x0 = NULL, Q = NULL, R = NULL, time_len = 100) {
  p <- nrow(C)
  k <- nrow(A)

  ## provide defaults
  if (is.null(Q)) {
    Q <- diag(k)
  }
  if (is.null(R)) {
    R <- diag(p)
  }

  ## initialize
  y <- matrix(nrow = time_len, ncol = p)
  x <- matrix(0, nrow = time_len, ncol = k)
  if (!is.null(x0)) {
    x[1, ] <- x0
  }

  ## sample
  y[1, ] <- C %*% x[1, ] + sqrtm(R) %*% rnorm(p)
  for (i in seq(2, time_len)) {
    x[i, ] <- A %*% x[i - 1, ] + sqrtm(Q) %*% rnorm(k)
    y[i, ] <- C %*% x[i, ] + sqrtm(R) %*% rnorm(p)
  }

  list("x" = x, "y" = y)
}
