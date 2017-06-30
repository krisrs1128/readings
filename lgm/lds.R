#! /usr/bin/env Rscript

#' File description -------------------------------------------------------------
#' An attempt to implement EM for the linear dynamic system as described in
#' section A.3.1 of
#'
#' Roweis, Sam, and Zoubin Ghahramani. "A unifying review of linear Gaussian
#' models." Neural computation 11.2 (1999): 305-345.
#'
#' author: sankaran.kris@gmail.com
#' date: 6/29/2017

library("expm")
#' Simulate from an LDS
#'
#' x[t] = Ax[t - 1] + v[t]
#' y[t] = C[t]x[t] + w[t]
#'
#' where v and w are iid normal with covariances Q and R, respectively
#' @examples
#' simulate()
#' lds <- simulate(0.99, 1, 0, diag(0.1, nrow = 1), diag(1, nrow = 1), time_len = 100)
#' plot(lds$y)
#' points(lds$x, col = "red")
simulate <- function(A, C, x0 = NULL, Q = NULL, R = NULL, time_len = 50) {
  ## in case function is given scalars
  A <- as.matrix(A)
  C <- as.matrix(C)

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
  sqrtR <- sqrtm(R)
  sqrtQ <- sqrtm(Q)
  y[1, ] <- C %*% x[1, ] + sqrtR %*% rnorm(p)
  for (i in seq(2, time_len)) {
    x[i, ] <- A %*% x[i - 1, ] + sqrtQ %*% rnorm(k)
    y[i, ] <- C %*% x[i, ] + sqrtR %*% rnorm(p)
  }

  list("x" = x, "y" = y)
}

lds_inference <- function(Y, A, C, Q, R, x0_1, V0_1) {
}

lds_learn <- function(Y, k, c) {
}
