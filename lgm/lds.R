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

#' ---- libraries ----
library("expm")

#' ---- utils ----
#' Simulate from an LDS
#'
#' x[t] = Ax[t - 1] + v[t]
#' y[t] = C[t]x[t] + w[t]
#'
#' where v and w are iid normal with covariances Q and R, respectively
#' @examples
#' A <- diag(0.99, nrow = 1)
#' C <- diag(1, nrow = 1)
#' Q <- diag(0.1, nrow = 1)
#' R <- diag(1, nrow = 1)
#' lds <- simulate(A, C, 0, Q, R, time_len = 100)
#' plot(lds$y)
#' points(lds$x, col = "red")
simulate <- function(A, C, x0 = NULL, Q = NULL, R = NULL, time_len = 50) {
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

#' @examples
#' A <- diag(0.99, nrow = 1)
#' C <- diag(1, nrow = 1)
#' Q <- diag(0.1, nrow = 1)
#' R <- diag(1, nrow = 1)
#' lds <- simulate(A, C, 0, Q, R, time_len = 100)
#' inf <- lds_inference(lds$y, A, C, Q, R, 0, 0.1)
#' plot(lds$y)
#' points(lds$x, col = "red")
#' points(inf$x_filter, col = "blue")
lds_inference <- function(y, A, C, Q, R, x01, v01) {
  time_len <- nrow(y)
  p <- nrow(C)
  k <- nrow(A)
  x_predict <- x01
  v_predict <- v01
  x_filter <- matrix(nrow = time_len, ncol = k)
  v_filter <- array(dim = c(k, k, time_len))

  #' forward pass
  for (i in seq_len(time_len)) {
    if (i > 1) {
      x_predict <- A %*% x_filter[i - 1, ]
      v_predict <- A %*% v_filter[,, i - 1] %*% t(A) + Q
    }

    K <- v_predict %*% t(C) %*% solve(C %*% v_predict %*% t(C) + R)
    x_filter[i, ] <- x_predict + K %*% (y[i, ] - C %*% x_predict)
    v_filter[,, i] <- v_predict - K %*% C %*% v_predict
  }

  list("x_filter" = x_filter, "v_filter" = v_filter)
}

lds_learn <- function(Y, k, c) {
}
