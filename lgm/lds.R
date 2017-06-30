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
#' lines(lds$x, col = "red")
#' lines(inf$x_filter, col = "blue")
#' lines(inf$x_smooth, col = "orange")
lds_inference <- function(y, A, C, Q, R, x01, v01) {
  time_len <- nrow(y)
  p <- nrow(C)
  k <- nrow(A)

  x_predict <- x01
  x_filter <- matrix(nrow = time_len, ncol = k)
  v_predict <- array(dim = c(k, k, time_len))
  v_predict[,, 1] <- v01
  v_filter <- array(dim = c(k, k, time_len))

  #' forward pass
  for (i in seq_len(time_len)) {
    if (i > 1) {
      x_predict <- A %*% x_filter[i - 1, ]
      v_predict[,, i] <- A %*% v_filter[,, i - 1] %*% t(A) + Q
    }

    K <- v_predict[,, i] %*% t(C) %*% solve(C %*% v_predict[,, i] %*% t(C) + R)
    x_filter[i, ] <- x_predict + K %*% (y[i, ] - C %*% x_predict)
    v_filter[,, i] <- v_predict[,, i] - K %*% C %*% v_predict[,, i]
  }

  #' backwards pass
  x_smooth <- matrix(nrow = time_len, ncol = k)
  x_smooth[time_len, ] <- x_filter[time_len, ]
  v_smooth <- array(dim = c(k, k, time_len))
  v_pair <- array(dim = c(k, k, time_len - 1))
  v_pair[,, time_len - 1] <- (diag(k) - K %*% C) %*% A %*% v_filter[,, time_len - 1]

  for (i in seq(time_len, 2)) {
    J_prev <- v_filter[,, i - 1] %*% t(A) %*% solve(v_predict[,, i])

    x_smooth[i - 1, ] <- x_filter[i - 1, ] +
      J_prev %*% (x_smooth[i, ] - A %*% x_filter[i - 1, ])
    v_smooth[,, i - 1] <- v_filter[,, i - 1] +
      J_prev %*% (v_smooth[,, i] - v_predict[,, i]) %*% t(J_prev)

    if (i < time_len) {
      v_pair[,, i - 1] <- v_filter[,, i] +
        J_next %*% (v_pair[,, i] - A %*% v_filter[,, i]) %*% t(J_prev)
    }
    J_next <- J_prev
  }

  list(
    "x_filter" = x_filter,
    "v_filter" = v_filter,
    "x_smooth" = x_smooth,
    "v_smooth" = v_smooth,
    "v_pair" = v_pair
  )
}

lds_learn <- function(Y, k, eps = 0.001) {
  time_len <- nrow(Y)
  p <- ncol(Y)

  #' initialize parameters
  A <- diag(0.5, nrow = k)
  C <- matrix(runif(p * k), p, k)
  Q <- diag(1, k)
  R <- diag(1, p)
  x01 <- rep(0, k)
  v01 <- diag(1, k)

  alpha <- t(y) %*% y

  #' not calculating loglik right now...
  for (iter in seq_len(100)) {

    #' E-step
    inf <- lds_inference(y, A, C, Q, R, x01, v01)
    x_smooth <- inf$x_smooth
    v_smooth <- inf$v_smooth
    v_pair <- inf$v_pair

    delta <- 0
    gamma <- 0
    beta <- 0

    for (i in seq_len(time_len)) {
      delta <- delta + t(y) %*% x_smooth[i, ]
      gamma <- gamma + x_smooth[i, ] %*% x_smooth[i, ] + v_smooth[,, i]

      if (i > 1) {
        beta <- beta + x_smooth[i, ] %*% t(x_smooth[i - 1, ]) + v_pair[,, i]
      }
    }

    gamma1 <- gamma -
      x_smooth[time_len, ] %*% t(x_smooth[time_len, ]) - v_smooth[,, time_len]
    gamma2 <- gamma -
      x_smooth[1, ] %*% t(x_smooth[1, ]) - v_smooth[,, 1]

    #' M-step
    C <- delta %*% solve(gamma)
    R <- (alpha - C %*% t(delta)) / time_len
    A <- beta %*% solve(gamma1)
    Q <- (gamma2 - A %*% t(beta)) / (time_len - 1)
    x01 <- x_smooth[1, ]
    v01 <- v_smooth[,, 1]
  }

  list(
    "A" = A,
    "C" = C,
    "Q" = Q,
    "R" = R,
    "v01" = x01,
    "v01" = V01
  )
}
