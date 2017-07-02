#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## An attempt to implement the algorithm described in
##
## Murphy, Kevin P. Switching kalman filters. technical report, UC Berkeley,
## 1998.
##
## author: sankaran.kris@gmail.com
## date: 07/01/2017

###############################################################################
## setup
###############################################################################

library("mvtnorm")

###############################################################################
## simulation
###############################################################################
#' @examples
#' ## perturbation like example
#' As <- list(diag(0.1, nrow = 1), diag(0.9, nrow = 1))
#' Cs <- list(diag(1, nrow = 1), diag(1, nrow = 1))
#' s <- c(rep(2, 10), rep(1, 10), rep(2, 25), rep(1, 10), rep(2, 10))
#' Qs <- list(diag(4, nrow = 1), diag(1, nrow = 1))
#' Rs <- list(diag(10, nrow = 1), diag(0.5, nrow = 1))
#' res <- simulate(As, Cs, s, 0, Qs, Rs)
#' plot(res$y)
#' plot(res$y, col = s)
#' lines(res$x, col = "blue")
#'
#' ## example from paper
#' As <- list(diag(0.99, nrow = 1), diag(0.9, nrow = 1))
#' Qs <- list(diag(1, nrow = 1), diag(10, nrow = 1))
#' Cs <- list(diag(1, nrow = 1), diag(1, nrow = 1))
#' s <- c(rep(2, 10), rep(1, 10), rep(2, 25), rep(1, 10), rep(2, 10))
#' Rs <- list(diag(.1, nrow = 1), diag(0.1, nrow = 1))
#' res <- simulate(As, Cs, s, 0, Qs, Rs)
#' plot(res$y)
#' plot(res$y, col = s)
#' lines(res$x, col = "blue")
simulate <- function(As ,Cs, s, x0 = NULL, Qs = NULL, Rs = NULL) {
  p <- nrow(Cs[[1]])
  k <- nrow(As[[1]])
  m <- length(unique(s))
  time_len <- length(s)

  ## provide defaults
  if (is.null(Qs)) {
    Qs <- rep(diag(k), m)
  }
  if (is.null(Rs)) {
    Rs <- rep(diag(p), m)
  }

  ## initialize
  y <- matrix(nrow = time_len, ncol = p)
  x <- matrix(0, nrow = time_len, ncol = k)
  if (!is.null(x0)) {
    x[1, ] <- x0
  }

  ## sample
  y[1, ] <- Cs[[s[1]]] %*% x[1, ] + sqrt(Rs[[s[1]]]) %*% rnorm(p)
  for (i in seq(2, time_len)) {
    x[i, ] <- As[[s[i]]] %*% x[i - 1, ] + sqrtm(Qs[[s[i]]]) %*% rnorm(k)
    y[i, ] <- Cs[[s[i]]] %*% x[i, ] + sqrtm(Rs[[s[i]]]) %*% rnorm(p)
  }

  list("x" = x, "y" = y)
}

###############################################################################
## inference
###############################################################################

#' Filter update
#'
#' See section A.1
filter <- function(x_prev, v_prev, y_cur, theta) {
  x_pred <- theta$A %*% x_prev
  v_pred <- theta$A %*% v_prev %*% t(theta$A) + theta$Q
  S <- theta$C %*% v_pred %*% t(C) + theta$R
  K <- v_pred %*% t(theta$C) %*% solve(S)
  e <- y_cur - theta$C %*% x_pred
  log_lik <- dmvrnorm(e, 0, S)

  list(
   "x_cur" = x_pred + K %*% e,
   "v_cur" = v_pred - K %*% S %*% t(K),
   "v_pair" = (diag(length(x)) - K %*% theta$C) %*% theta$A %*% v_prev
  )
}

###############################################################################
## learning
###############################################################################
