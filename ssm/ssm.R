#! /usr/bin/env Rscript

#' File description -------------------------------------------------------------
#' This is an attempt to perform switching state space model estimation using variational methods, as described in
#'
#' Ghahramani, Zoubin, and Geoffrey E. Hinton. Switching state-space models.
#' Technical Report CRG-TR-96-3 DRAFT, Dept. of Computer Science, University of
#' Toronto, 1996.
#'
#' author: sankaran.kris@gmail.com
#' date: 06/30/2017

#' ---- libraries ----
library("expm")

#' ---- utils ----
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
  if (is.null(R)) {
    Rs <- rep(diag(p), m)
  }

  ## initialize
  y <- matrix(nrow = time_len, ncol = p)
  x <- matrix(0, nrow = time_len, ncol = k)
  if (!is.null(x0)) {
    x[1, ] <- x0
  }

  ## sample
  y[1, ] <- Cs[[s[1]]] %*% x[1, ] + sqrtR %*% rnorm(p)
  for (i in seq(2, time_len)) {
    x[i, ] <- As[[s[i]]] %*% x[i - 1, ] + sqrtm(Qs[[s[i]]]) %*% rnorm(k)
    y[i, ] <- Cs[[s[i]]] %*% x[i, ] + sqrtm(Rs[[s[i]]]) %*% rnorm(p)
  }

  list("x" = x, "y" = y)
}

#' ---- qt-update ----
trace <- function(A) {
  sum(diag(A))
}

#' mf_x is T x P x M array of mean field estimates for x
#' mf_xx is T x P x P x M array of mean field estimates for xx^t
qt_update <- function(y, mf_x, C, R, tau = 1) {
  time_len <- nrow(y)
  M <- dim(X)[3]

  log_q <- matrix(nrow = time_len, ncol = m)
  for (m in seq_len(M)) {
    for (i in seq_len(i)) {
      log_q[i, m] <- t(y[i, ]) %*% solve(R) %*% y[i, ] +
        t(y[i, ]) %*% solve(R) %*% C[,, m] %*% mf_x[i,, m] +
        trace(t(C[,, m]) %*% solve(R) %*% t(C[,, m]) %*% mf_xx[i,,, m])
    }
  }

  - 0.5 / tau * log_q
}

#' ---- ht-update ----

#' ---- smoothing-estimates ----
#' ---- state-space-learning ----
#' ---- hmm-learning ----
