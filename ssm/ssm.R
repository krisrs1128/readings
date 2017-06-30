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
ht_update <- function(qt, pi, phi, tau) {
  log_alpha <- forwards(phi, qt, pi)
  log_beta <- backwards(phi, qt)
  log_ht <- log_alpha + log_beta - log(tau)
}

#' @examples
#' phi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
#' lik <- matrix(runif(20), 10, 2)
#' forwards(phi, lik, c(0.5, 0.5))
forwards <- function(phi, qt, pi) {
  K <- nrow(phi)
  time_len <- nrow(qt)
  log_alpha <- matrix(0, time_len, K)

  ## base case
  log_alpha[1, ] <- log(p0) + qt[1, ]
  log_alpha[1, ] <- normalize_log_space(log_alpha[1, ])

  for (i in seq(2, time_len)) {
    log_alpha[i, ] <- qt[i, ] + log(t(phi) %*% exp(log_alpha[i - 1, ]))
    log_alpha[i, ] <- normalize_log_space(log_alpha[i, ])
  }

  log_alpha
}

#' @examples
#' phi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
#' lik <- matrix(runif(20), 10, 2)
#' backwards(phi, lik)
backwards <- function(phi, qt) {
  K <- nrow(phi)
  time_len <- nrow(qt)

  log_beta <- matrix(0, time_len, K)
  for (i in seq(time_len - 1, 1)) {
    for (k in seq_len(K)) {
      log_beta[i, k] <- lse(log_beta[i + 1, ] + qt[i + 1, ] + log(phi[k, ]))
    }
  }

  log_beta
}

#' ---- smoothing-estimates ----
#' ---- state-space-learning ----
#' ---- hmm-learning ----
