#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## This is an attempt to perform switching state space model estimation using variational methods, as described in
##
## Ghahramani, Zoubin, and Geoffrey E. Hinton. Switching state-space models.
## Technical Report CRG-TR-96-3 DRAFT, Dept. of Computer Science, University of
## Toronto, 1996.
##
## author: sankaran.kris@gmail.com
## date: 06/30/2017

## ---- libraries ----
library("expm")
library("Rcpp")
library("tidyverse")
library("MCMCpack")
sourceCpp("utils.cpp")

## ---- utils ----
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

## ---- qt-update ----
#' @examples
#' As <- list(diag(1, nrow = 1), diag(1, nrow = 1))
#' Cs <- list(diag(1, nrow = 1), diag(1, nrow = 1))
#' s <- c(rep(1, 50), rep(2, 50))
#' Qs <- list(diag(0.1, nrow = 1), diag(0.1, nrow = 1))
#' Rs <- list(diag(5, nrow = 1), diag(0.5, nrow = 1))
#' res <- simulate(As, Cs, s, 0, Qs, Rs)
#' y <- res$y
#'
#' test <- ssm_em(y, 2, 1, n_iter = 25)
#' test$lds_param[[1]]$R
#' test$lds_param[[2]]$R
#' plot(y)
#' points(test$lds_infer[[1]]$x_smooth * test$lds_param[[1]]$C[1, 1], col = "blue")
#' points(test$lds_infer[[2]]$x_smooth * test$lds_param[[2]]$C[1, 1], col = "red")
ssm_em <- function(y, M = 2, K = 1, n_iter = 10) {
  time_len <- nrow(y)
  p <- ncol(y)

  ## initialize parameters
  lds_param <- initialize_lds(M, K, p)
  lds_infer <- lds_inference_multi(y, lds_param, matrix(1, time_len, M))
  hmm_param <- list(
    "pi" = rep(1 / M, M),
    "phi" = rdirichlet(M, rep(1, M))
  )

  tau <- rep(1, n_iter)
  for (iter in seq_len(n_iter)) {
    cat(sprintf("iteration %s\n", iter))

    ## E-step
    log_q <- matrix(nrow = time_len, ncol = M)
    log_ht <- matrix(nrow = time_len, ncol = M)
    for (m in seq_len(M)) {
      log_q[, m] <- qt_update(
        y,
        lds_infer[[m]]$x_smooth,
        lds_infer[[m]]$v_smooth,
        lds_param[[m]]$C,
        lds_param[[m]]$R,
        tau[iter]
      )
    }

    log_q <- normalize_log(log_q)

    log_alpha <- forwards(hmm_param$phi, log_q, hmm_param$pi)
    log_beta <- backwards(hmm_param$phi, log_q)
    log_xi <- two_step_marginal(hmm_param$phi, log_q, log_alpha, log_beta)
    log_ht <- normalize_log(tau[iter] * (log_alpha + log_beta))

    lds_infer <- lds_inference_multi(y, lds_param, exp(log_ht))

    ## M-step
    for (m in seq_len(M)) {
      lds_param[[m]] <- lds_learn(
        y,
        lds_infer[[m]]$x_smooth,
        lds_infer[[m]]$v_smooth,
        lds_infer[[m]]$v_pair,
        exp(log_ht)[, m]
      )
      print(lds_param[[m]])
    }
    hmm_param <- hmm_learn(log_xi, log_ht)

  }
  list(
    "lds_param" = lds_param,
    "hmm_param" = hmm_param,
    "log_ht" = log_ht,
    "log_q" = log_q,
    "lds_infer" = lds_infer
  )
}

## mf_x is T x P x M array of mean field estimates for x
## mf_xx is T x P x P x M array of mean field estimates for xx^t
qt_update <- function(y, x, v, C, R, tau = 1) {
  time_len <- nrow(y)

  log_q <- vector(length = time_len)
  for (i in seq_len(time_len)) {
    xx <- x[i, ] %*% t(x[i, ]) + v[,, i]
    log_q[i] <- - (1 / 2) * t(y[i, ]) %*% solve(R) %*% y[i, ] +
      t(y[i, ]) %*% solve(R) %*% C %*% x[i, ] -
      (1 / 2) * trace(t(C) %*% solve(R) %*% t(C) %*% xx)
  }

  1 / (tau) * log_q
}

## @examples
## phi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
## lik <- matrix(runif(20), 10, 2)
## forwards(phi, lik, c(0.5, 0.5))
forwards <- function(phi, log_q, pi) {
  K <- nrow(phi)
  time_len <- nrow(log_q)
  log_alpha <- matrix(0, time_len, K)

  ## base case
  log_alpha[1, ] <- log(pi) + log_q[1, ]
  log_alpha[1, ] <- log_alpha[1, ] - lse(log_alpha[1, ])

  for (i in seq(2, time_len)) {
    log_alpha[i, ] <- log_q[i, ] + log(t(phi) %*% exp(log_alpha[i - 1, ]))
    log_alpha[i, ] <- log_alpha[i, ] - lse(log_alpha[i, ])
  }

  log_alpha
}

## @examples
## phi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
## lik <- matrix(runif(20), 10, 2)
## backwards(phi, lik)
backwards <- function(phi, log_q) {
  K <- nrow(phi)
  time_len <- nrow(log_q)

  log_beta <- matrix(0, time_len, K)
  for (i in seq(time_len - 1, 1)) {
    for (k in seq_len(K)) {
      log_beta[i, k] <- lse(log_beta[i + 1, ] + log_q[i + 1, ] + log(phi[k, ]))
    }
  }

  log_beta
}

## @examples
## pi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
## lik <- matrix(runif(20), 10, 2)
## beta <- exp(backwards(pi, lik))
## alpha <- forwards(pi, lik)$alpha
## two_step_marginal(pi, lik, alpha, beta)
## @references Section 17.4.3.2 of "Machine Learning" by Murpy
two_step_marginal <- function(phi, log_q, log_alpha, log_beta) {
  K <- nrow(phi)
  time_len <- nrow(log_q)
  log_xi <- array(0, dim = c(time_len - 1, K, K))
  for (i in seq_len(time_len - 1)) {
    for (j in seq_len(K)) {
      for (k in seq_len(K)) {
        log_xi[i, j, k] <- log_alpha[i, k] + log_q[i + 1, j] + log_beta[i + 1, j] + log(phi[k, j])
      }
    }
    log_xi[i,, ] <- log_xi[i,, ] - lse(log_xi[i,, ])
  }

  log_xi
}

## ---- smoothing-estimates ----
## @examples
## A <- diag(0.99, nrow = 1)
## C <- diag(1, nrow = 1)
## Q <- diag(0.1, nrow = 1)
## R <- diag(1, nrow = 1)
## lds <- simulate(A, C, 0, Q, R, time_len = 100)
## v_weights <- c(rep(1, 40), rep(0.001, 20), rep(1, 40))
## inf <- lds_inference(lds$y, A, C, Q, R, 0, 0.1, v_weights)
## inf2 <- lds_inference(lds$y, A, C, Q, R, 0, 0.1)
## plot(lds$y)
## lines(lds$x, col = "red")
## lines(inf$x_smooth, col = "orange")
## lines(inf2$x_smooth, col = "purple")
lds_inference <- function(y, A, C, Q, R, x01, v01, v_weights = NULL) {
  time_len <- nrow(y)
  p <- nrow(C)
  k <- nrow(A)

  if (is.null(v_weights)) {
    v_weights <- rep(1, time_len)
  }
  v_weights[v_weights < 1e-4] <- 1e-4

  ## initialize for the forwards pass
  x_predict <- x01
  x_filter <- matrix(nrow = time_len, ncol = k)
  v_predict <- array(dim = c(k, k, time_len))
  v_predict[,, 1] <- v01
  v_filter <- array(dim = c(k, k, time_len))

  ## make the forward pass
  for (i in seq_len(time_len)) {
    if (i > 1) {
      x_predict <- A %*% x_filter[i - 1, ]
      v_predict[,, i] <- A %*% v_filter[,, i - 1] %*% t(A) + Q
    }

    K <- v_predict[,, i] %*% t(C) %*% solve(C %*% v_predict[,, i] %*% t(C) + R / v_weights[i])
    x_filter[i, ] <- x_predict + K %*% (y[i, ] - C %*% x_predict)
    v_filter[,, i] <- v_predict[,, i] - K %*% C %*% v_predict[,, i]
  }

  ## initialize for backwards pass
  x_smooth <- matrix(nrow = time_len, ncol = k)
  x_smooth[time_len, ] <- x_filter[time_len, ]
  v_smooth <- array(dim = c(k, k, time_len))
  v_smooth[,, time_len] <- v_filter[,, time_len]
  v_pair <- array(dim = c(k, k, time_len - 1))
  v_pair[,, time_len - 1] <- (diag(k) - K %*% C) %*% A %*% v_filter[,, time_len - 1]

  ## make the backwards pass
  for (i in seq(time_len, 2)) {
    J_prev <- v_filter[,, i - 1] %*% t(A) %*% solve(v_predict[,, i])

    x_smooth[i - 1, ] <- x_filter[i - 1, ] +
      J_prev %*% (x_smooth[i, ] - A %*% x_filter[i - 1, ])
    v_smooth[,, i - 1] <- v_filter[,, i - 1] +
      J_prev %*% (v_smooth[,, i] - v_predict[,, i]) %*% t(J_prev)

    if (i < time_len) {
      v_pair[,, i - 1] <- v_filter[,, i] %*% t(J_prev) +
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

lds_inference_multi <- function(y, lds_param, weights) {
  lds_infer <- list()

  M <- length(lds_param)
  for (m in seq_len(M)) {
    lds_infer[[m]] <- lds_inference(
      y,
      lds_param[[m]]$A,
      lds_param[[m]]$C,
      lds_param[[m]]$Q,
      lds_param[[m]]$R,
      lds_param[[m]]$x01,
      lds_param[[m]]$v01,
      weights[, m]
    )
  }
  lds_infer
}

initialize_lds <- function(M, K, p) {
  lds_param <- vector(mode = "list", length = M)
  for (m in seq_len(M)) {
    lds_param[[m]] <- list()
    lds_param[[m]]$A <- diag(1, K)
    lds_param[[m]]$C <- matrix(1, p, K)
    lds_param[[m]]$Q <- diag(1, K)
    lds_param[[m]]$R <- diag(1, p)
    lds_param[[m]]$x01 <- rep(runif(K), K)
    lds_param[[m]]$v01 <- diag(1, K)
  }

  lds_param
}

## ---- state-space-learning ----
row_mult <- function(A, weights) {
  apply(sweep(A, 1, weights, "*"), 3, sum)
}

lds_learn <- function(y, x_smooth, v_smooth, v_pair, weights = NULL) {
  time_len <- nrow(y)
  if (is.null(weights)) {
    weights <- rep(1, time_len)
  }
  weights[weights < 1e-4] <- 1e-4
  K <- nrow(v_smooth)

  xx <- array(dim = c(time_len, K, K))
  xx_cross <- array(dim = c(time_len - 1, K, K))

  for (i in seq_len(time_len)) {
    xx[i,, ] <- v_smooth[,, i] + x_smooth[i, ] %*% t(x_smooth[i, ])
    if (i > 1) {
      xx_cross[i - 1,, ] <- v_pair[,, i - 1] +
        t(x_smooth[i - 1, ]) %*% t(x_smooth[i, ])
    }
  }

  C <- (t(y) %*% diag(weights) %*% x_smooth) %*% solve(row_mult(xx, weights))

  R <- 0
  for (i in seq_len(time_len)) {
    R <- R + weights[i] * (y[i, ] %*% t(y[i, ]) - C %*% x_smooth[i, ] %*% t(y[i, ]))
  }
  R <- R / sum(weights)

  A <- row_mult(xx_cross, weights[-1]) %*% solve(row_mult(xx[-1,,, drop = F], weights[-1]))
  Q <- (1 / sum(weights)) * (row_mult(xx, weights) - A %*% row_mult(xx_cross, weights[-1]))

  x01 <- x_smooth[1, ]
  v01 <- v_smooth[,, 1]

  list("C" = C, "R" = R, "A" = A, "Q" = Q, "x01" = x01, "v01" = v01)
}

## ---- hmm-learning ----
expected_njk <- function(log_xi) {
  K <- ncol(log_xi)
  log_njk <- matrix(nrow = K, ncol = K)
  for (j in seq_len(K)) {
    for (k in seq_len(K)) {
      log_njk[j, k] <- lse(log_xi[, j, k])
    }
  }

  exp(log_njk)
}

## gamma is time_len x K matrix of smoothing probabilities
## xi is the (time_len  - 1) x K x K two step marginal array
hmm_learn <- function(log_xi, log_gamma) {
  phi <- expected_njk(log_xi)
  phi <- phi / rowSums(phi)
  pi <- exp(log_gamma[1, ])
  list("phi" = phi, "pi" = pi)
}

## ---- elementary ----
trace <- function(A) {
  sum(diag(A))
}

normalize_log <- function(x) {
  t(apply(x, 1, function(x) { x - lse(x)} ))
}
