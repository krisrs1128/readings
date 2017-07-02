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

forwards <- function(y, thetas, pi, Phi) {
  time_len <- nrow(y)
  J <- nrow(Phi)
  K <- nrow(thetas$A)

  x <- array(dim = c(time_len, K, J)) ## filtered means
  v <- array(dim = c(time_len, K, K, J)) ## filtered covariances
  m <- matrix(nrow = time_len, ncol = K) ## filtered state probabilities (unnormalized)

  ## initialize
  for (j in seq_len(J)) {
    x[1,, j] <- thetas[[j]]$mu
    v[1,,, j] <- thetas[[j]]$Sigma
    m[1, ] <- pi
  }

  ## make the forwards pass
  for (cur_time in seq(2, time_len)) {
    for (j in seq_len(J)) {
      x_step <- vector(mode = "list", length = J)
      v_step <- vector(mode = "list", length = J)
      v_pairs_step <- vector(mode = "list", length = J)
      log_lik <- vector(length = J)
      m_pairs <- vector(length = J)

      ## compute filter statistics for each possible previous i
      for (i in seq_len(J)) {
        filt_res <- filter(
          x[cur_time - 1,, i],
          v[cur_time - 1,,, i],
          y[cur_time, ],
          thetas[[j]]
        )
        x_step[i] <- filt_res$x_cur
        v_step[i] <- filt_res$v_cur
        v_pairs_step[i] <- filt_res$v_pair
        log_lik[i] <- filt_res$log_lik
      }

      for (i in seq_len(J)) {
        m_pairs[i] <- exp(log_lik[j]) * Phi[i, j] * m[cur_time - 1, i]
      }
      m[cur_time, j] <- sum(m_pairs)
      collapse_res <- collapse(x_step, v_step, m_pairs / m[cur_time, j])
      x[cur_time,, j] <- collapse_res$mu_x
      v[cur_time,,, j] <- collapse_res$v_xy
    }
  }

  list("x" = x, "v" = v, "m" = m)
}


backwards <- function(x, v, m, theta, phi) {
  time_len <- nrow(y)
  J <- nrow(Phi)
  K <- nrow(thetas$A)

  xs <- array(dim = c(time_len, K, J)) ## smoothed means
  vs <- array(dim = c(time_len, K, K, J)) ## smoothed covariances
  ms <- matrix(nrow = time_len, ncol = K) ## smoothed state probabilities (unnormalized)

  ## initialize
  for (j in seq_len(J)) {
    xs[time_len,, j] <- x[time_len,, j]
    vs[time_len,,, j] <- v[time_len,, j]
    ms[time_len, ] <- m[time_len, ]
  }

}

#' Filter update
#'
#' See section A.1
filter <- function(x_prev, v_prev, y_cur, theta) {
  x_pred <- theta$A %*% x_prev
  v_pred <- theta$A %*% v_prev %*% t(theta$A) + theta$Q
  S <- theta$C %*% v_pred %*% t(C) + theta$R
  K <- v_pred %*% t(theta$C) %*% solve(S)
  e <- y_cur - theta$C %*% x_pred

  list(
   "x_cur" = x_pred + K %*% e,
   "v_cur" = v_pred - K %*% S %*% t(K),
   "v_pair" = (diag(length(x)) - K %*% theta$C) %*% theta$A %*% v_prev,
   "log_lik" = dmvrnorm(e, 0, S)
  )
}

#' Smoother update
#'
#' See section A.2
#'
#' xs and vs refer to smoothed data and covariances. x and v are filtered
#' analogs.
smooth <- function(xs_next, vs_next, x_cur, v_cur, v_next, v_pair, theta) {
  x_pred <- theta$A %*% x_cur
  v_pred <- theta$A %*% v_cur %*% t(theta$A) + theta$Q
  J <- v_cur %*% t(theta$A) %*% solve(v_pred)

  list(
    "xs_cur" = x_cur + J %*% (xs_next - x_pred),
    "vs_cur" = v_cur + J %*% (vs_next - v_pred) %*% t(J),
    "vs_pair" = v_pair + (vs_next - v_next) %*% solve(v_next) %*% v_pair
  )
}

#' Collapse-crossing operation
#'
#' See section A.3
#'
#' mu_xs, mu_ys, and v_xys are lists whose elements are the j^th mean vector /
#' covariance matrix. Maybe we can rewrite it as arrays, but for now I'm fine
#' with this.
collapse_cross <- function(mu_xs, mu_ys, v_xys, ps) {
  K <- length(ps)
  mu_x <- 0
  for (j in seq_along(ps)) {
    mu_x <- mu_x + ps[j] * mu_xs[[j]]
    mu_y <- mu_y + ps[j] * mu_ys[[j]]
  }

  v_xy <- matrix(0, ncol(v_xys[[1]]), nrow(v_xys[[2]]))
  for (j in seq_along(ps)) {
    v_xy <- v_xy +
      ps[j] * (v_xys[[j]] + (mu_xs[[j]] - mu_x) %*% t(mu_ys[[j]] - mu_y))
  }

  list(
    "mu_x" = mu_x,
    "mu_y" = mu_y,
    "v_xy" = v_xy
  )
}
  
#' Collapsing operation
#'
#' See section A.3
  collapse <- function(mu_xs, v_xs, ps) {
  collapse_cross(mu_xs, mu_xs, v_xs, ps)
}

###############################################################################
## learning
###############################################################################
