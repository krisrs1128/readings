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
#' Q <- diag(2, nrow = 1)
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


#' @examples
#' A <- diag(0.9, nrow = 1)
#' C <- diag(1, nrow = 1)
#' Q <- diag(2, nrow = 1)
#' R <- diag(1, nrow = 1)
#' res <- simulate(A, C, 0, Q, R)
#' filt <- kalman_filter(res$y, A, C, R, Q)
#' points(filt$mu, col = "red")
#' points(filt$mu_pred, col = "red")
kalman_filter <- function(y, A, C, R, Q) {
  time_len <- length(y)
  p <- nrow(C)
  k <- nrow(A)

  sigma <- array(0, dim = c(k, k, time_len))
  sigma_pred <- array(0, dim = c(k, k, time_len))
  K <- array(0, dim = c(k, k, time_len))
  mu <- matrix(0, time_len, k)
  mu_pred <- matrix(0, time_len, k)

  for (i in seq(2, time_len)) {
    ## predict
    mu_pred[i, ] <- A %*% mu[i - 1, ]
    sigma_pred[,, i] <- A %*% sigma[,, i - 1] %*% t(A) + Q
    K[,, i] <- solve(solve(sigma_pred[,, i]) + t(C) %*% R %*% t(C)) %*% t(C) %*% solve(R)

    ## update
    sigma[,, i] <- (diag(k) - K[,, i]) %*% C %*% sigma_pred[,, i]
    mu[i, ] <- mu_pred[i, ] + K[,, i] %*% (y[i, ] - C %*% mu_pred[i, ])
  }

  list(
    "mu" = mu,
    "sigma" = sigma,
    "mu_pred" = mu_pred,
    "sigma_pred" = sigma_pred,
    "K" = K
  )

}

#' @examples
#' A <- diag(0.9, nrow = 1)
#' C <- diag(1, nrow = 1)
#' Q <- diag(2, nrow = 1)
#' R <- diag(1, nrow = 1)
#' res <- simulate(A, C, 0, Q, R)
#' filt <- kalman_filter(res$y, A, C, R, Q)
#' backwards_pass(res$y, filt$mu_pred, filt$sigma_pred, filt$K, A, C, R)
backwards_pass <- function(y, mu_pred, sigma_pred, K, A, C, R) {
  time_len <- length(y)
  k <- nrow(A)
  p <- nrow(C)
  r <- matrix(0, time_len, k)
  N <- array(0, dim = c(k, k, time_len))

  ## derivable from other variables, but stored for convenience
  u <- matrix(0, time_len, k)
  M_inv <- array(0, dim = c(k, k, time_len))
  Ct <- array(0, dim = c(k, k, time_len))

  for (i in seq(time_len, 1)) {
    v_pred_inv <- solve(C %*% sigma_pred[,, i] %*% t(C) + R)
    M <- v_pred_inv + t(K[,, i]) %*% N[,, i] %*% K[,, i]
    M_inv[,, i] <- solve(M)
    u[i] <- v_pred_inv %*% (y[i, ] - mu_pred[i, ]) - t(K[,, i]) %*% r[i, ]
    Ct[,, i] <- M %*% C - t(K[,, i]) %*% N[,, i] %*% A

    if (i > 1) {
      r[i - 1] <- t(C) %*% u[i] + t(A) %*% r[i]
      N[,, i - 1] <- t(C) %*% v_pred_inv %*% C + t(A - K[,, i] %*% C) %*% N[,, i] %*% (A - K[,, i] %*% C)
    }

  }

  list(
    "r" = r,
    "N" = N,
    "u" = u,
    "M_inv" = M_inv,
    "Ct" = Ct
  )
}

#' @examples
#' A <- diag(0.9, nrow = 1)
#' C <- diag(1, nrow = 1)
#' Q <- diag(2, nrow = 1)
#' R <- diag(1, nrow = 1)
#' res <- simulate(A, C, 0, Q, R)
#' filt <- kalman_filter(res$y, A, C, R, Q)
#' back <- backwards_pass(res$y, filt$mu_pred, filt$sigma_pred, filt$K, A, C, R)
#' y_prime <- forward_scan(res$y, back$u, back$Ct, back$M_inv, filt$K, C)
#' plot(res$y)
#' points(y_prime, col = "green")
forward_scan <- function(y, u, Ct, M_inv, K, C) {
  k <- ncol(u)
  time_len <- nrow(u)
  b <- matrix(0, time_len, k)

  for (i in seq_len(time_len)) {
    u[i, ] <- u[i, ] - Ct[,, i] %*% b[i, ]
    delta <- M_inv[,, i] %*% u[i, ] + sqrtm(as.matrix(M_inv[,, i])) %*% rnorm(k)
    u[i, ] <- u[i, ] - solve(M_inv[,, i]) %*% delta
    b[i + 1] <- t(A - K[,, i] %*% C) %*% b[i, ] - K[,, i] %*% delta
    y[i, ] <- y[i, ] - delta
  }

  list("y" = y, "u" = u)
}

#' @examples
#' A <- diag(0.9, nrow = 1)
#' C <- diag(1, nrow = 1)
#' Q <- diag(2, nrow = 1)
#' R <- diag(1, nrow = 1)
#' res <- simulate(A, C, 0, Q, R)
#' filt <- kalman_filter(res$y, A, C, R, Q)
#' back <- backwards_pass(res$y, filt$mu_pred, filt$sigma_pred, filt$K, A, C, R)
#' fw <- forward_scan(res$y, back$u, back$Ct, back$M_inv, filt$K, C)
#' bw <- backward_scan(fw$y, fw$u, back$Ct, back$M_inv, filt$K, C)
#' plot(res$y)
#' points(fw$y, col = "green")
#' points(bw$y, col = "purple")
backward_scan <- function(y, u, Ct, M_inv, K, C) {
  k <- ncol(u)
  time_len <- nrow(u)
  b <- matrix(0, time_len, k)

  for (i in seq(time_len, 1)) {
    u[i, ] <- u[i, ] - K[,, i] %*% b[i, ]
    delta <- M_inv[,, i] %*% u[i, ] + sqrtm(as.matrix(M_inv[,, i])) %*% rnorm(k)
    u[i, ] <- u[i, ] - solve(M_inv[,, i]) %*% delta
    b[i - 1] <- t(A - K[,, i] %*% C) %*% b[i, ] - Ct[,, i] %*% delta
    y[i, ] <- y[i, ] - delta
  }

  list("y" = y, "u" = u)
}
