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
As <- list(diag(0.999, nrow = 1), diag(0.0, nrow = 1))
Cs <- list(diag(1, nrow = 1), diag(2, nrow = 1))
s <- c(rep(1, 15), rep(2, 25), rep(1, 10))
Qs <- list(diag(1, nrow = 1), diag(0.1, nrow = 1))
res <- simulate(As, Cs, s, 0, Qs)
plot(res$y, col = s)
points(res$x, col = "blue")

simulate <- function(As ,Cs, s, x0 = NULL, Qs = NULL, R = NULL) {
  p <- nrow(Cs[[1]])
  k <- nrow(As[[1]])
  time_len <- length(s)

  ## provide defaults
  if (is.null(Qs)) {
    Q <- rep(diag(k), length(unique(s)))
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
  y[1, ] <- Cs[[s[1]]] %*% x[1, ] + sqrtR %*% rnorm(p)
  for (i in seq(2, time_len)) {
    print(As[[s[i]]])
    x[i, ] <- As[[s[i]]] %*% x[i - 1, ] + sqrtm(Qs[[s[i]]]) %*% rnorm(k)
    y[i, ] <- Cs[[s[i]]] %*% x[i, ] + sqrtR %*% rnorm(p)
  }

  list("x" = x, "y" = y)
}
