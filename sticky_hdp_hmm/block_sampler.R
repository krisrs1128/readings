#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Block sampler for the HDP-HMM. Based on Algorithm 10 of ``Bayesian
## Nonparametric Learning of Complex Dynamical Phenomena''
##
## author: sankaran.kris@gmail.com
## date: 6/21/2017

## ---- libraries ----
library("mvtnorm")
library("jsonlite")
source("auxiliary.R")
source("utils.R")
source("simulate.R")

## ---- functions ----
#' @examples
#' K <- 20
#' gamma <- 2
#' alpha <- 0.9
#' kappa <- 3
#'
#' z <- c(rep(1, 10), rep(2, 5), rep(1, 10), rep(3, 3), rep(2, 10))
#' lambda <- list(zeta = .02, theta = c(0, 0), nu = 10, Delta = diag(c(0.1, 0.1)))
#' theta <- emission_parameters(K, lambda)
#' y <- emissions(z, theta)
#' plot(y[, 1], col = z)
block_sampler <- function(y, L = 20) {
  ## initialize state space
  Pi = matrix(1 / L, L, L)
  z <- kmeans(y, L)
}

messages <- function(Pi, y, emission) {
  L <- nrow(Pi)
  time_len <- nrow(y)
  m <- matrix(1, time_len, L)

  for (i in seq(time_len - 1, 1)) {
    y_dens <- multi_dmvnorm(y[i, ], emission)
    m[i, ] <- Pi %*% (y_dens * m[i + 1, ])
  }
  m
}

multi_dmvnorm <- function(yt, emission) {
  L <- length(emission)
  y_dens <- vector(length = L)
  for (l in seq_len(L)) {
    y_dens[l] <- dmvnorm(emission[[l]]$mu, emission[[l]]$sigma)
  }
  y_dens
}

sample_z <- function(Pi, y, emission, m) {
  time_len <- nrow(y)
  z <- vector(length = time_len)
  for (i in seq(2, time_len)) {
    y_dens <- multi_dvmnorm(y[i, ], emission)
    f <- Pi[z[i -1], ] * y_dens * m[i, ]
    z[i] <- sample(seq_along(f), f / sum(f))
  }
  z
}

sample_beta <- function(gamma, m_bar) {
  L <- ncol(m_bar)
  rdirichlet(1, gamma / L + colSums(m_bar))[1, ]
}

sample_pi <- function(alpha, beta, kappa, n) {
  L <- length(beta)
  Pi <- matrix(L, L)
  for (l in seq_len(L)) {
    u <- alpha * beta + n[l, ]
    u[l] <- mu[l] + kappa
    Pi[l, ] <- rdirichlet(u)[1, ]
  }
  Pi
}

sample_mu <- function(mu0, sigma0, y, sigma) {
  sigma_bar <- solve(solve(sigma0) + nrow(y) * solve(sigma))
  mu_bar <- sigma_bar %*% (solve(sigma0) %*% mu0 + sigma %*% rowSums(y))
  rnorm(1, mu_bar, sigma_bar)
}

sample_sigma <- function(nu, delta, y, mu) {
  nu_bar <- nu + nrow(y)
  e <- y - rep(1, nrow(y)) %*% matrix(mu, nrow = 1)
  nu_delta_bar <- nu * delta + t(e) %*% e
  riwishart(nu_delta_bar, nu_bar)
}

update_emission <- function(y, z, emission, nu, delta, sigma0, mu0) {
  L <- length(emission)
  for (i in seq_len(n_iter)) {
    for (l in seq_along(L)) {
      emission[[l]]$mu <- sample_mu(mu0, sigma0, y[z == l, ], emission[[l]]$sigma)
      emission[[l]]$sigma <- sample_sigma(nu, delta, y[z == l, ], emission[[l]]$mu)
    }
  }

  emission
}
