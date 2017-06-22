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
    nu_delta_coef <- (emission[[l]]$zeta + 1) / (emission[[l]]$zeta * (emission[[l]]$nu - d - 1))
    y_dens[l] <- dmvnorm(emission[[l]]$theta, nu_delta_coef * emission[[l]]$nu_delta)
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
    Pi[l, ] <- rdirichlet(u)
  }
  Pi
}
