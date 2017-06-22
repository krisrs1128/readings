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
#' z_clust <- kmeans(y, 20)$cluster
#' plot(y[, 1], col = 'white', asp = 1)
#' text(y[, 1], labels = z_clust)
#'
#' block_sampler(y)
block_sampler <- function(y, hyper = list(), lambda = list()) {
  ## merge default opts
  hyper <- merge_default_hyper(hyper)
  lambda <- merge_default_lambda(lambda)

  ## initialize state space
  Pi <- matrix(1 / L, L, L, dimnames = list(1:L, 1:L))
  init_clust <- kmeans(y, L)
  z <- init_clust$cluster
  beta <- setNames(rep(0.1, L), 1:L)

  theta <- list()
  for (l in as.character(1:L)) {
    theta[[l]] <- list("mu" = init_clust$centers[l, ], "sigma" = 0.75 * cov(y))
  }

  for (i in seq_len(hyper$n_iter)) {
    cat(sprintf("iteration %s\n", i))
    msg <- messages(Pi, y, theta)
    z <- sample_z(Pi, y, theta, msg)

    m <- sample_m(z, hyper$alpha, beta, hyper$kappa)
    w <- sample_override(diag(m), hyper$kappa / (hyper$kappa + hyper$alpha), beta)
    beta <- sample_beta(m, w, gamma)

    Pi <- sample_pi(z, hyper$alpha, beta, hyper$kappa)
    theta <- sample_theta(y, z, theta, lambda, hyper$theta_iter)
    state <- list(z = z, beta = beta, theta = theta)
  }

  state
}

messages <- function(Pi, y, theta) {
  time_len <- nrow(y)
  modes <- colnames(Pi)
  m <- matrix(1, time_len, nrow(Pi),
              dimnames = list(1:time_len, modes))

  for (i in seq(time_len - 1, 1)) {
    y_dens <- multi_dmvnorm(y[i, ], theta)
    m[i, ] <- Pi %*% (y_dens * m[i + 1, ])
  }
  m
}

sample_z <- function(Pi, y, theta, msg) {
  time_len <- nrow(y)
  z <- rep(1, time_len)
  for (i in seq(2, time_len)) {
    y_dens <- multi_dmvnorm(y[i, ], theta)
    f <- Pi[z[i -1], ] * y_dens * msg[i, ]
    z[i] <- sample(seq_along(f), 1, prob = f / sum(f))
  }
  z
}

sample_beta <- function(m, w, gamma) {
  m_bar <- m
  diag(m_bar) <- diag(m_bar) - w

  setNames(
    rdirichlet(1, gamma / ncol(m_bar) + colSums(m_bar))[1, ],
    colnames(m)
  )
}

sample_pi <- function(z, alpha, beta, kappa) {
  modes <- names(beta)
  n <- transition_counts(z, modes)
  Pi <- matrix(0, length(modes), length(modes),
               dimnames = list(modes, modes))
  for (l in modes) {
    u <- alpha * beta + n[l, ]
    u[l] <- u[l] + kappa
    Pi[l, ] <- rdirichlet(1, u)[1, ]
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

sample_theta <- function(y, z, theta, lambda, n_iter) {
  modes <- names(theta)
  for (i in seq_len(n_iter)) {
    for (l in modes) {
      theta[[l]]$mu <- sample_mu(lambda$mu0, lambda$sigma0, y[z == l, ], theta[[l]]$sigma)
      theta[[l]]$sigma <- sample_sigma(lambda$nu, lambda$delta, y[z == l, ], theta[[l]]$mu)
    }
  }

  theta
}
