#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Functions used to sample the auxiliary variables in the HDP-HMM. These are in
## their own file because this step is shared in both the direct assignment and
## the block samplers.
##
## author: sankaran.kris@gmail.com
## date: 6/21/2017

#' Random CRT distributed variable
#'
#' This is the distribution of the number of tables associated with a CRP after
#' m people have been seated.
#'
#' See https://dukespace.lib.duke.edu/dspace/handle/10161/7204
rcrt <- function(gamma, m) {
  probs <- gamma / (seq_len(m) - 1 + gamma)
  sum(rbinom(m, 1, probs))
}

sample_m <- function(z, alpha, beta, kappa) {
  modes <- setdiff(names(beta), "new")
  m <- matrix(0, length(modes), length(modes),
              dimnames = list(modes, modes))
  n <- transition_counts(z)

  for (j in modes) {
    for (k in modes) {
      m[j, k] <- rcrt(alpha * beta[k] + kappa * (j == k), n[j, k])
    }
  }

  m
}

sample_override <- function(m_diag, rho, beta) {
  K <- length(m_diag)
  w <- vector(length = K)
  for (k in seq_len(K)) {
    w[k] <- rbinom(1, m_diag[k], rho / (rho + beta[k] * (1 - rho)))
  }

  w
}
