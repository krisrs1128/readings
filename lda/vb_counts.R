###############################################################################
# Variational Bayesian inference for counts version of LDA
#
###############################################################################

###############################################################################
# Helper functions
###############################################################################

#' @title Difference of Digammas
#' @description This appears in the entropy of Dirichlet variables:
#' Psi(alpha_{k}) - Psi(sum alpha_{k})
#' where Psi are digamma functions
xi <- function(alpha) {
  digamma(alpha) - digamma(sum(alpha))
}

#' @param x A length K vector of dirichlet probabilities
#' @param alpha A length K vector of dirichlet hyperparametrs
dirichlet_entropy <- function(x, alpha) {
  sum(x * xi(alpha))
}

#' @title Compute ExpectedEntropy across many Dirichlets
#' @param x An N x K matrix of dirichlet probabilities
#' @param
dirichlet_entropies <- function(x, alpha) {
  N <- nrow(x)
  entropies <- vector(length = N)
  for (i in seq_len(N)) {
    entropies[i] <- dirichlet_entropy(x[i, ], alpha[i, ])
  }
}

multinomial_entropy <- function(N, p) {
  N * sum(p * log(p))
}

#' @param nkv_tilde [3-d array] An array of document x topic x word
#' variational parameters.
#' @param theta_tilde [matrix] An array of document x topic variational
#' parameters.
#' @param beta_tilde [matrix] An array of topic x vocabulary variational
#' parameters
#' @param alpha [scalar] The dirichlet hyperparameter for the theta (document
#' topic) mixture proportions
#' @param eta [scalar] The analog for alpha on the beta (vocabulary)
#' @return elbo scalar The evidence lower bound, which should increase after
#' every variational update.
#' proportions.
lower_bound <- function(nkv_tilde,
                        theta_tilde,
                        beta_tilde,
                        alpha,
                        eta) {
  N <- sum(nkv)
  xi_theta_tilde <- t(apply(theta_tilde, 1, xi))
  xi_beta_tilde <- t(apply(beta_tilde, 1, xi))

  # complete data likelihood terms
  N * nkv_xi_theta(xi_theta_tilde, nkv_tilde) +
  N * nkv_xi_beta(xi_beta_tilde, nkv_tilde) +
  alpha * sum(xi_theta_tilde) +
  eta * sum(xi_beta_tilde) +

  #  start entropy terms
  N * multinomial_entropy(N, as.numeric(nkv_tilde)) +
  sum(dirichlet_entropies(theta, theta)) +
  sum(dirichlet_entropies(beta, beta))
}

nkv_xi_theta <- function(nkv_tilde, xi_theta_tilde) {
  V <- dim(nkv_tilde)[3]
  products  <- vector(length = V)
  for (v in seq_len(V)) {
    products[v] <- sum(nkv_tilde[,, v] * xi_tilde)
  }
  sum(products)
}

nkv_xi_beta <- function(nkv_tilde, xi_beta_tilde) {
  D <- dim(nkv_tilde)[1]
  products  <- vector(length = D)
  for (d in seq_len(D)) {
    products[d] <- sum(nkv_tilde[d,,] * xi_beta_tilde)
  }
  sum(products)
}
