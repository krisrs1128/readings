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
  sum(x * xi(alpha), na.rm = TRUE) # sometimes probabilities are too small for digamma
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

###############################################################################
# Variational Lower Bound
###############################################################################

ndkv_xi_theta <- function(ndkv_tilde, xi_theta_tilde) {
  V <- dim(ndkv_tilde)[3]
  products  <- vector(length = V)
  for (v in seq_len(V)) {
    products[v] <- sum(ndkv_tilde[,, v] * xi_theta_tilde)
  }
  sum(products)
}

ndkv_xi_beta <- function(ndkv_tilde, xi_beta_tilde) {
  D <- dim(ndkv_tilde)[1]
  products  <- vector(length = D)
  for (d in seq_len(D)) {
    products[d] <- sum(ndkv_tilde[d,, ] * xi_beta_tilde)
  }
  sum(products)
}

#' @param ndkv_tilde [3-d array] An array of document x topic x word
#' variational parameters.
#' @param theta_tilde [matrix] An array of document x topic variational
#' parameters.
#' @param beta_tilde [matrix] An array of topic x vocabulary variational
#' parameters
#' @param alpha [vector] The dirichlet hyperparameter for the theta (document
#' topic) mixture proportions
#' @param eta [vector] The analog for alpha on the beta (vocabulary)
#' @return elbo scalar The evidence lower bound, which should increase after
#' every variational update.
#' proportions.
lower_bound <- function(ndkv_tilde,
                        theta_tilde,
                        beta_tilde,
                        alpha,
                        eta) {
  N <- dim(ndkv_tilde)[1]
  xi_theta_tilde <- t(apply(theta_tilde, 1, xi))
  xi_beta_tilde <- t(apply(beta_tilde, 1, xi))

  # complete data likelihood terms
  N * ndkv_xi_theta(ndkv_tilde, xi_theta_tilde)
  N * ndkv_xi_beta(ndkv_tilde, xi_beta_tilde) +
  sum(alpha * xi_theta_tilde) +
  sum(eta * xi_beta_tilde) +

  #  start entropy terms
  N * multinomial_entropy(N, as.numeric(ndkv_tilde)) +
  sum(dirichlet_entropies(theta_tilde, theta_tilde)) +
  sum(dirichlet_entropies(beta_tilde, beta_tilde))
}

###############################################################################
# E and M steps
###############################################################################

#' @param nv [vector] A vector of word counts associated with the i^th document
#' @param beta [matrix] A topics x vocabulary matrix of word probabilities
#' across topics.
#' @param alpha [scalar] The dirichlet hyperparameter for the theta (document
#' topic) mixture proportions.
e_step <- function(nv, beta_tilde, alpha, n_iter = 10) {
  theta_tilde <- alpha
  xi_beta_tilde <- apply(beta_tilde, 2, xi)
  K <- nrow(beta_tilde)
  V <- ncol(beta_tilde)

  for (iter in seq_len(n_iter)) {
    nkv_tilde <- exp(xi_beta_tilde + xi(theta_tilde) %*% matrix(1, 1, V))
    nkv_tilde[is.na(nkv_tilde)] <- 0

    for (v in seq_len(V)) {
      nkv_tilde[, v] <- nkv_tilde[, v] / sum(nkv_tilde[, v])
    }

    theta_tilde <- rowSums((matrix(1, K, 1) %*% nv) * nkv_tilde)
  }

  list(
    "theta_tilde" = theta_tilde,
    "nkv_tilde" = nkv_tilde
  )
}

#' @param eta [scalar] The analog for alpha on the beta (vocabulary)
#' @param S [matrix] A topics x vocabulary matrix of sufficient statistics
m_step <- function(eta, S) {
  eta + S
}

###############################################################################
# Overall variational bayes
###############################################################################

#' @param ndv [matrix] A document by vocabulary matrix of word counts.
#' @example
#' ndv <- matrix(sample(1:10, 200 * 100, replace = T), 200, 100)
#' vb_counts(ndv, rep(1, 3), rep(1, 100))
vb_counts <- function(ndv, alpha, eta, beta_tilde_init = NULL, n_iter = 100) {
  K <- length(alpha)
  D <- nrow(ndv)
  V <- ncol(ndv)

  if (is.null(beta_tilde_init)) {
    beta_tilde_init <- t(rdirichlet(V, alpha))
  }
  beta_tilde <- beta_tilde_init
  elbo <- vector(length = n_iter)

  latent_data <- list(
    "theta_tilde" = matrix(0, D, K),
    "ndkv_tilde" = array(0, c(D, K, V))
  )

  for (iter in seq_len(n_iter)) {
    S <- matrix(0, K, V)

    for (d in seq_len(D)) {
      cur_latent_data <- e_step(ndv[d, ], beta_tilde, alpha)
      latent_data$theta_tilde[d, ] <- cur_latent_data$theta_tilde
      latent_data$ndkv_tilde[d,, ] <- cur_latent_data$nkv_tilde
      S <- S + (matrix(1, K, 1) %*% ndv[d, ]) * latent_data$ndkv_tilde[d,, ]
    }

    elbo[iter] <- lower_bound(
      latent_data$ndkv_tilde,
      latent_data$theta_tilde,
      beta_tilde,
      alpha,
      eta
    )
    cat(sprintf("iter %d | elbo %f \n ", iter, elbo[iter]))

    beta_tilde <- m_step(eta, S)
  }

  list (
    "elbo" = elbo,
    "theta_tilde" = latent_data$theta_tilde,
    "ndkv_tilde" = latent_data$ndkv_tilde,
    "beta_tilde" = beta_tilde
  )
}
