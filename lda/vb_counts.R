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
psi <- function(alpha) {
  digamma(alpha) - digamma(sum(alpha))
}

elbo <- function(ndv, vb_params, alpha, eta) {
  expected_complete_density(ndv, vb_params, alpha, eta) +
    dirichlet_entropy(vb_params$beta_tilde) +
    dirichlet_entropy(vb_params$theta_tilde) +
    cat_entropy(ndv, vb_params$ndvk_tilde)
}


#' @description E_{q}{log p(x, z)} for LDA
#' @examples
#' ndv <- matrix(sample(1:10, size = 200 * 15, replace = TRUE), 200, 15)
#' ndvk <- array(sample(1:10, size = 200 * 15 * 3, replace = TRUE), c(200, 3, 15))
#' expected_complete_density(ndv,
#'                           list(beta_tilde = matrix(runif(3 * 15), 3, 15),
#'                                ndvk_tilde = ndvk,
#'                                theta_tilde = matrix(runif(200 * 3), 200, 3)),
#'                           rep(1, 3), rep(1, 15))
expected_complete_density <- function(ndv, vb_params, alpha, eta) {
  V <- ncol(vb_params$beta_tilde)
  K <- nrow(vb_params$beta_tilde)
  D <- nrow(vb_params$theta_tilde)

  result <- 0
  for (k in seq_len(K)) {
    psi_beta_k <- psi(vb_params$beta_tilde[k, ])
    for (v in seq_len(V)) {
      result  <- result +
        (sum(ndv[, v] * vb_params$ndvk_tilde[, k, v]) + eta[v] - 1) * psi_beta_k[v]
    }
  }

  for (d in seq_len(D)) {
    psi_theta_d <- psi(vb_params$theta_tilde[d, ])
    for (k in seq_len(K)) {
      result <- result +
        (sum(ndv[d, ] * vb_params$ndvk_tilde[d, k, ]) + alpha[k] - 1) * psi_theta_d[k]
    }
  }

  result
}

#' @param alpha A p x k matrix of dirichlet parameters. Each row is considered
#' the parameter to a dirichlet distribution.
#' @return
#' sum_{j, k} log Gamma(alpha_{jk}) +
#' sum_{j} log Gamma(sum_{k} alpha_{jk}) -
#' sum_{j, k} (alpha_{jk} - 1) * psi_{v}(alpha_{j})
#'
#' @examples
#' alpha <- rdirichlet(20, c(1, 1, 1, 1))
#' dirichlet_entropy(alpha)
dirichlet_entropy <- function(alpha) {
  psi_ks  <- t(apply(alpha, 1, psi))
  sum(lgamma(alpha)) +
    sum(lgamma(rowSums(alpha))) +
    sum((alpha - 1) * psi_ks)
}

multinomial_entropy <- function(N, p) {
  N * sum(p * log(p))
}

###############################################################################
# Variational Lower Bound
###############################################################################

ndvk_psi_theta <- function(ndvk_tilde, psi_theta_tilde) {
  V <- dim(ndvk_tilde)[3]
  products  <- vector(length = V)
  for (v in seq_len(V)) {
    products[v] <- sum(ndvk_tilde[,, v] * psi_theta_tilde)
  }
  sum(products)
}

ndvk_psi_beta <- function(ndvk_tilde, psi_beta_tilde) {
  D <- dim(ndvk_tilde)[1]
  products  <- vector(length = D)
  for (d in seq_len(D)) {
    products[d] <- sum(ndvk_tilde[d,, ] * psi_beta_tilde)
  }
  sum(products)
}

#' @param ndvk_tilde [3-d array] An array of document x topic x word
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
lower_bound <- function(ndvk_tilde,
                        theta_tilde,
                        beta_tilde,
                        alpha,
                        eta) {
  N <- dim(ndvk_tilde)[1]
  psi_theta_tilde <- t(apply(theta_tilde, 1, psi))
  psi_beta_tilde <- t(apply(beta_tilde, 1, psi))

  # complete data likelihood terms
  N * ndvk_psi_theta(ndvk_tilde, psi_theta_tilde)
  N * ndvk_psi_beta(ndvk_tilde, psi_beta_tilde) +
  sum(alpha * psi_theta_tilde) +
  sum(eta * psi_beta_tilde) +

  #  start entropy terms
  N * multinomial_entropy(N, as.numeric(ndvk_tilde)) +
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
e_step <- function(nv, beta_tilde, alpha, n_iter = 100) {
  K <- nrow(beta_tilde)
  V <- ncol(beta_tilde)

  psi_beta_tilde <- matrix(0, K, V)
  for (k in seq_len(K)) {
    psi_beta_tilde[k, ] <- psi(beta_tilde[k, ])
  }

  nkv_tilde <- matrix(0, K, V)
  theta_tilde <- alpha
  for (iter in seq_len(n_iter)) {
    theta_tilde_old <- theta_tilde
    theta_tilde <- alpha

    for (v in seq_len(V)) {
      for (k in seq_len(K)) {
        nkv_tilde[k, v] <- exp(psi_beta_tilde[k, v] + psi(theta_tilde_old)[k])
      }
      nkv_tilde[, v] <- nkv_tilde[, v] / sum(nkv_tilde[, v])
      for (k in seq_len(K)) {
        theta_tilde[k] <- theta_tilde[k] + nv[v] * nkv_tilde[k, v]
      }
    }
  }

  list(
    "theta_tilde" = theta_tilde,
    "nkv_tilde" = nkv_tilde
  )
}

#' @param eta [scalar] The analog for alpha on the beta (vocabulary)
#' @param S [matrix] A topics x vocabulary matrix of sufficient statistics
m_step <- function(eta, S) {
  K <- nrow(S)
  V <- ncol(S)

  beta_tilde <- matrix(nrow = K, ncol = V)
  for (v in seq_len(V)) {
    beta_tilde[, v] <- eta[v] + S[, v]
  }
  beta_tilde
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
    beta_tilde_init <- matrix(rgamma(K * V, shape = 100, scale = 1 / 100), K, V)
  }
  beta_tilde <- beta_tilde_init
  elbo <- vector(length = n_iter)

  latent_data <- list(
    "theta_tilde" = matrix(0, D, K),
    "ndvk_tilde" = array(0, c(D, K, V))
  )

  for (iter in seq_len(n_iter)) {
    S <- matrix(0, K, V)

    for (d in seq_len(D)) {
      cur_latent_data <- e_step(ndv[d, ], beta_tilde, alpha)
      latent_data$theta_tilde[d, ] <- cur_latent_data$theta_tilde
      latent_data$ndvk_tilde[d,, ] <- cur_latent_data$nkv_tilde
      S <- S + (matrix(1, K, 1) %*% ndv[d, ]) * latent_data$ndvk_tilde[d,, ]
    }

    elbo[iter] <- lower_bound(
      latent_data$ndvk_tilde,
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
    "ndvk_tilde" = latent_data$ndvk_tilde,
    "beta_tilde" = beta_tilde
  )
}
