
###############################################################################
# Dynamic unigram model
#
# This is a simple model used to describe the essential ideas in Blei and
# Lafferty's Dynamic Topic Models.
###############################################################################

###############################################################################
# Generate data from dynamic unigram model
###############################################################################

#' @title Generate topic parameters before putting in logit
#' @examples
#' beta <- topic_params(15, 1000)
#' plot(beta[, 2])
topic_params <- function(V, n_times, sigma = 1, beta0 = NULL) {
  beta <- matrix(nrow = n_times, ncol = V)
  if (is.null(beta0)) {
    beta0 <- rnorm(V, sd = sigma)
  }

  beta[1, ] <- beta0
  for (i in seq_len(n_times - 1)) {
    beta[i + 1, ] <- beta[i, ] + rnorm(V, 0, sigma)
  }

  beta
}

softmax <- function(mu) {
  exp(mu) / sum(exp(mu))
}

#' @examples
#' beta <- topic_params(10, 200, sigma = 0.1)
#' counts <- word_counts(beta)
#' library("reshape2")
#' library("ggplot2")
#' ggplot(melt(counts)) +
#'   geom_line(aes(x = Var1, y = value, group = Var2))
word_counts <- function(beta, document_lengths = NULL) {
  if (is.null(document_lengths)) {
    document_lengths <- rep(1000, nrow(beta))
  }

  n_times <- nrow(beta)
  V <- ncol(beta)
  counts <- matrix(nrow = n_times, ncol = V)

  for (i in seq_len(n_times)) {
    counts[i, ] <- rmultinom(1, document_lengths[i], softmax(beta[i, ]))
  }

  counts
}

###############################################################################
# calculate ELBO
###############################################################################

#' @title Evidence lower bound in dynamic unigram model
#' @examples
#' # simulate 200 timepoints, 100 words
#' beta <- topic_params(100, 200, sigma = 0.1)
#' ntv <- word_counts(beta)
#' mtv <- beta
#' Vt <- matrix(1, 200, 100)
#' zeta <- runif(200)
#' sigma <- 1
#' evidence_lower_bound(ntv, mtv, Vt, zeta, sigma)
evidence_lower_bound <- function(ntv, mtv, Vt, zeta, sigma) {
  n_times <- nrow(ntv)
  n_words <- ncol(ntv)
  res <- 0

  # contribution from expected likelihood
  res <- res + sum(ntv * mtv)
  for(i in seq_len(n_times)) {
    nt <- sum(ntv[i, ])
    for (v in seq_len(n_words)) {
      res <- res +
        nt * (exp(mtv[i, v] + .5 * Vt[i, v]) / zeta[i] - log(zeta[i]) + 1)
    }
  }

  # contribution from expected prior
  for (i in seq_len(n_times - 1)) {
    res <- res +
      1 / (2 * sigma ^ 2) * sum((mtv[i + 1, ] - mtv[i]) ^ 2) -
      1 / sigma ^ 2 * sum(Vt[i, ])
  }
  res <- res - 1 / (2 * sigma ^ 2) * (sum(Vt[1, ]) + sum(Vt[n_times, ]))

  # contribution from entropy
  res + (1 / 2) * sum(log(Vt))
}

###############################################################################
# Calculate required derivatives
###############################################################################

#' @title mean derivative for a single word
#' @description dm_tilde / dbeta in the appendix of Blei and Lafferty's Dynamic
#' Topic Models. i^{th} row are derivatives of mean at time i. j^{th} column are
#' derivatives with respect to time j.
#' @param nu [length T vector] A length t vector of variational parameters for
#' this one word (across) all times.
#' @param sigma [scalar] The assumed variance in the diffusion of the underlying
#' latent process.
#' @param Vt [length T vector] The smoothed variances for this one word across
#' all times.
mean_derivative <- function(nu, sigma, Vt) {
  dm_dbeta <- mean_derivative_forwards(nu, sigma, Vt)
  mean_derivative_backwards(dm_dbeta, sigma, Vt)
}

# forwards pass in derivatives for optimizing elbo wrt pseudo-observations
mean_derivative_forwards <- function(nu, sigma, Vt) {
  n_times <- length(Vt)
  dm_dbeta <- matrix(0, nrow = n_times, ncol = n_times)

  for (j in seq_len(n_times)) {
    for (i in seq_len(n_times - 1)) {
      gamma <- nu[i + 1] / (Vt[i] + sigma ^ 2 + nu[i + 1])
      dm_beta[i + 1, j] <- gamma * dm_debata[i, j]
      if (i == j) {
        dm_beta[i, j] <- dm_beta[i, j] + (1 - gamma)
      }
    }
  }

  dm_debta
}

# backwards pass in derivatives for optimizing elbo wrt pseudo-observations
mean_derivative_backwards <- function(dm_dbeta, sigma, Vt) {
  n_times <- length(Vt)
  dmtilde_dbeta <- matrix(0, nrow = n_times, ncol = n_times)
  dmtilde_dbeta[n_times, ] <- dm_dbeta[n_times, ]

  for (j in seq_len(n_times)) {
    for (i in (n_times - 1):1) {
      gamma <- sigma ^ 2 / (Vt[i] + sigma ^ 2)
      dmtilde_dbeta[i, j] <- gamma * dm_dbeta[i, j] +
        (1 - gamma) * dmtilde_dbeta[i + 1, j]
    }
  }

  dmtilde_dbeta
}

#' @param Vt [T x V matrix] a matrix of variances for word v at time t.
update_zeta <- function(mtv, Vt) {
  rowSums(exp(mtv + (1 / 2) * Vt))
}

gradient_beta_v <- function(nv, nt, m_tilde_v, dmtilde_dbeta, Vt_tilde, zeta,
                            sigma) {
  n_times <- length(nv)
  dl_dbeta <- rep(0, n_times)

  for (j in seq_len(n_times)) {
    for (i in seq_len(n_times - 1)) {
      dl_beta[j]  <- dl_dbeta[j] -
        1 / (sigma ^ 2) * (
          m_tilde_v[i + 1, s] - m_tilde_v[i, s]
        ) +
        dmtilde_dbeta[i, j] * (
          nv[i + 1] - (nt[i] * exp(m_tilde_v[i + 1] + Vt_tilde[i + 1] / 2))
        )
    }
  }
  dl_dbeta
}

unigram_optimize <- function(ntv, sigma, beta_hat_0, xi_hat_0, n_iter, stepsize,
                             n_steps) {
  stop("VB for unigram model is not implemented yet")
  n_times <- nrow(ntv)
  V <- ncol(ntv)

  mtv <- matrix(nrow = n_times, ncol = V)
  mtv_tilde <- matrix(nrow = n_times, ncol = V)
  Vt <- matrix(nrow = n_times, ncol = V)
  Vt_tilde <- matrix(nrow = n_times, ncol = V)
  elbo <- matrix(nrow = n_iter, ncol = V)

  zeta <- zeta0
  beta_hat <- beta_hat_0

  for (iter in seq_len(n_iter)) {
    for (v in seq_len(V)) {
      #kalman_wrapper
      zeta <- update_zeta(mtv, Vt)
      beta_hat_0[, v] <- update_beta(
        ntv[, v],
        rowSums(ntv),
        m_tilde_v,
        Vt_tilde_v,
        zeta,
        sigma
      )

      elbo[iter, v] <- evidence_lower_bound(
        ntv,
        mtv_tilde,
        Vt_tilde,
        zeta,
        sigma
      )
    }
  }
}
