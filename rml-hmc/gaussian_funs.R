
###############################################################################
# Functions for sampling normal model using (riemannian )langevin monte carlo
###############################################################################

grad_posterior <- function(theta, x) {
  mu <- theta["mu"]
  sigma <- theta["sigma"]

  grad <- c((1 / sigma ^ 2) * sum(x - mu),
            - length(x) / sigma + (1 / sigma ^ 3) * sum( (x - mu) ^ 2 ))
  names(grad) <- c("mu", "sigma")
  grad
}

proposal <- function(theta, x, eps) {
  proposal_mean(theta, x, eps) + rnorm(2, 0, eps)
}

proposal_mean <- function(theta, x, eps) {
  grad <- grad_posterior(theta, x)
  theta + (eps ^ 2 / 2) * grad
}

log_likelihood_ratio <- function(x, theta_new, theta_cur) {
  n <- length(x)

  sigma_cur <- theta_cur["sigma"]
  sigma_new <- theta_new["sigma"]

  mu_cur <- theta_cur["mu"]
  mu_new <- theta_new["mu"]

  n * (log(sigma_cur) - log(sigma_new)) -
    1 / (2 * sigma_new ^ 2) * sum( (x - mu_new) ^ 2) +
      1 / (2 * sigma_cur ^ 2) * sum( (x - mu_cur) ^ 2)
}

log_transition_ratio <- function(x, theta_new, theta_cur, eps) {
  mu_forwards <- proposal_mean(theta_cur, x, eps)
  mu_reversed <- proposal_mean(theta_new, x, eps)

  - 1 / (2 * eps ^ 2) * (
    - sum( (theta_cur - mu_reversed) ^ 2) + sum( (theta_new - mu_forwards) ^ 2)
  )
}

accept_proposal <- function(theta_new, theta_cur, x, eps) {
  log_acceptance_ratio <- log_transition_ratio(x, theta_new, theta_cur, eps) +
    log_likelihood_ratio(x, theta_new, theta_cur)

  if (log_acceptance_ratio > 0) {
    return (TRUE)
  } else {
    u <- runif(1)
    if (u > exp(log_acceptance_ratio))
      return (TRUE)
  }
  return (FALSE)
}

mcmc <- function(x, theta0, eps, n_iter) {
  thetas <- matrix(NA, nrow = n_iter, ncol = length(theta0))
  colnames(thetas) <- names(theta0)

  acceptances <- rep(0, n_iter)
  thetas[1, ] <- theta0

  for (i in 2:n_iter) {
    theta_new <- proposal(thetas[i - 1, ], x, eps)

    if (accept_proposal(theta_new, thetas[i - 1, ], x, eps)) {
      acceptances[i] <- 1
      thetas[i, ] <- theta_new
    } else {
      thetas[i, ] <- thetas[i - 1, ]
    }

  }

  list(thetas = thetas, acceptances = acceptances)
}
