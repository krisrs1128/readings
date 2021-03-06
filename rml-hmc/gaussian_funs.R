
###############################################################################
# Functions for sampling normal model using (riemannian )langevin monte carlo
###############################################################################

###############################################################################
# helpers
###############################################################################

information <- function(theta) {
  (1 / theta["sigma"] ^ 2) * diag(c(1, 2))
}

mat_pow <- function(X, pow = 1) {
  eigen_X <- eigen(X)
  eigen_X$vectors %*% diag(eigen_X$values ^ pow) %*% t(eigen_X$vectors)
}

grad_U <- function(theta, x) {
  mu <- theta["mu"]
  sigma <- theta["sigma"]

  grad <- c( (-1 / sigma ^ 2) * sum(x - mu),
            length(x) / sigma - (1 / sigma ^ 3) * sum( (x - mu) ^ 2 ))
  names(grad) <- c("mu", "sigma")
  grad
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

#' @title Mahalanobis distances between two vectors
maha_dist <- function(x, y, mu, Sigma) {
  t(x - mu) %*% solve(Sigma) %*% (y - mu)
}

#' @Calculate log of ratio of gaussian densities
#' @description This is the log of the ratio of densities
#' N(x_1 | mu_1, Sigma_1) / N(x_2 | mu_2, Sigma_2).
log_gaussian_ratio <- function(x_1, x_2, mu_1, mu_2, Sigma_1, Sigma_2) {
  0.5 * (log(det(Sigma_2)) -
           log(det(Sigma_1)) -
           maha_dist(x_1, x_1, mu_1, Sigma_1) +
           maha_dist(x_2, x_2, mu_2, Sigma_2))
}

check_ratio <- function(log_acceptance_ratio) {
  if (log_acceptance_ratio > 0) {
    return (TRUE)
  } else {
    u <- runif(1)
    if (u < exp(log_acceptance_ratio))
      return (TRUE)
  }
  return (FALSE)
}

#' Generic MCMC
mcmc <- function(x, theta0, n_iter, proposal_fun, accept_fun) {
  thetas <- matrix(NA, nrow = n_iter, ncol = length(theta0))
  colnames(thetas) <- names(theta0)

  acceptances <- rep(0, n_iter)
  thetas[1, ] <- theta0

  for (i in 2:n_iter) {
    theta_new <- proposal_fun(thetas[i - 1, ], x)

    accept <- accept_fun(theta_new, thetas[i - 1, ], x)

    if (accept) {
      acceptances[i] <- 1
      thetas[i, ] <- theta_new
    } else {
      thetas[i, ] <- thetas[i - 1, ]
    }

  }

  list(thetas = thetas, acceptances = acceptances)
}

###############################################################################
# usual langevin
###############################################################################

langevin_proposal <- function(eps) {
  function(theta, x) {
    proposal_mean(theta, x, eps) + rnorm(2, 0, sqrt(2 * eps))
  }
}

proposal_mean <- function(theta, x, eps) {
  grad <- grad_U(theta, x)
  theta - (eps / 2) * grad
}

log_langevin_ratio <- function(x, theta_new, theta_cur, eps) {
  mu_forwards <- proposal_mean(theta_cur, x, eps)
  mu_reversed <- proposal_mean(theta_new, x, eps)
  Sigma <- 2 * diag(rep(eps, 2))

  log_gaussian_ratio(
    theta_cur,
    theta_new,
    mu_reversed,
    mu_forwards,
    Sigma,
    Sigma
  )
}

langevin_acceptance <- function(eps) {
  function(theta_new, theta_cur, x) {
    log_p_ratio <- log_likelihood_ratio(x, theta_new, theta_cur) # could be replaced by log guassian ratio
    log_q_ratio <- log_langevin_ratio(x, theta_new, theta_cur, eps)
    check_ratio(log_p_ratio + log_q_ratio)
  }
}

langevin_mcmc <- function(x, theta0, n_iter, eps) {
  mcmc(x, theta0, n_iter, langevin_proposal(eps), langevin_acceptance(eps))
}

###############################################################################
## rmc langevin
###############################################################################

rmc_proposal <- function(eps) {
  function(theta, x) {
    G <- length(x) * information(theta)
    theta_star <- rmc_proposal_mean(theta, x, eps) +
      mat_pow(G, -0.5) %*% rnorm(2, 0, sqrt(2 * eps))
    setNames(as.numeric(theta_star), c("mu", "sigma"))
  }
}

rmc_proposal_mean <- function(theta, x, eps) {
  n <- length(x)
  grad <- grad_U(theta, x)

  G <- n * information(theta)
  theta - (eps / 2) * (
    solve(G) %*% grad + theta["sigma"] / n
  )
}

log_rmc_ratio <- function(x, theta_new, theta_cur, eps) {
  n <- length(x)
  mu_forwards <- rmc_proposal_mean(theta_cur, x, eps)
  mu_reversed <- rmc_proposal_mean(theta_new, x, eps)

  Sigma_reversed <- 2 * eps * solve(n * information(theta_cur))
  Sigma_forwards <- 2 * eps * solve(n * information(theta_new))

  log_gaussian_ratio(
    theta_cur,
    theta_new,
    mu_reversed,
    mu_forwards,
    Sigma_reversed,
    Sigma_forwards
  )
}

rmc_acceptance <- function(eps) {
  function(theta_new, theta_cur, x) {
    log_p_ratio <- log_likelihood_ratio(x, theta_new, theta_cur) # could be replaced by log guassian ratio
    log_q_ratio <- log_rmc_ratio(x, theta_new, theta_cur, eps)
    check_ratio(log_p_ratio + log_q_ratio)
  }
}

rmc_mcmc <- function(x, theta0, n_iter, eps) {
  mcmc(x, theta0, n_iter, rmc_proposal(eps), rmc_acceptance(eps))
}
