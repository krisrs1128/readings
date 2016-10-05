sigmoid <- function(z) {
  1 / (1 + exp(-z))
}

get_eps <- function(n_iter, gamma, eps_init, eps_final) {
  xi <- (eps_init / eps_final) ^ (1 / gamma)
  b <- (xi - n_iter) / (1 - xi)
  a <- eps_init * (b + 1) ^ gamma

  a / (b + seq_len(n_iter)) ^ gamma
}

grad_log_prior <- function(beta) {
  - sign(beta)
}

log_prior <- function(beta) {

}

log_likelihood <- function(X, y, beta) {
  sum(log(sigmoid(y * X %*% beta)))
}

grad_log_likelihood <- function(X, y, beta) {
  sig <- sigmoid(-y * X %*% beta)
  sweep(X, 1, y * sig, FUN = "*")
}

sgld <- function(X, y, beta0, grad_log_prior, grad_log_likelihood, eps = NULL,
                 n_iter = 1000, batch_size = 10) {

  betas <- matrix(NA, nrow = n_iter, ncol = length(beta0))
  lls <- vector(length = n_iter)
  betas[1, ] <- beta0
  lls[1] <- log_likelihood(X, y, beta0)
  N <- nrow(X)
  batches <- rep(1:batch_size, length.out = N)

  for (i in 2:n_iter) {
    if (i %% 50 == 0) {
      cat(sprintf("iter %d \n", i))
    }
    cur_ix <- batches == (i %% batch_size + 1)
    betas[i, ] <- betas[i - 1, ] + (eps[i] / 2) * (
      grad_log_prior(betas[i - 1, ]) +
        N / batch_size * colSums(
          grad_log_likelihood(X[cur_ix, ], y[cur_ix], betas[i - 1, ])
        )
    ) + rnorm(1, 0, sqrt(eps[i]))

    lls[i] <- log_likelihood(X, y, betas[i, ])
  }

  list(betas = betas, lls = lls)
}
