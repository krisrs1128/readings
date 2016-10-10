###############################################################################
# Logistic Regrsesion SGLD experiment
###############################################################################

source("sgld_logit.R")
library("distr")
library("reshape2")
library("ggplot2")

# simulate logistic regression data
n <- 1000
p <- 2
true_beta <- DExp(rate = .1)@r(p)
X <- matrix(rnorm(n * p, 0, 1), n, p)
py <- sigmoid(X %*% true_beta)
y <- ifelse(runif(n) < py, 1, -1)
train_ix <- 1:(0.9 * n)


# paramters for the SGLD
n_iter <- 5000
batch_size <- 5

n_search <- 1000
sgd_search <- data.frame(
  "gamma" = numeric(n_search),
  "eps_init" = numeric(n_search),
  "eps_final" = numeric(n_search),
  "error_beta" = numeric(n_search),
  "error_p_train" = numeric(n_search),
  "error_p_test" = numeric(n_search)
)

# random search loop
set.seed(100)
beta0 <- rep(0, length(true_beta))
for (i in seq_len(n_search)) {
  eps_init <- exp(runif(1, -7, -1))
  eps_final <- exp(runif(1, -7, -1)) * eps_init
  gamma <- runif(1, 0.5, 1)
  cat(sprintf("search: %d, eps_init: %f | eps_final: %f | gamma: %f \n",
              i, eps_init, eps_final, gamma))

  res <- sgld(X[train_ix, ],
              y[train_ix],
              beta0,
              grad_log_prior,
              grad_log_likelihood,
              eps = get_eps(n_iter, gamma, eps_init, eps_final),
              n_iter = n_iter,
              batch_size = batch_size)

  beta_hat <- colMeans(res$betas[(0.9 * n_iter) : n_iter, ])
  p_test <- sigmoid(X[-train_ix, ] %*% true_beta)
  p_hat_test <- sigmoid(X[-train_ix, ] %*% beta_hat)
  p_train <- sigmoid(X[train_ix, ] %*% true_beta)
  p_hat_train <- sigmoid(X[train_ix, ] %*% beta_hat)

  sgd_search[i, "eps_init"] <- eps_init
  sgd_search[i, "eps_final"] <- eps_final
  sgd_search[i, "gamma"] <- gamma
  sgd_search[i, "error_beta"] <- sqrt(mean((beta_hat - true_beta) ^ 2))
  sgd_search[i, "error_p_train"] <- sqrt(mean( (p_hat_test - p_test) ^ 2))
  sgd_search[i, "error_p_test"] <- sqrt(mean( (p_hat_train - p_train) ^ 2))
}

# evaluate results
p <- list()
for (error_type in c("log(error_beta)", "error_p_train", "error_p_test")) {
  p[[error_type]] <- list(
    ggplot(sgd_search) +
      geom_point(aes_string(x = "log(eps_init)", y = error_type, col = "log(eps_final)")),
    ggplot(sgd_search) +
      geom_point(aes_string(x = "log(eps_final)", y = error_type, col = "log(eps_init)")),
    ggplot(sgd_search) +
      geom_point(aes_string(x = "log(eps_init)", y = "log(eps_final)", col = error_type)),
    ggplot(sgd_search) +
      geom_point(aes_string(x = "gamma", y = error_type, col = "log(eps_init)"))
  )
}

# refit using smallest eps within 1-sd of min
error_sd <- sd(sgd_search$error_p_test)
opt_params <- sgd_search %>%
  filter(error_p_test < min(error_p_test) + error_sd) %>%
  filter(eps_init == min(eps_init))

eps <- get_eps(n_iter, opt_params$gamma, opt_params$eps_init, opt_params$eps_final)
res <- sgld(X, y, beta0, grad_log_prior,
            grad_log_likelihood, eps = eps,
            n_iter = n_iter, batch_size = batch_size)
plot(sigmoid(X %*% true_beta), sigmoid(X %*% res$betas[.1 * n_iter, ]))
plot(sigmoid(X %*% true_beta), sigmoid(X %*% res$betas[n_iter, ]))

ggplot(data.frame(iter = 1:n_iter, llp = log(-res$log_posterior))) +
  geom_point(aes(x = iter, y = llp), size = 0.1)

ggplot(data.frame(iter = seq_len(n_iter), beta = res$betas)) +
  geom_point(aes(x = iter, y = beta.1), size = 0.1)

ggplot(data.frame(iter = seq_len(n_iter), beta = res$betas) %>%
         filter(iter > 0.8 * n_iter)) +
  geom_point(aes(x = iter, y = beta.1), size = 0.1)

ggplot(data.frame(iter = seq_len(n_iter), beta = res$betas) %>%
         filter(iter > 0.8 * n_iter)) +
  geom_point(aes(x = iter, y = beta.2), size = 0.1)

# just for comparison
glm(y ~ ., family = "binomial", data = data.frame(y = y == 1, X))

# visualize the posterior, after thinning and removing burn-in
full_eps <- get_eps(100 * n_iter, opt_params$gamma, opt_params$eps_init, opt_params$eps_final)
full_res <- sgld(X, y, beta0, grad_log_prior,
            grad_log_likelihood, eps = full_eps,
            n_iter = 100 * n_iter, batch_size = batch_size)

betas_posterior <- full_res$betas[(0.1 * 100 * n_iter) : (100 * n_iter), ]
betas_posterior <- betas_posterior[rep(1:100, length.out = nrow(betas_posterior)) == 1, ]
ggplot(melt(betas_posterior)) +
  geom_histogram(aes(x = value), bins = 100) +
  facet_grid(~ Var2, scale = "free_x")

ggplot(data.frame(betas_posterior)) +
  geom_point(aes(x = X1, y = X2), size = .2, alpha = 0.3)
