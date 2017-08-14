#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Functions to accompany the study of nonparametric ANOVA.
##
## author: sankaran.kris@gmail.com
## date: 06/19/2017

factor_x <- function(n_rows, p_levels) {
  factors <- vector("list", length = length(p_levels))
  for (j in seq_along(factors)) {
    K <- length(p_levels[[j]])
    factors[[j]] <- sample(
      LETTERS[1:K],
      n_rows,
      prob = p_levels[[j]],
      replace = TRUE
    )
  }
  names(factors) <- paste0("factor_", seq_along(factors))
  as_data_frame(factors)
}

gaussian_null <- function(opts, ...) {
  Y <- matrix(
    rnorm(opts$n * opts$p, opts$mu_x, opts$sigma_x),
    opts$n, opts$p
  )
  X <- factor_x(opts$n, opts$p_levels)
  adonis(Y ~ ., data = X, method = opts$distance, ...)
}

low_rank_null <- function(opts, ...) {
  Y <- matrix(
    rnorm(opts$n * opts$p, opts$mu_x, opts$sigma_x),
    opts$n, opts$p
  )
  svd_y <- svd(Y)
  d <- svd_y$d
  d[-c(1:opts$K)] <- 0
  Y_hat <- svd_y$u %*% diag(d) %*% t(svd_y$v)

  X <- factor_x(opts$n, opts$p_levels)
  adonis(Y_hat ~ ., data = X, ...)
}

nb_null <- function(opts, ...) {
  Y <- matrix(
    rnbinom(opts$n * opts$p, opts$size, opts$prob),
    opts$n, opts$p
  )
  X <- factor_x(opts$n, opts$p_levels)
  adonis(Y ~ ., data = X, method = opts$distance, ...)
}

plot_perms <- function(mod, fname, output_dir = "../doc/figure/") {
  p <- ggplot(as_data_frame(mod$f.perms)) +
    geom_histogram(aes(x = V1), bins = 500)
  ggsave(file.path(output_dir, fname), p)
}

perm_histo <- function(f_stats) {
  f_stats <- do.call(bind_rows, f_stats)
  ggplot(f_stats) +
    geom_histogram(aes(x = f_perm), bins = 100)
}

save_fig <- function(p, fname, output_dir = "../doc/figure/") {
  ggsave(file.path(output_dir, fname), p)
}

perm_data <- function(mod, opts) {
  f_data <- mod$f.perms %>%
    as_data_frame() %>%
    rownames_to_column("iter") %>%
    gather(dimension, f_perm, -iter)
  opts$p2 <- length(opts$p_levels)
  opts$p_levels <- NULL
  opts_df <- as_data_frame(opts)
  bind_cols(opts_df[rep(1, nrow(f_data)), ], f_data)
}

logit <- function(x) {
  1 / (1 + exp(-x))
}

print_iter <- function(i, m = 10) {
  if (i %% m == 0) {
    cat(sprintf("iteration %s\n", i))
  }
}

#' Sample from Mixture of Gaussians
#'
#' @examples
#' mu_mat <- matrix(rnorm(5 * p2, 0, 4), 5, p2)
#' X <- rmix_gauss(1000, mu_mat)
#' hist(X[, 2], breaks = 100)
rmix_gauss <- function(n, mu_mat, Sigma = NULL, probs = NULL) {
  K <- nrow(mu_mat)
  p <- ncol(mu_mat)

  if (is.null(probs)) {
    probs <- rep(1 / K, K)
  }
  if (is.null(Sigma)) {
    Sigma <- diag(nrow = p)
  }

  z <- sample(1:K, size = n, prob = probs, replace = TRUE)
  x <- matrix(nrow = n, ncol = p)
  for (i in seq_len(n)) {
    x[i, ] <- rmvnorm(1, mu_mat[z[i], ], Sigma)
  }

  list("z" = z, "x" = x)
}

#' Simulate from an LDA model
#'
#' @examples
#' rlda(100, 10)
rlda <- function(n, V, lambda = 500, K = 4, alpha0 = 1) {
  beta <- rdirichlet(K, alpha0 * rep(1, V))
  x <- matrix(nrow = n, ncol = V)
  for (i in seq_len(n)) {
    x[i, ] <- rmultinom(
      1,
      rpois(1, lambda),
      prob = t(beta) %*% theta[i, ]
    )
  }

  list("x" = x, "theta" = theta)
}

sample_probs <- function(probs) {
  y <- vector(length = length(probs))
  for (i in seq_along(probs)) {
    y[i] <- sample(0:1, 1, replace = TRUE, prob = c(1 - probs[i], probs[i]))
  }

  y
}

plot_species_counts <- function(x, u, y) {
  plot_df <- data.frame("x" = x, "u" = u, "y" = y) %>%
    melt(
      id.vars = c("u.1", "u.2", "y"),
      variable = "species",
      value.name = "count"
    )

  ggplot(plot_df) +
    geom_point(
      aes(x = u.1, y = count, size = u.2, col = as.factor(y)),
      alpha = 0.4
    ) +
    scale_y_continuous(breaks = trans_breaks(identity, identity, n = 3)) +
    scale_size_continuous(range = c(0.005, 1.5)) +
    scale_color_brewer(palette = "Set2") +
    facet_wrap(~species, scales = "free") +
    labs(col = "y") +
    theme(
      strip.text.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 7)
    )
}
