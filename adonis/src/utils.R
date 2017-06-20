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
  adonis(Y ~ ., data = X, ...)
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
  adonis(Y ~ ., data = X, ...)
}

plot_perms <- function(mod, fname, output_dir = "../doc/figure/") {
  p <- ggplot(as_data_frame(mod$f.perms)) +
    geom_histogram(aes(x = V1), bins = 500)
  ggsave(file.path(output_dir, fname), p)
}
