#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Some experiments with adonis under other settings. Focused more on mixture
## and latent factor models. These generally maintain exchangeability between
## sites, so the permutation p-values are in fact valid.
##
## author: sankaran.kris@gmail.com
## date: 08/14/2017

###############################################################################
## Libraries and setup
###############################################################################

library("tidyverse")
library("scales")
library("reshape2")
library("kernlab")
library("mvtnorm")
library("vegan")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)
## variables used throughout experiment
n_sim <- 2000
p1 <- 20
p2 <- 2
n <- 100
u <- matrix(runif(n * p2), n, p2)
K <- kernelMatrix(rbfdot(sigma = 5), u)

###############################################################################
## utility functions
###############################################################################

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
  X <- matrix(nrow = n, ncol = p)
  for (i in seq_len(n)) {
    X[i, ] <- rmvnorm(1, mu_mat[z[i], ], Sigma)
  }

  X
}


###############################################################################
## simulation experiment
###############################################################################

## generate matrix y according to logit
probs <- t(logit(rmvnorm(n_sim, sigma = K)))
y <- matrix(0, n, n_sim)
for (i in seq_len(n)) {
  for (j in seq_len(n_sim)) {
    y[i, j] <- sample(
      0:1, 1,
      prob = c(1 - probs[i, j], probs[i, j])
    )
  }
}

## setup for factor analysis and mixture of gaussians
k <- 1
sigma <- 1
W <- matrix(runif(p2 * k), p2, k)
Sigma <- W %*% t(W) + sigma * diag(nrow = p2)
mu_mat <- matrix(rnorm(4 * p2, 0, 2), 4, p2)

models <- list()
for (i in seq_len(n_sim)) {
  print_iter(i)

  ## factor analysis
  x <- rmvnorm(n, sigma = W %*% t(W) + sigma * diag(nrow = p2))
  models$factor$adonis[[i]] <- adonis(x ~ y[, i] + u, method = "euclidean")
  models$factor$logistic[[i]] <- glm(y[, i] ~ x + u, family = binomial())
  models$factor$lm[[i]] <- lm(x[, 1] ~ y[, i] + u)

  ## mixture of gaussians
  x <- rmix_gauss(n, mu_mat)
  models$mix_gauss$adonis[[i]] <- adonis(x ~ y[, i] + u, method = "euclidean")
  models$mix_gauss$logistic[[i]] <- glm(y[, i] ~ x + u, family = binomial())
  models$mix_gauss$lm[[i]] <- lm(x[, 1] ~ y[, i] + u)
}

pvals <- list(
  "factor" = list(
    "adonis" = sapply(models$factor$adonis, function(x) x$aov.tab["y[, i]", "Pr(>F)"]),
    "logistic" = as.numeric(sapply(models$factor$logistic, function(x) coef(summary(x))[2:(p2 + 1), "Pr(>|z|)"])),
    "lm" = sapply(models$factor$lm, function(x) coef(summary(x))["y[, i]", "Pr(>|t|)"])
  ),
  "mix_gauss" = list(
    "adonis" = sapply(models$mix_gauss$adonis, function(x) x$aov.tab["y[, i]", "Pr(>F)"]),
    "logistic" = as.numeric(sapply(models$mix_gauss$logistic, function(x) coef(summary(x))[2:(p2 + 1), "Pr(>|z|)"])),
    "lm" = sapply(models$mix_gauss$lm, function(x) coef(summary(x))["y[, i]", "Pr(>|t|)"])
  )
)

m_pvals <- melt(pvals)
colnames(m_pvals) <- c("pval", "method", "mechanism")

m_pvals$method <- factor(m_pvals$method, levels = c("adonis", "logistic", "lm"))
p <- ggplot(m_pvals) +
  geom_histogram(aes(x = pval), binwidth = 0.02) +
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2)) +
  facet_grid(method ~ mechanism, scales = "free_y")
p
