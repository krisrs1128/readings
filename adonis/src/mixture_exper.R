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

library("MCMCpack")
library("tidyverse")
library("scales")
library("reshape2")
library("kernlab")
library("mvtnorm")
library("vegan")
source("utils.R")

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
n_sim <- 100
p1 <- 20
p2 <- 2
n <- 100
u <- matrix(runif(n * p2), n, p2)
K <- kernelMatrix(rbfdot(sigma = 5), u)

###############################################################################
## simulation experiment
###############################################################################

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
  y <- sample_probs(rep(0.5, n))
  models$factor$adonis[[i]] <- adonis(x ~ y, method = "euclidean", perm = 99)
  models$factor$logistic[[i]] <- glm(y ~ x, family = binomial())
  models$factor$lm[[i]] <- lm(x[, 1] ~ y)

  ## mixture of gaussians
  x <- rmix_gauss(n, mu_mat)
  y <- sample_probs(runif(nrow(mu_mat)[x$z]))
  models$mix_gauss$adonis[[i]] <- adonis(x$x ~ y, method = "euclidean", perm = 99)
  models$mix_gauss$logistic[[i]] <- glm(y ~ x$x, family = binomial())
  models$mix_gauss$lm[[i]] <- lm(x$x[, 1] ~ y)

  ## latent dirichlet allocation
  x <- rlda(n, p2)
  y <- sample_probs(logit(x$theta %*% rnorm(4, 0, 2.5))[, 1])
  models$lda$adonis[[i]] <- adonis(x$x ~ y + x$theta[, -1], method = "euclidean", perm = 99)
  models$lda$logistic[[i]] <- glm(y ~ x$x + x$theta[, -1], family = binomial())
  models$lda$lm[[i]] <- lm(x$x[, 1] ~ y + x$theta[, -1])
}

pvals <- list(
  "factor" = list(
    "adonis" = sapply(models$factor$adonis, function(x) x$aov.tab["y", "Pr(>F)"]),
    "logistic" = as.numeric(sapply(models$factor$logistic, function(x) coef(summary(x))[2:(p2 + 1), "Pr(>|z|)"])),
    "lm" = sapply(models$factor$lm, function(x) coef(summary(x))["y", "Pr(>|t|)"])
  ),
  "mix_gauss" = list(
    "adonis" = sapply(models$mix_gauss$adonis, function(x) x$aov.tab["y", "Pr(>F)"]),
    "logistic" = as.numeric(sapply(models$mix_gauss$logistic, function(x) coef(summary(x))[2:(p2 + 1), "Pr(>|z|)"])),
    "lm" = sapply(models$mix_gauss$lm, function(x) coef(summary(x))["y", "Pr(>|t|)"])
  ),
  "lda" = list(
    "adonis" = sapply(models$lda$adonis, function(x) x$aov.tab["y", "Pr(>F)"]),
    "logistic" = as.numeric(sapply(models$lda$logistic, function(x) coef(summary(x))[2:(p2 + 1), "Pr(>|z|)"])),
    "lm" = sapply(models$lda$lm, function(x) coef(summary(x))["y", "Pr(>|t|)"])
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
