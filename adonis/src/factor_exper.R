#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Some experiments with adonis under other settings. Focused more on mixture
## and latent factor models. These generally maintain exchangeability between
## sites, so the permutation p-values are in fact valid.
##
## author: sankaran.kris@gmail.com
## date: 08/14/2017

library("tidyverse")
library("scales")
library("reshape2")
library("kernlab")
library("mvtnorm")
library("vegan")

## variables used throughout experiment
n_sim <- 100
p1 <- 20
p2 <- 2
n <- 100
u <- matrix(runif(n * p2), n, p2)
K <- kernelMatrix(rbfdot(sigma = 5), u)

logit <- function(x) {
  1 / (1 + exp(-x))
}

print_iter <- function(i, m = 10) {
  if (i %% m == 0) {
    cat(sprintf("iteration %s\n", i))
  }
}


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

## write Sigma for factor analysis
k <- 1
sigma <- 1
W <- matrix(runif(p2 * k), p2, k)
Sigma <- W %*% t(W) + sigma * diag(nrow = p2) # factor analysis Sigma

## Factor setup
models <- list()
for (i in seq_len(n_sim)) {
  print_iter(i)
  x <- rmvnorm(n, sigma = W %*% t(W) + sigma * diag(nrow = p2))
  models$factor$adonis[[i]] <- adonis(x ~ y[, i] + u, method = "euclidean")
  models$factor$logistic[[i]] <- glm(y[, i] ~ x + u, family = binomial())
  models$factor$lm[[i]] <- lm(x[, 1] ~ y[, i] + u)
}

pvals <- list(
  "factor" = list(
    "adonis" = sapply(models$factor$adonis, function(x) x$aov.tab["y[, i]", "Pr(>F)"]),
    "logistic" = as.numeric(sapply(models$factor$logistic, function(x) coef(summary(x))[2:(p2 + 1), "Pr(>|z|)"])),
    "lm" = sapply(models$factor$lm, function(x) coef(summary(x))["y[, i]", "Pr(>|t|)"])
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
