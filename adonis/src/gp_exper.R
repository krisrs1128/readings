#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## An experiment with Adonis inspired by
##
## Dismantling the Mantel tests (Guillot & Rousset)
##
## author: sankaran.kris@gmail.com
## date: 08/01/2017

###############################################################################
## libraries and setup
###############################################################################
library("kernlab")
library("expm")
library("vegan")
library("tidyverse")
library("reshape2")
library("mvtnorm")

print_iter <- function(i) {
  if (i %% 10 == 0) {
    cat(sprintf("iteration %s\n", i))
  }
}

logit <- function(x) {
  1 / (1 + exp(-x))
}

###############################################################################
## simulate underlying category data
###############################################################################

## variables used throughout experiment
n_sim <- 40
p1 <- 20
p2 <- 2
n <- 100
u <- matrix(runif(n * p2), n, p2)
K <- kernelMatrix(rbfdot(sigma = 5), u)

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

###############################################################################
## perform poisson and gaussian adonis / lm experiments
###############################################################################

models <- list(
  "poisson" = list(
    "adonis" = list(),
    "logistic" = list(),
    "lm" = list()
  ),
  "gaussian" = list(
    "adonis" = list(),
    "logistic" = list(),
    "lm" = list()
  )
)

## Poisson setup
for (i in seq_len(n_sim)) {
  print_iter(i)
  f <- t(rmvnorm(p1, sigma = K))
  x <- matrix(rpois(n * p1, lambda = exp(f)), n, p1)
  models$poisson$adonis[[i]] <- adonis(x ~ y[, i] + u, method = "bray")
  models$poisson$logistic[[i]] <- glm(y[, i] ~ x + u)
  models$poisson$lm[[i]] <- lm(x[, 1] ~ y[, i] + u)
}

## Gaussian setup
for (i in seq_len(n_sim)) {
  print_iter(i)
  x <- t(rmvnorm(p1, sigma = K))
  models$gaussian$adonis[[i]] <- adonis(x ~ y[, i] + u, method = "euclidean")
  models$gaussian$logistic[[i]] <- glm(y[, i] ~ x + u)
  models$gaussian$lm[[i]] <- lm(x[, 1] ~ y[, i] + u)
}

###############################################################################
## visualize the resulting p-values
###############################################################################
pvals <- list(
  "poisson" = list(
    "adonis" = sapply(models$poisson$adonis, function(x) x$aov.tab["y[, i]", "Pr(>F)"]),
    "logistic" = sapply(models$poisson$logistic, function(x) coef(summary(x))[2:(2 + p1), "Pr(>|t|)"]),
    "lm" = sapply(models$poisson$lm, function(x) coef(summary(x))["y[, i]", "Pr(>|t|)"])
  ),
  "gaussian" = list(
    "adonis" = sapply(models$gaussian$adonis, function(x) x$aov.tab["y[, i]", "Pr(>F)"]),
    "logistic" = sapply(models$gaussian$logistic, function(x) coef(summary(x))[2:(2 + p1), "Pr(>|t|)"]),
    "lm" = sapply(models$gaussian$lm, function(x) coef(summary(x))["y[, i]", "Pr(>|t|)"])
  )
)

m_pvals <- melt(pvals)
colnames(m_pvals) <- c("pval", "rsv", "sim", "method", "mechanism")
head(m_pvals)

ggplot(m_pvals) +
  geom_histogram(
    aes(x = pval)
  ) +
  facet_grid(method ~ mechanism, scales = "free_y")

## save to file
dir.create("data")
save(p_vals, models, file = "data/pois_exper.rda")
