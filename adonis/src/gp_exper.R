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
## actual simulation experiment
###############################################################################

## variables used throughout experiment
n_sim <- 1000
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

## Poisson setup
models <- list()
for (i in seq_len(n_sim)) {
  print_iter(i)
  f <- t(rmvnorm(p1, sigma = K))
  x <- matrix(rpois(n * p1, lambda = exp(f)), n, p1)
  models[[i]] <- adonis(x ~ y[, i] + u, method = "bray")
}

p_vals <- sapply(models, function(x) { x$aov.tab$"Pr(>F)"[1] })

## Okay, so this is rejecting a lot
hist(p_vals, breaks = 40)
dir.create("data")
save(p_vals, models, file = "data/pois_exper.rda")

## Gaussian setup
models <- list()
for (i in seq_len(n_sim)) {
  print_iter(i)
  x <- t(rmvnorm(p1, sigma = K))
  models[[i]] <- adonis(x ~ y[, i] + u, method = "euclidean")
}

p_vals <- sapply(models, function(x) { x$aov.tab$"Pr(>F)"[1] })
hist(p_vals, breaks = 40)
save(p_vals, models, file = "data/gauss_exper.rda")

## What happens if we're just doing logistic regression? Things seem fine.
models <- list()
for (i in seq_len(n_sim)) {
  print_iter(i)
  x <- t(rmvnorm(p1, sigma = K))
  models[[i]] <- glm(y[, i] ~ x + u)
}

p_vals <- sapply(models, function(x) {
  coef(summary(x))[, "Pr(>|t|)"]
})
hist(p_vals[2, ], breaks = 20)
hist(p_vals[3, ], breaks = 20)
save(p_vals, models, file = "data/logit_exper.rda")
