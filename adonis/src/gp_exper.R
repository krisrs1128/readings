#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## An experiment with Adonis inspired by
##
## Dismantling the Mantel tests (Guillot & Rousset)
##
## author: sankaran.kris@gmail.com
## date: 07/31/2017

###############################################################################
## libraries and setup
###############################################################################
library("kernlab")
library("expm")
library("vegan")
library("tidyverse")

print_iter <- function(i) {
  if (i %% 10 == 0) {
    cat(sprintf("iteration %s\n", i))
  }
}

## variables used throughout experiment
n_sim <- 1000
p1 <- 10
p2 <- 2
n <- 100
u <- matrix(runif(n * p, 0, 20), n, p2)

K <- kernelMatrix(rbf, u)
sqrtK <- Re(sqrtm(K))
y <- sqrtK %*% matrix(rnorm(n * n_sim), n, n_sim)
y <- y > 0

## Poisson setup
models <- list()
p_vals <- vector(length = n_sim)

for (i in seq_len(n_sim)) {
  print_iter(i)
  f <- sqrtK %*% matrix(rnorm(n * p1), n, p1)
  x <- matrix(rpois(n * p1, lambda = exp(f)), n, p1)
  models[[i]] <- adonis(x ~ y[, i], method = "bray")
}

p_vals <- lapply(models, function(x) { x$aov.tab$"Pr(>F)"[1] })

## this actuall seems conservative in this case?
hist(p_vals, breaks = 20)
save(p_vals, models, file = "data/pois_exper.rda")

## Gaussian setup
models <- list()
p_vals <- vector(length = n_sim)
for (i in seq_len(n_sim)) {
  print_iter(i)
  x <- sqrtK %*% matrix(rnorm(n * p1), n, p1)
  y <- as.integer(sqrtK %*% rnorm(n) > 0)
  models[[i]] <- adonis(x ~ y[, i], method = "euclidean")
}

p_vals <- sapply(models, function(x) { x$aov.tab$"Pr(>F)"[1] })
hist(p_vals, breaks = 20)
save(p_vals, models, file = "data/gauss_exper.rda")
