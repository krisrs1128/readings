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

p1 <- 10
p2 <- 2
n <- 100
u <- matrix(runif(n * p, 0, 10), n, p2)

K <- kernelMatrix(rbf, u)
sqrtK <- Re(sqrtm(K))

n_sim <- 1000
models <- list()
p_vals <- vector(length = n_sim)
for (i in seq_len(n_sim)) {
  if (i %% 10 == 0) {
    cat(sprintf("iteration %s\n", i))
  }
  f <- sqrtK %*% matrix(rnorm(n * p1), n, p1)
  x <- matrix(rpois(n * p1, lambda = exp(f)), n, p1)
  y <- as.integer(sqrtK %*% rnorm(n) > 0)
  models[[i]] <- adonis(x ~ y, data = df, method = "bray")
  p_vals[i] <- models[[i]]$aov.tab$"Pr(>F)"[1]
}

hist(p_vals, breaks = 20)
