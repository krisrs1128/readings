#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some experiments on nonparametric ANOVA under a few different settings.
##
## author: sankaran.kris@gmail.com
## date: 06/19/2017

## ---- libraries ----
library("tidyverse")
library("vegan")
source("utils.R")
theme_set(ggscaffold::min_theme())

## ---- gaussian ----
n_perm <- 1000
opts <- list(
  "mu_x" = 0,
  "sigma_x" = 1,
  "p_levels" = list(c(0.5, 0.5), c(0.5, 0.5), c(0.8, 0.2)),
  "distribution" = "gaussian",
  "distance" = "euclidean"
)

f_stats <- list()
i <- 1
for (n in seq(40, 1000, length.out = 5)) {
  for (p in seq(2, 2000, length.out = 5)) {
    cat(sprintf("regime n: %f \t p: %f\n", n, p))
    opts <- modifyList(opts, list("n" = n, "p" = p))
    mod <- gaussian_null(opts, permutations = n_perm)
    f_stats[[i]] <- perm_data(mod, opts)
    i <- i + 1
  }
}

ggplot(f_stats) +
  geom_histogram(aes(x = f_perm)) +
  facet_grid(n ~ p)

## ---- low-rank ----
f_list <- list()
i <- 1
opts$n <- 280
opts$p <- 2000
for (k in seq(2, 20, length.out = 5)) {
  cat(sprintf("regime k: %f\n", k))
  opts <- modifyList(opts, list("K" = k))
  mod <- gaussian_null(opts, permutations = n_perm)
  f_stats[[i]] <- perm_data(mod, opts)
  i <- i + 1
}

ggplot(f_stats) +
  geom_histogram(aes(x = f_perm)) +
  facet_grid(~ k)

## ---- num-perms ----
f_list <- list()
i <- 1
for (nperm in 10 ^ (2:5)) {
  cat(sprintf("regime n_perm: %f\n", n_perm))
  opts <- modifyList(opts, list("K" = k))
  mod <- gaussian_null(opts, permutations = n_perm)
  f_stats[[i]] <- perm_data(mod, opts)
  i <- i + 1
}

ggplot(f_stats) +
  geom_histogram(aes(x = f_perm)) +
  facet_grid(~ n_perm)

## ---- negative-binomial ----
f_list <- list()
i <- 1
for (p in seq(0.01, 0.9, length.out = 5)) {
  for (size in seq(1, 1000, length.out = 5)) {
    opts <- modifyList(opts, list("prob" = 0.01, "size" = 100))
    mod <- nb_null(opts, permutations = n_perm)
    f_stats[[i]] <- perm_data(mod, opts)
    i <- i + 1
  }
}

ggplot(f_stats) +
  geom_histogram(aes(x = f_perm)) +
  facet_grid(~ p)

## ---- many-factors ----
i <- 1
f_list <- list()
for (p2 in seq(2, 500, length.out = 4)) {
  cat(sprintf("regime p2: %f\n", p2))
  opts <- modifyList(opts, list("p_levels" = rep(list(c(0.5, 0.5)), p2)))
  mod <- gaussian_null(opts, permutations = n_perm)
  f_stats[[i]] <- perm_data(mod, opts)
  i <- i + 1
}

ggplot(f_stats) +
  geom_histogram(aes(x = f_perm)) +
  facet_grid(~ p2)
