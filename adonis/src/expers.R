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
for (n in seq(40, 250, length.out = 5)) {
  for (p in seq(2, 100, length.out = 5)) {
    cat(sprintf("regime n: %f \t p: %f\n", n, p))
    opts <- modifyList(opts, list("n" = n, "p" = p))
    mod <- gaussian_null(opts, permutations = n_perm)
    f_stats[[i]] <- perm_data(mod, opts)
    i <- i + 1
  }
}

p <- perm_histo(f_stats) +
  facet_grid(n ~ p, scales = "free")
save_fig(p, "vary_n_p.png")

## ---- low-rank ----
f_stats <- list()
i <- 1
opts$n <- 100
opts$p <- 1000
for (k in seq(2, 75, length.out = 6)) {
  cat(sprintf("regime k: %f\n", k))
  opts <- modifyList(opts, list("K" = k))
  mod <- low_rank_null(opts, permutations = n_perm)
  f_stats[[i]] <- perm_data(mod, opts)
  i <- i + 1
}

p <- perm_histo(f_stats) +
  facet_grid(K ~ .)
save_fig(p, "vary_rank.png")

## ---- num-perms ----
f_stats <- list()
i <- 1
for (cur_perm in 10 ^ (seq(2, 4, length.out = 4))) {
  cat(sprintf("regime cur_perm: %f\n", cur_perm))
  mod <- gaussian_null(opts, permutations = cur_perm)
  opts$n_perm <- cur_perm
  f_stats[[i]] <- perm_data(mod, opts)
  i <- i + 1
}

p <- perm_histo(f_stats) +
  facet_grid(n_perm ~ ., scales = "free")
save_fig(p, "vary_perms.png")

## ---- negative-binomial ----
f_stats <- list()
i <- 1
opts$distance <- "bray"
for (prob in seq(0.01, 0.99, length.out = 5)) {
  for (size in seq(1, 500, length.out = 5)) {
    cat(sprintf("regime p: %f \t size: %f \n", prob, size))
    opts <- modifyList(opts, list("prob" = prob, "size" = size))
    mod <- nb_null(opts, permutations = n_perm)
    f_stats[[i]] <- perm_data(mod, opts)
    i <- i + 1
  }
}

p <- perm_histo(f_stats) +
  facet_grid(prob ~ size, scales = "free")
save_fig(p, "vary_nb_params.png")

## ---- many-factors ----
i <- 1
f_stats <- list()
opts$distance <- "euclidean"
opts$n <- 500
for (p2 in seq(2, 450, length.out = 4)) {
  cat(sprintf("regime p2: %f\n", p2))
  opts$p_levels <- rep(list(c(0.5, 0.5)), p2)
  mod <- gaussian_null(opts, permutations = n_perm)
  f_stats[[i]] <- perm_data(mod, opts)
  i <- i + 1
}

p <- perm_histo(f_stats) +
  facet_wrap(~ p2, scales = "free")
save_fig(p, "vary_factors.png")
