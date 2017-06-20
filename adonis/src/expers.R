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
    opts <- modifyList(opts, list("n" = 100, "p" = 10))
    mod <- gaussian_null(opts, permutations = n_perm)
    f_stats[[i]] <- perm_data(mod, opts)
    i <- i + 1
  }
}

## ---- low-rank ----
opts <- modifyList(opts, list("n" = 40, "p" = 2000, "K" = 10))
mod <- low_rank_null(opts, permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

## ---- negative-binomial ----
opts <- modifyList(
  opts,
  list("prob" = 0.01, "size" = 1, distribution = "nb", distance = "bray")
)
mod <- nb_null(opts, permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("prob" = 0.01, "size" = 1))
mod <- nb_null(opts, permutations = n_perm / 100)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("prob" = 0.01, "size" = 100))
mod <- nb_null(opts, permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("prob" = 0.01, "size" = 100))
mod <- nb_null(opts, permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

## ---- many-factors ----
opts <- modifyList(opts, list("p_levels" = rep(list(c(0.5, 0.5)), 500)))
mod <- gaussian_null(opts, permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

## ---- plot-all ----
f_stats <- do.call(bind_rows, f_stats)

ggplot(
  f_stats %>%
  filter(distribution == "gaussian", is.na(K))
) +
  geom_histogram(
    aes(x = f_perm)
  ) +
  facet_wrap(n~p, scales = "free")
