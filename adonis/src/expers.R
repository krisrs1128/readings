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
n_perm <- 100000
opts <- list(
  "n" = 20,
  "p" = 2,
  "mu_x" = 0,
  "sigma_x" = 1,
  "p_levels" = list(c(0.5, 0.5), c(0.5, 0.5), c(0.9, 0.1))
)
f_stats <- list()

mod <- gaussian_null(opts, method = "euclidean", permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("n" = 100, "p" = 10))
mod <- gaussian_null(opts, method = "euclidean", permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("n" = 20, "p" = 200))
mod <- gaussian_null(opts, method = "euclidean", permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("n" = 20, "p" = 2000))
mod <- gaussian_null(opts, method = "euclidean", permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

## ---- low-rank ----
opts <- modifyList(opts, list("n" = 20, "p" = 2000, "K" = 10))
mod <- low_rank_null(opts, method = "euclidean", permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

## ---- negative-binomial ----
opts <- modifyList(opts, list("prob" = 0.01, "size" = 1))
mod <- nb_null(opts, permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("prob" = 0.01, "size" = 1))
mod <- nb_null(opts, permutations = n_perm / 100)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("prob" = 0.01, "size" = 100))
mod <- nb_null(opts, permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))

opts <- modifyList(opts, list("prob" = 0.01, "size" = 100))
mod <- nb_null(opts, permutations = n_perm, method = "euclidean")
f_stats <- c(f_stats, list(perm_data(mod, opts)))

## ---- many-factors ----
opts <- modifyList(opts, list("p_levels" = rep(list(c(0.5, 0.5)), 500)))
mod <- gaussian_null(opts, method = "euclidean", permutations = n_perm)
f_stats <- c(f_stats, list(perm_data(mod, opts)))
