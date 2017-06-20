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
opts <- list(
  "n" = 20,
  "p" = 2,
  "mu" = 0,
  "sigma" = 1,
  "p_levels" = list(c(0.5, 0.5), c(0.5, 0.5), c(0.9, 0.1))
)

mod <- gaussian_null(opts, method = "euclidean", permutations = 10000)
plot_perms(mod)

opts <- modifyList(opts, list("n" = 100, "p" = 10))
mod <- gaussian_null(opts, method = "euclidean", permutations = 10000)
plot_perms(mod)


opts <- modifyList(opts, list("n" = 20, "p" = 200))
mod <- gaussian_null(opts, method = "euclidean", permutations = 10000)
plot_perms(mod)

opts <- modifyList(opts, list("n" = 20, "p" = 2000))
mod <- gaussian_null(opts, method = "euclidean", permutations = 10000)
plot_perms(mod)

## ---- low-rank ----
opts <- modifyList(opts, list("n" = 20, "p" = 2000, "K" = 10))
mod <- low_rank_null(opts, method = "euclidean", permutations = 10000)
plot_perms(mod)
