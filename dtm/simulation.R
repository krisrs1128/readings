#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This describes a basic simulation using the dynamic unigram model

## ---- setup ----
# List of packages for session
.packages  <-  c(
  "data.table",
  "plyr",
  "dplyr",
  "knitr",
  "reshape2",
  "ggplot2",
  "rstan"
)

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Load packages into session
sapply(.packages, require, character.only = TRUE)
set.seed(05122016)
cat("\014")  # Clear console
rm(list = ls()) # Delete all existing variables
graphics.off() # Close all open plots

# minimal theme for ggplots
theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

source("./unigram.R")

# Code Block -------------------------------------------------------------------

# generate data
N <- 100
V <- 15
T <- N
times <- 1:N
times_mapping <- times
sigma <- 0.1

beta <- topic_params(V, N, sigma = sigma)
X <- word_counts(beta, rep(500, N))
head(X)

stan_data <- list(
  N = N,
  V = V,
  T = T,
  sigma = sigma,
  times = times,
  times_mapping = times_mapping,
  X = X
)

m <- stan_model("unigram.stan")
stan_fit <- vb(m, data = stan_data)
samples <- extract(stan_fit)

beta_hat <- samples$beta %>%
  apply(c(2, 3), mean)
