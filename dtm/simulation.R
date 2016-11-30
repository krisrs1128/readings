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

## ---- generate_data ----
N <- 2000
V <- 15
T <- N / 20
times <- 1:T
times_mapping <- rep(1:T, each = 20)
sigma <- 1

beta <- topic_params(V, T, sigma = sigma)
X <- word_counts(beta[times_mapping, ], rep(500, N))
head(X)

## ---- run_model ----
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

## ---- visualize_beta ----
beta_hat <- samples$beta %>%
  apply(c(2, 3), mean) %>%
  melt(varnames = c("sample", "word"))
beta_sd <- samples$beta %>%
  apply(c(2, 3), sd) %>%
  melt(varnames = c("sample", "word"), value.name = "sd")
beta_hat <- beta_hat %>%
  left_join(beta_sd) %>%
  mutate(type = "estimate")

mbeta <- beta %>%
  melt(varnames = c("sample", "word")) %>%
  mutate(type = "truth", sd = NA)
mbeta <- rbind(mbeta, beta_hat)

ggplot(mbeta) +
  geom_line(aes(x = sample, y = value, col = type)) +
  geom_line(aes(x = sample, y = value - 1.96 * sd, col = type), linetype = 2) +
  geom_line(aes(x = sample, y = value + 1.96 * sd, col = type), linetype = 2) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~word)
