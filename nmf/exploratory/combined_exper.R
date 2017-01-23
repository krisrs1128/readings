#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Visualize a combined experiment across various NMF / fitting
## parameterizations.
##
## author: kriss1@stanford.edu

## ---- libraries ----
## assumed running from NMF directory
library("data.table")
source(file.path("src", "nmf_utils.R"))

## ---- theta-reshape ----
## extract theta information from the fits
fits <- list.files("fits", "fit-*", full.names = TRUE)
expers <- fromJSON(
  file.path("batch", "config.json"),
  simplifyVector = FALSE
)

theta_fits <- reshape_all_samples(
  fits,
  file.path("batch", "config.json"),
  "theta",
  c("i", "k")
)
theta_fits$method <- basename(as.character(theta_fits$method))

## ---- visualize-thetas -----
## Visualize the fitted thetas, according to a few different simulation properties
plot_opts <- list(
  "x" = "value_1",
  "y" = "value_2",
  "fill" = "log(..level..)",
  "fill_type" = "gradient",
  "facet_terms" = c("N", "inference", "P"),
  "group" = "i",
  "alpha" = 0.05,
  "h" = 0.1,
  "mean_col" = "#e34a33",
  "x_lim" = c(0, 3.5),
  "y_lim" = c(0, 5.4),
  "text_size" = 2,
  "panel_border" = 0.2
)

## first, visualization in the non-zero-inflated case
gamma_pois_data <- theta_fits %>%
  filter(zero_inf_prob == 0, method == "nmf_gamma_poisson.stan")

theta_plots <- scores_contours(gamma_pois_data, plot_opts)
##theta_plots$grouped

zinf_data <- theta_fits %>%
  filter(zero_inf_prob != 0, P == 75, N == 100)
plot_opts$facet_terms <- c("zero_inf_prob", "inference", "method")
theta_plots <- scores_contours(zinf_data, plot_opts)
##theta_plots$grouped

## ---- visualize-betas ----
beta_fits <- reshape_all_samples(
  fits,
  file.path("batch", "config.json"),
  "beta",
  c("v", "k")
)
beta_fits$method <- basename(as.character(beta_fits$method))

plot_opts$facet_terms <- c("N", "inference", "P")
plot_opts$group <- "v"

gamma_pois_data <- beta_fits %>%
  filter(zero_inf_prob == 0, method == "nmf_gamma_poisson.stan")

beta_plots <- scores_contours(gamma_pois_data, plot_opts)
ggsave("~/test.png", beta_plots$grouped)

zinf_data <- beta_fits %>%
  filter(zero_inf_prob != 0, P == 75, N == 100)
plot_opts$facet_terms <- c("zero_inf_prob", "inference", "method")
theta_plots <- scores_contours(zinf_data, plot_opts)
##theta_plots$grouped
