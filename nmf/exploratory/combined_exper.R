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
fit_ids <- stringr::str_extract(fits, "[0-9]+")
exper_ids <- sapply(expers, function(x) { x$id })

theta_fits <- list()
for (i in seq_along(fits)) {
  ## retrieved simulated thetas
  set.seed(01112017)
  cur_exper <- expers[[which(exper_ids == fit_ids[[i]])]]
  cur_data <- nmf_sim(cur_exper$sim_opts)

  ## retrieve fitted thetas
  theta_fits[[i]] <- reshape_samples(
    get(load(fits[[i]]))$theta,
    cur_data$theta,
    c("i", "k")
  )

  ## join in simulation parameters
  cur_config <- data.frame(c(expers[[i]]$sim_opts, expers[[i]]$model_opts))
  theta_fits[[i]] <- cbind(theta_fits[[i]], cur_config)
}

theta_fits <- rbindlist(theta_fits)
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

gamma_pois_data <- theta_fits %>%
  filter(zero_inf_prob == 0, method == "nmf_gamma_poisson.stan")

theta_plots <- scores_contours(gamma_pois_data, plot_opts)
theta_plots$grouped
