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
## get original thetas
set.seed(01112017)
cur_data <- nmf_sim(expers[[i]]$sim_opts)

## extract theta information from the fits
fits <- list.files(pattern = "fits/fit*.rda", full.names = TRUE)
theta_fits <- list()
for (i in seq_along(fits)) {
  theta_fits[[i]] <- reshape_samples(
    get(load(fits[[i]]))$theta,
    cur_data$theta,
    c("i", "k")
  )
}

head(theta_fits[[10]])
