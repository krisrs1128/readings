#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## We don't want to run all the experiment configurations in parallel or series,
## we need something intermediate. This script runs a subset of the experiments
## and can be called by a SLURM submitter.
##
## author: kriss1@stanford.edu

## ---- libraries ----
library("jsonlite")
source("/scratch/users/kriss1/programming/readings/nmf/src/nmf_utils.R")

## ---- parse-args ----
args <- commandArgs(trailingOnly = TRUE)
expers <- fromJSON(args[[1]])
subset_ix <- as.integer(args[[2]])

for (i in seq_along(expers)) {
  set.seed(01112017)
  if (expers[[i]]$batch != i) next
  cur_data <- sim_data(expers[[i]]$sim_opts)
  cur_fit <- fit_model(cur_data, expers[[i]]$model_opts)
  save(cur_fit, file = expers[[i]]$output_dir)
}
