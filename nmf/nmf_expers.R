#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Experiments with NMF, under different parameterizations and fitting
## techniques. Here, we consider the parameterizations in the parameterizations
## table and then save the resulting fits in a fits/ directory. These are
## visualized separately.
##
## author: kriss1@stanford.edu

## ---- libraries ----
library("SLURMHelpers")
source("./nmf_utils.R")

## ---- configuration ----
## create the configuration JSON file

## loop over unique values in the "batch" field of the json file

## within each value of the loop, send of a SLURM job that runs the required
## stan jobs
