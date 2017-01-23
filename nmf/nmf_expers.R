#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Experiments with NMF, under different parameterizations and fitting
## techniques. Here, we consider the parameterizations in the parameterizations
## table and then save the resulting fits in a fits/ directory. These are
## visualized separately.
##
## author: kriss1@stanford.edu

## ---- libraries ----
library("plyr")
library("dplyr")
library("data.table")
library("rstan")
source("./nmf_utils.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(01112017)

## ---- configuration ----
