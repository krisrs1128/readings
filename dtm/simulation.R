#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This describes a basic simulation using the dynamic unigram model

## ---- setup ----
# List of packages for session
.packages  <-  c("data.table", "plyr", "dplyr", "knitr", "reshape2", "ggplot2")

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
beta <- topic_params(15, 1000, sigma = .1)
X <- word_counts(beta, rep(500, 1000))
head(X)
