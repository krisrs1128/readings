#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## A simple bootstrap experiment in variable transformation, as described in "An
## Introduction to the boostrap". The idea is to automatically recover the usual
## VST for the correlation in a bivariate normal model automatically using the
## bootstrap.

# Setup packages ---------------------------------------------------------------
library("ggplot2")

## Load packages into session
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  text = element_text(family = "Ubuntu Regular", color = "#22211d"),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

## ---- utils ----
generate_correlated <- function(N, rho) {
  z <- rnorm(N)
  cbind(X1 = z, X2 = rho * z + rnorm(N, 0, sqrt(1 - rho ^ 2)))
}

## ---- simulate ----
N <- 10
rho <- .4
R <- 15000

cors <- vector(length = R)
for (i in seq_len(R)) {
  X <- generate_correlated(N, rho)
  cors[i] <- cor(X)[1, 2]
}

## ---- usual-transformation ----
ggplot(data.frame(cors)) +
  geom_histogram(aes(x = cors), binwidth = 0.01)
ggplot(data.frame(cors)) +
  geom_histogram(aes(x = atanh(cors)), binwidth = 0.01)
