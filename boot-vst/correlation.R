#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## A simple bootstrap experiment in variable transformation, as described in "An
## Introduction to the boostrap". The idea is to automatically recover the usual
## VST for the correlation in a bivariate normal model automatically using the
## bootstrap.

# Setup packages ---------------------------------------------------------------
library("ggplot2")

## Load packages into session
sapply(.packages, require, character.only = TRUE)
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
