#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Plot results of mix gp on simulated data
##
## author: sankaran.kris@gmail.com
## date: 09/13/2017

library("tidyverse")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

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

post <- read_tsv("data/post_mix.csv", col_names = FALSE)
data <- read_tsv("data/mix_data.csv", col_names = FALSE)

ggplot() +
  geom_line(
    data = post,
    aes(x = X2, y = X3, group = X1),
  ) +
  geom_point(
    data = data,
    aes(x = X1, y = X2)
  ) +
  ylim(-30, 10)

bump <- read_tsv("data/bump_data.csv", col_names = FALSE)
ggplot(bump) +
  geom_point(aes(x = X1, y = X2))
