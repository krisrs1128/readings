#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Visualize results of HMC sampling for GP hyperparameters
##
## author: sankaran.kris@gmail.com
## date: 08/29/2017

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

## Load data
train_data <- read_tsv("data/train_data.tsv", col_names = F) %>%
  rename(x = X1, y = X2)
post <- read_tsv("data/post_data.tsv", col_names = F) %>%
  rename(iter = X1, x_new = X2, mu = X3, var = X4)

## Plot
ggplot() +
  geom_line(
    data = post,
    aes(x = x_new, y = mu, group = iter),
    alpha = 0.05
  ) +
  geom_line(
    data = post,
    aes(x = x_new, y = mu + 1.96 * sqrt(var), group = iter),
    alpha = 0.05, col = "red"
  ) +
  geom_line(
    data = post,
    aes(x = x_new, y = mu - 1.96 * sqrt(var), group = iter),
    alpha = 0.05, col = "red"
  ) +
  geom_point(
    data = train_data,
    aes(x = x, y = y)
  )
