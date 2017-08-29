#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Visualize results of HMC sampling for GP hyperparameters
##
## author: sankaran.kris@gmail.com
## date: 08/29/2017

library("tidyverse")

## Load data
train_data <- read_tsv("train_data.csv", col_names = F) %>%
  rename(x = X, y = X2)
post <- read_tsv("post_data.csv", col_names = F) %>%
  rename(iter = X1, x_new = X2, mu = X3, var = X4)

## Plot
ggplot() +
  geom_line(
    data = post,
    aes(x = x_new, y = mu, alpha = iter)
  ) +
  geom_line(
    data = post,
    aes(x = x_new, y = mu + 1.96 * var, alpha = iter),
    col = "red"
  ) +
  geom_line(
    data = post,
    aes(x = x_new, y = mu - 1.96 * var, alpha = iter),
    col = "red"
  ) +
  geom_point(
    data = train_data,
    aes(x = x, y = y)
  ) +
  scale_alpha(range = c(0, 0.03))
