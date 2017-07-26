#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some experiments using the scan sampler on simulated and real data
##
## author: sankaran.kris@gmail.com
## date: 07/18/2017

library("reshape2")
library("tidyverse")
library("abind")
theme_set(ggscaffold::min_theme())

## A simulated examples
A <- diag(0.7, nrow = 1)
C <- diag(1, nrow = 1)
Q <- diag(1, nrow = 1)
R <- diag(1, nrow = 1)
set.seed(1000)
res <- simulate(A, C, 0, Q, R)
y <- res$y
y[y < 0] <- 0
y_samples <- sampler(y, A, C, Q, R, n_iter = 100)
y_samples <- abind(res$y, y_samples, along = 3)
my <- melt(
  y_samples,
  varnames = c("i", "p", "iter"),
  value.name = "y"
)

p <- ggplot() +
  geom_point(
    data = my %>% filter(y < 0),
    aes(x = i, y = y),
    size = .5, alpha = 0.05
  ) +
  geom_point(
    data = my %>% filter(iter == 2),
    aes(x = i, y = y),
    col = "red", size = 1, alpha = 1
  ) +
  geom_point(
    data = my %>% filter(iter == 1),
    aes(x = i, y = y),
    col = "red", size = 1, alpha = 0.5
  )

p + ylim(-10, 5)

## antibiotics data
library("treelapse")
library("phyloseq")
data(abt)
abt <- abt %>%
  filter_taxa(function(x) var(x) > 4, prune = TRUE) %>%
  subset_samples(ind == "F") %>%
  transform_sample_counts(asinh)
y <- get_taxa(abt)[3, ]

y_samples <- sampler(y, A, C, Q, R, n_iter = 100)

y_samples <- abind(y, y_samples, along = 3)
my <- melt(
  y_samples,
  varnames = c("i", "p", "iter"),
  value.name = "y"
)

p <- ggplot() +
  geom_point(
    data = my %>% filter(y < 0),
    aes(x = i, y = y),
    size = .5, alpha = 0.05
  ) +
  geom_point(
    data = my %>% filter(iter == 2),
    aes(x = i, y = y),
    col = "red", size = 1, alpha = 1
  ) +
  geom_point(
    data = my %>% filter(iter == 1),
    aes(x = i, y = y),
    col = "red", size = 1, alpha = 0.5
  )

p + ylim(-15, 5)