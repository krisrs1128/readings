#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Experiment running gamma-poisson factorization model, on simulated data.
##
## reference:https://arxiv.org/pdf/1506.03431.pdf

## ---- libraries ----
library("plyr")
library("dplyr")
library("data.table")
library("reshape2")
library("ggplot2")
library("ggscaffold")
library("rstan")
source("./nmf_utils.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(01082017)

## ---- simulate ----
sim_data <- nmf_sim(list(zero_inf_prob = 0.8))
y <- sim_data$y
theta <- sim_data$theta
beta <- sim_data$beta

## ---- heatmap ----
y_df <- y %>%
  melt(
    varnames = c("row", "col"),
    value.name = "fill"
  )
plot_opts <- list(
  "x" = "row",
  "y" = "col",
  "x_order" = order(theta[, 1]),
  "y_order" = order(beta[, 1])
)

ggheatmap(y_df, plot_opts)

## ---- pca ----
compare_data <- data.frame(
  theta,
  princomp(scale(y))$scores
)

theme_set(min_theme())
ggplot(compare_data) +
  geom_point(aes(x = Comp.1, Comp.2, size = X1, col = X2))

## ---- overdispersion ----
yy <- sort(rpois(N * P, mean(y)))
qq_df <- data.frame(
  y = c(sort(y), yy),
  label = c(rep("zinf-gamma-poisson", N * P), rep("poisson", N * P))
)

ggplot(qq_df) +
  geom_histogram(aes(x = y), binwidth = .5) +
  facet_grid(label ~ .)

ggplot(data.frame(
  mu = rowMeans(y),
  sigma = apply(y, 1, sd)
)) +
  geom_point(
    aes(x = mu, y = sigma)
  ) +
  coord_fixed() +
  geom_abline(slope = 1)

## ---- stan-fit ----
fit <- extract(
  stan(file = "nmf_gamma_poisson_zero.stan", data = stan_data, chains = 1)
)
save(fit, file = "nmf_zero.rda")

## ---- examine ----
theta_fit <- reshape_samples(fit$theta, theta, c("i", "k"))

plot_opts <- list(
  "x" = "value_1",
  "y" = "value_2",
  "fill" = "log(..level..)",
  "fill_type" = "gradient",
  "group" = "i",
  "alpha" = 0.05,
  "h" = 0.3,
  "mean_col" = "#e34a33",
  "x_lim" = c(0, 6),
  "y_lim" = c(0, 8),
  "text_size" = 2,
  "panel_border" = 0.2
)

theta_plots <- scores_contours(theta_fit, plot_opts)
theta_plots$grouped

## ---- plot-beta ----
beta_fit <- reshape_samples(fit$beta, beta, c("v", "k"))

plot_opts$group <- "v"
beta_plots <- scores_contours(beta_fit, plot_opts)
beta_plots$grouped
