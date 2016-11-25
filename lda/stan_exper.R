#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# Generate data for multinomial mixture model, and run rstan code for it
# Based on
# https://github.com/stan-dev/stan/releases/download/v2.12.0/stan-reference-2.12.0.pdf
# page 157

## ---- setup ----
library("rstan")
library("data.table")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
source("./lda_counts.R")
set.seed(11242016)

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

# Code Block -------------------------------------------------------------------

# generate data
D <- 50
K <- 5
V <- 300
words_per_doc <- rep(2500, D)
alpha <- rep(1, K)
gamma <- rep(1, V)
lda_data <- generate_data(D, words_per_doc, alpha, gamma)

Nv <- lda_data$Nv %>%
  dcast(document ~ word, fill = 0) %>%
  select(-document) %>%
  as.matrix()

heatmap(Nv)

# run STAN model
stan_data <- list(
  K = K,
  V = V,
  D = D,
  n = Nv,
  alpha = alpha
)

m <- stan_model(file = "lda_counts.stan")
stan_fit <- vb(m, stan_data)

# get samples
samples <- extract(stan_fit)

# beta
samples_beta <- melt(samples$beta)

beta_hat <- samples_beta %>%
  group_by(Var2, Var3) %>%
  summarise(mean = mean(value)) %>%
  dcast(Var2 ~ Var3)

# which estimated topics line up with which true ones?
lda_data$beta
pairs(cbind(t(beta_hat[, -1]), t(lda_data$beta)))
match_ix <- apply(cor(t(lda_data$beta), t(beta_hat[, -1])), 1, which.max)

# recovered cluster memberships
samples_theta <- melt(
  samples$theta[,, match_ix],
  varnames= c("iteration", "n", "k")
)

theta_hat <- samples_theta %>%
  group_by(n, k) %>%
  summarise(mean = mean(value, na.rm = TRUE)) %>%
  mutate(type = "estimate") %>%
  data.table()

theta <- lda_data$theta %>%
  melt(varnames = c("n", "k"), value.name = "mean") %>%
  mutate(type = "truth")

theta_hat <- rbind(theta, theta_hat) %>%
  dcast(n + k  ~ type, value.var ="mean")

ggplot(theta_hat) +
  geom_point(aes(x = truth, y = estimate, col = as.factor(k)), size = .5) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed()
