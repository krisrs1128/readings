#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This is an application of the dynamic unigram model to the antibiotics data

## ---- setup ----
library("rstan")
library("data.table")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("phyloseq")
library("treelapse")
library("feather")
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
## ---- get_data ----
data(abt)
n_samples <- ncol(otu_table(abt))
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .4 * n_samples, prune = TRUE) %>%
  subset_samples(ind == "D")
hist(colSums(otu_table(abt)), 15)

## ---- vis_times ----
raw_times <- sample_data(abt)$time
X <- log(1 + t(otu_table(abt)@.Data))
X[] <- as.integer(round(X, 2) * 100)

times <- 4 * round(raw_times / 4)
times_mapping <- match(times, unique(times))
times <- unique(times)

## ----  run_model ----
N <- nrow(X)
V <- ncol(X)
T <- length(times)
sigma <- 0.1
delta <- 0.1

stan_data <- list(
  N = N,
  V = V,
  T = T,
  K = 4,
  sigma = sigma,
  delta = delta,
  times = times,
  times_mapping = times_mapping,
  X = X
)

m <- stan_model("dtm.stan")
stan_fit <- vb(m, data = stan_data)
samples <- rstan::extract(stan_fit)

## ---- visualize_theta ----
alpha_hat <- apply(samples$alpha, c(2, 3), mean)
theta_hat <- t(apply(alpha_hat, 1, softmax))

theta_hat <- cbind(
  sample = colnames(otu_table(abt)),
  theta_hat[times_mapping, ],
  sample_data(abt)
) %>%
  melt(
    id.vars = c("sample", "ind", "condition", "time"),
    variable.name = "cluster"
  )

ggplot(theta_hat) +
  geom_tile(aes(x = sample, y = cluster, fill = value)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#5BBABA", limits = c(0, 1))

ggplot(theta_hat) +
  geom_line(aes(x = time, y = value, col = cluster))
