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

source("./unigram.R")

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

times <- 5 * round(raw_times / 5)
times_mapping <- match(times, unique(times))
times <- unique(times)

## ----  run_model ----
m <- stan_model(file = "lda_counts.stan")

N <- nrow(X)
V <- ncol(X)
T <- length(times)
sigma <- .1

stan_data <- list(
  N = N,
  V = V,
  T = T,
  sigma = sigma,
  times = times,
  times_mapping = times_mapping,
  X = X
)

m <- stan_model("unigram.stan")
stan_fit <- vb(m, data = stan_data)
samples <- rstan::extract(stan_fit)

## ---- visualize_beta ----
beta_hat <- samples$beta %>%
  apply(c(2, 3), mean) %>%
  melt(varnames = c("time", "i"))
beta_sd <- samples$beta %>%
  apply(c(2, 3), sd) %>%
  melt(varnames = c("time", "i"), value.name = "sd")
beta_hat <- beta_hat %>%
  left_join(beta_sd) %>%
  mutate(type = "estimate")
beta_hat$time <- times[beta_hat$time]
beta_hat$rsv <- rownames(otu_table(abt))[beta_hat$i]

ggplot(beta_hat %>% filter(word < 20)) +
  geom_line(aes(x = time, y = value, col = type)) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~rsv)

ggplot(beta_hat) +
  geom_tile(aes(x = rsv, y = time, fill = value))

## ---- view_taxa ----
taxa <- tax_table(abt)
beta_hat$group <- taxa@.Data[beta_hat$i, "Taxon_5"]
beta_hat$rsv <- factor(beta_hat$rsv, levels = rownames(taxa@.Data))

ggplot(beta_hat) +
  geom_line(aes(x = time, y = value, group = rsv), alpha = 0.1) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~group)

ggplot(beta_hat) +
  geom_tile(aes(x = time, y = rsv, fill = value)) +
  facet_wrap(~group, scale = "free_y")
