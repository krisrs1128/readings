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
X <- asinh(t(otu_table(abt)@.Data))
X[] <- as.integer(round(X, 2) * 100)

times <- 4 * round(raw_times / 4)
times_mapping <- match(times, unique(times))
times <- unique(times)

## ----  run_model ----
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

## ---- prepare-beta ----
taxa <- data.table(tax_table(abt)@.Data)
taxa <- data.table(rsv = rownames(tax_table(abt)), tax_table(abt)@.Data)

beta_hat <- samples$beta %>%
  melt(
    varnames = c("iteration", "time", "rsv_ix"),
    value.name = "beta"
  )
beta_hat$rsv <- rownames(otu_table(abt))[beta_hat$rsv_ix]
beta_hat$time <- times[beta_hat$time]
beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  group_by(time) %>%
  mutate(prob = softmax(beta)) 

group_order <- sort(table(taxa$Taxon_5), decreasing = TRUE)
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = names(group_order))
beta_hat$rsv <- factor(taxa[beta_hat$rsv_ix]$rsv, levels = rownames(tax_table(abt)))

## ---- unigram-series ----
plot_opts <- list(
  "x" = "time",
  "y" = "mean_prob",
  "col" = "Taxon_5",
  "alpha" = 0.4,
  "group" = "rsv",
  "facet_terms" = c("Taxon_5", ".")
)
gglines(
  beta_hat %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:8]) %>%
  group_by(rsv, time) %>%
  summarise(mean_prob = mean(prob), Taxon_5 = Taxon_5[1]) %>%
  as.data.frame(),
  plot_opts
) +
  theme(
    strip.text.y = element_blank(),
    legend.position = "bottom"
  )

## save, to avoid recomputing in the future
write_feather(beta_hat, "beta_unigram.feather")

## ---- unigram-beta-boxplots ----
plot_opts <- list(
  "x" = "rsv",
  "y" = "prob",
  "fill" = "Taxon_5",
  "col" = "Taxon_5",
  "outlier.size" = 0.01,
  "alpha" = 0.4,
  "facet_terms" = c("time", "Taxon_5"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x"
)
ggboxplot(
  beta_hat %>%
  filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:8]) %>%
  as.data.frame(),
  plot_opts
) +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )
