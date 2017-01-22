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
## ---- get-data ----
data(abt)
n_samples <- ncol(otu_table(abt))
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .4 * n_samples, prune = TRUE) %>%
  subset_samples(ind == "D")
hist(colSums(otu_table(abt)), 15)

## ---- heatmaps ----
heatmap(otu_table(abt))
heatmap(asinh(otu_table(abt)))

## ----  lda ----
m <- stan_model(file = "lda_counts.stan")
X <- asinh(t(otu_table(abt)@.Data))
X[] <- as.integer(round(X, 2) * 100)

stan_data <- list(
  K = 4,
  V = ncol(X),
  D = nrow(X),
  n = X,
  alpha = rep(1e-10, 4)
)

stan_fit <- vb(m, stan_data, eta = .1, adapt_engaged = FALSE, grad_samples = 1, iter = 2000)
samples <- rstan::extract(stan_fit)

## ---- extract_beta ----
# underlying RSV distributions
beta_hat <- samples$beta %>%
  melt() %>%
  setnames(c("iterations", "cluster", "rsv_ix", "beta"))
beta_hat$rsv <- rownames(tax_table(abt))[beta_hat$rsv_ix]

taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

beta_hat <- beta_hat %>%
  left_join(taxa)

sorted_taxa <- names(sort(table(beta_hat$Taxon_5), decreasing = TRUE))
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = sorted_taxa)
beta_hat$rsv <- factor(beta_hat$rsv, levels = taxa$rsv)

## ---- visualize_beta ----
# might want to set prior for more extreme decay

## ---- extrac_theta ----
samples_theta <- melt(samples$theta)
theta_hat <- samples_theta %>%
  group_by(Var2, Var3) %>%
  summarise(mean = mean(value))
theta_hat$Var2 <- rownames(X)[theta_hat$Var2]
colnames(theta_hat) <- c("sample", "cluster", "theta")

sample_info <- sample_data(abt)
sample_info$sample <- rownames(sample_info)
theta_hat$cluster <- paste("Cluster", theta_hat$cluster)

theta_hat <- theta_hat %>%
  left_join(sample_info, by = "sample")

## ---- visualize_theta ----
plot_opts <- list(
  "x" = "time",
  "y" = "cluster",
  "fill" = "theta"
)
p1 <- ggheatmap(theta_hat, plot_opts)

plot_opts <- list(
  "x" = "as.factor(time)",
  "y" = "theta",
  "facet_terms" = c("cluster", "time")
)
p1 <- ggboxplot(data.frame(theta_hat), plot_opts)

## ---- save_results ----
dir.create("results")
write_feather(theta_hat, path = "theta.feather")
write_feather(beta_hat, path = "beta.feather")
