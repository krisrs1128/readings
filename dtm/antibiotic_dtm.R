#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This is an application of the dynamic unigram model to the antibiotics data

## ---- setup ----
args <- commandArgs(trailingOnly = TRUE)
cur_ix <- args[1]

library("rstan")
library("data.table")
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("feather")
library("phyloseq")
library("treelapse")
library("ggscaffold")
set.seed(11242016)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Code Block -------------------------------------------------------------------
## ---- get_data ----
data(abt)
n_samples <- ncol(otu_table(abt))
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .4 * n_samples, prune = TRUE) %>%
  subset_samples(ind == "D")

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

param_grid <- expand.grid(
  sigma = c(.001, .005, .01, .05, .1, .5),
  delta = c(.001, .005, .01, .05, .1, .5),
  K = c(2, 3)
)

m <- stan_model("dtm.stan")
timestamp <- gsub("[^0-9]", "", Sys.time())
stan_data <- list(
  N = N,
  V = V,
  T = T,
  K = param_grid[cur_ix, "K"],
  sigma = param_grid[cur_ix, "sigma"],
  delta = param_grid[cur_ix, "delta"],
  times = times,
  times_mapping = times_mapping,
  X = X
)
print(timestamp)
stan_fit <- vb(m, data = stan_data)
dir.create("fits")
save(stan_fit, file = sprintf("fits/dtm_fit-%s.rda", cur_ix))

################################################################################
## Once we've found a saved model to visualize, we can use the code below. The
## model comparisons can be done by looking at some convergence diagnostics.
################################################################################

retrieve_ix <- 24 # fit that seems reasonably interpretable
stan_fit <- get(load(sprintf("fits/dtm_fit-%s.rda", retrieve_ix)))
samples <- rstan::extract(stan_fit)

## ---- visualize_theta ----
softmax <- function(mu) {
  exp(mu) / sum(exp(mu))
}

theta_hat <- apply(samples$alpha, c(1, 2), softmax) %>%
  melt(
    varnames = c("cluster", "iteration", "time"),
    value.name = "theta"
  ) %>%
  filter(cluster < param_grid[retrieve_ix, "K"])
theta_hat$time <- times[theta_hat$time]

cur_samples <- data.frame(sample_data(abt))
cur_samples$time <- 4 * round(cur_samples$time / 4)

theta_hat <- cur_samples %>%
  right_join(theta_hat)

plot_opts <- list(
  "x" = "as.factor(time)",
  "y" = "theta",
  "fill" = "as.factor(cluster)",
  "col" = "as.factor(cluster)",
  "col_colors" = brewer.pal(param_grid[retrieve_ix, "K"] - 1, "Set2"),
  "fill_colors" = brewer.pal(param_grid[retrieve_ix, "K"] - 1, "Set2"),
  "facet_terms" = c(".", "condition"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x"
)
p <- ggboxplot(theta_hat, plot_opts) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    fill = "Cluster",
    x = "time"
  )
ggsave(p, file = sprintf("figure/%s%s-theta.png", cur_ix, timestamp))
write_feather(theta_hat, "../lda/results/dtm_theta_hat.feather")

## ---- visualize_beta ----
beta_hat <- apply(samples$beta, c(1, 2, 3), softmax) %>%
  melt(
    varnames = c("rsv_ix", "iteration", "time", "cluster")
  )
beta_hat$time <- times[beta_hat$time]

## merge in taxonomic information (for labeling evolutionary families)
taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- factor(rownames(tax_table(abt)))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]
sorted_taxa <- names(sort(table(taxa$Taxon_5), decreasing = TRUE))
taxa$Taxon_5 <- factor(taxa$Taxon_5, levels = sorted_taxa)

beta_hat$rsv <- taxa$rsv[beta_hat$rsv_ix]
beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  filter(
    Taxon_5 %in% sorted_taxa[1:8],
    cluster < param_grid[retrieve_ix, "K"]
  )

plot_opts <- list(
  "x" = "rsv",
  "y" = "value",
  "col" = "Taxon_5",
  "fill" = "Taxon_5",
  "facet_terms" = c("time", "cluster", "Taxon_5"),
  "facet_scales" = "free_x",
  "facet_space" = "free_x"
)

p <- ggboxplot(beta_hat, plot_opts) +
  scale_y_continuous(limits = c(1 / 421 * 0.98, 1 / 421 * 1.03), expand = c(0, 0)) +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_blank()
  )
ggsave(p, file = sprintf("figure/%s%s-beta.png", cur_ix, timestamp))
write_feather(beta_hat, "../lda/results/dtm_abt_beta.feather")
