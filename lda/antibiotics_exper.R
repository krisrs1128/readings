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
X <- log(1 + t(otu_table(abt)@.Data))
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
samples_beta <- melt(samples$beta)

beta_hat <- samples_beta %>%
  group_by(Var2, Var3) %>%
  summarise(mean = mean(value)) %>%
  dcast(Var2 ~ Var3)

taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

colnames(beta_hat) <- c("cluster", colnames(X))
beta_hat <- beta_hat %>%
  melt(id.vars = "cluster", variable.name = "rsv") %>%
  left_join(taxa, by = "rsv")

sorted_taxa <- names(sort(table(beta_hat$Taxon_5), decreasing = TRUE))
beta_hat$Taxon_5 <- factor(
  beta_hat$Taxon_5,
  levels = sorted_taxa
)
beta_hat$cluster <- paste("Cluster", beta_hat$cluster)

## ---- visualize_beta ----
# might want to set prior for more extreme decay
ggplot(beta_hat) +
  geom_bar(aes(x = reorder(rsv, -value, mean), y = value, fill = as.factor(cluster)),
           stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = -90, size = 3))

ggplot(beta_hat %>%
         filter(Taxon_5 %in% levels(beta_hat$Taxon_5)[1:8])) +
  geom_bar(aes(x = reorder(rsv, -value, min), y = value, fill = Taxon_5),
           stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  facet_grid(cluster~Taxon_5, scale = "free_x", space = "free_x") +
  theme(
    axis.text.x = element_text(angle = -90, size = 3),
    strip.text.x = element_blank()
  )

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
ggplot(theta_hat) +
  geom_point(aes(x = time, y = theta, col = condition), size = .4) +
  geom_line(aes(x = time, y = theta, col = condition), size = 0.5) +
  facet_wrap(~cluster) +
  scale_color_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = -90))

ggplot(theta_hat) +
  geom_tile(aes(x = time, y = cluster, fill = theta)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#5BBABA", limits = c(0, 1)) +
  facet_grid(~condition, scale = "free_x", space = "free_x") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    panel.border = element_rect(fill = "transparent", size = .4, color = "#C9C9C9"),
    panel.spacing = unit(0, "line")
  )

## ---- save_results ----
dir.create("results")
write_feather(theta_hat, path = "results/theta.feather")
write_feather(beta_hat, path = "results/beta.feather")