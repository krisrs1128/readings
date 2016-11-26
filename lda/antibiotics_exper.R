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

# get samples
samples <- rstan::extract(stan_fit)

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

# might want to set prior for more extreme decay
ggplot(beta_hat) +
  geom_bar(aes(x = reorder(rsv, -value, mean), y = value, fill = as.factor(cluster)),
           stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = -90, size = 3))

ggplot(beta_hat %>%
         filter(Taxon_5 %in% sorted_taxa[1:5])) +
  geom_bar(aes(x = reorder(rsv, -value, min), y = value),
           stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  facet_grid(cluster~Taxon_5, scale = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = -90, size = 3))

# study sample cluster memberships
samples_theta <- melt(samples$theta)
theta_hat <- samples_theta %>%
  group_by(Var2, Var3) %>%
  summarise(mean = mean(value))
theta_hat$Var2 <- rownames(X)[theta_hat$Var2]
colnames(theta_hat) <- c("sample", "cluster", "theta")

sample_info <- sample_data(abt)
sample_info$sample <- rownames(sample_info)

theta_hat <- theta_hat %>%
  left_join(sample_info, by = "sample")

ggplot(theta_hat) +
  geom_line(aes(x = time, y = theta, col = as.factor(cluster))) +
  facet_wrap(~cluster) +
  scale_color_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = -90))

dir.create("results")
write.csv(samples_theta, file = "results/theta.csv", row.names = FALSE)
write.csv(samples_beta, file = "results/beta.csv", row.names = FALSE)
