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

## ---- visualize_beta ----
mu_hat <- apply(samples$mu, c(2, 3, 4), mean)
beta_hat <- apply(samples$beta, c(2, 3, 4), mean)

taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

mbeta_hat <- melt(beta_hat, varnames = c("time_ix", "cluster", "rsv_ix"))
mbeta_hat$rsv <- rownames(otu_table(abt))[mbeta_hat$rsv_ix]
mbeta_hat$time <- times[mbeta_hat$time_ix]
mbeta_hat <- mbeta_hat %>%
  left_join(taxa)

sorted_taxa <- names(sort(table(mbeta_hat$Taxon_5), decreasing = TRUE))
mbeta_hat$Taxon_5 <- factor(
  mbeta_hat$Taxon_5,
  levels = sorted_taxa
)



mbeta_hat$rsv <- factor(mbeta_hat$rsv, levels = taxa$rsv)

levels(mbeta_hat$Taxon_5)

mbeta_hat2 <- mbeta_hat %>%
  right_join(data.frame(Taxon_5 = levels(mbeta_hat$Taxon_5)))

mbeta_hat2$rsv <- factor(mbeta_hat2$rsv, levels = unique(mbeta_hat2$rsv))
mbeta_hat2$Taxon_5 <- factor(mbeta_hat2$Taxon_5, levels(mbeta_hat$Taxon_5))


ggplot(mbeta_hat2 %>%
         filter(Taxon_5 %in% levels(mbeta_hat$Taxon_5)[1:8],
                cluster %in% 1:4)
       ) +
  geom_bar(aes(x = rsv, y = value, fill = Taxon_5), stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(limits = c(0, .03), breaks = c(0, .02), oob = rescale_none) +
  facet_grid(time~cluster) +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.05),
    panel.spacing = unit(0, "line"),
    axis.text.x = element_blank(),
  )
