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
## ---- get_data ----
data(abt)
n_samples <- ncol(otu_table(abt))
abt <- abt %>%
  filter_taxa(function(x) sum(x != 0) > .4 * n_samples, prune = TRUE) %>%
  subset_samples(ind == "D")

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

param_grid <- expand.grid(
  sigma = c(.001, .005, .01, .05, .1, .5),
  delta = c(.001, .005, .01, .05, .1, .5),
  K = c(2, 3, 4)
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
samples <- rstan::extract(stan_fit)

## ---- visualize_theta ----
softmax <- function(mu) {
  exp(mu) / sum(exp(mu))
}

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
levels(theta_hat$cluster) <- rev(levels(theta_hat$cluster))

feather::write_feather(theta_hat, "dtm_theta_hat.feather")
ggplot(theta_hat) +
  geom_tile(aes(x = reorder(sample, time), y = cluster, fill = value)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#5BBABA", limits = c(0, 1))

p1 <- ggplot(theta_hat) +
  geom_tile(aes(x = time, y = cluster, fill = value)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#5BBABA", limits = c(0, 1)) +
  guides(fill = guide_legend(keywidth = .5, keyheight = .5)) +
  labs(x = "time") +
  coord_fixed(4)

ggplot(theta_hat) +
  geom_line(aes(x = time, y = value, col = cluster))

## ---- visualize_beta ----
mu_hat <- apply(samples$mu, c(2, 3, 4), mean)
beta_hat <- apply(samples$beta, c(2, 3, 4), mean)

taxa <- data.table(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]
sorted_taxa <- names(sort(table(taxa$Taxon_5), decreasing = TRUE))

mbeta_hat <- melt(beta_hat, varnames = c("time_ix", "cluster", "rsv_ix"))
mbeta_hat$rsv <- rownames(otu_table(abt))[mbeta_hat$rsv_ix]
mbeta_hat$time <- times[mbeta_hat$time_ix]
mbeta_hat <- mbeta_hat %>%
  left_join(taxa) %>%
  right_join(data.frame(Taxon_5 = factor(sorted_taxa)))
mbeta_hat$rsv <- factor(mbeta_hat$rsv, levels = unique(mbeta_hat$rsv))
mbeta_hat$Taxon_5 <- factor(mbeta_hat$Taxon_5, sorted_taxa)

p2 <- ggplot(mbeta_hat %>%
         filter(Taxon_5 %in% levels(mbeta_hat$Taxon_5)[1:8])
       ) +
  geom_bar(aes(x = rsv, y = value, fill = Taxon_5), stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(breaks = c(0, .02), limits = c(0, .05), oob = scales::rescale_none) +
  facet_grid(time ~ cluster) +
  guides(fill = guide_legend(keywidth = .5, keyheight = .5)) +
  theme(
    panel.border = element_rect(fill = "transparent", size = 0.3),
    panel.spacing = unit(0, "line"),
    axis.text.x = element_blank(),
    legend.text = element_text(size = 5),
    legend.title = element_blank()
  )

feather::write_feather(mbeta_hat, "dtm_mbeta_hat.feather")

## ---- multiplot ----
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
png(paste0("figure/", cur_ix, timestamp, "figure.png"))
multiplot(p1, p2, layout = matrix(c(1, 2, 2, 2), ncol = 1))
dev.off()
