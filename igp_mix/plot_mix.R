#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Plot results of mix gp on simulated data
##
## author: sankaran.kris@gmail.com
## date: 09/13/2017

library("tidyverse")
library("forcats")
library("reshape2")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

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

#' Preprocessing class memberships
#'
#' Useful column names / types for gibbs samples for c (the length n class
#' membership variable).
preprocess_c <- function(c_samples, data) {
  colnames(c_samples) <- c("iter", "sample", "class")
  c_samples$truth <- as.factor(data$class[c_samples$sample])
  c_samples$class <- factor(c_samples$class, levels = names(sort(table(c_samples$class), decreasing = TRUE)))
  c_samples$class_group <- fct_lump(c_samples$class, n = 7)
  c_samples$sample <- factor(c_samples$sample, levels = order(data$class))
  c_samples
}

#' Count Cooccurrences
#'
#' We can study how many times pairs of samples get placed in the same class
#' across gibbs sampling iterations.
cooccurrence_counts <- function(x) {
  samples <- unique(x$sample)
  n <- length(samples)
  counts <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(samples, samples)
  )

  for (i in unique(x$iter)) {
    message("Processing ", i)
    cur_data <- x %>%
      filter(iter == i) %>%
      select(sample, class)
    cur_data <- setNames(cur_data$class, cur_data$sample)

    for (j in samples) {
      for (k in samples) {
        if (cur_data[j] == cur_data[k]) {
          counts[j, k] <- 1 + counts[j, k]
        }
      }
    }

  }
  counts
}

#' Plot Mix-GP Fits
#'
#' Plot the raw data and the estimated fits using the final sampled
#' hyperparameters. A better strategy of course would look at the posterior
#' means of the hyperparameters after accounting for label switching (see the
#' plot_c function for that though).
plot_fits <- function(data, post) {
  ggplot() +
    geom_line(
      data = post,
      aes(x = x, y = y, group = class),
      ) +
    geom_point(
      data = data,
      aes(x = x, y = y, col = as.factor(class))
    )
}

#' Plot Classes
#'
#' Plot the class memberships for individual samples over many iterations, after
#' arranging the samples according to their true class memberships.
plot_c <- function(c_samples) {
  ggplot(c_samples) +
    geom_tile(
      aes(x = sample, y = iter, fill = class_group)
    ) +
    geom_text(
      aes(x = sample, y = iter, label = class),
      size = 0.75,
      alpha = 0.5
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    facet_grid(. ~ truth, scale = "free", space = "free") +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      panel.border = element_rect(size = 0.5, fill = "transparent"),
      panel.spacing = unit(0, "cm")
    )
}

#' Plot Cooccurrences
#'
#' Plot the count of the number of times samples are placed into the same
#' classes, after arranging samples according to their true class memberships.
plot_cooccurrence <- function(counts, data) {
  classes <- min(data$class):max(data$class)
  diag(counts) <- 0
  mcounts <- counts %>%
    melt(
      varnames = c("i", "j"),
      value.name = "count"
    ) %>%
    mutate(
      class_group1 = factor(data$class[i], levels = rev(classes)),
      class_group2 = factor(data$class[j], levels = classes),
      i = factor(i, levels = order(data$class)),
      j = factor(j, levels = order(data$class))
    )

  ggplot(mcounts) +
    geom_tile(
      aes(
        x = i,
        y = j,
        fill = count)
    ) +
    scale_fill_gradient(low = "white", high = "black") +
    facet_grid(class_group2 ~ class_group1, scales = "free", space = "free")
}

###############################################################################
## Data generated from true mixture of GPs model
###############################################################################
post <- read_csv("data/mix_post.csv", col_names = FALSE)
colnames(post) <- c("class", "x", "y")
data <- read_csv("data/mix_data.csv", col_names = FALSE)
colnames(data) <- c("class", "x", "y")
plot_fits(data, post)
dir.create("figure")
ggsave("figure/gp_fits.png")

## read in the cluster membership data
c_samples <- read_csv("data/samples/sim0914/c.csv", col_names = FALSE) %>%
  preprocess_c(data)
plot_c(c_samples)
ggsave("figure/gp_c_samples.png")

## plot co-occurrence of classes between samples
counts <- cooccurrence_counts(c_samples)
plot_cooccurrence(counts, data)
ggsave("figure/gp_cooccurrence.png")

###############################################################################
## Data generated from toy bump function
###############################################################################
post <- read_csv("data/bump_posteriors.csv", col_names = FALSE)
colnames(post) <- c("class", "x", "y")
bump <- read_csv("data/bump_data.csv", col_names = FALSE)
colnames(bump) <- c("class", "x", "y")
plot_fits(bump, post)
ggsave("figure/bump_fits.png")

c_samples <- read_csv("data/samples/bump0914/c.csv", col_names = FALSE) %>%
  preprocess_c(bump)
plot_c(c_samples)
ggsave("figure/bump_c_samples.png")

counts <- cooccurrence_counts(c_samples)
plot_cooccurrence(counts)
ggsave("figure/bump_cooccurrence.png")
