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

preprocess_c <- function(c_samples) {
  colnames(c_samples) <- c("iter", "sample", "class")
  c_samples$truth <- as.factor(data$class[c_samples$sample])
  c_samples$class <- factor(c_samples$class, levels = names(sort(table(c_samples$class), decreasing = TRUE)))
  c_samples$class_group <- fct_lump(c_samples$class, n = 7)
  c_samples$sample <- factor(c_samples$sample, levels = order(data$class))
  c_samples
}

coocurrence_counts <- function(x) {
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


###############################################################################
## Data generated from true mixture of GPs model
###############################################################################
post <- read_csv("data/mix_post.csv", col_names = FALSE)
colnames(post) <- c("class", "x", "y")
data <- read_csv("data/mix_data.csv", col_names = FALSE)
colnames(data) <- c("class", "x", "y")

## plot the raw data and the final estimated fits
ggplot() +
  geom_line(
    data = post,
    aes(x = x, y = y, group = class),
  ) +
  geom_point(
    data = data,
    aes(x = x, y = y, col = as.factor(class))
  )

## read in the cluster membership data
c_samples <- read_csv("data/samples/sim0914/c.csv", col_names = FALSE) %>%
  preprocess_c()

## plot the class memberships over time
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

## plot co-occurrence of classes between samples
counts <- coocurrence_counts(c_samples)
mcounts <- melt(
  counts,
  varnames = c("i", "j"),
  value.name = "count"
)

ggplot(mcounts) +
  geom_tile(
    aes(x = i, y = j, fill = count)
  ) +
  scale_fill_gradient(low = "white", high = "black")

## see how often samples were placed in the same cluster
bump <- read_csv("data/bump_data.csv", col_names = FALSE)
post <- read_csv("data/bump_posteriors.csv", col_names = FALSE)
ggplot() +
  geom_point(
    data = bump,
    aes(x = X1, y = X2)
  ) +
  geom_line(
    data = post,
    aes(x = X2, y = X3, group = X1)
  )
