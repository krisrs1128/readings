#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# Writing believable models of microbiome data can be hard, since
# there tends to be some non-standard structure -- log-tails,
# zero-inflation, and unequal sampling depths, for example. That said,
# we would like to simulate data in order to evaluate the quality of
# different methods (we don't have any ground truth otherwise). A
# compromise between a full parametric simulation and an application
# on real data is a semiparametric simulation, where we resample from
# the real data but introduce known structure. This is the idea
# pursued here, to artificially introduce clusters and gradients in
# samples from the preterm birth study.
#
# author: kriss1@stanford.edu

## ---- setup ----
library("phyloseq")
library("plyr")
library("dplyr")
library("ggplot2")
library("reshape2")
source("./bootstrap_mixture.R")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 8),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

## ---- utils ----
#' Simulate a Nonparametric Mixture with Dirichlet Weights
#'
#' This wraps bootstrap_mixture in the case that the weights are
#' sampled from a Dirichlet(q).
#'
#' @param x_list [list of matrices or vectors] The K data sets to use
#'   for generating nonparametric CDFs, on which we mix.
#' @param B [integer] The number of bootstrap replications to use
#' @param gamma [numeric] The Dirichlet rate parameter used in
#'   simulating the mixture proportions.
simulate_gradient <- function(x_list, B, gamma) {
  K <- length(x_list)
  theta <- matrix(rgamma(B * K, gamma), B, K)
  theta <- (1 / rowSums(theta)) * theta
  list(theta = theta, x_star = bootstrap_mixture(x_list, theta))
}

#' Reshape PCA data for plotting
#'
#' This melts the values of theta along rows, and identifies the
#' mixture component with the maximal contribution to a bootstrap
#' samples.
#'
#' @param theta [B x K matrix] A matrix giving mixing proportions
#'   between the \hat{F}_{k}, for each of the B bootstrap samples.
#' @param scores [matrix] A matrix of pca scores (columns give
#'   components).
#' @return pca_data [data.frame] A data.frame supplementing the raw
#'   scores with information from theta.
merge_pca <- function(theta, scores) {
  K <- ncol(theta)
  data.frame(B = 1:nrow(theta), theta = theta, scores) %>%
    melt(measure.vars = paste0("theta.", 1:K)) %>%
    group_by(B) %>%
    mutate(max_theta = which.max(value))
}

#' Plot PCA on Nonparametric Mixture
#'
#' @param X [data.frame] The reshaped PCA data, as output by merge_pca
#' @param evals [numeric] A vector of eigenvalues from the PCA, used
#'   to adjust the aspect of the resulting figure.
plot_gradient <- function(X, evals) {
  ggplot(X %>% filter(variable == "theta.1")) +
    geom_point(aes(x = Comp.1, y = Comp.2, col = sqrt(value), size = Comp.3),
               alpha = 0.8) +
    scale_size(range = c(.01, 3)) +
    scale_color_gradient(low = "#9EABC8", high = "#D09DB3", limits = c(0, 1)) +
    facet_grid(max_theta ~ .) +
    theme(
      panel.border = element_rect(fill = "transparent", size = 0.3),
      panel.spacing = unit(0, "line")
    ) +
    labs(col = "sqrt(theta(1))") +
    coord_fixed(sqrt(evals[2] / evals[1]))
}

## ---- get-groups ----
# Download microbiome data
pregnancy_path <- "http://statweb.stanford.edu/~susan/papers/Pregnancy/PregnancyClosed15.Rdata"
tmp <- tempfile()
download.file(pregnancy_path, tmp)
load(tmp)

site <- "Vaginal_Swab"
ps <- PSPreg[[site]] %>%
  filter_taxa(function(x) sum(x > 1) > 0.05 * length(x), TRUE)

X <- data.frame(
  outcome = factor(sample_data(ps)$Outcome, levels = c("Term", "Marginal", "Preterm", "VeryPreterm")),
  otu_table(ps)
)

## ---- vis-groups ----
ggplot(melt(X, id.vars = "outcome")) +
  geom_histogram(aes(x = log(1 + value)), bins = 100) +
  facet_grid(outcome ~ ., scale = "free_y")

ggplot(melt(X[, 1:10], id.vars = "outcome")) +
  geom_histogram(aes(x = log(1 + value)), bins = 100) +
  facet_grid(outcome ~ variable, scale = "free_y")

## ---- simulate-gradient ----
x_list <- X %>%
  dlply(.(outcome), function(x) {
    as.matrix(x[, -1])
  })

x_list <- x_list[c(1, 3)]
B <- 1000
K <- length(x_list)
gammas <- c(0.01, 1, 10)

gradients <- list()
for (i in seq_along(gammas)) {
  gradients[[i]] <- simulate_gradient(x_list, B, gammas[i])
}

## ---- plot-gradients ----
for (i in seq_along(gradients)) {
  cur_pca <- princomp(gradients[[i]]$x_star)
  plot_gradient(
    merge_pca(gradients[[i]]$theta, cur_pca$scores),
    cur_pca$sdev
  ) %>%
    print()

  cur_pca <- princomp(scale(asinh(gradients[[i]]$x_star)))
  plot_gradient(
    merge_pca(gradients[[i]]$theta, cur_pca$scores),
    cur_pca$sdev
  ) %>%
    print()
}
