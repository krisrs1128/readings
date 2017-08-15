#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## An experiment with Adonis inspired by
##
## Dismantling the Mantel tests (Guillot & Rousset)
##
## author: sankaran.kris@gmail.com
## date: 08/01/2017

###############################################################################
## libraries and setup
###############################################################################
library("kernlab")
library("vegan")
library("nlme")
library("tidyverse")
library("reshape2")
library("mvtnorm")
library("viridis")
library("scales")
source("utils.R")
theme_set(
  ggscaffold::min_theme(list(text_size = 7, subtitle_size = 9))
)

###############################################################################
## simulate underlying categorical data
###############################################################################

## variables used throughout experiment
n_sim <- 1000
p1 <- 20
p2 <- 2
n <- 100
u <- matrix(runif(n * p2), n, p2)
K <- kernelMatrix(rbfdot(sigma = 5), u)
probs <- t(logit(rmvnorm(n_sim, sigma = K)))

###############################################################################
## perform poisson and gaussian adonis / lm experiments
###############################################################################

models <- list(
  "poisson" = list(
    "adonis" = list(),
    "logistic" = list(),
    "lm" = list()
  ),
  "gaussian" = list(
    "adonis" = list(),
    "logistic" = list(),
    "lm" = list(),
    "gls" = list()
  )
)

## Poisson setup
for (i in seq_len(n_sim)) {
  print_iter(i)
  f <- t(rmvnorm(p1, sigma = K))
  x <- matrix(rpois(n * p1, lambda = exp(f)), n, p1)
  y <- sample_probs(probs[, i])
  models$poisson$adonis[[i]] <- adonis(x ~ y + u, method = "bray", perm = 20)
  models$poisson$logistic[[i]] <- glm(y ~ x + u, family = binomial())
  models$poisson$lm[[i]] <- glm(x[, 1] ~ y + u, family = "poisson")
  models$poisson$gls[[i]] <- gls(x[, 1] ~ y + u)
}

## Gaussian setup
for (i in seq_len(n_sim)) {
  print_iter(i)
  y <- sample_probs(probs[, i])
  x <- t(rmvnorm(p1, sigma = K))
  models$gaussian$adonis[[i]] <- adonis(x ~ y + u, method = "euclidean", perm = 19)
  models$gaussian$logistic[[i]] <- glm(y ~ x + u, family = binomial())
  models$gaussian$lm[[i]] <- lm(x[, 1] ~ y + u)
  models$gaussian$gls[[i]] <- gls(x[, 1] ~ y + u)
}

###############################################################################
## visualize the resulting p-values
###############################################################################
pvals <- list(
  "poisson" = list(
    "adonis" = sapply(models$poisson$adonis, function(x) x$aov.tab["y", "Pr(>F)"]),
    "logistic" = sapply(models$poisson$logistic, function(x) coef(summary(x))[2:(p1 + 1), "Pr(>|z|)"]),
    "lm" = sapply(models$poisson$lm, function(x) coef(summary(x))["y", "Pr(>|z|)"]),
    "gls" = sapply(models$poisson$gls, function(x) coef(summary(x))["y", "p-value"])
  ),
  "gaussian" = list(
    "adonis" = sapply(models$gaussian$adonis, function(x) x$aov.tab["y", "Pr(>F)"]),
    "logistic" = sapply(models$gaussian$logistic, function(x) coef(summary(x))[2:(p1 + 1), "Pr(>|z|)"]),
    "lm" = sapply(models$gaussian$lm, function(x) coef(summary(x))["y", "Pr(>|t|)"]),
    "gls" = sapply(models$gaussian$gls, function(x) coef(summary(x))["y", "p-value"])
  )
)

m_pvals <- melt(pvals)
colnames(m_pvals) <- c("pval", "rsv", "sim", "method", "mechanism")

m_pvals$method <- factor(m_pvals$method, levels = c("adonis", "logistic", "lm", "gls"))
p <- ggplot(m_pvals) +
  geom_histogram(aes(x = pval), binwidth = 0.02) +
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2)) +
  facet_grid(method ~ mechanism, scales = "free_y")
ggsave("../doc/gp_exper/figure/pvals_comparison.png", width = 4, height = 2)

## save to file
dir.create("data")
save(pvals, models, file = "data/exper.rda")

###############################################################################
## illustrate some of the simulation mechanism
###############################################################################
f <- t(rmvnorm(p1, sigma = K))
x <- matrix(rpois(n * p1, lambda = exp(f)), n, p1)

p <- plot_species_counts(x, u, y)
ggsave("../doc/gp_exper/figure/x_poisson.png", p, width = 4.5, height = 1.7)

x <- t(rmvnorm(p1, sigma = K))
p <- plot_species_counts(x, u, y)
ggsave("../doc/gp_exper/figure/x_gaussian.png", p, width = 4, height = 1.7)

plot_df <- data.frame("y" = y, "p" = probs, "u" = u)
p <- ggplot(plot_df) +
  geom_point(
    aes(x = u.1, y = p, size = u.2),
    alpha = 0.4
  ) +
  geom_point(
    aes(x = u.1, y = y, size = u.2, col = as.factor(y)),
    shape = 3, alpha = 0.9
  ) +
  scale_color_brewer("y", palette = "Set2") +
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 3)) +
  scale_size_continuous(range = c(0.1, 3)) +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7)
  )
ggsave("../doc/gp_exper/figure/p_gp.png", p, width = 4, height = 2)

###############################################################################
## How might we be able to salvage the linear model approach?
###############################################################################

## the model is not valid as is because the errors are correlated
plot(u[, 1], resid(models$gaussian$lm[[1]]))
plot(u[, 2], resid(models$gaussian$lm[[1]]))


###############################################################################
## Can we still try to leverage interesting distances?
###############################################################################
