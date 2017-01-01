#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## A simple bootstrap experiment in variable transformation, as described in "An
## Introduction to the boostrap". The idea is to automatically recover the usual
## VST for the correlation in a bivariate normal model automatically using the
## bootstrap.

# Setup packages ---------------------------------------------------------------
library("ggplot2")

## Load packages into session
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  text = element_text(family = "Ubuntu Regular", color = "#22211d"),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

## ---- utils ----
generate_correlated <- function(N, rho) {
  z <- rnorm(N)
  cbind(X1 = z, X2 = rho * z + rnorm(N, 0, sqrt(1 - rho ^ 2)))
}

## ---- simulate ----
N <- 10
rho <- .4
R <- 15000
B <- c(500, 500)

cors <- vector(length = R)
for (i in seq_len(R)) {
  X <- generate_correlated(N, rho)
  cors[i] <- cor(X)[1, 2]
}

## ---- usual-transformation ----
ggplot(data.frame(cors)) +
  geom_histogram(aes(x = cors), binwidth = 0.01)
ggplot(data.frame(cors)) +
  geom_histogram(aes(x = atanh(cors)), binwidth = 0.01)

## ---- bootstrap-simulation----

## Use the current version of X to estimate this transformation
vst_ests <- list(
  fits = vector(length = B[1]),
  variances = vector(length = B[1])
)

for (i in seq_len(B[1])) {
  if (i %% 20 == 0) {
    cat(sprintf("outer bootstrap iteration %d\n", i))
  }

  cur_x <- X[sample(N, replace = TRUE), ]
  vst_ests$fits[i] <- cor(cur_x)[1, 2]
  cur_rhos <- vector(length = B[2])
  for (j in seq_len(B[2])) {
    cur_rhos[j] <- cor(cur_x[sample(N, replace = TRUE), ])[1, 2]
  }
  vst_ests$variances[i] <- var(na.omit(cur_rhos))
}

## ---- bootstrap-transformation ----
s_smooth <- ksmooth(
  vst_ests$fits,
  sqrt(vst_ests$variances),
  kernel = "box",
  bandwidth = 0.05
)

ggplot() +
  geom_point(
    data = data.frame(vst_ests),
    aes(x = fits, y = sqrt(variances))
  ) +
  geom_line(
    data = data.frame(s_smooth),
    aes(x = x, y = y)
  )
 
