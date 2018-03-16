#! /usr/bin/env Rscript
##
## Wavlets study of the antibiotics data
##
## author: krissankaran@stanford.edu
## date: 03/16/2018

library("phyloseq")
library("wavethresh")
library("caret")
source("ZDaub.r")

abt <- get(load("abt.rda")) %>%
  subset_samples(ind == "D")

x <- get_taxa(abt) %>%
  t() %>%
  asinh() %>%
  scale()
times <- sample_data(abt)$time

model_ix <- abt %>%
  taxa_sums() %>%
  which.max()
z <- ZDaub(times, resolution = 2 ^ 18, filterNumber = 7)
df <- data.frame(x = x[, model_ix], z)
fit <- train(x ~ ., data = df, method = "rf")
fit

plot(times, x[, model_ix])
points(times, predict(fit), col = "red")

plot(x[, model_ix], predict(fit))
