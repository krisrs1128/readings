#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulate data from a mixture of gaussian processes, as in
## http://mlg.eng.cam.ac.uk/zoubin/papers/iMGPE.pdf
## but with finite K.
##
## author: sankaran.kris@gmail.com
## date: 07/24/2017

set.seed(724207)
library("tidyverse")
library("kernlab")
library("expm")

## Parameters
K <- 3
time_len <- 50
changepoints <- sample(1:50, K - 1)
z <- vector(length = time_len)

k <- 1
for (i in seq_len(time_len)) {
  z[i] <- k
  if (i %in% changepoints) {
    k <- k + 1
  }
}

kernels <- list()
gp <- matrix(nrow = time_len, ncol = K)
for (k in seq_len(K)) {
  kernels[[k]] <- kernelMatrix(rbfdot(sigma = 20 * rexp(1)), seq_len(time_len) / 50)
  gp[, k] <- Re(sqrtm(kernels[[k]])) %*% rnorm(time_len)
}

y <- vector(length = time_len)
for (i in seq_len(time_len)) {
  y[i] <- gp[i, z[i]]
}
plot(y)

## a different one, continuous
y <- rnorm(time_len, 0, 0.1)
cur_ix <- (0.5 * time_len):time_len
y[cur_ix] <- gp[cur_ix, 1] - gp[cur_ix[1], 1] + rnorm(length(cur_ix), 0, 0.1)
plot(y)

write.csv(y, "data/y.csv", row.names = FALSE, col.names = FALSE)
write.csv(seq_len(time_len) / 50, "data/x.csv", row.names = FALSE, .names = FALSE)
