#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Code for the direct assignment collapsed gibbs sampler for the sticky
## HDP-HMM, as described in Algorithm 8 of Emily Fox's thesis
##
## author: Kris Sankaran
## date: 06/19/2017

#' @examples
#' z <- c(1, 1, 2, 1, 1, 3, 3, 3, 1)
#' transition_counts(z)
transition_counts <- function(z) {
  K <- max(z)
  n <- matrix(0, nrow = K, ncol = K)
  time_len <- length(z)

  for (i in seq_len(time_len - 1)) {
    n[z[i], z[i + 1]] <- n[z[i], z[i + 1]] + 1
  }
  n
}
