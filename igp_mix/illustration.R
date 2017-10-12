#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Some plots of GPs just for illustration purposes.
##
## author: sankaran.kris@gmail.com
## date: 10/08/2017

library("tidyverse")
library("simData")

bandwidths <- c(rep(5, 10), rep(50, 10), rep(500, 10))
sim <- matrix(nrow = 30, ncol = 50)
for (i in seq_along(bandwidths)) {
  sim[i, ] <- gp_data(1:50, bandwidths[i])
}

msim <- melt(sim, varnames = c("series", "time"))
msim$bandwidth <- as.factor(bandwidths[msim$series])

ggplot(msim) +
  geom_line(
    aes(x = time, y = value, group = series),
    alpha = 0.5,
    size = 0.5
  ) +
  facet_wrap(~bandwidth)
ggsave("figure/gp_bandwidths.png", width = 4, height = 1.5)
