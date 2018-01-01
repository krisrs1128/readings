#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Beta regression example
##
## https://cran.r-project.org/web/packages/rstanarm/vignettes/betareg.html
##
## author: sankaran.kris@gmail.com
## date: 12/20/2017

library("rstanarm")
library("tidyverse")

###############################################################################
## Simulating some data
###############################################################################
eta <- c(1, -0.2)
gamma <- c(1.8, 0.4)
N <- 200
x <- rnorm(N, 2, 2)
z <- rnorm(N, 0, 2)
mu <- binomial(link = logit)$linkinv(eta[1] + eta[2]*x)
phi <- binomial(link = log)$linkinv(gamma[1] + gamma[2]*z)
y <- rbeta(N, mu * phi, (1 - mu) * phi)
dat <- data.frame(cbind(y, x, z))

###############################################################################
## Fitting a beta regression model
###############################################################################

## doesn't seem to work...
bfit <- stan_betareg(
  y ~ x | z,
  data = dat,
  link = "logit",
  link.phi = "log"
)

summary(bfit)
