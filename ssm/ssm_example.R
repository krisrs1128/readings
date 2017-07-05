#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## An example applying the state space switching model to toy data.
##
## author: sankaran.kris@gmail.com
## date: 07/05/2017

## ---- libraries ----
library("expm")
library("Rcpp")
library("tidyverse")
library("forcats")
library("MCMCpack")
sourceCpp("utils.cpp")
source("ssm.R")
theme_set(ggscaffold::min_theme(list(
                        "legend_position" = "right",
                        "border_size" = 0.2
                      ))
          )

###############################################################################
## simulate data
###############################################################################
set.seed(0701)
As <- list(diag(0.99, nrow = 1), diag(0.05, nrow = 1))
Cs <- list(diag(1, nrow = 1), diag(1, nrow = 1))
s <- c(rep(1, 50,), rep(2, 50))
Qs <- list(diag(0.5, nrow = 1), diag(0.1, nrow = 1))
Rs <- list(diag(0.1, nrow = 1), diag(4, nrow = 1))
res <- simulate(As, Cs, s, 1, Qs, Rs)
y <- res$y

y_df <- data_frame(
  "i" = seq_len(nrow(y)),
  "y" = y[, 1],
  "state" = as_factor(as.character(s)),
  "x1" = res$x[, 1]
)
p <- ggplot(y_df) +
  geom_point(aes(x = i, y = y)) +
  geom_line(aes(x = i, y = x), alpha = 0.5) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))
ggsave("~/Desktop/lab_meetings/20170705/figure/ssm_truth.pdf", p)

###############################################################################
## fit an SSM
###############################################################################
n_iter <- 50
tau <- c(seq(3, 1, length.out = n_iter))
test <- ssm_em(y, 2, 1, n_iter = n_iter, tau)
test$lds_param[[1]]$R
test$lds_param[[2]]$R
y_df$ssm1 <- as.numeric(test$lds_infer[[1]]$x_smooth * test$lds_param[[1]]$C[1, 1])
y_df$ssm2 <- as.numeric(test$lds_infer[[2]]$x_smooth * test$lds_param[[2]]$C[1, 1])
y_df$p <- exp(test$log_ht[, 2])

p <- ggplot(y_df) +
  geom_point(aes(x = i, y = y)) +
  geom_line(aes(x = i, y = x), alpha = 0.5) +
  geom_line(aes(x = i, y = ssm2, col = p), size = 1) +
  geom_line(aes(x = i, y = ssm1, col = p), size = 1) +
  scale_color_gradient2(low ="#008071", midpoint = 0.5, high = "#710080") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))
ggsave("~/Desktop/lab_meetings/20170705/figure/ssm_fit.pdf", p)
