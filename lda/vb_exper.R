source("lda_counts.R")
source("vb_counts.R")

lda_data <- generate_data(200, rep(1000, 200), rep(1, 3), rep(1, 15))
ndv <- dcast(lda_data$Nv, document ~ word, fill = 0) %>%
  select(-document) %>%
  as.matrix()

res <- vb_counts(ndv, rep(1, 3), rep(1, 15), n_iter = 10)
theta_tilde <- res$theta_tilde
for(i in 1:nrow(theta_tilde)) {
  theta_tilde[i, ] <- theta_tilde[i, ] / sum(theta_tilde[i, ])
}

plot(sqrt(lda_data$theta[, 3]), sqrt(theta_tilde[, 2]), asp =1 )
abline(b = 1, a = 0, col = 'red')

beta_tilde <- res$beta_tilde
for (i in 1:nrow(res$beta_tilde)) {
  beta_tilde[i, ] <- beta_tilde[i, ] / sum(beta_tilde[i, ])
}

lda_data$beta
beta_tilde

plot(lda_data$beta[2, ], beta_tilde[2, ], asp = 1)
abline(b = 1, a = 0)

plot(lda_data$beta[3, ], beta_tilde[1, ], asp = 1)
abline(b = 1, a = 0)

plot(lda_data$beta[1, ], beta_tilde[3, ], asp = 1)
abline(b = 1, a = 0)
