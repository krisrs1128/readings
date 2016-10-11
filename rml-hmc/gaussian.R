
###############################################################################
# example one, normal posterior using langevin monte carlo
###############################################################################

source("gaussian_funs.R")

x <- rnorm(300, 0, 10)
res <- langevin_mcmc(x, c("mu" = 5, "sigma" = 40), 200, 0.1)
mean(res$acceptances)
colnames(res$thetas)
plot(res$thetas[, 1])
plot(res$thetas[, 2])
plot(res$thetas, col = "white", asp = 1)
lines(res$thetas)

x <- rnorm(300, 0, 10)
res <- langevin_mcmc(x, c("mu" = 5, "sigma" = 40), 2000, 0.1)
mean(res$acceptances)
colnames(res$thetas)
plot(res$thetas[, 1])
plot(res$thetas[, 2])
plot(res$thetas, col = "white", asp = 1)
lines(res$thetas)


samples_df <- data.frame(
  res$thetas,
  accept = res$acceptances,
  iter = seq_len(nrow(res$thetas))
)

#save.image("~/Desktop/langevin_exper.Rdata")
