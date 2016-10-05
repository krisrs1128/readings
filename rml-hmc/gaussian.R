
###############################################################################
# example one, normal posterior using langevin monte carlo
###############################################################################

source("gaussian_funs.R")

x <- rnorm(30, 0, 10)
res <- mcmc(x, c("mu" = 5, "sigma" = 40), .75, 200)
mean(res$acceptances)
colnames(res$thetas)
plot(res$thetas[, 1])
plot(res$thetas[, 2])
plot(res$thetas, col = "white")
lines(res$thetas)

x <- rnorm(300, 0, 10)
res <- mcmc(x, c("mu" = 5, "sigma" = 40), 0.45, 1000)
mean(res$acceptances)
colnames(res$thetas)
plot(res$thetas[, 1])
plot(res$thetas[, 2])
plot(res$thetas, col = "white")
lines(res$thetas)
