
###############################################################################
# example one, normal posterior using langevin monte carlo
###############################################################################

source("gaussian_funs.R")
library("ggplot2")

x <- rnorm(300, 0, 10)
langevin_res <- langevin_mcmc(x, c("mu" = 5, "sigma" = 100), 800, 0.1)
rmc_res <- rmc_mcmc(x, c("mu" = 5, "sigma" = 100), 800, 0.1)
samples <- data.frame(
  method = c(rep("langevin", 800), rep("riemann", 800)),
  rbind(langevin_res$thetas, rmc_res$thetas),
  accept = as.factor(c(langevin_res$acceptances, rmc_res$acceptances)),
  iter = c(1:800, 1:800)
)

## ---- rml_vis ----
ggplot(samples) +
  geom_point(aes(x = sigma, y = mu, col = accept), size = 0.5) +
  geom_path(aes(x = sigma, y = mu), alpha = 0.1) +
  facet_grid(method~.) +
  scale_color_manual(values = c("#F09465", "#5EBD90")) +
  coord_fixed()
