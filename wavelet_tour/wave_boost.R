## Example from Wavelet-Based Gradient Boosting [Dubossarsky et. al. 2016]
##
##
## 03/14/2018
## sankaran.kris@gmail.com

library("wavethresh")
source("ZDaub.r")
x1 <- runif(1000)
x2 <- runif(1000)
x3 <- runif(1000)
y <- (2.2*(4*sin(4*pi*x1)-sign(x1-0.3)-sign(0.72-x1)) +10*sign(x2-(1/3))-7*sign(x2-(2/3))+rnorm(1000))
Z1 <- ZDaub(x1,range.x=c(0,1),numLevels=7)
Z2 <- ZDaub(x2,range.x=c(0,1),numLevels=7)
Z3 <- ZDaub(x3,range.x=c(0,1),numLevels=7)
Z <- cbind(Z1,Z2,Z3)

plot(x1, y)
plot(x1, Z1[, 1])
plot(x1, Z1[, 2])
