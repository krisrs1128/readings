#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## HSAUR rstanarm vignette, from
##
## https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
##
## author: sankaran.kris@gmail.com
## date: 12/20/2017

###############################################################################
## Libraries and setup
###############################################################################
library("rstanarm")

###############################################################################
## Data and basic (frequentistic) glm
###############################################################################
data("womensrole", package = "HSAUR3")
womensrole$total <- womensrole$agree + womensrole$disagree
womensrole_glm_1 <- glm(cbind(agree, disagree) ~ education + gender,
                        data = womensrole, family = binomial(link = "logit"))
round(coef(summary(womensrole_glm_1)), 3)

###############################################################################
## Bayesian GLM isntead
###############################################################################
womensrole_bglm_1 <- stan_glm(
  cbind(agree, disagree) ~ education + gender,
  data = womensrole,
  family = binomial(link = "logit"), 
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7)
)

womensrole_bglm_1
