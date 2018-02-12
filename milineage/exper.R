#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Experiment using the milineage package
##
## author: sankaran.kris@gmail.com
## date: 02/12/2018

library("miLineage")
example(miLineage)
browseVignettes("miLineage")

## x4freq is not an intuitive name. Why note covariates?
## why not default the x4 matrices to null as well?
## Also, why not give some simple default to case? Like, all default to all the
## variables being of interest?
example(ZIGDM)
test <- ZIGDM(OTU.toy.reorder, NULL, NULL, case, test.type = "Disp", 1)

## would be nice if this function returned more than a p-value
## Like, what if we want to study the fit of the model, convergence of the EM,
## etc.?
test
