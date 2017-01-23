#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Experiments with NMF, under different parameterizations and fitting
## techniques. Here, we consider the parameterizations in the parameterizations
## table and then save the resulting fits in a fits/ directory. These are
## visualized separately.
##
## author: kriss1@stanford.edu

## ---- libraries ----
library("jsonlite")
library("SLURMHelpers")
base_dir <- "/scratch/users/kriss1/programming/readings/nmf"
source(file.path(base_dir, "src", "nmf_utils.R"))

## ---- configuration ----
## create the configuration JSON file
output_dir <- file.path(base_dir, "fits")
batch_dir <- file.path(base_dir, "batch")
dir.create(batch_dir)
config_path <- file.path(batch_dir, "config.json")


sim_factors <- list(
  "N" = c(100),
  "P" = c(75, 125),
  "zero_inf_prob" = c(0, 0.2)
)
model_factors <- list(
  "inference" = c("gibbs", "vb"),
  "method" = c(
    file.path(base_dir, "src", "nmf_gamma_poisson.stan"),
    file.path(base_dir, "src", "nmf_gamma_poisson_zero.stan")
  )
)

write_configs(
  sim_factors,
  model_factors,
  n_batches = 3,
  config_path = config_path
)

## ---- submit-jobs ----
## loop over unique values in the "batch" field of the json file
configs <- fromJSON(config_path, simplifyVector = FALSE)
batches <- sapply(configs, function(x) { x$batch })
batch_opts <- list("mem_alloc" = 4000)

for (i in seq_along(unique(batches))) {
    batch_script <- file.path(batch_dir, paste0("batch-", i, ".sbatch"))
    rscript_file <- file.path(base_dir, "src", "nmf_script.R")
    rscript_cmd <- paste("Rscript", rscript_file, config_path, i)

    create_job(
      batch_script,
      paste0("nmf-", i),
      rscript_cmd,
      batch_opts
    )
   system(paste("sbatch", batch_script))
}
