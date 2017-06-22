#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Non sampler-specific utilities that we use in the HDP-HMM samplers.
##
## author: sankaran.kris@gmail.com
## date: 6/21/2017

append_to_file <- function(name, x) {
  write.table(
    x, name,
    append = TRUE,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(name)
  )
}

write_state <- function(out_dir, state, iter) {
  dir.create(out_dir, recursive = TRUE)
  z_mat <- cbind(iter, t(as.matrix(state$z)))
  colnames(z_mat) <- c("iter", paste0("time_", seq_along(state$z)))
  append_to_file(file.path(out_dir, "z.csv"), z_mat)

  beta_mat <- cbind(iter, t(as.matrix(state$beta)))
  colnames(z_mat) <- c("iter", paste0("k_", seq_along(state$beta)))
  append_to_file(file.path(out_dir, "beta.csv"), beta_mat)

  w_mat <- cbind(iter, t(as.matrix(state$w)))
  colnames(w_mat) <- c("iter", paste0("k_", seq_along(state$w)))
  append_to_file(file.path(out_dir, "w.csv"), w_mat)

  em <- state$emission
  em$iter <- iter
  cat(toJSON(em), file = file.path(out_dir, "emissions.json"), append = TRUE)

  m_mat <- cbind(iter, t(as.matrix(state$m)))
  colnames(m_mat) <- c("iter", paste0("k_", seq_len(ncol(m_mat))))
  append_to_file(file.path(out_dir, "m.csv"), m_mat)
}

## ---- block-sampler-utils ----
merge_default_hyper <- function(opts = list()) {
  default_opts <- list(
    "L" = 20,
    "n_iter" = 1000,
    "theta_iter" = 2,
    "kappa" = 1,
    "alpha" = 1
  )
  modifyList(default_opts, opts)
}

merge_default_lambda <- function(opts = list()) {
  default_opts <- list(
    "mu0" = c(0, 0),
    "sigma0" = diag(2),
    "nu" = 3,
    "delta" = matrix(c(1, 0, 0, 1), nrow = 2)
  )
  modifyList(default_opts, opts)
}

multi_dmvnorm <- function(yt, theta) {
  modes <- names(theta)
  y_dens <- setNames(seq_along(modes), modes)
  for (l in modes) {
    y_dens[l] <- dmvnorm(yt, theta[[l]]$mu, theta[[l]]$sigma, log = TRUE)
  }
  y_dens
}
