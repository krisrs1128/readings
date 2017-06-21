
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
  ## colnames(z_mat) <- c("iter", paste0("time_", seq_along(state$z)))
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
