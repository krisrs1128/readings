
###############################################################################
# kalman filtering
###############################################################################

#' @title Kalman filtering
#' @param A Evolution matrix across timepoints
#' @param C Observation matrix across timepoints
#' @examples
#' p <- 10
#' n_times <- 1000
#' y <- topic_params(p, n_times)
#' mu0 <- rnorm(p)
#' sigma0 <- diag(p)
#' A <- array(dim = c(p, p, n_times))
#' C <- array(dim = c(p, p, n_times))
#' Q <- array(dim = c(p, p, n_times))
#' R <- array(dim = c(p, p, n_times))
#' for (i in seq_len(n_times)) {
#'   A[,,i] <- diag(p)
#'   R[,,i] <- diag(p)
#'   Q[,,i] <- diag(p)
#'   C[,,i] <- diag(p)
#' }
#'
#' y_tilde <- kalman_filter(y, mu0, sigma0, A, C, Q, R)
#' plot(y[, 1], type = "l")
#' lines(y_tilde$filtered$mu[, 1], col = "red")
kalman_filter <- function(y, mu0, sigma0, A, C, Q, R) {
  n_times <- nrow(y)
  p <- ncol(y)

  filtered_values <- list(
    "mu" = matrix(nrow = n_times, ncol = p),
    "sigma" = array(dim = c(p, p, n_times))
  )

  pred_values <- list(
    "mu_pred" = matrix(nrow = n_times, ncol = p),
    "sigma_pred" = array(dim = c(p, p, n_times))
  )

  filtered_values$mu[1, ] <- mu0
  filtered_values$sigma[,, 1] <- sigma0
  for (i in seq_len(n_times - 1)) {
    pred_step <- measurement_update(
      filtered_values$mu[i, ],
      filtered_values$sigma[,, i],
      A[,, i + 1],
      Q[,, i + 1]
    )

    cur_updates <- prediction_update(
      y[i, ],
      C[,, i + 1],
      pred_step$mu_pred,
      pred_step$sigma_pred,
      R[,, i]
    )

    filtered_values$mu[i + 1, ] <- cur_updates$mu_cur
    filtered_values$sigma[,, i + 1] <- cur_updates$sigma_cur
    pred_values$mu[i + 1, ] <- pred_step$mu_pred
    pred_values$sigma[,, i + 1] <- pred_step$sigma_pred

  }

  list(
    "filtered" = filtered_values,
    "predicted" = pred_values
  )
}

measurement_update <- function(mu_prev, sigma_prev, A_cur, Q_cur) {
  list(
    "mu_pred" = A_cur %*% mu_prev,
    "sigma_pred" = A_cur %*% sigma_prev %*% t(A_cur) + Q_cur
  )
}

prediction_update <- function(y_cur, C_cur, mu_pred, sigma_pred, R_cur) {
  y_hat <- as.numeric(C_cur %*% mu_pred)
  K_cur <- solve(solve(sigma_pred) + t(C_cur) %*% R_cur %*% C_cur) %*% t(C_cur) %*% solve(R_cur)

  list(
    "mu_cur" = mu_pred + K_cur %*% (y_cur - y_hat),
    "sigma_cur" = (diag(nrow(sigma_pred)) - K_cur %*% C_cur) %*% sigma_pred
  )
}

smooth_filtered_estimates <- function(mu_filtered, sigma_filtered, mu_pred, sigma_pred, A) {
  n_times <- nrow(mu_filtered)
  p <- ncol(mu_filtered)

  mu_smoothed <- matrix(nrow = n_times, ncol = p)
  sigma_smoothed <- array(dim = dim(sigma_filtered))

  mu_smoothed[n_times, ] <- mu_filtered[n_times, ]
  sigma_smoothed[,, n_times] <- sigma_filtered[,, n_times]

  for (i in (n_times - 1):1) {
    J <- sigma_filtered[,, i] %*% t(A[,, i + 1]) %*% solve(sigma_pred[,, i + 1])

    mu_smoothed[i, ] <- mu_filtered[i, ] +
      J %*% (mu_smoothed[i + 1, ] - mu_pred[i + 1, ])
    sigma_smoothed[,, i] <- sigma_filtered[,, i] +
      J %*% (sigma_smoothed[,, i + 1] - sigma_pred[,, i + 1]) %*% t(J)
  }

  list(
    "mu_smoothed" = mu_smoothed,
    "sigma_smoothed" = sigma_smoothed
  )
}
