#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This gives a simple R interface to exponential family PCA. Before
# running this function, make sure that matlab is avialable on the
# path specified (e.g., run ml matlab; which matlab).


#' Wrapper for Exponential Family PCA
#'
#' This calls the matlab function expca_wrapper.m from R, which is
#' itself a wrapper of https://github.com/lydiatliu/epca/. Not all
#' parameter options are supported here.
#'
#' @param X [data.frame or matrix] The matrix on which to perform epca
#' @param exp_fam_type [string] The type of epca to run.
#' @param r [integer] The number of components to return
#' @param path_to_epca [string] The full path to the exponential pca
#'   repo.
#' @param path_to_matlab [string] The full path to the matlab
#'   executable.
#' @return A list with the following components
#'    $U: The ePCA scores
#'    $V: The ePCA directions
#'    $evals: The ePCA eigenvalues
epca <- function(X, exp_fam_type = "poisson", r = 3,
                 path_to_epca = "~/epca/software/",
                 path_to_matlab = "/share/sw/licensed/MATLAB-R2016b/bin/matlab") {
  options(matlab.path = path_to_matlab)
  tmp_input <- tempfile()
  tmp_output <- tempdir()

  # write data to run pca on
  write.table(
    X,
    file = tmp_input,
    sep = ",",
    col.names = FALSE,
    row.names = FALSE
  )

  # run through matlab
  script <- sprintf(
    "epca_wrapper('%s', '%s', '%s', '%s', %d, %f, %f)",
    tmp_input,
    tmp_output,
    path_to_epca,
    exp_fam_type,
    r, 1, 1
  )
  matlabr::run_matlab_code(script)

  # return
  list(
    "U" = read.csv(file.path(tmp_output, "U.txt"), header = FALSE),
    "V" = read.csv(file.path(tmp_output, "V.txt"), header = FALSE),
    "evals" = read.csv(file.path(tmp_output, "evals.txt"), header = FALSE)
  )
}
