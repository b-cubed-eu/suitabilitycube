#' Estimate a global bandwidth for hypervolume KDE
#'
#' Uses hypervolume::estimate_bandwidth on the predictor matrix.
#' Returns 1 if estimation fails or produces non-finite values.
#'
#' @param env_df A data.frame of (numeric) standardized predictors.
#' @return A numeric bandwidth (scalar or vector as returned by the estimator).
#' @export
compute_global_bw <- function(env_df) {
  bw <- try(hypervolume::estimate_bandwidth(as.matrix(env_df)), silent = TRUE)
  if (inherits(bw, "try-error") || any(!is.finite(bw))) 1 else bw
}
