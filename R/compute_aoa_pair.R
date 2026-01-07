#' Compute AOA and DI for present and future rasters
#'
#' Builds a training table from occurrence predictors (dropping zero-variance
#' numeric columns), then computes CAST::aoa for a present and future raster.
#'
#' @param env_train_df Training data.frame of predictors at occurrences.
#' @param new_present terra::SpatRaster of predictors for the present.
#' @param new_future terra::SpatRaster of predictors for the future.
#'
#' @return A list with elements "present" and "future", each being CAST::aoa output.
#' @export
compute_aoa_pair <- function(env_train_df, new_present, new_future) {
  keep <- vapply(env_train_df, function(x) is.numeric(x) && stats::sd(x, na.rm = TRUE) > 0, logical(1))
  trn  <- env_train_df[, keep, drop = FALSE]
  vars <- colnames(trn)

  list(
    present = CAST::aoa(newdata = new_present, train = trn, variables = vars),
    future  = CAST::aoa(newdata = new_future,  train = trn, variables = vars)
  )
}
