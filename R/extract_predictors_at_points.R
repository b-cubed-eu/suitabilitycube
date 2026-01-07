#' Extract environmental predictors at occurrence points
#'
#' Extracts raster values at point locations, keeps only selected predictor
#' variables, removes rows with missing values, and drops numeric predictors
#' with zero variance.
#'
#' @param rst A terra::SpatRaster containing environmental predictors.
#' @param occ_sf An sf object with point geometries (occurrence locations).
#' @param pred_vars Character vector of predictor variable names to keep.
#'
#' @return A data.frame of predictor values at occurrences.
extract_predictors_at_points <- function(rst, occ_sf, pred_vars) {
  if (!is.na(sf::st_crs(occ_sf)) && !is.na(terra::crs(rst)) &&
      sf::st_crs(occ_sf)$wkt != terra::crs(rst)) {
    occ_sf <- sf::st_transform(occ_sf, crs = terra::crs(rst))
  }

  vals <- terra::extract(rst, occ_sf)
  df <- dplyr::select(as.data.frame(vals), dplyr::all_of(pred_vars)) |> stats::na.omit()

  keep <- vapply(df, function(x) is.numeric(x) && stats::sd(x, na.rm = TRUE) > 0, logical(1))
  df[, keep, drop = FALSE]
}
