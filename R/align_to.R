#' Align a raster stack to a target grid (bilinear)
#'
#' @param src A terra SpatRaster.
#' @param target A terra SpatRaster used as target geometry.
#' @return A SpatRaster aligned to target.
#' @export
align_to <- function(src, target) {
  terra::resample(src, target, method = "bilinear")
}
