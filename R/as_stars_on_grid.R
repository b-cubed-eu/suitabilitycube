#' Aggregate a single-layer raster/stars to a polygon grid and return a stars object
#'
#' Converts input to stars (if needed), aligns CRS (warp), aggregates values over
#' polygon geometries, and sets the geometry dimension name to "cell".
#'
#' Special case: if `name` contains "AOA" (case-insensitive), aggregation uses
#' a mode function (useful for binary/categorical AOA masks).
#'
#' @param sr A terra::SpatRaster or a stars object (single-layer preferred).
#' @param grid An sf object of polygon grid cells (from make_country_grid()).
#' @param fun Aggregation function (e.g., mean). Ignored for AOA where mode is used.
#' @param name Optional name for the resulting stars attribute (e.g., "AOA", "DI").
#' @param na.rm Logical; remove NA values in aggregation.
#'
#' @return A stars object aggregated to the grid, with dimension "cell".
#' @export
as_stars_on_grid <- function(sr, grid, fun = mean, name = NULL, na.rm = TRUE) {
  stopifnot(inherits(grid, "sf"))

  mode_fun <- function(x, na.rm = TRUE) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  if (!is.null(name) && grepl("AOA", name, ignore.case = TRUE)) {
    fun <- mode_fun
  }

  # 1) convert to stars if needed
  sr_st <- if (inherits(sr, "SpatRaster")) stars::st_as_stars(sr) else sr
  if (!inherits(sr_st, "stars")) stop("`sr` must be a terra::SpatRaster or a stars object.")

  # keep one layer if multiple
  if (length(names(sr_st)) != 1L) sr_st <- sr_st[1]

  # 2) align CRS
  crs_sr   <- sf::st_crs(sr_st)
  crs_grid <- sf::st_crs(grid)
  if (!is.na(crs_sr) && !is.na(crs_grid) && crs_sr != crs_grid) {
    sr_st <- stars::st_warp(sr_st, crs = crs_grid)
  }

  # 3) aggregate over polygons
  agg <- suppressWarnings(
    stats::aggregate(sr_st, by = sf::st_geometry(grid), FUN = fun, na.rm = na.rm)
  )

  # 4) rename geometry dimension to "cell" and set values to grid geometry
  if ("geometry" %in% names(stars::st_dimensions(agg))) {
    agg <- stars::st_set_dimensions(agg, "geometry", names = "cell")
  }
  agg <- stars::st_set_dimensions(agg, "cell", values = sf::st_geometry(grid))

  # 5) optional rename
  if (!is.null(name)) names(agg) <- name

  agg
}
