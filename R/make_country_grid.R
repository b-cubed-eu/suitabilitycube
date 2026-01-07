#' Build an analysis grid over a country polygon
#'
#' Creates a regular grid (square or hex) over the country extent and assigns a
#' sequential cell ID. Output CRS is set to EPSG:4326 unless the input is already
#' in another CRS (you can reproject beforehand if needed).
#'
#' @param country_sf An sf polygon/multipolygon of the country boundary.
#' @param cellsize_deg Grid cell size in degrees (e.g., 0.25 ~ 25 km at equator).
#' @param square Logical; TRUE for square cells, FALSE for hexagonal.
#'
#' @return An sf object of grid polygons with a `cell` column.
#' @export
make_country_grid <- function(country_sf, cellsize_deg = 0.25, square = FALSE) {
  stopifnot(inherits(country_sf, "sf"))
  stopifnot(is.numeric(cellsize_deg), length(cellsize_deg) == 1, cellsize_deg > 0)

  # buffer(0) to fix potential invalid geometries
  country_sf <- sf::st_buffer(country_sf, 0)

  grid_cells <- sf::st_make_grid(
    country_sf,
    cellsize = cellsize_deg,
    what = "polygons",
    square = isTRUE(square)
  ) |>
    sf::st_as_sf() |>
    dplyr::mutate(cell = seq_len(dplyr::n()))

  # You forced EPSG:4326 in the script. Keep that behavior.
  sf::st_crs(grid_cells) <- 4326
  grid_cells
}
