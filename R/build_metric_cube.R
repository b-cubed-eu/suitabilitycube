#' Build a metric cube (AOA or DI) for present and future
#'
#' Aggregates AOA/DI rasters to a grid per species and stacks them into a stars cube
#' with dimensions: cell × taxon × time.
#'
#' @param aoa_di_by_species Named list by species. Each element must have
#'   `$present` and `$future`, with sub-elements "AOA" and/or "DI" (as in CAST::aoa output).
#' @param species_vec Character vector of species names to include (order matters).
#' @param grid_cells sf polygon grid (from make_country_grid()).
#' @param metric One of "AOA" or "DI".
#' @param fun Aggregation function used for DI (mean by default). AOA uses mode automatically.
#'
#' @return A stars object with one attribute named as `metric`.
#' @export
build_metric_cube <- function(aoa_di_by_species, species_vec, grid_cells,
                              metric = c("AOA", "DI"), fun = mean) {
  metric <- match.arg(metric)
  stopifnot(is.list(aoa_di_by_species))
  stopifnot(is.character(species_vec), length(species_vec) > 0)
  stopifnot(inherits(grid_cells, "sf"))

  # keep only species with non-null data
  species_vec <- species_vec[
    species_vec %in% names(aoa_di_by_species) &
      !vapply(aoa_di_by_species[species_vec], is.null, logical(1))
  ]
  if (length(species_vec) == 0) stop("No valid species found in aoa_di_by_species.")

  pres_list <- lapply(
    species_vec,
    \(sp) as_stars_on_grid(aoa_di_by_species[[sp]]$present[[metric]], grid_cells, fun = fun, name = metric)
  )
  fut_list <- lapply(
    species_vec,
    \(sp) as_stars_on_grid(aoa_di_by_species[[sp]]$future[[metric]],  grid_cells, fun = fun, name = metric)
  )

  pres <- do.call(c, pres_list) |>
    stars::st_redimension() |>
    stars::st_set_dimensions(2, values = species_vec, names = "taxon")

  fut <- do.call(c, fut_list) |>
    stars::st_redimension() |>
    stars::st_set_dimensions(2, values = species_vec, names = "taxon")

  cube <- c(pres, fut, along = list(time = c("present", "future"))) |>
    stars::st_set_dimensions(1, names = "cell")

  names(cube) <- metric
  cube
}
