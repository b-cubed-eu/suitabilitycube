#' Merge AOA, DI and HV cubes into a single multi-attribute stars object
#'
#' @param AOA_cube stars cube with attribute "AOA".
#' @param DI_cube  stars cube with attribute "DI".
#' @param HV_cube  stars cube with attribute "HV".
#' @param species_vec Character vector of species names to set on the taxon dimension.
#' @param time_values Character vector for time dimension (default: c("present","future")).
#'
#' @return A multi-attribute stars object containing AOA, DI, HV.
#' @export
merge_suitability_cube <- function(AOA_cube, DI_cube, HV_cube,
                                   species_vec, time_values = c("present", "future")) {
  stopifnot(inherits(AOA_cube, "stars"), inherits(DI_cube, "stars"), inherits(HV_cube, "stars"))
  stopifnot(is.character(species_vec))
  stopifnot(is.character(time_values), length(time_values) == 2)

  data_cube <- c(AOA_cube, DI_cube, HV_cube)
  data_cube <- stars::st_set_dimensions(data_cube, "taxon", values = species_vec)
  data_cube <- stars::st_set_dimensions(data_cube, "time", values = time_values)
  data_cube
}
