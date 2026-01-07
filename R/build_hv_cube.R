#' Build an HV cube aligned to an existing metric cube
#'
#' Creates a stars cube with the same dimensions as `template_cube` (cell × taxon × time).
#' HV values are scalars per species and are placed only in the "present" time slice.
#' The "future" slice is NA.
#'
#' @param hv_by_species Named list of hypervolume volumes by species (numeric/NA).
#' @param template_cube A stars cube (typically the AOA cube) providing dimensions.
#' @param species_vec Character vector of species names/order to use.
#'
#' @return A stars object with attribute "HV".
#' @export
build_hv_cube <- function(hv_by_species, template_cube, species_vec) {
  stopifnot(is.list(hv_by_species))
  stopifnot(inherits(template_cube, "stars"))
  stopifnot(is.character(species_vec), length(species_vec) > 0)

  dims  <- stars::st_dimensions(template_cube)
  shape <- dim(template_cube[[1]])  # first attribute array shape

  # gather HV scalars in species order (missing species -> NA)
  hv_vals <- vapply(species_vec, function(sp) as.numeric(hv_by_species[[sp]]), numeric(1))

  hv_arr <- array(
    NA_real_,
    shape,
    dimnames = list(NULL, dims$taxon$values, dims$time$values)
  )

  i_present <- match("present", dims$time$values)
  if (is.na(i_present)) stop("template_cube must have time dimension including 'present'.")

  for (j in seq_along(species_vec)) {
    hv_arr[, j, i_present] <- hv_vals[j]
  }

  stars::st_as_stars(list(HV = hv_arr), dimensions = dims)
}
