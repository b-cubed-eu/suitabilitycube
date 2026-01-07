#' Download GBIF occurrences for multiple species
#'
#' @param species_vec Character vector of species scientific names.
#' @param country_iso ISO2 country code.
#' @param years Length-2 vector.
#' @param limit Max records per species.
#' @return Named list of sf objects (or NULL entries).
#' @export
gbif_occ_list <- function(species_vec, country_iso, years, limit) {
  setNames(lapply(species_vec, \(sp) gbif_occ_sf(sp, country_iso, years, limit)), species_vec)
}
