#' Default parameters (edit this list in your workflow)
#'
#' This function returns the same `params` list used in the tutorial scripts.
#' Users should assign it and edit values, e.g.:
#' `params <- params_template(); params$country_name <- "Italy"`.
#'
#' @return A named list of parameters.
#' @export
params_template <- function() {
  list(
    species      = c("Bufo bufo", "Bufotes viridis", "Bombina variegata"),
    country_name = "Italy",
    country_iso  = "IT",
    res_arcmin   = 2.5,
    ssp_code     = "245",
    gcm_model    = "BCC-CSM2-MR",
    period       = "2041-2060",
    outdir       = tempdir(),      # change to a persistent path if desired
    gbif_years   = c(2010, 2020),
    gbif_limit   = 20000,
    cor_thr      = 0.7,
    cor_frac     = 0.10,
    grid_cellsize_deg = 0.25,      # ~25 km
    grid_square      = FALSE       # FALSE = hex, TRUE = square
  )
}
