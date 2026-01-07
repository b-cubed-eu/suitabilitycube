#' Download GBIF occurrences for one species and return sf points
#'
#' @param scientific_name Species scientific name.
#' @param country_iso ISO2 country code (e.g., "IT").
#' @param years Length-2 vector: c(start, end).
#' @param limit Max records returned by GBIF call.
#' @return An sf object (EPSG:4326) or NULL if none found.
#' @export
gbif_occ_sf <- function(scientific_name, country_iso = "IT", years = c(2000, 2020), limit = 20000) {
  key <- rgbif::name_backbone(name = scientific_name)$usageKey
  if (is.null(key) || is.na(key)) return(NULL)

  dat <- rgbif::occ_search(
    taxonKey = key,
    country  = country_iso,
    year     = paste0(years[1], ",", years[2]),
    basisOfRecord = "HUMAN_OBSERVATION",
    occurrenceStatus = "PRESENT",
    hasCoordinate   = TRUE,
    limit = limit
  )$data

  if (is.null(dat) || nrow(dat) == 0) return(NULL)

  dat <- dat[!is.na(dat$decimalLongitude) & !is.na(dat$decimalLatitude), ]
  dat <- dat[!duplicated(dat[, c("decimalLongitude","decimalLatitude")]), ]

  sf::st_as_sf(dat, coords = c("decimalLongitude","decimalLatitude"), crs = 4326, remove = FALSE)
}
