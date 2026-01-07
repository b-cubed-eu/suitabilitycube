#' Compute Gaussian hypervolume volume
#'
#' Computes a Gaussian hypervolume using hypervolume::hypervolume_gaussian and
#' returns the volume. Returns NA if there are insufficient unique points or if
#' the computation fails.
#'
#' @param env_df A data.frame (typically standardized) with only numeric predictors.
#' @param bw Bandwidth for KDE (as from compute_global_bw()).
#' @param spp Samples per point passed to hypervolume_gaussian.
#' @return Numeric hypervolume volume, or NA_real_ on failure.
hyp_calc <- function(env_df, bw, spp = 50L) {
  X <- as.matrix(env_df)

  if (nrow(unique(X)) < (ncol(X) + 1L)) return(NA_real_)

  hv <- try(
    hypervolume::hypervolume_gaussian(
      X,
      kde.bandwidth = bw,
      samples.per.point = spp,
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (inherits(hv, "try-error")) return(NA_real_)
  hv@Volume
}
