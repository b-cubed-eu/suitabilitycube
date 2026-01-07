#' Drop highly correlated raster layers using a greedy scheme
#'
#' @param rst A terra SpatRaster with multiple layers.
#' @param thr Correlation threshold.
#' @param frac Fraction of cells to sample for correlation calculation.
#' @param seed Random seed.
#' @return A list with elements: cor (matrix), selected (character), dropped (character).
#' @export
drop_high_corr <- function(rst, thr = 0.7, frac = 0.1, seed = 42) {
  set.seed(seed)

  sz  <- max(1000, round(terra::ncell(rst) * frac))
  smp <- terra::spatSample(rst, size = sz, method = "random", na.rm = TRUE, as.points = FALSE)
  mat <- as.matrix(smp)

  keep_cols <- which(colSums(!is.na(mat)) > 0)
  mat <- mat[, keep_cols, drop = FALSE]

  cm  <- suppressWarnings(stats::cor(mat, use = "pairwise.complete.obs", method = "pearson"))

  to_drop <- character(0)
  vars <- colnames(cm)

  avg_abs <- sort(colMeans(abs(cm), na.rm = TRUE), decreasing = TRUE)
  for (v in names(avg_abs)) {
    if (v %in% to_drop) next
    high <- setdiff(names(which(abs(cm[v, ]) > thr)), v)
    to_drop <- union(to_drop, high)
  }

  list(cor = cm, selected = setdiff(vars, to_drop), dropped = to_drop)
}
