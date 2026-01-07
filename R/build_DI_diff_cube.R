#' Build a pairwise DI-difference cube
#'
#' Creates a stars cube with DI differences for each species pair:
#' \deqn{DI\_diff = DI(species_i) - DI(species_j)}
#'
#' Output dimensions are: cell × comparison × time, where "comparison" encodes
#' the ordered pair (sp1 - sp2).
#'
#' @param data_cube A stars object containing a "DI" attribute, with dimensions
#'   cell × taxon × time (as created by suitabilitycube cube-building functions).
#' @param pairs Optional. If NULL (default), uses all pairwise combinations i < j.
#'   Otherwise, provide either:
#'   - a character matrix/data.frame with two columns (species labels), or
#'   - a numeric matrix with two columns (1-based taxon indices).
#'
#' @return A stars object with attribute "DI_diff" and dimensions
#'   cell × comparison × time.
#' @export
build_DI_diff_cube <- function(data_cube, pairs = NULL) {
  stopifnot(inherits(data_cube, "stars"))

  dims <- stars::st_dimensions(data_cube)
  if (is.null(dims$taxon) || is.null(dims$time)) {
    stop("data_cube must have 'taxon' and 'time' dimensions.")
  }

  tax <- trimws(as.character(dims$taxon$values))

  di <- data_cube["DI"]
  arr <- di[[1]]  # array [cell, taxon, time]

  stopifnot(length(tax) >= 2)
  ncell <- dim(arr)[1]
  ntax  <- dim(arr)[2]
  ntime <- dim(arr)[3]

  # Resolve pairs (default: all i<j)
  if (is.null(pairs)) {
    pairs_idx <- t(utils::combn(ntax, 2))
  } else if (is.character(pairs)) {
    pairs <- as.matrix(pairs)
    if (ncol(pairs) != 2) stop("`pairs` must have two columns (sp1, sp2).")

    pairs_idx <- cbind(
      match(trimws(pairs[, 1]), tax),
      match(trimws(pairs[, 2]), tax)
    )
  } else {
    pairs_idx <- as.matrix(pairs)
    if (ncol(pairs_idx) != 2) stop("`pairs` must have two columns (sp1, sp2).")
  }

  if (any(is.na(pairs_idx))) stop("Pairs include unknown species labels or invalid indices.")

  K <- nrow(pairs_idx)
  arr_diff <- array(NA_real_, dim = c(ncell, K, ntime))
  comp_labels <- character(K)

  for (k in seq_len(K)) {
    i <- pairs_idx[k, 1]
    j <- pairs_idx[k, 2]
    arr_diff[, k, ] <- arr[, i, ] - arr[, j, ]
    comp_labels[k]  <- paste0(tax[i], " - ", tax[j])
  }

  # Create new dimensions object (cell and time from original, new comparison)
  dims_new <- dims
  names(dims_new)[2] <- "comparison"
  dims_new$comparison$values <- comp_labels

  stars::st_as_stars(list(DI_diff = arr_diff), dimensions = dims_new)
}
