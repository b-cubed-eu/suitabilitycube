#' Z-standardize predictors with mean imputation for missing values
#'
#' Numeric columns are standardized (mean 0, sd 1). Missing values are replaced
#' with the column mean prior to standardization. Columns with zero variance
#' become all zeros.
#'
#' @param df A data.frame of predictors.
#' @return A data.frame with transformed numeric columns.
#' @export
z_transform <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      nas <- is.na(col)
      m <- mean(col, na.rm = TRUE)
      col[nas] <- m
      s <- stats::sd(col, na.rm = TRUE)

      if (is.na(s) || s == 0) return(rep(0, length(col)))
      (col - m) / s
    } else {
      col
    }
  })
  df
}
