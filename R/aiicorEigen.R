#' First eigenvalue from average inter-item association (Kaiser's formula)
#'
#' @description
#' Computes the average inter-item association, the first eigenvalue using
#' Kaiser's formula, and the percentage of variance explained by the first
#' component.
#'
#' @param data A data.frame containing item responses.
#' @param rmethod Type of association: \code{"pearson"}, \code{"spearman"},
#'   or \code{"poly"}.
#' @param nd Number of digits to round all numeric results (default = 3).
#'
#' @return A data.frame with:
#' \itemize{
#'   \item \code{n.items}: Number of items (\eqn{k}).
#'   \item \code{rmean}: Average inter-item association
#'         (mean of the lower triangular elements of the association matrix,
#'         excluding the diagonal).
#'   \item \code{lambda1.Kaiser}: First eigenvalue according to Kaiser's
#'         formula, \eqn{1 + (k - 1)\,\bar{r}}.
#'   \item \code{pct.var}: Percentage of total variance explained by the first
#'         component.
#' }
#'
#' @references
#' Cureton, E. E. (1971). A Measure of the Average Intercorrelation.
#' Educational and Psychological Measurement, 31(3), 627-628.
#' \doi{10.1177/001316447103100303}
#'
#' Kaiser, H. F. (1968). A Measure of the Average Intercorrelation.
#' Educational and Psychological Measurement, 28(2), 245-247.
#' \doi{10.1177/001316446802800203}
#'
#' @examples
#' \donttest{
#' # Example with continuous items (Pearson correlation)
#' set.seed(123)
#' x1 <- rnorm(200)
#' x2 <- 0.6 * x1 + rnorm(200, sd = 0.8)
#' x3 <- 0.6 * x1 + rnorm(200, sd = 0.8)
#' dat <- data.frame(x1, x2, x3)
#'
#' aiicorEigen(dat, rmethod = "pearson")
#' aiicorEigen(dat, rmethod = "spearman")
#' }
#'
#' \donttest{
#' # Example with ordinal items (polychoric correlation)
#' # Uncomment if the 'psych' package is installed
#' dat_ord <- data.frame(
#'   item1 = cut(x1, breaks = 4, labels = FALSE),
#'   item2 = cut(x2, breaks = 4, labels = FALSE),
#'   item3 = cut(x3, breaks = 4, labels = FALSE))
#'
#' aiicorEigen(dat_ord, rmethod = "poly")
#'}
#'
#' @export
aiicorEigen <- function(data,
                        rmethod = c("pearson", "spearman", "poly"),
                        nd = 3) {
  rmethod <- match.arg(rmethod)

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  k <- ncol(data)
  if (k < 2L) {
    stop("'data' must contain at least two item columns.")
  }

  if (rmethod %in% c("pearson", "spearman")) {
    is_num <- vapply(data, is.numeric, logical(1))
    if (!all(is_num)) {
      stop("For rmethod = 'pearson' or 'spearman', all columns in 'data' must be numeric.")
    }
  }

  # ---- Association matrix ----
  R <- switch(
    rmethod,
    "pearson" = stats::cor(data, use = "pairwise.complete.obs", method = "pearson"),

    "spearman" = stats::cor(data, use = "pairwise.complete.obs", method = "spearman"),

    "poly" = {
      if (!requireNamespace("psych", quietly = TRUE)) {
        stop("Package 'psych' is required for polychoric correlations.")
      }
      psych::polychoric(data)$rho
    }
  )

  # ---- Average inter-item association ----
  mean_assoc <- mean(R[lower.tri(R)], na.rm = TRUE)

  # ---- Kaiser's formula ----
  lambda1 <- 1 + (k - 1) * mean_assoc
  pct_var_percent <- (lambda1 / k) * 100

  # ---- Output with rounding ----
  out <- data.frame(
    n.items = k,
    rmean = round(mean_assoc, nd),
    lambda1.Kaiser = round(lambda1, nd),
    pct.var = round(pct_var_percent, nd)
  )

  attr(out, "association_matrix") <- R
  attr(out, "cor.type") <- rmethod

  out
}
