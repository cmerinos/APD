#' First eigenvalue from average inter-item association (Kaiser's formula)
#'
#' @description
#' Computes the average inter-item association, the first eigenvalue using
#' Kaiser's formula, and the percentage of variance explained by the first
#' component.
#'
#' @param data A data.frame containing item responses.
#' @param cor Type of association: "pearson", "poly", or "gamma".
#' @param nd Number of digits to round all numeric results (default = 3).
#'
#' @return A data.frame with:
#' \itemize{
#'   \item \code{n_items}: Number of items (\eqn{k}).
#'   \item \code{mean_assoc}: Average inter-item association
#'         (mean of the lower triangular elements of the association matrix,
#'         excluding the diagonal).
#'   \item \code{lambda1_kaiser}: First eigenvalue according to Kaiser's
#'         formula, \eqn{1 + (k - 1)\,\bar{r}}.
#'   \item \code{pct_var_percent}: Same quantity expressed as a percentage
#'         (\code{pct_var * 100}).
#' }
#'
#'@references
#'Cureton, E. E. (1971). A Measure of the Average Intercorrelation. 
#'Educational and Psychological Measurement, 31(3), 627-628. 
#'https://doi.org/10.1177/001316447103100303
#'
#'Kaiser, H. F. (1968). A Measure of the Average Intercorrelation. 
#'Educational and Psychological Measurement, 28(2), 245-247. 
#'https://doi.org/10.1177/001316446802800203
#'
#'
#' @examples
#' # Simple example with continuous items (Pearson correlation)
#' set.seed(123)
#' x1 <- rnorm(200)
#' x2 <- 0.6 * x1 + rnorm(200, sd = 0.8)
#' x3 <- 0.6 * x1 + rnorm(200, sd = 0.8)
#' dat <- data.frame(x1, x2, x3)
#'
#' iiacorEigen(dat, cor = "pearson")
#'
#' Example with ordinal items (polychoric correlation)
#' (uncomment if the 'psych' package is installed)
#' library(psych)
#' dat_ord <- data.frame(
#'   item1 = cut(x1, breaks = 4, labels = FALSE),
#'   item2 = cut(x2, breaks = 4, labels = FALSE),
#'   item3 = cut(x3, breaks = 4, labels = FALSE)
#' )
#' iiacorEigen(dat_ord, cor = "poly")
#'
#' @export
iiacorEigen <- function(data, cor = c("pearson", "poly", "gamma"), nd = 3) {
  cor <- match.arg(cor)
  
  if (!is.data.frame(data))
    data <- as.data.frame(data)
  
  k <- ncol(data)
  if (k < 2L)
    stop("'data' must contain at least two item columns.")
  
  # ---- Association matrix ----
  R <- switch(
    cor,
    "pearson" = stats::cor(data, use = "pairwise.complete.obs"),
    
    "poly" = {
      if (!requireNamespace("psych", quietly = TRUE))
        stop("Package 'psych' is required for polychoric correlations.")
      psych::polychoric(data)$rho
    },
    
    "gamma" = {
      if (!requireNamespace("vcd", quietly = TRUE))
        stop("Package 'vcd' is required for gamma.")
      G <- matrix(1, nrow = k, ncol = k)
      colnames(G) <- rownames(G) <- colnames(data)
      for (i in 1:(k-1)) {
        for (j in (i+1):k) {
          tab_ij <- table(data[[i]], data[[j]])
          g_ij <- suppressWarnings(vcd::GKgamma(tab_ij)$gamma)
          G[i, j] <- G[j, i] <- g_ij
        }
      }
      G
    }
  )
  
  # ---- Average inter-item association ----
  mean_assoc <- mean(R[lower.tri(R)])
  
  # ---- Kaiser's formula ----
  lambda1 <- 1 + (k - 1) * mean_assoc
  pct_var_percent <- (lambda1 / k) * 100
  
  # ---- Output with rounding ----
  out <- data.frame(
    n_items         = k,
    mean_assoc      = round(mean_assoc, nd),
    lambda1_kaiser  = round(lambda1, nd),
    pct_var_percent = round(pct_var_percent, nd)
  )
  
  attr(out, "association_matrix") <- R
  attr(out, "cor.type") <- cor
  
  out
}
