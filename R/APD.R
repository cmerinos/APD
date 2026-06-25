#' Average Proportional Distance (APD)
#'
#' @description
#' Average Proportional Distance (APD) is a measure of a test's internal
#' consistency that focuses on the average difference between item scores.
#' To calculate it, the absolute difference is obtained for every pair of
#' item scores, these differences are averaged, and the result is divided
#' by the number of response options minus one to obtain a proportional
#' score ranging from 0 to 1.
#'
#' @details
#' The APD is computed in three main steps:
#' \enumerate{
#'   \item Compute the absolute difference between every possible pair of item
#'   scores for each respondent.
#'   \item Average all absolute differences across respondents and item pairs to
#'   obtain the average difference \code{AD}.
#'   \item Divide \code{AD} by \code{ncat - 1} to rescale it to the
#'   proportional metric:
#'   \deqn{APD = AD / (n_{\mathrm{cat}} - 1).}
#' }
#'
#' The APD ranges from 0 (perfect internal consistency) to 1 (maximum possible
#' inconsistency given the response scale). Lower values reflect greater similarity
#' among item scores, indicating stronger internal consistency.
#'
#' When \code{ci = TRUE}, APD confidence intervals are obtained by nonparametric
#' bootstrap resampling. Three interval types are available:
#' \itemize{
#'   \item \code{"bca"} - bias-corrected and accelerated interval.
#'   \item \code{"perc"} - percentile interval.
#'   \item \code{"norm"} - normal-approximation interval.
#' }
#'
#' @section Interpretation:
#' APD quantifies the average disagreement between item scores on a 0 to 1 scale.
#' The following informal guidelines may help interpretation:
#' \itemize{
#'   \item \strong{APD < 0.20}: typically indicates very good internal consistency.
#'   \item \strong{0.20 <= APD <= 0.25}: often acceptable depending on the construct.
#'   \item \strong{APD > 0.25}: may signal weaker internal consistency or heterogeneous item content.
#' }
#' These values were suggested by Sturman et al. (2009), but they are not
#' strict cutoffs and should be interpreted alongside other reliability
#' evidence (e.g., alpha, omega) and substantive test characteristics.
#'
#' @references
#' Sturman, D., Cribbie, R. A., & Flett, G. L. (2009).
#' The average distance between item values: A novel approach for estimating
#' internal consistency. \emph{Educational and Psychological Measurement},
#' 69(6), 913-932.
#'
#' @param data A \code{data.frame} or matrix containing item responses. Each
#' column represents an item scored on the same categorical scale.
#'
#' @param ncat Integer. Number of response categories for the items
#' (e.g., 5 for a 1-5 scale). Used to convert the average difference into
#' the proportional metric.
#'
#' @param ci Logical. Should a bootstrap confidence interval for APD be computed?
#'
#' @param conf.level Numeric (0 to 1). Confidence level
#' (e.g., \code{0.95}) when \code{ci = TRUE}.
#'
#' @param B Integer. Number of bootstrap resamples used when computing confidence intervals.
#'
#' @param cimethod Character. Type of confidence interval:
#' \code{"bca"}, \code{"perc"}, or \code{"norm"}.
#'
#' @param nd Integer. Number of digits to round the results. Default is 3.
#'
#' @return
#' A one-row \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{Estimate}: character string identifying the reported coefficient.
#'   \item \code{AD}: average absolute difference between all pairs of item scores.
#'   \item \code{APD}: average proportional distance.
#'   \item \code{lwr.ci}: lower confidence bound for APD when \code{ci = TRUE};
#'   otherwise \code{NA}.
#'   \item \code{upr.ci}: upper confidence bound for APD when \code{ci = TRUE};
#'   otherwise \code{NA}.
#' }
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- data.frame(matrix(sample(1:5, 200 * 5, replace = TRUE), ncol = 5))
#'
#' # APD without confidence interval
#' APD(x, ncat = 5, ci = FALSE, conf.level = 0.95, B = 500,
#'     cimethod = "perc", nd = 3)
#'
#' # APD with 95% percentile confidence interval
#' APD(x, ncat = 5, ci = TRUE, conf.level = 0.95, B = 300,
#'     cimethod = "perc", nd = 3)
#'
#' # APD with BCa confidence interval
#' APD(x, ncat = 5, ci = TRUE, conf.level = 0.95, B = 300,
#'     cimethod = "bca", nd = 3)
#' }
#'
#' @export
APD <- function(data, ncat, ci, conf.level, B, cimethod, nd = 3) {

  # Generate all possible combinations of item columns
  column_combinations <- combn(ncol(data), 2)

  # Initialize a list to store the absolute differences for each item pair
  diff_abs_list <- list()

  # Compute the absolute differences for each pair of columns
  for (i in 1:ncol(column_combinations)) {

    # Extract the indices of the two columns to be compared
    col1 <- column_combinations[1, i]
    col2 <- column_combinations[2, i]

    # Compute the absolute difference between the two item scores
    diff_abs_list[[i]] <- abs(data[[col1]] - data[[col2]])
  }

  # Combine all absolute differences into a single vector
  all_diff_abs <- unlist(diff_abs_list)

  # Compute the average absolute difference
  AD <- mean(all_diff_abs, na.rm = TRUE)

  # Compute the Average Proportional Distance
  APD <- round(AD / (ncat - 1), nd)

  # Initialize confidence interval limits
  LCI <- NA
  UCI <- NA

  # Compute confidence intervals if requested
  if (ci) {

    # Internal function to compute APD in a bootstrap sample
    boot_apd <- function(data, ncat) {

      # Draw a bootstrap sample of rows with replacement
      sample_data <- data[sample(1:nrow(data), replace = TRUE), ]

      # Generate all possible item-pair combinations in the bootstrap sample
      column_combinations <- combn(ncol(sample_data), 2)

      # Initialize a list to store bootstrap absolute differences
      diff_abs_list <- list()

      # Compute the absolute differences for each pair of columns
      for (i in 1:ncol(column_combinations)) {
        col1 <- column_combinations[1, i]
        col2 <- column_combinations[2, i]
        diff_abs_list[[i]] <- abs(sample_data[[col1]] - sample_data[[col2]])
      }

      # Combine all bootstrap absolute differences
      all_diff_abs <- unlist(diff_abs_list)

      # Compute bootstrap AD
      AD_boot <- mean(all_diff_abs, na.rm = TRUE)

      # Compute bootstrap APD
      APD_boot <- AD_boot / (ncat - 1)

      return(APD_boot)
    }

    # Generate bootstrap APD estimates
    boot_samples <- lapply(1:B, function(x) boot_apd(data, ncat))
    boot_samples <- unlist(boot_samples)

    # Compute the confidence interval according to the selected method
    if (cimethod == "bca") {

      # Bias-corrected and accelerated confidence interval
      boot_obj <- boot::boot(
        data = data,
        statistic = function(data, idx) boot_apd(data[idx, ], ncat),
        R = B
      )

      # See: https://stats.stackexchange.com/questions/37918/why-is-the-error-estimated-adjustment-a-is-na-generated-from-r-boot-package
      ci_bca <- boot::boot.ci(
        boot_obj,
        type = "bca",
        conf = conf.level,
        L = boot::empinf(boot_obj, index = 1L, type = "jack")
      )

      # Extract BCa lower and upper bounds
      LCI <- round(ci_bca$bca[4], nd)
      UCI <- round(ci_bca$bca[5], nd)

    } else if (cimethod == "perc") {

      # Percentile confidence interval
      LCI <- round(as.numeric(stats::quantile(
        boot_samples,
        probs = (1 - conf.level) / 2
      )), nd)

      UCI <- round(as.numeric(stats::quantile(
        boot_samples,
        probs = 1 - (1 - conf.level) / 2
      )), nd)

    } else if (cimethod == "norm") {

      # Normal-approximation confidence interval based on bootstrap mean and SD
      se <- stats::sd(boot_samples)
      mean_boot <- mean(boot_samples)
      z <- stats::qnorm((1 + conf.level) / 2)

      LCI <- round(mean_boot - z * se, nd)
      UCI <- round(mean_boot + z * se, nd)
    }
  }

  # Create the output data frame in compact format
  resultados <- data.frame(
    Estimate = "APD",
    AD = round(AD, nd),
    APD = APD,
    lwr.ci = LCI,
    upr.ci = UCI
  )

  # Return results
  return(resultados)
}
