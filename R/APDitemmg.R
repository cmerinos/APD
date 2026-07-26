#' Item-Level Average Proportional Distance for Multiple Groups (APDitemmg)
#'
#' @description
#' Computes the item-level Average Proportional Distance (APD) separately
#' within multiple groups defined by a grouping variable.
#'
#' @details
#' \code{APDitemmg()} extends the logic of \code{APDitem()} to multigroup
#' settings by estimating item-level APD indices independently within each group.
#'
#' For each item, the procedure computes its absolute discrepancy with all
#' remaining items across respondents within the same group, averages those
#' discrepancies, and rescales the result to the proportional metric:
#'
#' \deqn{
#' APD_j = \frac{AD_j}{K - 1}
#' }
#'
#' where \eqn{AD_j} is the average absolute discrepancy involving item
#' \eqn{j}, and \eqn{K} is the number of response categories.
#'
#' Conceptually, this means that for each item, its absolute discrepancy
#' with all other items is computed across respondents and averaged.
#' This extension is congruent with the spirit of the original APD method,
#' because Sturman et al. (2009) emphasized that APD provides direct
#' information about differences among item scores and penalizes large
#' average discrepancies between items.
#'
#' The resulting item-level APD values provide a fine-grained diagnostic
#' of the contribution of each item to the overall response inconsistency
#' within each group.
#'
#' When \code{ci = TRUE}, bootstrap confidence intervals are computed
#' independently for each item within each group.
#'
#' @param data A \code{data.frame} or matrix containing item responses and,
#' optionally, the grouping variable.
#'
#' @param group A grouping variable. It can be:
#' \itemize{
#'   \item a character string indicating the name of the grouping column
#'   in \code{data}, or
#'   \item a vector of group memberships with length equal to
#'   \code{nrow(data)}.
#' }
#'
#' @param ncat Integer. Number of response categories for the items
#' (e.g., \code{5} for a 1--5 response scale).
#'
#' @param ci Logical. Should bootstrap confidence intervals be computed?
#' Default is \code{FALSE}.
#'
#' @param conf.level Numeric (0--1). Confidence level for bootstrap
#' confidence intervals. Default is \code{0.95}.
#'
#' @param B Integer. Number of bootstrap resamples used when
#' \code{ci = TRUE}. Default is \code{500}.
#'
#' @param cimethod Character. Type of bootstrap confidence interval.
#' Available options are:
#' \itemize{
#'   \item \code{"perc"} — percentile interval.
#'   \item \code{"bca"} — bias-corrected and accelerated interval.
#'   \item \code{"norm"} — normal approximation interval.
#' }
#'
#' @return
#' A named list in which each element corresponds to a group.
#' Each group contains a \code{data.frame} with the following columns:
#'
#' \itemize{
#'   \item \code{item}: item name.
#'   \item \code{AvDiff}: average absolute discrepancy involving the item.
#'   \item \code{APD}: proportional item-level APD value.
#'   \item \code{\%cont}: percentage contribution of the item to the sum
#'   of item-level APD values within the group.
#'   \item \code{lwr.ci}: lower confidence interval bound (if \code{ci = TRUE}).
#'   \item \code{upr.ci}: upper confidence interval bound (if \code{ci = TRUE}).
#' }
#'
#' @references
#' Sturman, D., Cribbie, R. A., & Flett, G. L. (2009).
#' The average distance between item values: A novel approach for estimating
#' internal consistency. \emph{Journal of Psychoeducational Assessment}, 27(5), 409-420.
#' \doi{10.1177/0734282908330937}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' dat <- data.frame(
#'   group = rep(c("A", "B"), each = 100),
#'   Item1 = sample(1:5, 200, replace = TRUE),
#'   Item2 = sample(1:5, 200, replace = TRUE),
#'   Item3 = sample(1:5, 200, replace = TRUE),
#'   Item4 = sample(1:5, 200, replace = TRUE)
#' )
#'
#' APDitemmg(
#'   data = dat,
#'   group = "group",
#'   ncat = 5,
#'   ci = FALSE
#' )
#' }
#'
#' @export
APDitemmg <- function(data,
                      group,
                      ncat,
                      ci = FALSE,
                      conf.level = 0.95,
                      B = 500,
                      cimethod = c("perc", "bca", "norm")) {

  cimethod <- match.arg(cimethod)

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data.frame or matrix.")
  }

  data <- as.data.frame(data)

  if (missing(group)) {
    stop("'group' must be provided.")
  }

  if (is.character(group) && length(group) == 1) {
    if (!group %in% names(data)) {
      stop("Grouping column not found in 'data'.")
    }
    g <- data[[group]]
    data_items <- data[, setdiff(names(data), group), drop = FALSE]
  } else {
    if (length(group) != nrow(data)) {
      stop("When 'group' is a vector, its length must equal nrow(data).")
    }
    g <- group
    data_items <- data
  }

  if (!all(vapply(data_items, is.numeric, logical(1)))) {
    stop("All item columns must be numeric.")
  }

  g <- as.factor(g)
  split_data <- split(data_items, g)

  out <- lapply(split_data, function(dsub) {
    APDitem(
      data = dsub,
      ncat = ncat,
      ci = ci,
      conf.level = conf.level,
      B = B,
      cimethod = cimethod
    )
  })

  out
}
