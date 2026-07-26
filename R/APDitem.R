#' Item-Level Average Proportional Distance (APDitem)
#'
#' @description
#' Computes the item-level Average Proportional Distance (APD) for a single
#' data frame or matrix of item responses. For each item, the function
#' calculates its average absolute discrepancy against all remaining items
#' across respondents, and rescales it to the proportional metric.
#' For a scale with p items and K response categories, the total APD is defined as
#' the average absolute discrepancy across all within-person item pairs, divided
#' by K-1. Extending this logic, the item-level APD for item j is the average absolute
#' discrepancy between that item and all remaining items across respondents, also
#' divided by K-1. Under this definition, the total APD equals the mean of the
#' item-level APD values.
#'
#' @details
#' For each item \eqn{j}, the average distance is computed as:
#'
#' \deqn{
#' AvDiff_j = \frac{1}{n(p-1)} \sum_{i=1}^{n} \sum_{k \neq j} |x_{ij} - x_{ik}|
#' }
#'
#' The proportional index is:
#'
#' \deqn{
#' APD_j = \frac{AvDiff_j}{n_{cat} - 1}
#' }
#'
#' Under this definition, the mean of the item-level APD values equals the
#' total APD computed from the same data matrix.
#'
#' When \code{ci = TRUE}, bootstrap confidence intervals are computed for the
#' item-level APD values.
#'
#' @param data A \code{data.frame} or matrix containing item responses.
#' Each column represents an item scored on the same categorical scale.
#' @param ncat Integer. Number of response categories for the items
#' (e.g., \code{5} for a 1--5 scale).
#' @param ci Logical. Should bootstrap confidence intervals be computed?
#' Default is \code{FALSE}.
#' @param conf.level Numeric. Confidence level for the intervals.
#' Default is \code{0.95}.
#' @param B Integer. Number of bootstrap resamples. Default is \code{500}.
#' @param cimethod Character. Confidence interval method:
#' \code{"perc"}, \code{"bca"}, or \code{"norm"}.
#'
#' @return
#' A \code{data.frame} with one row per item and the following columns:
#' \itemize{
#'   \item \code{item}: item name.
#'   \item \code{AvDiff}: average absolute discrepancy of the item against all remaining items.
#'   \item \code{APD}: proportional item-level APD.
#'   \item \code{\%cont}: percentage contribution of the item to the sum of item-level APD values.
#'   \item \code{lwr.ci}: lower confidence bound (if \code{ci = TRUE}).
#'   \item \code{upr.ci}: upper confidence bound (if \code{ci = TRUE}).
#' }
#'
#' @export
APDitem <- function(data,
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

  if (ncol(data) < 2) {
    stop("'data' must contain at least 2 item columns.")
  }

  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop("All columns in 'data' must be numeric.")
  }

  if (!is.numeric(ncat) || length(ncat) != 1 || is.na(ncat) || ncat <= 1) {
    stop("'ncat' must be a single number greater than 1.")
  }

  if (!is.logical(ci) || length(ci) != 1 || is.na(ci)) {
    stop("'ci' must be TRUE or FALSE.")
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }

  if (!is.numeric(B) || length(B) != 1 || is.na(B) || B < 10) {
    stop("'B' must be a single number >= 10.")
  }

  if (is.null(colnames(data))) {
    colnames(data) <- paste0("Item", seq_len(ncol(data)))
  }

  .apditem_core <- function(dat, ncat) {
    p <- ncol(dat)
    avdiff <- numeric(p)

    for (j in seq_len(p)) {
      others <- setdiff(seq_len(p), j)
      diffs  <- lapply(others, function(k) abs(dat[[j]] - dat[[k]]))
      avdiff[j] <- mean(unlist(diffs), na.rm = TRUE)
    }

    apd <- avdiff / (ncat - 1)
    pct <- (apd / sum(apd, na.rm = TRUE)) * 100

    out <- data.frame(
      item   = colnames(dat),
      AvDiff = avdiff,
      APD    = apd,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    out[["%cont"]] <- pct
    out
  }

  res <- .apditem_core(data, ncat = ncat)

  if (!ci) {
    res$AvDiff   <- round(res$AvDiff, 3)
    res$APD      <- round(res$APD, 3)
    res[["%cont"]] <- round(res[["%cont"]], 2)
    return(res)
  }

  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required when ci = TRUE.")
  }

  boot_fun <- function(dat, idx) {
    dat_boot <- dat[idx, , drop = FALSE]
    .apditem_core(dat_boot, ncat = ncat)$APD
  }

  boot_obj <- boot::boot(data = data, statistic = boot_fun, R = B)

  lwr <- upr <- rep(NA_real_, ncol(data))

  if (cimethod == "perc") {
    alpha <- (1 - conf.level) / 2
    lwr <- apply(boot_obj$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
    upr <- apply(boot_obj$t, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  }

  if (cimethod == "norm") {
    z <- stats::qnorm((1 + conf.level) / 2)
    m <- colMeans(boot_obj$t, na.rm = TRUE)
    s <- apply(boot_obj$t, 2, stats::sd, na.rm = TRUE)
    lwr <- m - z * s
    upr <- m + z * s
  }

  if (cimethod == "bca") {
    for (j in seq_len(ncol(data))) {
      ci_j <- try(
        boot::boot.ci(boot_obj, type = "bca", conf = conf.level, index = j),
        silent = TRUE
      )

      if (!inherits(ci_j, "try-error") && !is.null(ci_j$bca)) {
        lwr[j] <- ci_j$bca[4]
        upr[j] <- ci_j$bca[5]
      }
    }
  }

  res$lwr.ci <- lwr
  res$upr.ci <- upr

  res$AvDiff   <- round(res$AvDiff, 3)
  res$APD      <- round(res$APD, 3)
  res[["%cont"]] <- round(res[["%cont"]], 2)
  res$lwr.ci   <- round(res$lwr.ci, 3)
  res$upr.ci   <- round(res$upr.ci, 3)

  res
}


#' Item-Level Average Proportional Distance for Multiple Groups (APDitemmg)
#'
#' @description
#' Computes item-level Average Proportional Distance (APD) separately for each
#' group defined by a grouping variable.
#'
#' @details
#' The function splits the data according to \code{group}, applies
#' \code{APDitem()} to each subgroup, and returns a named list of results.
#'
#' @param data A \code{data.frame} or matrix containing item responses and,
#' if needed, the grouping variable.
#' @param group A grouping variable. It can be:
#' \itemize{
#'   \item a character string naming the grouping column in \code{data}, or
#'   \item a vector of group memberships of length equal to \code{nrow(data)}.
#' }
#' @param ncat Integer. Number of response categories for the items.
#' @param ci Logical. Should bootstrap confidence intervals be computed?
#' Default is \code{FALSE}.
#' @param conf.level Numeric. Confidence level for the intervals.
#' Default is \code{0.95}.
#' @param B Integer. Number of bootstrap resamples. Default is \code{500}.
#' @param cimethod Character. Confidence interval method:
#' \code{"perc"}, \code{"bca"}, or \code{"norm"}.
#'
#' @return
#' A named list. Each element corresponds to one group and contains a
#' \code{data.frame} with columns:
#' \code{item}, \code{AvDiff}, \code{APD}, \code{\%cont}, \code{lwr.ci}, and \code{upr.ci}
#' (the last two only if \code{ci = TRUE}).
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
