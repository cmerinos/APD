#' Root-Mean-Square Inter-Item Correlation (Meyer, 1975)
#'
#' @description
#' Computes the root-mean-square inter-item correlation (RMSIIC), introduced by
#' Meyer (1975) as a measure of the average inter-item association based on the
#' squared correlation matrix. The function can also compute group-specific
#' estimates and pairwise group comparisons using the MOVER method.
#'
#' @param data A numeric data frame or matrix with items in columns.
#'   Only complete cases are used automatically.
#' @param method Type of correlation: \code{"pearson"} (default) or
#'   \code{"spearman"}.
#' @param nd Number of digits for rounding. Default is 3.
#' @param ci Logical. If \code{TRUE}, bootstrap confidence intervals are computed.
#' @param B Number of bootstrap replications. Default is 1000.
#' @param conf.level Confidence level for the interval. Default is 0.95.
#' @param method.ci Type of bootstrap interval: \code{"perc"} (percentile),
#'   \code{"bca"} (bias-corrected and accelerated), or \code{"norm"}.
#' @param group Optional vector indicating group membership (same length as
#'   \code{nrow(data)}). If supplied, RMSIIC is computed separately for each
#'   group and pairwise differences are compared using MOVER intervals.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{group_results}: A data frame with one row per group (or a single
#'         row for the total sample), containing:
#'         \itemize{
#'           \item \code{group}: Group label.
#'           \item \code{n.items}: Number of items.
#'           \item \code{rms.iic}: RMSIIC estimate.
#'           \item \code{lwr.ci}, \code{upr.ci}: Bootstrap confidence limits
#'                 (if \code{ci=TRUE}).
#'         }
#'   \item \code{comparisons}: A data frame with MOVER-based pairwise group
#'         comparisons:
#'         \itemize{
#'           \item \code{group1}, \code{group2}: Groups being compared.
#'           \item \code{diff}: Difference in RMSIIC (group1 - group2).
#'           \item \code{lwr.ci}, \code{upr.ci}: MOVER confidence limits.
#'         }
#' }
#'
#' @details
#' Let \eqn{R} be the \eqn{p \times p} correlation matrix of a set of items.
#' Meyer (1975) defined the root-mean-square inter-item correlation as:
#'
#' \deqn{
#'   \mathrm{RMSIIC} =
#'   \sqrt{
#'     \frac{
#'       \mathrm{tr}(R^{2}) - p
#'     }{
#'       p(p - 1)
#'     }
#'   } ,
#' }
#'
#' where \eqn{\mathrm{tr}(R^{2})} denotes the trace of the squared correlation
#' matrix. This equals the square root of the average squared inter-item
#' correlation, excluding the diagonal.
#'
#' Meyer argued that RMSIIC is preferable to the simple mean correlation because:
#' \itemize{
#'   \item it treats negative correlations as positive contributions (via
#'         squaring), avoiding sign cancellation,
#'   \item it gives proportionally more weight to strong correlations,
#'   \item it reflects overall item homogeneity more sensitively.
#' }
#'
#' When \code{group} is supplied, the statistic is computed independently within
#' each group. Differences between groups are evaluated using the MOVER method,
#' combining group-specific confidence intervals into an interval for the
#' difference without assuming asymptotic normality.
#'
#' Bootstrap confidence intervals may be computed using percentile, BCa, or
#' normal-based methods.
#'
#' @references
#' Meyer, E. P. (1975). A measure of the average intercorrelation.
#' \emph{Educational and Psychological Measurement, 35}(1), 67–72.
#' https://doi.org/10.1177/001316447503500107
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(200), ncol = 5)
#'
#' # Single-sample RMSIIC
#' rmsiic(X)
#'
#' # Bootstrap CI
#' rmsiic(X, ci = TRUE, B = 500)
#'
#' # Group comparison using MOVER
#' g <- rep(c("A","B"), each = 20)
#' rmsiic(X[1:40, ], group = g, ci = TRUE, B = 300)
#'
#' @export
rmsiic <- function(data,
                   method     = c("pearson", "spearman"),
                   nd         = 3,
                   ci         = FALSE,
                   B          = 1000,
                   conf.level = 0.95,
                   method.ci  = c("perc", "bca", "norm"),
                   group      = NULL) {
  # ---- ARGUMENT CHECKS ----
  if (!is.data.frame(data)) data <- as.data.frame(data)

  method    <- match.arg(method)
  method.ci <- match.arg(method.ci)

  if (!is.numeric(nd) || nd < 0 || nd %% 1 != 0)
    stop("Argument 'nd' must be a non-negative integer.")

  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("Argument 'conf.level' must be a number between 0 and 1.")

  if (!is.numeric(B) || B <= 0 || B %% 1 != 0)
    stop("Argument 'B' must be a positive integer.")

  if (!is.null(group) && length(group) != nrow(data))
    stop("Length of 'group' must match the number of rows in 'data'.")

  if (ci && method.ci == "bca" && !requireNamespace("boot", quietly = TRUE))
    stop("Package 'boot' is required for BCa intervals.", call. = FALSE)

  # ---- HANDLE MISSING VALUES ----
  if (is.null(group)) {
    cc <- stats::complete.cases(data)
  } else {
    cc <- stats::complete.cases(data) & !is.na(group)
  }

  if (!all(cc)) {
    data  <- data[cc, , drop = FALSE]
    if (!is.null(group)) group <- group[cc]
    message("Rows with missing values were removed.")
  }

  if (nrow(data) <= 3)
    stop("Not enough data after removing missing values.")

  # ---- NUMERIC COLUMNS ONLY ----
  num_cols <- vapply(data, is.numeric, logical(1))
  if (!all(num_cols)) {
    data <- data[, num_cols, drop = FALSE]
    message("Non-numeric columns were removed.")
  }

  if (ncol(data) < 2)
    stop("At least two numeric items are required.")

  p <- ncol(data)

  # ---- AUXILIARY FUNCTION: RMSIIC FOR ONE GROUP ----
  compute_rmsiic <- function(sub_data, group_name = "total") {

    if (nrow(sub_data) <= 3 || ncol(sub_data) < 2) {
      stop(paste0("Group '", group_name,
                  "' does not have enough data (need >= 4 rows and >= 2 items)."))
    }

    R <- stats::cor(sub_data, method = method)

    # Aviso si hay correlaciones negativas
    if (any(R[upper.tri(R, diag = FALSE)] < 0)) {
      warning("Some inter-item correlations are negative in the correlation matrix.")
    }

    sum_sq_total <- sum(R^2)
    sum_sq_diag  <- p
    sum_sq_off   <- sum_sq_total - sum_sq_diag
    rms_val      <- sqrt(sum_sq_off / (p * (p - 1)))

    res_row <- data.frame(
      group   = group_name,
      n.items = p,
      rms.iic = round(rms_val, nd),
      lwr.ci  = NA_real_,
      upr.ci  = NA_real_,
      row.names = NULL
    )

    # ---- No CI requested ----
    if (!ci) return(res_row)

    # ---- Bootstrap CI ----
    n_g <- nrow(sub_data)

    stat_fun <- function(d, idx) {
      d_b <- d[idx, , drop = FALSE]
      Rb  <- stats::cor(d_b, method = method)
      sum_sq_total_b <- sum(Rb^2)
      sum_sq_off_b   <- sum_sq_total_b - p
      sqrt(sum_sq_off_b / (p * (p - 1)))
    }

    boot_obj <-
      if (requireNamespace("boot", quietly = TRUE)) {
        boot::boot(sub_data, statistic = stat_fun, R = B)
      } else {
        vals <- numeric(B)
        for (b in seq_len(B)) {
          idx_b <- sample.int(n_g, n_g, replace = TRUE)
          vals[b] <- stat_fun(sub_data, idx_b)
        }
        structure(list(t = matrix(vals, ncol = 1), t0 = rms_val), class = "boot")
      }

    alpha <- 1 - conf.level

    if (method.ci == "bca" && requireNamespace("boot", quietly = TRUE)) {
      ci_obj   <- boot::boot.ci(boot_obj, type = "bca", conf = conf.level)
      ci_limits <- ci_obj$bca[4:5]

    } else if (method.ci == "perc") {
      vals <- as.numeric(boot_obj$t)
      ci_limits <- stats::quantile(
        vals,
        probs = c(alpha/2, 1 - alpha/2),
        type = 6,
        na.rm = TRUE
      )

    } else { # "norm"
      vals <- as.numeric(boot_obj$t)
      m <- mean(vals, na.rm = TRUE)
      s <- stats::sd(vals,  na.rm = TRUE)
      z <- stats::qnorm(1 - alpha/2)
      ci_limits <- c(m - z*s, m + z*s)
    }

    res_row$lwr.ci <- round(ci_limits[1], nd)
    res_row$upr.ci <- round(ci_limits[2], nd)

    res_row
  }

  # ---- GROUP-WISE RESULTS ----
  if (is.null(group)) {
    group_results <- compute_rmsiic(data, group_name = "total")
  } else {
    g_levels <- unique(group)
    group_results <- do.call(
      rbind,
      lapply(g_levels, function(g) {
        compute_rmsiic(data[group == g, , drop = FALSE], group_name = g)
      })
    )
    rownames(group_results) <- NULL
  }

  # ---- MOVER COMPARISONS ----
  if (nrow(group_results) > 1) {
    gnames <- unique(group_results$group)
    combs  <- utils::combn(gnames, 2, simplify = FALSE)

    comparisons <- lapply(combs, function(gs) {
      g1 <- gs[1]; g2 <- gs[2]

      d1 <- group_results[group_results$group == g1, ]
      d2 <- group_results[group_results$group == g2, ]

      T1 <- d1$rms.iic; L1 <- d1$lwr.ci; U1 <- d1$upr.ci
      T2 <- d2$rms.iic; L2 <- d2$lwr.ci; U2 <- d2$upr.ci

      D <- T1 - T2

      lower_D <- D - sqrt((T1 - L1)^2 + (U2 - T2)^2)
      upper_D <- D + sqrt((U1 - T1)^2 + (T2 - L2)^2)

      data.frame(
        group1 = g1,
        group2 = g2,
        diff   = round(D,       nd),
        lwr.ci = round(lower_D, nd),
        upr.ci = round(upper_D, nd),
        row.names = NULL
      )
    })

    comparisons <- do.call(rbind, comparisons)

  } else {
    comparisons <- data.frame(
      group1 = character(0),
      group2 = character(0),
      diff   = numeric(0),
      lwr.ci = numeric(0),
      upr.ci = numeric(0)
    )
  }

  list(
    group_results = group_results,
    comparisons   = comparisons
  )
}
