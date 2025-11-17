#' @title Average Inter-Item Correlation with Bootstrap CI
#'
#' @description
#' Computes the average inter-item correlation (iia) for a set of items, with optional bootstrap
#' confidence intervals. It supports Pearson or Spearman correlations and allows grouped analysis.
#' Fisher's z-transformation is used to average correlations appropriately. When groups are provided,
#' the function also computes MOVER-type confidence intervals for the difference in iia between
#' all pairs of groups.
#'
#' @details
#' For each group, the function estimates the inter-item correlation matrix, extracts the lower-triangular
#' elements, and averages them after applying Fisher's \eqn{z = \mathrm{arctanh}(r)} transformation.
#' Bootstrap resampling is performed on the Fisher z values, and the resulting confidence intervals
#' are back-transformed to the correlation metric.
#'
#' When more than one group is present, pairwise differences in average inter-item correlations are
#' summarized using a simple MOVER-type (Method Of Variance Estimates Recovery) confidence interval.
#' The MOVER approach combines two confidence intervals into a confidence interval for the difference
#' between the corresponding estimates, using only the endpoints of the original intervals.
#'
#' @param data A \code{data.frame} of numeric item responses. Rows are observations, columns are items.
#' @param rmethod Character. Correlation method: \code{"pearson"} (default) or \code{"spearman"}.
#' @param conf.level Numeric. Confidence level for the interval (default = \code{0.95}).
#' @param nboot Integer. Number of bootstrap samples for CI computation (default = \code{1000}).
#' @param bmethod Character. Bootstrap CI method: \code{"bca"} (bias-corrected and accelerated),
#'   \code{"perc"} (percentile), or \code{"norm"} (normal approximation). Default = \code{"bca"}.
#' @param nd Integer. Number of decimal digits to round results (default = \code{3}).
#' @param group Optional vector indicating group membership (factor, character, or numeric).
#'   If supplied, group-wise iia will be calculated and pairwise differences between groups
#'   will be reported.
#'
#' @return
#' A list with two components:
#' \itemize{
#'   \item \code{group_results}: A \code{data.frame} containing, for each group:
#'   \itemize{
#'     \item \code{group}: Group label.
#'     \item \code{avg_r}: Average inter-item correlation (back-transformed from Fisher's z).
#'     \item \code{lwr.ci}, \code{upr.ci}: Lower and upper bounds of the bootstrap confidence interval.
#'     \item \code{min}, \code{max}: Minimum and maximum inter-item correlations in the matrix.
#'     \item \code{sd}: Standard deviation of inter-item correlations.
#'   }
#'
#'   \item \code{comparisons}: A \code{data.frame} containing all pairwise comparisons between groups:
#'   \itemize{
#'     \item \code{group1}, \code{group2}: The groups being compared.
#'     \item \code{diff}: Difference in average inter-item correlations
#'       (\code{avg_r_group1 - avg_r_group2}).
#'     \item \code{lwr.ci}, \code{upr.ci}: MOVER-type confidence interval limits for the difference.
#'   }
#'   If only one group is present (or \code{group = NULL}), \code{comparisons} is returned as an
#'   empty \code{data.frame}.
#' }
#'
#' @references
#' Briggs, S.R. and Cheek, J.M. (1986). The role of factor analysis in the development and evaluation
#' of personality scales. \emph{Journal of Personality, 54}, 106–148.
#' https://doi.org/10.1111/j.1467-6494.1986.tb00391.x
#'
#' Clark, L. A., & Watson, D. (1995). Constructing validity: Basic issues in objective scale development.
#' \emph{Psychological Assessment, 7}(3), 309–319. https://doi.org/10.1037/1040-3590.7.3.309
#'
#' Piedmont, R.L. (2014). Inter-item correlations. In A.C. Michalos (Ed.),
#' \emph{Encyclopedia of Quality of Life and Well-Being Research}. Springer, Dordrecht.
#' https://doi.org/10.1007/978-94-007-0753-5_1493
#'
#' Park, J., van den Broek, K. L., Bhullar, N., Ogunbode, C. A., Schermer, J. A., Doran, R., Ardi, R.,
#' Hanss, D., Maran, D. A., Albzour, M., Aquino, S. D., Ayanian, A. H., Chegeni, R., Chukwuorji,
#' J. B. C., Enea, V., Ghanbarian, E., Ghorayeb, J., Jiang, F., Kehinde, O. A., ... Yadav, R. (2022).
#' Comparison of the inter-item correlations of the Big Five Inventory-10 (BFI-10) between Western and
#' non-Western contexts. \emph{Personality and Individual Differences, 196}, 111751.
#' https://doi.org/10.1016/j.paid.2022.111751
#'
#' @examples
#' set.seed(123)
#'
#' data <- data.frame(
#'   item1 = rnorm(100),
#'   item2 = rnorm(100),
#'   item3 = rnorm(100)
#' )
#'
#' # Single-group average inter-item correlation
#' res1 <- iiacor(data)
#' res1$group_results
#'
#' # Two groups with MOVER comparison
#' grp <- rep(c("A", "B"), each = 50)
#' res2 <- iiacor(data, rmethod = "spearman", group = grp)
#' res2$group_results
#' res2$comparisons
#'
#' @export
iiacor <- function(data,
                   rmethod    = "pearson",
                   conf.level = 0.95,
                   nboot      = 1000,
                   bmethod    = "bca",
                   nd         = 3,
                   group      = NULL) {
  
  # Convert to data.frame if not already
  if (!is.data.frame(data)) data <- as.data.frame(data)
  
  # Argument checks
  if (!rmethod %in% c("pearson", "spearman"))
    stop("Argument 'rmethod' must be 'pearson' or 'spearman'.")
  
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("Argument 'conf.level' must be a number between 0 and 1.")
  
  if (!is.numeric(nboot) || nboot <= 0 || nboot %% 1 != 0)
    stop("Argument 'nboot' must be a positive integer.")
  
  if (!bmethod %in% c("bca", "perc", "norm"))
    stop("Argument 'bmethod' must be one of 'bca', 'perc', or 'norm'.")
  
  if (!is.numeric(nd) || nd < 0 || nd %% 1 != 0)
    stop("Argument 'nd' must be a non-negative integer.")
  
  if (!is.null(group) && length(group) != nrow(data))
    stop("Length of 'group' must match the number of rows in 'data'.")
  
  if (!requireNamespace("boot", quietly = TRUE))
    stop("Package 'boot' is required for bootstrap computation.", call. = FALSE)
  
  # --- HANDLE MISSING VALUES ---
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
    stop("Not enough data to compute inter-item correlations.")
  
  # --- KEEP ONLY NUMERIC COLUMNS ---
  num_cols <- vapply(data, is.numeric, logical(1))
  if (!all(num_cols)) {
    data <- data[, num_cols, drop = FALSE]
    message("Non-numeric columns were removed.")
  }
  
  if (ncol(data) < 2)
    stop("At least two numeric items are required to compute inter-item correlations.")
  
  # ----- AUXILIARY FUNCTION FOR ONE GROUP -----
  compute_metrics <- function(sub_data, group_name = "total") {
    
    if (nrow(sub_data) <= 3 || ncol(sub_data) < 2) {
      stop(paste0("Group '", group_name,
                  "' does not have enough data (need >= 4 rows and >= 2 items)."))
    }
    
    cor_matrix <- stats::cor(sub_data, method = rmethod)
    cor_values <- cor_matrix[lower.tri(cor_matrix)]
    
    fisher_z   <- atanh(cor_values)
    z_bar      <- mean(fisher_z)
    avg_r      <- tanh(z_bar)
    
    cor_sd  <- stats::sd(cor_values)
    cor_min <- min(cor_values)
    cor_max <- max(cor_values)
    
    # Bootstrap on Fisher's z-mean
    boot_fun <- function(d, idx) mean(d[idx])
    boot_obj <- boot::boot(fisher_z, boot_fun, R = nboot)
    
    z_samples <- as.numeric(boot_obj$t)
    r_samples <- tanh(z_samples)
    
    alpha <- 1 - conf.level
    
    if (bmethod == "bca") {
      ci_z <- boot::boot.ci(boot_obj, type = "bca", conf = conf.level)$bca[4:5]
      ci_r <- tanh(ci_z)
    } else if (bmethod == "perc") {
      ci_r <- stats::quantile(r_samples,
                              probs = c(alpha / 2, 1 - alpha / 2),
                              na.rm = TRUE)
    } else {  # "norm"
      se_r  <- stats::sd(r_samples, na.rm = TRUE)
      crit  <- stats::qnorm(1 - alpha / 2)
      ci_r  <- c(avg_r - crit * se_r, avg_r + crit * se_r)
    }
    
    data.frame(
      group  = group_name,
      avg_r  = round(avg_r,   nd),
      lwr.ci = round(ci_r[1], nd),
      upr.ci = round(ci_r[2], nd),
      min    = round(cor_min, nd),
      max    = round(cor_max, nd),
      sd     = round(cor_sd,  nd),
      row.names = NULL
    )
  }
  
  # ---- COMPUTE GROUP-WISE RESULTS ----
  if (is.null(group)) {
    group_results <- compute_metrics(data, group_name = "total")
  } else {
    groups <- unique(group)
    group_results <- do.call(
      rbind,
      lapply(groups, function(g) {
        sub_data <- data[group == g, , drop = FALSE]
        compute_metrics(sub_data, group_name = g)
      })
    )
    rownames(group_results) <- NULL
  }
  
  # ---- MOVER COMPARISONS BETWEEN GROUPS ----
  if (nrow(group_results) > 1) {
    
    # sanity check: one row per group
    if (any(table(group_results$group) != 1)) {
      stop("Internal error: 'group_results' must contain exactly one row per group.")
    }
    
    gnames <- unique(group_results$group)
    combs  <- utils::combn(gnames, 2, simplify = FALSE)
    
    comparisons <- lapply(combs, function(gs) {
      g1 <- gs[1]
      g2 <- gs[2]
      
      g1_dat <- group_results[group_results$group == g1, ]
      g2_dat <- group_results[group_results$group == g2, ]
      
      theta1 <- g1_dat$avg_r
      L1     <- g1_dat$lwr.ci
      U1     <- g1_dat$upr.ci
      
      theta2 <- g2_dat$avg_r
      L2     <- g2_dat$lwr.ci
      U2     <- g2_dat$upr.ci
      
      D <- theta1 - theta2
      
      lower_D <- D - sqrt((theta1 - L1)^2 + (U2 - theta2)^2)
      upper_D <- D + sqrt((U1 - theta1)^2 + (theta2 - L2)^2)
      
      data.frame(
        group1 = g1,
        group2 = g2,
        diff   = round(D,        nd),
        lwr.ci = round(lower_D,  nd),
        upr.ci = round(upper_D,  nd),
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
  
  # ---- RETURN (STYLE APDmg) ----
  list(
    group_results = group_results,
    comparisons   = comparisons
  )
}
