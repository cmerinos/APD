#' Average Proportional Distance for multiple groups
#'
#' @description
#' Computes the Average Proportional Distance (APD) for two or more groups.
#' APD is an index of internal consistency based on the average absolute
#' difference between item scores, scaled to the range \eqn{[0, 1]} by
#' dividing by \code{ncat - 1}.
#'
#' For each group, the function can optionally estimate bootstrap confidence
#' intervals. Pairwise differences in APD between groups are then quantified
#' using a MOVER-type (Method Of Variance Estimates Recovery) confidence
#' interval, which combines the uncertainty of the two APD estimates.
#' This MOVER approach is considered experimental but reasonable for
#' exploratory comparison of APD values across groups.
#'
#' @param data A \code{data.frame} of item responses. Rows are respondents and
#'   columns are items. All columns must be numeric or integer-encoded item
#'   scores with the same number of response categories.
#' @param ncat Integer. Number of response categories for the items. This is
#'   used to scale the average absolute difference to the proportional metric
#'   \eqn{[0, 1]} by dividing by \code{ncat - 1}.
#' @param group A vector (factor, character, or numeric) with the group
#'   membership for each row in \code{data}. Its length must be equal to
#'   \code{nrow(data)}.
#' @param ci Logical. If \code{TRUE} (default), bootstrap confidence intervals
#'   for APD are computed within each group.
#' @param conf.level Confidence level for the bootstrap confidence intervals.
#'   Default is \code{0.95}.
#' @param B Integer. Number of bootstrap resamples used to estimate confidence
#'   intervals within each group. Default is \code{1000}.
#' @param cimethod Character string indicating the bootstrap method passed
#'   to \code{\link[boot]{boot.ci}}. One of \code{"norm"}, \code{"basic"},
#'   \code{"stud"}, \code{"perc"}, or \code{"bca"} (default). The option
#'   \code{"stud"} is not currently implemented for this function.
#'
#'#' @details
#' The function applies a MOVER-type (Method Of Variance Estimates Recovery)
#' approach for constructing confidence intervals for the difference in APD
#' between groups. MOVER methods combine two confidence intervals into a
#' confidence interval for their difference, relying only on the endpoints of
#' the original intervals rather than on asymptotic distributional assumptions.
#' This family of methods has been shown to perform well in related problems,
#' such as confidence intervals for differences in proportions or risk
#' differences (Newcombe, 1998; Zou & Donner, 2008).
#'
#' The implementation used here adapts the simple MOVER formula for independent
#' estimates:
#'
#' \deqn{
#' \mathrm{CI}_{\Delta}
#' =
#' \left[
#' (\hat{\theta}_1 - \hat{\theta}_2)
#' \pm
#' \sqrt{(\hat{\theta}_1 - L_1)^2 + (U_2 - \hat{\theta}_2)^2}
#' \right],
#' }{}
#'
#' where \eqn{\hat{\theta}_g} is the APD estimate for group \eqn{g}, and
#' \eqn{[L_g, U_g]} is its confidence interval.
#'
#' This MOVER-based interval is considered experimental for APD because, unlike
#' proportions or risk measures, APD has no established sampling distribution or
#' variance expressions for combining confidence intervals. Nonetheless, MOVER
#' provides a transparent and distribution-free approximation that incorporates
#' the uncertainty from both groups' APD estimates. The resulting intervals are
#' intended for exploratory interpretation rather than formal inference.
#'
#' @return
#' A list with two elements:
#' \itemize{
#'   \item \code{group_results}: A \code{data.frame} containing, for each group:
#'     \itemize{
#'       \item \code{group}: Group label.
#'       \item \code{APD}: Estimated Average Proportional Distance.
#'       \item \code{LCI}: Lower confidence limit (or \code{NA} if \code{ci = FALSE}).
#'       \item \code{UCI}: Upper confidence limit (or \code{NA} if \code{ci = FALSE}).
#'     }
#'
#'   \item \code{comparisons}: A \code{data.frame} containing all pairwise
#'   comparisons between groups:
#'     \itemize{
#'       \item \code{group1}, \code{group2}: The groups compared.
#'       \item \code{Difference}: APD difference (\code{APD_group1 - APD_group2}).
#'       \item \code{Lower_CI}, \code{Upper_CI}: MOVER-based confidence interval
#'         for the difference. Empty if only one group is present.
#'     }
#' }
#'
#' @references
#'
#' Newcombe, R. G. (1998). Interval estimation for the difference between
#' independent proportions: comparison of eleven methods.
#' \emph{Statistics in Medicine, 17}, 873–890.
#' \doi{10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I}
#'
#' Sturman, D., Cribbie, R. A., & Flett, G. L. (2009).
#' The average distance between item values: A novel approach for estimating
#' internal consistency. \emph{Educational and Psychological Measurement},
#' 69(6), 913–932.
#'
#' Zou, G. Y., & Donner, A. (2008). Construction of confidence limits
#' about effect measures: A general approach.
#' \emph{Statistics in Medicine, 27}, 1693–1702.
#' \doi{10.1002/sim.3095}
#'
#' @examples
#' set.seed(123)
#'
#' # Simulated data: 30 respondents, 5 items, 2 groups
#' dat <- data.frame(
#'   item1 = sample(1:5, 30, replace = TRUE),
#'   item2 = sample(1:5, 30, replace = TRUE),
#'   item3 = sample(1:5, 30, replace = TRUE),
#'   item4 = sample(1:5, 30, replace = TRUE),
#'   item5 = sample(1:5, 30, replace = TRUE)
#' )
#'
#' grp <- rep(c("Group_A", "Group_B"), each = 15)
#'
#' res <- APDmg(data = dat, ncat = 5, group = grp,
#'              ci = TRUE, conf.level = 0.95, B = 200,
#'              cimethod = "bca")
#'
#' res$APD.group
#' res$comparisons
#'
#' @export
APDmg <- function(data, ncat, group, ci = TRUE, conf.level = 0.95, B = 1000, cimethod = "bca") {
  library(boot)

  # Validaciones
  if (!is.data.frame(data)) stop("El argumento 'data' debe ser un data.frame.")
  if (missing(group)) stop("El argumento 'group' es obligatorio para calcular APD por grupos.")
  if (length(group) != nrow(data)) stop("El vector 'group' debe tener la misma longitud que el número de filas de 'data'.")
  if (!cimethod %in% c("norm", "basic", "stud", "perc", "bca")) {
    stop("El método seleccionado no es válido. Usa uno de: 'norm', 'basic', 'stud', 'perc', 'bca'.")
  }

  # Función auxiliar para calcular APD para un grupo específico
  compute_apd <- function(sub_data, ncat) {
    column_combinations <- combn(ncol(sub_data), 2)
    diff_abs_list <- list()
    for (i in 1:ncol(column_combinations)) {
      col1 <- column_combinations[1, i]
      col2 <- column_combinations[2, i]
      diff_abs_list[[i]] <- abs(sub_data[[col1]] - sub_data[[col2]])
    }
    all_diff_abs <- unlist(diff_abs_list)
    AD <- mean(all_diff_abs, na.rm = TRUE)
    APD <- round(AD / (ncat - 1), 3)
    return(APD)
  }

  # Bootstrap para intervalos de confianza
  bootstrap_apd <- function(data, ncat) {
    boot_apd <- function(data, idx) {
      sample_data <- data[idx, , drop = FALSE]
      compute_apd(sample_data, ncat)
    }
    boot_obj <- boot(data, statistic = boot_apd, R = B)
    return(boot_obj)
  }

  # Calcular APD y CI por grupo
  unique_groups <- unique(group)
  group_results <- lapply(unique_groups, function(g) {
    sub_data <- data[group == g, , drop = FALSE]
    APD <- compute_apd(sub_data, ncat)

    # Intervalos de confianza
    LCI <- NA
    UCI <- NA
    if (ci) {
      boot_obj <- bootstrap_apd(sub_data, ncat)
      ci_obj <- boot.ci(boot_obj, type = cimethod, conf = conf.level)
      if (cimethod == "bca") {
        LCI <- round(ci_obj$bca[4], 3)
        UCI <- round(ci_obj$bca[5], 3)
      } else if (cimethod == "perc") {
        LCI <- round(ci_obj$percent[4], 3)
        UCI <- round(ci_obj$percent[5], 3)
      } else if (cimethod == "norm") {
        LCI <- round(ci_obj$normal[2], 3)
        UCI <- round(ci_obj$normal[3], 3)
      } else if (cimethod == "basic") {
        LCI <- round(ci_obj$basic[4], 3)
        UCI <- round(ci_obj$basic[5], 3)
      } else if (cimethod == "stud") {
        stop("El método 'stud' requiere estadísticas de error adicionales. Por favor, usa otro método.")
      }
    }

    return(data.frame(group = g, APD = APD, LCI = LCI, UCI = UCI))
  })

  # Combinar los resultados en un data.frame
  group_results <- do.call(rbind, group_results)

  # Comparaciones entre grupos
  if (nrow(group_results) > 1) {
    combinations <- combn(group_results$group, 2, simplify = FALSE)
    comparisons <- lapply(combinations, function(groups) {
      group1 <- groups[1]
      group2 <- groups[2]

      APD1 <- group_results$APD[group_results$group == group1]
      APD2 <- group_results$APD[group_results$group == group2]
      LCI1 <- group_results$LCI[group_results$group == group1]
      UCI1 <- group_results$UCI[group_results$group == group1]
      LCI2 <- group_results$LCI[group_results$group == group2]
      UCI2 <- group_results$UCI[group_results$group == group2]

      # Diferencia entre APDs
      difference <- APD1 - APD2

      # Método MOVER para intervalos de confianza
      lower_diff <- difference - sqrt((APD1 - LCI1)^2 + (UCI2 - APD2)^2)
      upper_diff <- difference + sqrt((UCI1 - APD1)^2 + (APD2 - LCI2)^2)

      data.frame(
        group1 = group1,
        group2 = group2,
        Difference = round(difference, 3),
        lwr.ci = round(lower_diff, 3),
        upr.ci = round(upper_diff, 3)
      )
    })
    comparisons <- do.call(rbind, comparisons)
  } else {
    comparisons <- data.frame(
      group1 = character(0),
      group2 = character(0),
      Difference = numeric(0),
      lwr.ci = numeric(0),
      upr.ci = numeric(0)
    )
  }

  # Retornar resultados
  return(list(
    APD.group = group_results,
    Comparisons = comparisons
  ))
}
