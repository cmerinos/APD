#' Average Proportional Distance (APD)
#'
#' @description
#' Average Proportional Distance (APD) is a measure of a test's internal
#' consistency that focuses on the average difference between item scores.
#' To calculate it, you find the absolute difference between every pair of
#' item scores, average those differences, and then divide by the number of
#' response options (minus one) to get a proportional score, which ranges
#' from 0 to 1.
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
#'   \item \code{"bca"} — bias-corrected and accelerated interval.
#'   \item \code{"perc"} — percentile interval.
#'   \item \code{"norm"} — normal-approximation interval.
#' }
#'
#' @section Interpretation:
#' APD quantifies the *average disagreement* between item scores on a 0–1 scale.
#' The following informal guidelines may help interpretation:
#' \itemize{
#'   \item \strong{APD < 0.20}: typically indicates very good internal consistency.
#'   \item \strong{0.20 ≤ APD ≤ 0.25}: often acceptable depending on the construct.
#'   \item \strong{APD > 0.25}: may signal weaker internal consistency or heterogeneous item content.
#' }
#' These values are recommended for Sturman et al. (2009), and they are not strict cutoffs and should
#' be interpreted alongside other reliability evidence (e.g., alpha, omega) and substantive test
#' characteristics.
#'
#' @references
#' Sturman, D., Cribbie, R. A., & Flett, G. L. (2009).
#' The average distance between item values: A novel approach for estimating
#' internal consistency. \emph{Educational and Psychological Measurement},
#' 69(6), 913–932.
#'
#' @param data A \code{data.frame} or matrix containing item responses. Each
#' column represents an item scored on the same categorical scale.
#'
#' @param ncat Integer. Number of response categories for the items (e.g., 5 for a 1–5 scale).
#' Used to convert the average difference into the proportional metric.
#'
#' @param ci Logical. Should a bootstrap confidence interval for APD be computed?
#'
#' @param conf.level Numeric (0–1). Confidence level (e.g., \code{0.95}) when \code{ci = TRUE}.
#'
#' @param B Integer. Number of bootstrap resamples used when computing CIs.
#'
#' @param cimethod Character. Type of confidence interval: \code{"bca"}, \code{"perc"}, or \code{"norm"}.
#'
#' @return
#' A \code{data.frame} with the following rows:
#' \itemize{
#'   \item \code{"Av. diff."}: the average absolute difference \code{AD}.
#'   \item \code{"Av. Prop Diff."}: the proportional difference \code{APD}.
#'   \item \code{"Low CI"}: bootstrap lower confidence bound (if \code{ci = TRUE}).
#'   \item \code{"Upper CI"}: bootstrap upper confidence bound (if \code{ci = TRUE}).
#' }
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- data.frame(matrix(sample(1:5, 200 * 5, replace = TRUE), ncol = 5))
#'
#' # APD without confidence interval
#' APD(x, ncat = 5, ci = FALSE, conf.level = 0.95, B = 500, cimethod = "perc")
#'
#' # APD with 95% percentile CI
#' APD(x, ncat = 5, ci = TRUE, conf.level = 0.95, B = 300, cimethod = "perc")
#'
#' # APD with BCa CI
#' APD(x, ncat = 5, ci = TRUE, conf.level = 0.95, B = 300, cimethod = "bca")
#' }
#'
#' @export
APD <- function(data, ncat, ci, conf.level, B, cimethod) {
  library(boot)
  library(parallel)

  # Generar todas las combinaciones de columnas (pares de items)
  column_combinations <- combn(ncol(data), 2)

  # Inicializar una lista para almacenar las diferencias absolutas entre pares de columnas
  diff_abs_list <- list()

  # Calcular las diferencias absolutas para cada par de columnas usando operaciones vectorizadas
  for (i in 1:ncol(column_combinations)) {
    # Obtener los índices de las columnas a comparar
    col1 <- column_combinations[1, i]
    col2 <- column_combinations[2, i]

    # Calcular la diferencia absoluta entre los valores de las dos columnas de forma vectorizada
    diff_abs_list[[i]] <- abs(data[[col1]] - data[[col2]])
  }

  # Unir todas las diferencias absolutas en un solo vector
  all_diff_abs <- unlist(diff_abs_list)

  # Calcular el promedio de las diferencias absolutas
  AD <- mean(all_diff_abs, na.rm = TRUE)

  # Calcular el APD total (Average Proportional Distance)
  APD <- round(AD / (ncat - 1), 3)

  # Inicializar valores para el límite inferior (LCI) y el límite superior (UCI) del intervalo de confianza
  LCI <- NA
  UCI <- NA

  # Si se desea calcular el intervalo de confianza
  if (ci) {
    # Función para calcular APD con una muestra bootstrap
    boot_apd <- function(data, ncat) {
      # Tomar una muestra bootstrap aleatoria con reemplazo
      sample_data <- data[sample(1:nrow(data), replace = TRUE), ]

      # Recalcular AD y APD para la muestra bootstrap
      column_combinations <- combn(ncol(sample_data), 2)
      diff_abs_list <- list()
      for (i in 1:ncol(column_combinations)) {
        col1 <- column_combinations[1, i]
        col2 <- column_combinations[2, i]
        diff_abs_list[[i]] <- abs(sample_data[[col1]] - sample_data[[col2]])
      }
      all_diff_abs <- unlist(diff_abs_list)
      AD_boot <- mean(all_diff_abs, na.rm = TRUE)
      APD_boot <- AD_boot / (ncat - 1)
      return(APD_boot)
    }

    # Generar las muestras bootstrap y calcular APD para cada una en paralelo
    library(boot)
    boot_samples <- lapply(1:B, function(x) boot_apd(data, ncat))
    boot_samples <- unlist(boot_samples)

    # Calcular el intervalo de confianza de acuerdo al método seleccionado
    if (cimethod == "bca") {
      # Método de corrección de sesgo y aceleración (BCa)
      boot_obj <- boot(data, statistic = function(data, idx) boot_apd(data[idx, ], ncat), R = B)
      # Look in https://stats.stackexchange.com/questions/37918/why-is-the-error-estimated-adjustment-a-is-na-generated-from-r-boot-package
      ci_bca <- boot.ci(boot_obj, type = "bca", conf = conf.level, L =empinf(boot_obj, index=1L, type="jack"))
      LCI <- round(ci_bca$bca[4], 3)  # Límite inferior del IC
      UCI <- round(ci_bca$bca[5], 3)  # Límite superior del IC
    } else if (cimethod == "perc") {
      # Método percentil
      LCI <- round(quantile(boot_samples, probs = (1 - conf.level) / 2), 3)
      UCI <- round(quantile(boot_samples, probs = 1 - (1 - conf.level) / 2), 3)
    } else if (cimethod == "norm") {
      # Método normal basado en la media y la desviación estándar
      se <- sd(boot_samples)  # Desviación estándar de las muestras bootstrap
      mean_boot <- mean(boot_samples)  # Media de las muestras bootstrap
      z <- qnorm((1 + conf.level) / 2)  # Valor crítico z para el nivel de confianza
      LCI <- round(mean_boot - z * se, 3)  # Límite inferior del IC
      UCI <- round(mean_boot + z * se, 3)  # Límite superior del IC
    }
  }

  # Crear el dataframe de resultados
  if (ci) {
    # Si se calculó el intervalo de confianza, incluir LCI y UCI en los resultados
    resultados <- data.frame(
      Parameters = c("Av. diff.", "Av. Prop Diff.", "Low CI", "Upper CI"),
      Value = c(round(AD, 3), APD, LCI, UCI)
    )
  } else {
    # Si no se calculó el intervalo de confianza, solo incluir AD y APD
    resultados <- data.frame(
      Parameters = c("Av. diff.", "Av. Prop Diff."),
      Value = c(round(AD, 3), APD)
    )
  }

  # Retornar el dataframe de resultados
  return(resultados)
}

