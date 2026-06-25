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
