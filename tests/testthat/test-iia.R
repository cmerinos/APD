test_that("iiacor() returns correct structure without groups", {
  set.seed(123)
  dat <- matrix(sample(1:5, 100, replace = TRUE), ncol = 5)

  res <- iiacor(dat)

  # Estructura general
  expect_type(res, "list")
  expect_true(all(c("group_results", "comparisons", "global_test") %in% names(res)))

  # ----- group_results -----
  gr <- res$group_results
  expect_s3_class(gr, "data.frame")
  expect_equal(nrow(gr), 1)
  expect_equal(names(gr), c("group", "avg_r", "lwr.ci", "upr.ci", "min", "max", "sd"))

  expect_equal(gr$group, "total")
  expect_true(is.numeric(gr$avg_r))
  expect_true(is.numeric(gr$lwr.ci))
  expect_true(is.numeric(gr$upr.ci))
  expect_true(is.numeric(gr$min))
  expect_true(is.numeric(gr$max))
  expect_true(is.numeric(gr$sd))

  # Rango razonable para correlaciones
  expect_true(gr$avg_r >= -1 && gr$avg_r <= 1)
  expect_true(gr$lwr.ci >= -1 && gr$lwr.ci <= 1)
  expect_true(gr$upr.ci >= -1 && gr$upr.ci <= 1)

  # Coherencia del IC
  expect_true(gr$lwr.ci <= gr$avg_r)
  expect_true(gr$avg_r <= gr$upr.ci)

  # Descriptivos
  expect_true(gr$min >= -1 && gr$min <= 1)
  expect_true(gr$max >= -1 && gr$max <= 1)
  expect_true(gr$sd >= 0)

  # ----- comparisons -----
  comp <- res$comparisons
  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 0)  # sin grupos → sin comparaciones

  # ----- global_test -----
  gt <- res$global_test
  expect_s3_class(gt, "data.frame")
  expect_equal(nrow(gt), 1)
  expect_equal(names(gt), c("Q", "df", "p.value", "I2", "k", "pooled"))

  # Con un solo grupo: Q, df, p.value, I2 son NA; k = 1
  expect_true(is.na(gt$Q))
  expect_true(is.na(gt$df))
  expect_true(is.na(gt$p.value))
  expect_true(is.na(gt$I2))
  expect_equal(gt$k, 1)

  # pooled debe coincidir (dentro de tolerancia) con avg_r
  expect_equal(gt$pooled, gr$avg_r, tolerance = 1e-12)
})

test_that("iiacor() with groups returns multiple group results and a comparison", {
  set.seed(123)
  dat <- matrix(sample(1:5, 200, replace = TRUE), ncol = 5)
  group <- rep(c("G1", "G2"), each = 20)

  res <- iiacor(dat, group = group)

  expect_type(res, "list")
  expect_true(all(c("group_results", "comparisons", "global_test") %in% names(res)))

  # ----- group_results -----
  gr <- res$group_results
  expect_s3_class(gr, "data.frame")
  expect_equal(nrow(gr), 2)
  expect_true(all(c("group", "avg_r", "lwr.ci", "upr.ci", "min", "max", "sd") %in% names(gr)))
  expect_true(all(gr$group %in% c("G1", "G2")))

  # Correlaciones promedio entre -1 y 1
  expect_true(all(gr$avg_r >= -1 & gr$avg_r <= 1))
  expect_true(all(gr$lwr.ci >= -1 & gr$lwr.ci <= 1))
  expect_true(all(gr$upr.ci >= -1 & gr$upr.ci <= 1))

  # Coherencia de IC
  expect_true(all(gr$lwr.ci <= gr$avg_r))
  expect_true(all(gr$avg_r <= gr$upr.ci))

  # ----- comparisons -----
  comp <- res$comparisons
  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 1)  # una comparación G1 vs G2
  expect_true(all(c("group1", "group2", "diff", "lwr.ci", "upr.ci", "z", "p.value") %in% names(comp)))
  expect_true(all(comp$group1 %in% c("G1", "G2")))
  expect_true(all(comp$group2 %in% c("G1", "G2")))
  expect_true(is.numeric(comp$diff))
  expect_true(is.numeric(comp$lwr.ci))
  expect_true(is.numeric(comp$upr.ci))

  # ----- global_test -----
  gt <- res$global_test
  expect_s3_class(gt, "data.frame")
  expect_equal(nrow(gt), 1)
  expect_equal(names(gt), c("Q", "df", "p.value", "I2", "k", "pooled"))

  # Con dos grupos debería haber k = 2
  expect_equal(gt$k, 2)
  # Q, df y p.value pueden ser numéricos no-NA (dependiendo de tu implementación),
  # así que solo verificamos que sean numéricos o NA válidos.
  expect_true(is.numeric(gt$Q) || is.na(gt$Q))
  expect_true(is.numeric(gt$df) || is.na(gt$df))
  expect_true(is.numeric(gt$p.value) || is.na(gt$p.value))
  expect_true(is.numeric(gt$I2) || is.na(gt$I2))
})

test_that("iiacor() is invariant to column order (no groups)", {
  set.seed(123)
  dat <- matrix(sample(1:5, 100, replace = TRUE), ncol = 5)
  dat_shuffled <- dat[, sample(1:5)]

  res1 <- iiacor(dat)
  res2 <- iiacor(dat_shuffled)

  gr1 <- res1$group_results
  gr2 <- res2$group_results

  # El promedio de correlación debe ser el mismo
  expect_equal(gr1$avg_r, gr2$avg_r, tolerance = 1e-12)
})
