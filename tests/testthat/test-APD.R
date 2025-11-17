test_that("APD() returns the correct structure", {
  set.seed(123)
  dat <- matrix(sample(1:5, 50, replace = TRUE), ncol = 5)

  res <- APD(dat, ncat = 5, ci = FALSE)

  # Debe ser un data.frame
  expect_s3_class(res, "data.frame")

  # Columnas correctas
  expect_equal(names(res), c("Parameters", "Value"))

  # Filas correctas
  expect_equal(res$Parameters, c("Av. diff.", "Av. Prop Diff."))

  # Valores numéricos
  expect_true(is.numeric(res$Value))
  expect_equal(length(res$Value), 2)
})

test_that("APD() is zero when all responses are identical", {
  dat_const <- matrix(3, nrow = 10, ncol = 5)

  res <- APD(dat_const, ncat = 5, ci = FALSE)

  expect_equal(res$Value, c(0, 0))
})

test_that("APD() is invariant to column order", {
  set.seed(123)
  dat <- matrix(sample(1:5, 50, replace = TRUE), ncol = 5)
  dat_shuffled <- dat[, sample(1:5)]

  res1 <- APD(dat, ncat = 5, ci = FALSE)
  res2 <- APD(dat_shuffled, ncat = 5, ci = FALSE)

  expect_equal(res1$Value, res2$Value, tolerance = 1e-12)
})
