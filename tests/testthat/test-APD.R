test_that("APD return a valid data frame", {

  set.seed(123)

  x <- data.frame(
    i1 = sample(1:5, 100, TRUE),
    i2 = sample(1:5, 100, TRUE),
    i3 = sample(1:5, 100, TRUE)
  )

  res <- APD(
    data = x,
    ncat = 5,
    ci = FALSE,
    conf.level = 0.95,
    B = 100,
    cimethod = "perc"
  )

  expect_s3_class(res, "data.frame")

  expect_named(
    res,
    c("Estimate", "AD", "APD", "lwr.ci", "upr.ci")
  )

  expect_equal(nrow(res), 1)

})


test_that("APD is zero when all responses are identical", {

  x <- data.frame(
    i1 = rep(3, 100),
    i2 = rep(3, 100),
    i3 = rep(3, 100)
  )

  res <- APD(
    data = x,
    ncat = 5,
    ci = FALSE,
    conf.level = 0.95,
    B = 100,
    cimethod = "perc"
  )

  expect_equal(res$APD, 0)

})
