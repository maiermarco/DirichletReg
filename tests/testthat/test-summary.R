test_that("summary expands dots using the fitted data", {
  predictors <- ArcticLake[, "depth", drop = FALSE]
  proportions <- ArcticLake[, 1L:3L] / rowSums(ArcticLake[, 1L:3L])
  predictors$Y <- DR_data(proportions)

  dot <- summary(DirichReg(Y ~ ., data = predictors))
  explicit <- summary(DirichReg(Y ~ depth, data = predictors))

  expect_identical(attr(dot$terms, "term.labels"), "depth")
  expect_equal(dot$coef.mat, explicit$coef.mat)
})
