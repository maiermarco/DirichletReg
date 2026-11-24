test_that("weighted models predict on a shorter set of new rows", {
  predictors <- ArcticLake["depth"]
  composition <- as.matrix(ArcticLake[, 1L:3L])
  predictors$Y <- DR_data(composition / rowSums(composition))
  predictors$case_weights <- rep(c(1, 2), length.out = nrow(predictors))
  rows <- c(18L, 3L, 27L)
  newdata <- predictors[rows, "depth", drop = FALSE]

  for(model in c("common", "alternative")){
    fit <- DirichReg(Y ~ depth, data = predictors, model = model,
                     weights = case_weights)
    expect_equal(
      unname(predict(fit, newdata = newdata)),
      unname(fitted(fit)[rows, , drop = FALSE])
    )
  }
})
