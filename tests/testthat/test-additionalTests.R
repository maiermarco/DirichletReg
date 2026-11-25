# library(DirichletReg); library(testthat); .sce()
if(as.numeric_version(packageDescription("DirichletReg")[["Version"]]) >= as.numeric_version("0.8-0")){

test_that("Prepare Data for Tests", {
  AL   <<- ArcticLake["depth"]
  expect_warning(AL$Y <<- DR_data(ArcticLake[, 1:3]))
})



# https://github.com/maiermarco/DirichletReg/pull/16
# Test if formula with ~ . works for the common param.:
test_that("Dot-Formula: Common Param.", {
  mod_dot       <- DirichReg(Y ~ ., data = AL)
  mod_explicit1 <- DirichReg(Y ~ depth, data = AL)
  mod_explicit2 <- DirichReg(Y ~ depth | depth | depth, data = AL)
  
  summary_dot       <- summary(mod_dot)
  summary_explicit1 <- summary(mod_explicit1)
  summary_explicit2 <- summary(mod_explicit2)
  
  expect_identical(attr(terms(summary_dot      , type = "full"), "term.labels"), "depth")
  expect_identical(attr(terms(summary_explicit1, type = "full"), "term.labels"), "depth")
  expect_identical(attr(terms(summary_explicit2, type = "full"), "term.labels"), "depth")
  
  for(this_type in c("full", "separate")){
    expect_identical(terms(mod_dot      , type = this_type), terms(summary_dot      , type = this_type))
    expect_identical(terms(mod_explicit1, type = this_type), terms(summary_explicit1, type = this_type))
    expect_identical(terms(mod_explicit2, type = this_type), terms(summary_explicit2, type = this_type))
  }
  
  expect_equal(summary_explicit1$coef.mat, summary_dot$coef.mat)
  expect_equal(summary_explicit2$coef.mat, summary_dot$coef.mat)
})



### TODO include tests for alt. param. and combinations (e.g., y ~ . | 1 etc.)



# https://github.com/maiermarco/DirichletReg/pull/15
# Discard weights when using predict.DirichletRegModel()

AL$case_weights <- rep(c(1.0, 2.0), length.out = nrow(AL))
rows <- c(18L, 3L, 27L, 1L)
new_data <- AL[rows, "depth", drop = FALSE]

test_that("weighted models predict on a shorter set of new rows", {
  for(model in c("common", "alternative")){
    fit <- DirichReg(Y ~ depth, data = AL, model = model, weights = case_weights)
    expect_equal(
      unname(predict(fit, newdata = new_data)),
      unname(fitted(fit)[rows, , drop = FALSE])
    )
  }
})

AL$case_weights <- NULL

}
