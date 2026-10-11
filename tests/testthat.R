library(testthat)     # nolint
library(DirichletReg) # nolint

old_options <- options("warn" = 1L) # nolint

test_check("DirichletReg")

options(old_options) # nolint
