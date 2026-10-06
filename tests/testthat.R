library(testthat)     # nolint
library(DirichletReg) # nolint

old_options <- options("warn" = 1L) # nolint

test_check("DirichletReg") # , reporter = default_reporter()

options(old_options) # nolint
