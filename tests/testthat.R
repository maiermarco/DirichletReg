library(testthat)
library(DirichletReg)

old_options <- options("warn" = 1L)

test_check("DirichletReg")

options(old_options)
