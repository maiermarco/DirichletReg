# some checks require DirichletReg >= 0.8-0
._DR_ver_0_8 <- as.numeric_version(packageDescription("DirichletReg")[["Version"]]) >= as.numeric_version("0.8-0")

skip_ext_cran <- function(){
  if(isFALSE(as.logical(Sys.getenv("NOT_CRAN", "true")))) testthat::skip("Skipping Extended Tests on CRAN")
}

#x <- function(ooo) writeClipboard(paste(collapse = "", deparse(ooo, width.cutoff = 500L, control = "exact"))) # TODO
#x(bd1) # TODO

################################################################################
### Basic Checks ###############################################################
################################################################################

compodat <- seq(0.1, 0.9, 0.1)

if(._DR_ver_0_8) test_that("Checks that should result in errors", {
  skip_ext_cran()
  expect_error(DR_data(compodat, trafo = FALSE), regexp = "\"trafo\" must be a small number > 0 or TRUE.", fixed = TRUE)
  expect_error(DR_data(compodat, trafo =   0.2), regexp = "\"trafo\" must be a small number > 0 or TRUE.", fixed = TRUE)
  expect_error(DR_data(compodat, trafo =  -1.0), regexp = "\"trafo\" must be a small number > 0 or TRUE.", fixed = TRUE)
  
  expect_error(DR_data(compodat, norm_tol = TRUE), regexp = "\"norm_tol\" must be a small number > 0.", fixed = TRUE)
  expect_error(DR_data(compodat, norm_tol =  0.2), regexp = "\"norm_tol\" must be a small number > 0.", fixed = TRUE)
  expect_error(DR_data(compodat, norm_tol = -1.0), regexp = "\"norm_tol\" must be a small number > 0.", fixed = TRUE)
})



################################################################################
### Test Vector Input ##########################################################
################################################################################

beta_data     <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_1NA <- c(NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_0   <- c(0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_1   <- c(1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_m1  <- c(-1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_p2  <- c(+2.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

bd1_check <- structure(c(0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data", "beta_data")), Y.original = if(._DR_ver_0_8) structure(c(0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data", "beta_data"))) else structure(list(Y.original = c(0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = c("1 - beta_data", "beta_data"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
bd2_check <- structure(c(NA, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_1NA", "beta_data_1NA")), Y.original = structure(c(NA, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_1NA", "beta_data_1NA"))), dims = 2L, dim.names = c("1 - beta_data_1NA", "beta_data_1NA"), obs = 9L, valid_obs = 8L, normalized = FALSE, transformed = FALSE, base = 1L, class = "DirichletRegData") # nolint
bd3_check <- structure(c(0x1.e38e38e38e38ep-1, 0x1.8888888888889p-1, 0x1.5b05b05b05b05p-1, 0x1.2d82d82d82d83p-1, 0x1p-1, 0x1.a4fa4fa4fa4fbp-2, 0x1.49f49f49f49f5p-2, 0x1.ddddddddddddcp-3, 0x1.27d27d27d27d2p-3, 0x1.c71c71c71c71cp-5, 0x1.ddddddddddddep-3, 0x1.49f49f49f49f4p-2, 0x1.a4fa4fa4fa4fbp-2, 0x1p-1, 0x1.2d82d82d82d83p-1, 0x1.5b05b05b05b05p-1, 0x1.8888888888889p-1, 0x1.b60b60b60b60cp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_0", "beta_data_0")), Y.original = if(._DR_ver_0_8) structure(c(0x1p+0, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x0p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_0", "beta_data_0"))) else structure(list(Y.original = c(0x0p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = c("1 - beta_data_0", "beta_data_0"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
bd4_check <- structure(c(0x1.c71c71c71c71cp-5, 0x1.8888888888889p-1, 0x1.5b05b05b05b05p-1, 0x1.2d82d82d82d83p-1, 0x1p-1, 0x1.a4fa4fa4fa4fbp-2, 0x1.49f49f49f49f5p-2, 0x1.ddddddddddddcp-3, 0x1.27d27d27d27d2p-3, 0x1.e38e38e38e38ep-1, 0x1.ddddddddddddep-3, 0x1.49f49f49f49f4p-2, 0x1.a4fa4fa4fa4fbp-2, 0x1p-1, 0x1.2d82d82d82d83p-1, 0x1.5b05b05b05b05p-1, 0x1.8888888888889p-1, 0x1.b60b60b60b60cp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_1", "beta_data_1")), Y.original = if(._DR_ver_0_8) structure(c(0x0p+0, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_1", "beta_data_1"))) else structure(list(Y.original = c(0x1p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = c("1 - beta_data_1", "beta_data_1"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint

test_that("Testing Well-Behaved Beta-Distributed Data", {
  skip_ext_cran()
  expect_message(bd1 <- DR_data(beta_data), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_identical(bd1, bd1_check)
  expect_s3_class(bd1, "DirichletRegData")
  expect_no_error(print(bd1)) |> capture.output() # capture.output is used, so print methods do not clutter the output
  expect_no_error(print(summary(bd1))) |> capture.output()
})

if(._DR_ver_0_8){ # versions before 0.8-0 could not handle this scenario
  test_that("Testing a Beta-Distributed Data with Missings", {
    skip_ext_cran()
    expect_message(bd2 <- DR_data(beta_data_1NA), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
    expect_identical(bd2, bd2_check)
    expect_s3_class(bd2, "DirichletRegData")
    expect_no_error(print(bd2)) |> capture.output()
    expect_no_error(summary(bd2)) |> capture.output()
  })
}

test_that("Testing a Beta-Distributed Data with a Zero", {
  skip_ext_cran()
  expect_warning(bd3 <- DR_data(beta_data_0), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_identical(bd3, bd3_check)
  expect_s3_class(bd3, "DirichletRegData")
  expect_no_error(print(bd3)) |> capture.output()
  expect_no_error(summary(bd3)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data with a One", {
  skip_ext_cran()
  expect_warning(bd4 <- DR_data(beta_data_1), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_identical(bd4, bd4_check)
  expect_s3_class(bd4, "DirichletRegData")
  expect_no_error(print(bd4)) |> capture.output()
  expect_no_error(summary(bd4)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data LESS THAN ZERO or GREATER THAN ONE", {
  skip_ext_cran()
  expect_error(DR_data(beta_data_m1), regexp = if(._DR_ver_0_8) "beta distribution cannot safely be assumed" else '"Y" contains values < 0', fixed = TRUE)
  expect_error(DR_data(beta_data_p2), regexp = "beta distribution cannot safely be assumed", fixed = TRUE)
})

rm(list = c("bd1_check", "bd2_check", "bd3_check", "bd4_check"))



################################################################################
### Test Single-Column Matrix Input ############################################
################################################################################

beta_matrix     <- cbind(beta_data    )
beta_matrix_1NA <- cbind(beta_data_1NA)
beta_matrix_0   <- cbind(beta_data_0  )
beta_matrix_1   <- cbind(beta_data_1  )
beta_matrix_m1  <- cbind(beta_data_m1 )
beta_matrix_p2  <- cbind(beta_data_p2 )

# component naming changes in 0.8-0: variable names used if present
mbd1_check <- structure(c(0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_matrix", "beta_matrix")), Y.original = if(._DR_ver_0_8) structure(c(0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_matrix", "beta_matrix"))) else structure(list(beta_data = c(0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_matrix", "beta_matrix"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
mbd2_check <- structure(c(NA, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_matrix_1NA", "beta_matrix_1NA")), Y.original = if(._DR_ver_0_8) structure(c(NA, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_matrix_1NA", "beta_matrix_1NA"))) else structure(list(beta_data_1NA = c(NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_matrix_1NA", "beta_matrix_1NA"), obs = 9L, valid_obs = 8L, normalized = FALSE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
mbd3_check <- structure(c(0x1.e38e38e38e38ep-1, 0x1.8888888888889p-1, 0x1.5b05b05b05b05p-1, 0x1.2d82d82d82d83p-1, 0x1p-1, 0x1.a4fa4fa4fa4fbp-2, 0x1.49f49f49f49f5p-2, 0x1.ddddddddddddcp-3, 0x1.27d27d27d27d2p-3, 0x1.c71c71c71c71cp-5, 0x1.ddddddddddddep-3, 0x1.49f49f49f49f4p-2, 0x1.a4fa4fa4fa4fbp-2, 0x1p-1, 0x1.2d82d82d82d83p-1, 0x1.5b05b05b05b05p-1, 0x1.8888888888889p-1, 0x1.b60b60b60b60cp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_matrix_0", "beta_matrix_0")), Y.original = if(._DR_ver_0_8) structure(c(0x1p+0, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x0p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_matrix_0", "beta_matrix_0"))) else structure(list(beta_data_0 = c(0x0p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_matrix_0", "beta_matrix_0"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
mbd4_check <- structure(c(0x1.c71c71c71c71cp-5, 0x1.8888888888889p-1, 0x1.5b05b05b05b05p-1, 0x1.2d82d82d82d83p-1, 0x1p-1, 0x1.a4fa4fa4fa4fbp-2, 0x1.49f49f49f49f5p-2, 0x1.ddddddddddddcp-3, 0x1.27d27d27d27d2p-3, 0x1.e38e38e38e38ep-1, 0x1.ddddddddddddep-3, 0x1.49f49f49f49f4p-2, 0x1.a4fa4fa4fa4fbp-2, 0x1p-1, 0x1.2d82d82d82d83p-1, 0x1.5b05b05b05b05p-1, 0x1.8888888888889p-1, 0x1.b60b60b60b60cp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_matrix_1", "beta_matrix_1")), Y.original = if(._DR_ver_0_8) structure(c(0x0p+0, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_matrix_1", "beta_matrix_1"))) else structure(list(beta_data_1 = c(0x1p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_matrix_1", "beta_matrix_1"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint

test_that("Testing Well-Behaved Beta-Distributed Data (as a Matrix)", {
  skip_ext_cran()
  expect_message(mbd1 <- DR_data(beta_matrix), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_identical(mbd1, mbd1_check)
  expect_s3_class(mbd1, "DirichletRegData")
  expect_no_error(print(mbd1)) |> capture.output()
  expect_no_error(summary(mbd1)) |> capture.output()
})

# works for versions before 0.8-0 too
test_that("Testing a Beta-Distributed Data with Missings (as a Matrix)", {
  skip_ext_cran()
  expect_message(mbd2 <- DR_data(beta_matrix_1NA), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_identical(mbd2, mbd2_check)
  expect_s3_class(mbd2, "DirichletRegData")
  expect_no_error(print(mbd2)) |> capture.output()
  expect_no_error(summary(mbd2)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data with a Zero (as a Matrix)", {
  skip_ext_cran()
  expect_warning(mbd3 <- DR_data(beta_matrix_0), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_identical(mbd3, mbd3_check)
  expect_s3_class(mbd3, "DirichletRegData")
  expect_no_error(print(mbd3)) |> capture.output()
  expect_no_error(summary(mbd3)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data with a One (as a Matrix)", {
  skip_ext_cran()
  expect_warning(mbd4 <- DR_data(beta_matrix_1), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_identical(mbd4, mbd4_check)
  expect_s3_class(mbd4, "DirichletRegData")
  expect_no_error(print(mbd4)) |> capture.output()
  expect_no_error(summary(mbd4)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data LESS THAN ZERO or GREATER THAN ONE (as a Matrix)", {
  skip_ext_cran()
  expect_error(DR_data(beta_matrix_m1), regexp = if(._DR_ver_0_8) "beta distribution cannot safely be assumed" else '"Y" contains values < 0', fixed = TRUE)
  expect_error(DR_data(beta_matrix_p2), regexp = "beta distribution cannot safely be assumed", fixed = TRUE)
})

rm(list = c("mbd1_check", "mbd2_check", "mbd3_check", "mbd4_check",
  "beta_matrix", "beta_matrix_0", "beta_matrix_1", "beta_matrix_1NA", "beta_matrix_m1", "beta_matrix_p2"))



################################################################################
### Test Single-Column Data Frame Input ########################################
################################################################################

beta_df     <- data.frame(beta_data    )
beta_df_1NA <- data.frame(beta_data_1NA)
beta_df_0   <- data.frame(beta_data_0  )
beta_df_1   <- data.frame(beta_data_1  )
beta_df_m1  <- data.frame(beta_data_m1 )
beta_df_p2  <- data.frame(beta_data_p2 )

# component naming changes in 0.8-0: variable names used if present
dbd1_check <- structure(c(0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_df", "beta_df")), Y.original = if(._DR_ver_0_8) structure(c(0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_df", "beta_df"))) else structure(list(beta_data = c(0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_df", "beta_df"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
dbd2_check <- structure(c(NA, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_df_1NA", "beta_df_1NA")), Y.original = if(._DR_ver_0_8) structure(c(NA, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_df_1NA", "beta_df_1NA"))) else structure(list(beta_data_1NA = c(NA, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_df_1NA", "beta_df_1NA"), obs = 9L, valid_obs = 8L, normalized = FALSE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
dbd3_check <- structure(c(0x1.e38e38e38e38ep-1, 0x1.8888888888889p-1, 0x1.5b05b05b05b05p-1, 0x1.2d82d82d82d83p-1, 0x1p-1, 0x1.a4fa4fa4fa4fbp-2, 0x1.49f49f49f49f5p-2, 0x1.ddddddddddddcp-3, 0x1.27d27d27d27d2p-3, 0x1.c71c71c71c71cp-5, 0x1.ddddddddddddep-3, 0x1.49f49f49f49f4p-2, 0x1.a4fa4fa4fa4fbp-2, 0x1p-1, 0x1.2d82d82d82d83p-1, 0x1.5b05b05b05b05p-1, 0x1.8888888888889p-1, 0x1.b60b60b60b60cp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_df_0", "beta_df_0")), Y.original = if(._DR_ver_0_8) structure(c(0x1p+0, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x0p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_df_0", "beta_df_0"))) else structure(list(beta_data_0 = c(0x0p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_df_0", "beta_df_0"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
dbd4_check <- structure(c(0x1.c71c71c71c71cp-5, 0x1.8888888888889p-1, 0x1.5b05b05b05b05p-1, 0x1.2d82d82d82d83p-1, 0x1p-1, 0x1.a4fa4fa4fa4fbp-2, 0x1.49f49f49f49f5p-2, 0x1.ddddddddddddcp-3, 0x1.27d27d27d27d2p-3, 0x1.e38e38e38e38ep-1, 0x1.ddddddddddddep-3, 0x1.49f49f49f49f4p-2, 0x1.a4fa4fa4fa4fbp-2, 0x1p-1, 0x1.2d82d82d82d83p-1, 0x1.5b05b05b05b05p-1, 0x1.8888888888889p-1, 0x1.b60b60b60b60cp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_df_1", "beta_df_1")), Y.original = if(._DR_ver_0_8) structure(c(0x0p+0, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333334p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_df_1", "beta_df_1"))) else structure(list(beta_data_1 = c(0x1p+0, 0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_df_1", "beta_df_1"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint

test_that("Testing Well-Behaved Beta-Distributed Data (as a Data Frame)", {
  skip_ext_cran()
  expect_message(dbd1 <- DR_data(beta_df), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_identical(dbd1, dbd1_check)
  expect_s3_class(dbd1, "DirichletRegData")
  expect_no_error(print(dbd1)) |> capture.output()
  expect_no_error(summary(dbd1)) |> capture.output()
})

# works for versions before 0.8-0 too
test_that("Testing a Beta-Distributed Data with Missings (as a Data Frame)", {
  skip_ext_cran()
  expect_message(dbd2 <- DR_data(beta_df_1NA), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_identical(dbd2, dbd2_check)
  expect_s3_class(dbd2, "DirichletRegData")
  expect_no_error(print(dbd2)) |> capture.output()
  expect_no_error(summary(dbd2)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data with a Zero (as a Data Frame)", {
  skip_ext_cran()
  expect_warning(dbd3 <- DR_data(beta_df_0), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_identical(dbd3, dbd3_check)
  expect_s3_class(dbd3, "DirichletRegData")
  expect_no_error(print(dbd3)) |> capture.output()
  expect_no_error(summary(dbd3)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data with a One (as a Data Frame)", {
  skip_ext_cran()
  expect_warning(dbd4 <- DR_data(beta_df_1), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_identical(dbd4, dbd4_check)
  expect_s3_class(dbd4, "DirichletRegData")
  expect_no_error(print(dbd4)) |> capture.output()
  expect_no_error(summary(dbd4)) |> capture.output()
})

test_that("Testing a Beta-Distributed Data LESS THAN ZERO or GREATER THAN ONE (as a Data Frame)", {
  skip_ext_cran()
  expect_error(DR_data(beta_df_m1), regexp = if(._DR_ver_0_8) "beta distribution cannot safely be assumed" else '"Y" contains values < 0', fixed = TRUE)
  expect_error(DR_data(beta_df_p2), regexp = "beta distribution cannot safely be assumed", fixed = TRUE)
})

rm(list = c("dbd1_check", "dbd2_check", "dbd3_check", "dbd4_check",
  "beta_df", "beta_df_0", "beta_df_1", "beta_df_1NA", "beta_df_m1", "beta_df_p2"))

rm(list = c("beta_data", "beta_data_1NA", "beta_data_0", "beta_data_1", "beta_data_m1", "beta_data_p2"))










DR_AL_check <- structure(c(0x1.8cccccccccccdp-1, 0x1.7020c49ba5e36p-1, 0x1.0395810624dd3p-1, 0x1.0c1175db89e2cp-1, 0x1.6666666666666p-1, 0x1.547ae147ae148p-1, 0x1.b95810624dd3p-2, 0x1.116872b020c4ap-1, 0x1.3d70a3d70a3d7p-3, 0x1.449ba5e353f7dp-2, 0x1.50624dd2f1aap-1, 0x1.6872b020c49bbp-1, 0x1.645a1cac08312p-3, 0x1.b22d0e5604189p-4, 0x1.872b020c49ba6p-2, 0x1.ba5e353f7cedap-4, 0x1.78d4fdf3b645ap-3, 0x1.78d4fdf3b645ap-5, 0x1.3f7ced916872bp-3, 0x1.46a7ef9db22d1p-2, 0x1.851eb851eb852p-4, 0x1.5e353f7ced917p-3, 0x1.ae147ae147ae1p-4, 0x1.87427bcc092b8p-5, 0x1.a9fbe76c8b43ap-6, 0x1.d2f1a9fbe76c9p-4, 0x1.126e978d4fdf4p-4, 0x1.1a9fbe76c8b44p-4, 0x1.47ae147ae147bp-5, 0x1.2f684bda12f69p-4, 0x1.89374bc6a7efap-5, 0x1.70a3d70a3d70ap-5, 0x1.0e5604189374cp-4, 0x1.12b4eabb112b6p-4, 0x1.2f684bda12f69p-4, 0x1.eb851eb851eb8p-5, 0x1.020c49ba5e355p-4, 0x1.999999999999ap-6, 0x1.47ae147ae147bp-6, 0x1.8f5c28f5c28f6p-3, 0x1.fdf3b645a1cacp-3, 0x1.71a9fbe76c8b4p-2, 0x1.a41383b921a3p-2, 0x1.0f5c28f5c28f6p-2, 0x1.49ba5e353f7cfp-2, 0x1.1b22d0e560419p-1, 0x1.78d4fdf3b645ap-2, 0x1.16872b020c49cp-1, 0x1.a8f5c28f5c28fp-2, 0x1.1cac083126e98p-2, 0x1.28f5c28f5c28fp-2, 0x1.126e978d4fdf4p-1, 0x1.65604189374bcp-1, 0x1.b95810624dd2fp-2, 0x1.0dd2f1a9fbe77p-1, 0x1.0395810624dd3p-1, 0x1.e5604189374bcp-2, 0x1.020c49ba5e354p-1, 0x1.cdd2f1a9fbe77p-2, 0x1.11eb851eb851fp-1, 0x1.eb851eb851eb8p-2, 0x1.1ba5e353f7cedp-1, 0x1.16abaf81a7dd6p-1, 0x1.ced916872b021p-2, 0x1.0dd2f1a9fbe77p-1, 0x1.e04189374bc6ap-2, 0x1.fced916872b03p-2, 0x1.cbc6a7ef9db23p-2, 0x1.0874da5490875p-1, 0x1.fae147ae147aep-2, 0x1.f0a3d70a3d70ap-2, 0x1.0ac083126e979p-1, 0x1.e4d63af05e4d7p-2, 0x1.d36952421d36ap-2, 0x1.f4bc6a7ef9db2p-2, 0x1.1374bc6a7ef9ep-1, 0x1.eb851eb851eb8p-2, 0x1.e978d4fdf3b64p-2, 0x1.eb851eb851eb8p-6, 0x1.0624dd2f1a9fcp-5, 0x1.0e5604189374cp-3, 0x1.0f26423f2a5ep-4, 0x1.1eb851eb851ecp-5, 0x1.a9fbe76c8b43ap-7, 0x1.0624dd2f1a9fdp-6, 0x1.916872b020c4ap-4, 0x1.34395810624ddp-2, 0x1.126e978d4fdf4p-2, 0x1.0a3d70a3d70a4p-4, 0x1.89374bc6a7efap-8, 0x1.28f5c28f5c28fp-2, 0x1.916872b020c4ap-3, 0x1.7ef9db22d0e56p-3, 0x1.75c28f5c28f5cp-2, 0x1.3c6a7ef9db22dp-2, 0x1.eb851eb851eb8p-2, 0x1.5c28f5c28f5c3p-2, 0x1.d70a3d70a3d71p-3, 0x1.7ae147ae147aep-2, 0x1.65604189374bcp-2, 0x1.5d2f1a9fbe76dp-2, 0x1.a1c051832f1fcp-2, 0x1.0b4395810624ep-1, 0x1.6f9db22d0e56p-2, 0x1.db22d0e560418p-2, 0x1.bc6a7ef9db22dp-2, 0x1.05a1cac083127p-1, 0x1.a33c38605a33cp-2, 0x1.d3f7ced916873p-2, 0x1.e147ae147ae14p-2, 0x1.a6e978d4fdf3bp-2, 0x1.d67c8a60dd67dp-2, 0x1.e0bc9ac75e0bcp-2, 0x1.cdd2f1a9fbe77p-2, 0x1.989374bc6a7fp-2, 0x1.fae147ae147aep-2, 0x1.010624dd2f1aap-1), dim = c(39L, 3L), dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"), c("sand", "silt", "clay")), Y.original = if(._DR_ver_0_8) structure(c(0x1.8cccccccccccdp-1, 0x1.7020c49ba5e36p-1, 0x1.0395810624dd3p-1, 0x1.0b4395810624ep-1, 0x1.6666666666666p-1, 0x1.547ae147ae148p-1, 0x1.b95810624dd2fp-2, 0x1.116872b020c4ap-1, 0x1.3d70a3d70a3d7p-3, 0x1.449ba5e353f7dp-2, 0x1.50624dd2f1aap-1, 0x1.6872b020c49bbp-1, 0x1.645a1cac08312p-3, 0x1.b22d0e5604189p-4, 0x1.872b020c49ba6p-2, 0x1.ba5e353f7cedap-4, 0x1.78d4fdf3b645ap-3, 0x1.78d4fdf3b645ap-5, 0x1.3f7ced916872bp-3, 0x1.46a7ef9db22d1p-2, 0x1.851eb851eb852p-4, 0x1.5e353f7ced917p-3, 0x1.ae147ae147ae1p-4, 0x1.89374bc6a7efap-5, 0x1.a9fbe76c8b43ap-6, 0x1.d2f1a9fbe76c9p-4, 0x1.126e978d4fdf4p-4, 0x1.1a9fbe76c8b44p-4, 0x1.47ae147ae147bp-5, 0x1.2f1a9fbe76c8cp-4, 0x1.89374bc6a7efap-5, 0x1.70a3d70a3d70ap-5, 0x1.0e5604189374cp-4, 0x1.126e978d4fdf4p-4, 0x1.2f1a9fbe76c8cp-4, 0x1.eb851eb851eb8p-5, 0x1.020c49ba5e354p-4, 0x1.999999999999ap-6, 0x1.47ae147ae147bp-6, 0x1.8f5c28f5c28f6p-3, 0x1.fdf3b645a1cacp-3, 0x1.71a9fbe76c8b4p-2, 0x1.a2d0e56041893p-2, 0x1.0f5c28f5c28f6p-2, 0x1.49ba5e353f7cfp-2, 0x1.1b22d0e560418p-1, 0x1.78d4fdf3b645ap-2, 0x1.16872b020c49cp-1, 0x1.a8f5c28f5c28fp-2, 0x1.1cac083126e98p-2, 0x1.28f5c28f5c28fp-2, 0x1.126e978d4fdf4p-1, 0x1.65604189374bcp-1, 0x1.b95810624dd2fp-2, 0x1.0dd2f1a9fbe77p-1, 0x1.0395810624dd3p-1, 0x1.e5604189374bcp-2, 0x1.020c49ba5e354p-1, 0x1.cdd2f1a9fbe77p-2, 0x1.11eb851eb851fp-1, 0x1.eb851eb851eb8p-2, 0x1.1ba5e353f7cedp-1, 0x1.1810624dd2f1bp-1, 0x1.ced916872b021p-2, 0x1.0dd2f1a9fbe77p-1, 0x1.e04189374bc6ap-2, 0x1.fced916872b03p-2, 0x1.cbc6a7ef9db23p-2, 0x1.083126e978d5p-1, 0x1.fae147ae147aep-2, 0x1.f0a3d70a3d70ap-2, 0x1.0ac083126e979p-1, 0x1.e45a1cac08312p-2, 0x1.d2f1a9fbe76c9p-2, 0x1.f4bc6a7ef9db2p-2, 0x1.1374bc6a7ef9dp-1, 0x1.eb851eb851eb8p-2, 0x1.e978d4fdf3b64p-2, 0x1.eb851eb851eb8p-6, 0x1.0624dd2f1a9fcp-5, 0x1.0e5604189374cp-3, 0x1.0e5604189374cp-4, 0x1.1eb851eb851ecp-5, 0x1.a9fbe76c8b43ap-7, 0x1.0624dd2f1a9fcp-6, 0x1.916872b020c4ap-4, 0x1.34395810624ddp-2, 0x1.126e978d4fdf4p-2, 0x1.0a3d70a3d70a4p-4, 0x1.89374bc6a7efap-8, 0x1.28f5c28f5c28fp-2, 0x1.916872b020c4ap-3, 0x1.7ef9db22d0e56p-3, 0x1.75c28f5c28f5cp-2, 0x1.3c6a7ef9db22dp-2, 0x1.eb851eb851eb8p-2, 0x1.5c28f5c28f5c3p-2, 0x1.d70a3d70a3d71p-3, 0x1.7ae147ae147aep-2, 0x1.65604189374bcp-2, 0x1.5d2f1a9fbe76dp-2, 0x1.a3d70a3d70a3dp-2, 0x1.0b4395810624ep-1, 0x1.6f9db22d0e56p-2, 0x1.db22d0e560418p-2, 0x1.bc6a7ef9db22dp-2, 0x1.05a1cac083127p-1, 0x1.a2d0e56041893p-2, 0x1.d3f7ced916873p-2, 0x1.e147ae147ae14p-2, 0x1.a6e978d4fdf3bp-2, 0x1.d604189374bc6p-2, 0x1.e04189374bc6ap-2, 0x1.cdd2f1a9fbe77p-2, 0x1.989374bc6a7efp-2, 0x1.fae147ae147aep-2, 0x1.010624dd2f1aap-1), dim = c(39L, 3L), dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"), c("sand", "silt", "clay"))) else structure(list(sand = c(0x1.8cccccccccccdp-1, 0x1.7020c49ba5e36p-1, 0x1.0395810624dd3p-1, 0x1.0b4395810624ep-1, 0x1.6666666666666p-1, 0x1.547ae147ae148p-1, 0x1.b95810624dd2fp-2, 0x1.116872b020c4ap-1, 0x1.3d70a3d70a3d7p-3, 0x1.449ba5e353f7dp-2, 0x1.50624dd2f1aap-1, 0x1.6872b020c49bbp-1, 0x1.645a1cac08312p-3, 0x1.b22d0e5604189p-4, 0x1.872b020c49ba6p-2, 0x1.ba5e353f7cedap-4, 0x1.78d4fdf3b645ap-3, 0x1.78d4fdf3b645ap-5, 0x1.3f7ced916872bp-3, 0x1.46a7ef9db22d1p-2, 0x1.851eb851eb852p-4, 0x1.5e353f7ced917p-3, 0x1.ae147ae147ae1p-4, 0x1.89374bc6a7efap-5, 0x1.a9fbe76c8b43ap-6, 0x1.d2f1a9fbe76c9p-4, 0x1.126e978d4fdf4p-4, 0x1.1a9fbe76c8b44p-4, 0x1.47ae147ae147bp-5, 0x1.2f1a9fbe76c8cp-4, 0x1.89374bc6a7efap-5, 0x1.70a3d70a3d70ap-5, 0x1.0e5604189374cp-4, 0x1.126e978d4fdf4p-4, 0x1.2f1a9fbe76c8cp-4, 0x1.eb851eb851eb8p-5, 0x1.020c49ba5e354p-4, 0x1.999999999999ap-6, 0x1.47ae147ae147bp-6), silt = c(0x1.8f5c28f5c28f6p-3, 0x1.fdf3b645a1cacp-3, 0x1.71a9fbe76c8b4p-2, 0x1.a2d0e56041893p-2, 0x1.0f5c28f5c28f6p-2, 0x1.49ba5e353f7cfp-2, 0x1.1b22d0e560418p-1, 0x1.78d4fdf3b645ap-2, 0x1.16872b020c49cp-1, 0x1.a8f5c28f5c28fp-2, 0x1.1cac083126e98p-2, 0x1.28f5c28f5c28fp-2, 0x1.126e978d4fdf4p-1, 0x1.65604189374bcp-1, 0x1.b95810624dd2fp-2, 0x1.0dd2f1a9fbe77p-1, 0x1.0395810624dd3p-1, 0x1.e5604189374bcp-2, 0x1.020c49ba5e354p-1, 0x1.cdd2f1a9fbe77p-2, 0x1.11eb851eb851fp-1, 0x1.eb851eb851eb8p-2, 0x1.1ba5e353f7cedp-1, 0x1.1810624dd2f1bp-1, 0x1.ced916872b021p-2, 0x1.0dd2f1a9fbe77p-1, 0x1.e04189374bc6ap-2, 0x1.fced916872b03p-2, 0x1.cbc6a7ef9db23p-2, 0x1.083126e978d5p-1, 0x1.fae147ae147aep-2, 0x1.f0a3d70a3d70ap-2, 0x1.0ac083126e979p-1, 0x1.e45a1cac08312p-2, 0x1.d2f1a9fbe76c9p-2, 0x1.f4bc6a7ef9db2p-2, 0x1.1374bc6a7ef9dp-1, 0x1.eb851eb851eb8p-2, 0x1.e978d4fdf3b64p-2), clay = c(0x1.eb851eb851eb8p-6, 0x1.0624dd2f1a9fcp-5, 0x1.0e5604189374cp-3, 0x1.0e5604189374cp-4, 0x1.1eb851eb851ecp-5, 0x1.a9fbe76c8b43ap-7, 0x1.0624dd2f1a9fcp-6, 0x1.916872b020c4ap-4, 0x1.34395810624ddp-2, 0x1.126e978d4fdf4p-2, 0x1.0a3d70a3d70a4p-4, 0x1.89374bc6a7efap-8, 0x1.28f5c28f5c28fp-2, 0x1.916872b020c4ap-3, 0x1.7ef9db22d0e56p-3, 0x1.75c28f5c28f5cp-2, 0x1.3c6a7ef9db22dp-2, 0x1.eb851eb851eb8p-2, 0x1.5c28f5c28f5c3p-2, 0x1.d70a3d70a3d71p-3, 0x1.7ae147ae147aep-2, 0x1.65604189374bcp-2, 0x1.5d2f1a9fbe76dp-2, 0x1.a3d70a3d70a3dp-2, 0x1.0b4395810624ep-1, 0x1.6f9db22d0e56p-2, 0x1.db22d0e560418p-2, 0x1.bc6a7ef9db22dp-2, 0x1.05a1cac083127p-1, 0x1.a2d0e56041893p-2, 0x1.d3f7ced916873p-2, 0x1.e147ae147ae14p-2, 0x1.a6e978d4fdf3bp-2, 0x1.d604189374bc6p-2, 0x1.e04189374bc6ap-2, 0x1.cdd2f1a9fbe77p-2, 0x1.989374bc6a7efp-2, 0x1.fae147ae147aep-2, 0x1.010624dd2f1aap-1)), class = "data.frame", row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39")), dims = 3L, dim.names = c("sand", "silt", "clay"), obs = 39L, valid_obs = 39L, normalized = TRUE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
DR_AL_withNA_check <- structure(c(NA, NA, 0x1.0395810624dd3p-1, 0x1.0c1175db89e2cp-1, 0x1.6666666666666p-1, 0x1.547ae147ae148p-1, 0x1.b95810624dd3p-2, 0x1.116872b020c4ap-1, 0x1.3d70a3d70a3d7p-3, 0x1.449ba5e353f7dp-2, 0x1.50624dd2f1aap-1, 0x1.6872b020c49bbp-1, 0x1.645a1cac08312p-3, 0x1.b22d0e5604189p-4, 0x1.872b020c49ba6p-2, 0x1.ba5e353f7cedap-4, 0x1.78d4fdf3b645ap-3, 0x1.78d4fdf3b645ap-5, 0x1.3f7ced916872bp-3, 0x1.46a7ef9db22d1p-2, 0x1.851eb851eb852p-4, 0x1.5e353f7ced917p-3, 0x1.ae147ae147ae1p-4, 0x1.87427bcc092b8p-5, 0x1.a9fbe76c8b43ap-6, 0x1.d2f1a9fbe76c9p-4, 0x1.126e978d4fdf4p-4, 0x1.1a9fbe76c8b44p-4, 0x1.47ae147ae147bp-5, 0x1.2f684bda12f69p-4, 0x1.89374bc6a7efap-5, 0x1.70a3d70a3d70ap-5, 0x1.0e5604189374cp-4, 0x1.12b4eabb112b6p-4, 0x1.2f684bda12f69p-4, 0x1.eb851eb851eb8p-5, 0x1.020c49ba5e355p-4, 0x1.999999999999ap-6, 0x1.47ae147ae147bp-6, NA, NA, 0x1.71a9fbe76c8b4p-2, 0x1.a41383b921a3p-2, 0x1.0f5c28f5c28f6p-2, 0x1.49ba5e353f7cfp-2, 0x1.1b22d0e560419p-1, 0x1.78d4fdf3b645ap-2, 0x1.16872b020c49cp-1, 0x1.a8f5c28f5c28fp-2, 0x1.1cac083126e98p-2, 0x1.28f5c28f5c28fp-2, 0x1.126e978d4fdf4p-1, 0x1.65604189374bcp-1, 0x1.b95810624dd2fp-2, 0x1.0dd2f1a9fbe77p-1, 0x1.0395810624dd3p-1, 0x1.e5604189374bcp-2, 0x1.020c49ba5e354p-1, 0x1.cdd2f1a9fbe77p-2, 0x1.11eb851eb851fp-1, 0x1.eb851eb851eb8p-2, 0x1.1ba5e353f7cedp-1, 0x1.16abaf81a7dd6p-1, 0x1.ced916872b021p-2, 0x1.0dd2f1a9fbe77p-1, 0x1.e04189374bc6ap-2, 0x1.fced916872b03p-2, 0x1.cbc6a7ef9db23p-2, 0x1.0874da5490875p-1, 0x1.fae147ae147aep-2, 0x1.f0a3d70a3d70ap-2, 0x1.0ac083126e979p-1, 0x1.e4d63af05e4d7p-2, 0x1.d36952421d36ap-2, 0x1.f4bc6a7ef9db2p-2, 0x1.1374bc6a7ef9ep-1, 0x1.eb851eb851eb8p-2, 0x1.e978d4fdf3b64p-2, NA, NA, 0x1.0e5604189374cp-3, 0x1.0f26423f2a5ep-4, 0x1.1eb851eb851ecp-5, 0x1.a9fbe76c8b43ap-7, 0x1.0624dd2f1a9fdp-6, 0x1.916872b020c4ap-4, 0x1.34395810624ddp-2, 0x1.126e978d4fdf4p-2, 0x1.0a3d70a3d70a4p-4, 0x1.89374bc6a7efap-8, 0x1.28f5c28f5c28fp-2, 0x1.916872b020c4ap-3, 0x1.7ef9db22d0e56p-3, 0x1.75c28f5c28f5cp-2, 0x1.3c6a7ef9db22dp-2, 0x1.eb851eb851eb8p-2, 0x1.5c28f5c28f5c3p-2, 0x1.d70a3d70a3d71p-3, 0x1.7ae147ae147aep-2, 0x1.65604189374bcp-2, 0x1.5d2f1a9fbe76dp-2, 0x1.a1c051832f1fcp-2, 0x1.0b4395810624ep-1, 0x1.6f9db22d0e56p-2, 0x1.db22d0e560418p-2, 0x1.bc6a7ef9db22dp-2, 0x1.05a1cac083127p-1, 0x1.a33c38605a33cp-2, 0x1.d3f7ced916873p-2, 0x1.e147ae147ae14p-2, 0x1.a6e978d4fdf3bp-2, 0x1.d67c8a60dd67dp-2, 0x1.e0bc9ac75e0bcp-2, 0x1.cdd2f1a9fbe77p-2, 0x1.989374bc6a7fp-2, 0x1.fae147ae147aep-2, 0x1.010624dd2f1aap-1), dim = c(39L, 3L), dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"), c("sand", "silt", "clay")), Y.original = if(._DR_ver_0_8) structure(c(NA, NaN, 0x1.0395810624dd3p-1, 0x1.0b4395810624ep-1, 0x1.6666666666666p-1, 0x1.547ae147ae148p-1, 0x1.b95810624dd2fp-2, 0x1.116872b020c4ap-1, 0x1.3d70a3d70a3d7p-3, 0x1.449ba5e353f7dp-2, 0x1.50624dd2f1aap-1, 0x1.6872b020c49bbp-1, 0x1.645a1cac08312p-3, 0x1.b22d0e5604189p-4, 0x1.872b020c49ba6p-2, 0x1.ba5e353f7cedap-4, 0x1.78d4fdf3b645ap-3, 0x1.78d4fdf3b645ap-5, 0x1.3f7ced916872bp-3, 0x1.46a7ef9db22d1p-2, 0x1.851eb851eb852p-4, 0x1.5e353f7ced917p-3, 0x1.ae147ae147ae1p-4, 0x1.89374bc6a7efap-5, 0x1.a9fbe76c8b43ap-6, 0x1.d2f1a9fbe76c9p-4, 0x1.126e978d4fdf4p-4, 0x1.1a9fbe76c8b44p-4, 0x1.47ae147ae147bp-5, 0x1.2f1a9fbe76c8cp-4, 0x1.89374bc6a7efap-5, 0x1.70a3d70a3d70ap-5, 0x1.0e5604189374cp-4, 0x1.126e978d4fdf4p-4, 0x1.2f1a9fbe76c8cp-4, 0x1.eb851eb851eb8p-5, 0x1.020c49ba5e354p-4, 0x1.999999999999ap-6, 0x1.47ae147ae147bp-6, 0x1.8f5c28f5c28f6p-3, 0x1.fdf3b645a1cacp-3, 0x1.71a9fbe76c8b4p-2, 0x1.a2d0e56041893p-2, 0x1.0f5c28f5c28f6p-2, 0x1.49ba5e353f7cfp-2, 0x1.1b22d0e560418p-1, 0x1.78d4fdf3b645ap-2, 0x1.16872b020c49cp-1, 0x1.a8f5c28f5c28fp-2, 0x1.1cac083126e98p-2, 0x1.28f5c28f5c28fp-2, 0x1.126e978d4fdf4p-1, 0x1.65604189374bcp-1, 0x1.b95810624dd2fp-2, 0x1.0dd2f1a9fbe77p-1, 0x1.0395810624dd3p-1, 0x1.e5604189374bcp-2, 0x1.020c49ba5e354p-1, 0x1.cdd2f1a9fbe77p-2, 0x1.11eb851eb851fp-1, 0x1.eb851eb851eb8p-2, 0x1.1ba5e353f7cedp-1, 0x1.1810624dd2f1bp-1, 0x1.ced916872b021p-2, 0x1.0dd2f1a9fbe77p-1, 0x1.e04189374bc6ap-2, 0x1.fced916872b03p-2, 0x1.cbc6a7ef9db23p-2, 0x1.083126e978d5p-1, 0x1.fae147ae147aep-2, 0x1.f0a3d70a3d70ap-2, 0x1.0ac083126e979p-1, 0x1.e45a1cac08312p-2, 0x1.d2f1a9fbe76c9p-2, 0x1.f4bc6a7ef9db2p-2, 0x1.1374bc6a7ef9dp-1, 0x1.eb851eb851eb8p-2, 0x1.e978d4fdf3b64p-2, 0x1.eb851eb851eb8p-6, 0x1.0624dd2f1a9fcp-5, 0x1.0e5604189374cp-3, 0x1.0e5604189374cp-4, 0x1.1eb851eb851ecp-5, 0x1.a9fbe76c8b43ap-7, 0x1.0624dd2f1a9fcp-6, 0x1.916872b020c4ap-4, 0x1.34395810624ddp-2, 0x1.126e978d4fdf4p-2, 0x1.0a3d70a3d70a4p-4, 0x1.89374bc6a7efap-8, 0x1.28f5c28f5c28fp-2, 0x1.916872b020c4ap-3, 0x1.7ef9db22d0e56p-3, 0x1.75c28f5c28f5cp-2, 0x1.3c6a7ef9db22dp-2, 0x1.eb851eb851eb8p-2, 0x1.5c28f5c28f5c3p-2, 0x1.d70a3d70a3d71p-3, 0x1.7ae147ae147aep-2, 0x1.65604189374bcp-2, 0x1.5d2f1a9fbe76dp-2, 0x1.a3d70a3d70a3dp-2, 0x1.0b4395810624ep-1, 0x1.6f9db22d0e56p-2, 0x1.db22d0e560418p-2, 0x1.bc6a7ef9db22dp-2, 0x1.05a1cac083127p-1, 0x1.a2d0e56041893p-2, 0x1.d3f7ced916873p-2, 0x1.e147ae147ae14p-2, 0x1.a6e978d4fdf3bp-2, 0x1.d604189374bc6p-2, 0x1.e04189374bc6ap-2, 0x1.cdd2f1a9fbe77p-2, 0x1.989374bc6a7efp-2, 0x1.fae147ae147aep-2, 0x1.010624dd2f1aap-1), dim = c(39L, 3L), dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"), c("sand", "silt", "clay"))) else structure(list(sand = c(NA, NA, 0x1.0395810624dd3p-1, 0x1.0b4395810624ep-1, 0x1.6666666666666p-1, 0x1.547ae147ae148p-1, 0x1.b95810624dd2fp-2, 0x1.116872b020c4ap-1, 0x1.3d70a3d70a3d7p-3, 0x1.449ba5e353f7dp-2, 0x1.50624dd2f1aap-1, 0x1.6872b020c49bbp-1, 0x1.645a1cac08312p-3, 0x1.b22d0e5604189p-4, 0x1.872b020c49ba6p-2, 0x1.ba5e353f7cedap-4, 0x1.78d4fdf3b645ap-3, 0x1.78d4fdf3b645ap-5, 0x1.3f7ced916872bp-3, 0x1.46a7ef9db22d1p-2, 0x1.851eb851eb852p-4, 0x1.5e353f7ced917p-3, 0x1.ae147ae147ae1p-4, 0x1.89374bc6a7efap-5, 0x1.a9fbe76c8b43ap-6, 0x1.d2f1a9fbe76c9p-4, 0x1.126e978d4fdf4p-4, 0x1.1a9fbe76c8b44p-4, 0x1.47ae147ae147bp-5, 0x1.2f1a9fbe76c8cp-4, 0x1.89374bc6a7efap-5, 0x1.70a3d70a3d70ap-5, 0x1.0e5604189374cp-4, 0x1.126e978d4fdf4p-4, 0x1.2f1a9fbe76c8cp-4, 0x1.eb851eb851eb8p-5, 0x1.020c49ba5e354p-4, 0x1.999999999999ap-6, 0x1.47ae147ae147bp-6), silt = c(NA, NA, 0x1.71a9fbe76c8b4p-2, 0x1.a2d0e56041893p-2, 0x1.0f5c28f5c28f6p-2, 0x1.49ba5e353f7cfp-2, 0x1.1b22d0e560418p-1, 0x1.78d4fdf3b645ap-2, 0x1.16872b020c49cp-1, 0x1.a8f5c28f5c28fp-2, 0x1.1cac083126e98p-2, 0x1.28f5c28f5c28fp-2, 0x1.126e978d4fdf4p-1, 0x1.65604189374bcp-1, 0x1.b95810624dd2fp-2, 0x1.0dd2f1a9fbe77p-1, 0x1.0395810624dd3p-1, 0x1.e5604189374bcp-2, 0x1.020c49ba5e354p-1, 0x1.cdd2f1a9fbe77p-2, 0x1.11eb851eb851fp-1, 0x1.eb851eb851eb8p-2, 0x1.1ba5e353f7cedp-1, 0x1.1810624dd2f1bp-1, 0x1.ced916872b021p-2, 0x1.0dd2f1a9fbe77p-1, 0x1.e04189374bc6ap-2, 0x1.fced916872b03p-2, 0x1.cbc6a7ef9db23p-2, 0x1.083126e978d5p-1, 0x1.fae147ae147aep-2, 0x1.f0a3d70a3d70ap-2, 0x1.0ac083126e979p-1, 0x1.e45a1cac08312p-2, 0x1.d2f1a9fbe76c9p-2, 0x1.f4bc6a7ef9db2p-2, 0x1.1374bc6a7ef9dp-1, 0x1.eb851eb851eb8p-2, 0x1.e978d4fdf3b64p-2), clay = c(NA, NA, 0x1.0e5604189374cp-3, 0x1.0e5604189374cp-4, 0x1.1eb851eb851ecp-5, 0x1.a9fbe76c8b43ap-7, 0x1.0624dd2f1a9fcp-6, 0x1.916872b020c4ap-4, 0x1.34395810624ddp-2, 0x1.126e978d4fdf4p-2, 0x1.0a3d70a3d70a4p-4, 0x1.89374bc6a7efap-8, 0x1.28f5c28f5c28fp-2, 0x1.916872b020c4ap-3, 0x1.7ef9db22d0e56p-3, 0x1.75c28f5c28f5cp-2, 0x1.3c6a7ef9db22dp-2, 0x1.eb851eb851eb8p-2, 0x1.5c28f5c28f5c3p-2, 0x1.d70a3d70a3d71p-3, 0x1.7ae147ae147aep-2, 0x1.65604189374bcp-2, 0x1.5d2f1a9fbe76dp-2, 0x1.a3d70a3d70a3dp-2, 0x1.0b4395810624ep-1, 0x1.6f9db22d0e56p-2, 0x1.db22d0e560418p-2, 0x1.bc6a7ef9db22dp-2, 0x1.05a1cac083127p-1, 0x1.a2d0e56041893p-2, 0x1.d3f7ced916873p-2, 0x1.e147ae147ae14p-2, 0x1.a6e978d4fdf3bp-2, 0x1.d604189374bc6p-2, 0x1.e04189374bc6ap-2, 0x1.cdd2f1a9fbe77p-2, 0x1.989374bc6a7efp-2, 0x1.fae147ae147aep-2, 0x1.010624dd2f1aap-1)), row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"), class = "data.frame"), dims = 3L, dim.names = c("sand", "silt", "clay"), obs = 39L, valid_obs = 37L, normalized = TRUE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint

DA_PD_check <- structure(c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10, 0x1.0624dd2f1a9fcp-10, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1), dim = c(11L, 2L), dimnames = list(    NULL, c("y_2", "y_2")), Y.original = if(._DR_ver_0_8) structure(c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10, 0x1.0624dd2f1a9fcp-10, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2"))) else structure(list(y_2 = c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10), y_2 = c(0x1.0624dd2f1a9fcp-10, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = FALSE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
DA_PD_trafo_check <- structure(c(0x1.e843068430684p-1, 0x1.ba2e8ba2e8ba3p-1, 0x1.8ba2e8ba2e8bap-1, 0x1.5d1745d1745d1p-1, 0x1.2e8ba2e8ba2e9p-1, 0x1p-1, 0x1.a2e8ba2e8ba2fp-2, 0x1.45d1745d1745cp-2, 0x1.d1745d1745d16p-3, 0x1.1745d1745d174p-3, 0x1.7bcf97bcf97bdp-5, 0x1.7bcf97bcf97bdp-5, 0x1.1745d1745d174p-3, 0x1.d1745d1745d17p-3, 0x1.45d1745d1745ep-2, 0x1.a2e8ba2e8ba2fp-2, 0x1p-1, 0x1.2e8ba2e8ba2e9p-1, 0x1.5d1745d1745d2p-1, 0x1.8ba2e8ba2e8bap-1, 0x1.ba2e8ba2e8ba3p-1, 0x1.e843068430684p-1), dim = c(11L, 2L), dimnames = list(    NULL, c("y_2", "y_2")), Y.original = if(._DR_ver_0_8) structure(c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10, 0x1.0624dd2f1a9fcp-10, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2"))) else structure(list(y_2 = c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10), y_2 = c(0x1.0624dd2f1a9fcp-10, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = FALSE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint

DA_aPD_norm_check <- structure(c(0x1.fff2e188c0121p-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10, 0x1.a3cee7fdbec39p-14, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1), dim = c(11L, 2L), dimnames = list(    NULL, c("y_2", "y_2")), Y.original = if(._DR_ver_0_8) structure(c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10, 0x1.a36e2eb1c432dp-14, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2"))) else structure(list(y_2 = c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10), y_2 = c(0x1.a36e2eb1c432dp-14, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = TRUE, transformed = FALSE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint
DA_aPD_norm_trafo_check <- structure(c(0x1.e8ae41650bb36p-1, 0x1.ba2e8ba2e8ba3p-1, 0x1.8ba2e8ba2e8bap-1, 0x1.5d1745d1745d1p-1, 0x1.2e8ba2e8ba2e9p-1, 0x1p-1, 0x1.a2e8ba2e8ba2fp-2, 0x1.45d1745d1745cp-2, 0x1.d1745d1745d16p-3, 0x1.1745d1745d174p-3, 0x1.7bcf97bcf97bdp-5, 0x1.751be9af44cb1p-5, 0x1.1745d1745d174p-3, 0x1.d1745d1745d17p-3, 0x1.45d1745d1745ep-2, 0x1.a2e8ba2e8ba2fp-2, 0x1p-1, 0x1.2e8ba2e8ba2e9p-1, 0x1.5d1745d1745d2p-1, 0x1.8ba2e8ba2e8bap-1, 0x1.ba2e8ba2e8ba3p-1, 0x1.e843068430684p-1), dim = c(11L, 2L), dimnames = list(    NULL, c("y_2", "y_2")), Y.original = if(._DR_ver_0_8) structure(c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10, 0x1.a36e2eb1c432dp-14, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2"))) else structure(list(y_2 = c(0x1.ff7ced916872bp-1, 0x1.ccccccccccccdp-1, 0x1.999999999999ap-1, 0x1.6666666666666p-1, 0x1.3333333333333p-1, 0x1p-1, 0x1.999999999999ap-2, 0x1.3333333333332p-2, 0x1.9999999999998p-3, 0x1.9999999999998p-4, 0x1.0624dd2f1aap-10), y_2 = c(0x1.a36e2eb1c432dp-14, 0x1.999999999999ap-4, 0x1.999999999999ap-3, 0x1.3333333333334p-2, 0x1.999999999999ap-2, 0x1p-1, 0x1.3333333333333p-1, 0x1.6666666666667p-1, 0x1.999999999999ap-1, 0x1.ccccccccccccdp-1, 0x1.ff7ced916872bp-1)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = TRUE, transformed = TRUE, base = if(._DR_ver_0_8) 1L else 0x1p+0, class = "DirichletRegData") # nolint

AL <- ArcticLake[, 1L:3L]
mat_AL <- as.matrix(AL)

AL_withNA <- AL
AL_withNA[1L, 1L] <-  NA
AL_withNA[2L, 1L] <-  NaN
#AL_withNA[3L, 1L] <-  Inf
#AL_withNA[4L, 1L] <- -Inf

Perfect_Data <- cbind("y_2" = c(0.001, seq(0.1, 0.9, 0.1), 0.999))
Perfect_Data <- cbind("y_1" = 1.0 - Perfect_Data, Perfect_Data)

almost_Perfect_Data <- Perfect_Data
almost_Perfect_Data[1L, 2L] <- 0.0001

test_that("Test Various Scenarios with the Arctic Lake Dataset", {
  skip_ext_cran()
  expect_warning(DR_AL <- DR_data(AL), regexp = "normalization forced", fixed = TRUE) # nolint
  expect_identical(DR_AL, DR_AL_check)
  expect_equivalent(DR_AL, suppressWarnings(DR_data(mat_AL)))
  expect_s3_class(DR_AL, "DirichletRegData")
  expect_no_error(print(DR_AL)) |> capture.output()
  expect_no_error(summary(DR_AL)) |> capture.output()
  
  expect_warning(DR_AL_withNA <- DR_data(AL_withNA), regexp = "normalization forced", fixed = TRUE) # nolint
  expect_identical(DR_AL_withNA[, ], DR_AL_withNA_check[, ])
  expect_s3_class(DR_AL_withNA, "DirichletRegData")
  expect_no_error(print(DR_AL_withNA)) |> capture.output()
  expect_no_error(summary(DR_AL_withNA)) |> capture.output()
  
  expect_warning(DR_AL_withNA_base3 <- DR_data(AL_withNA, base = 3), regexp = "normalization forced", fixed = TRUE) # nolint
  expect_identical(attr(DR_AL_withNA_base3, "base"), if(._DR_ver_0_8) 3L else 0x1.8p+1)

  # process "perfect" data
  expect_no_condition(DA_PD <- DR_data(Perfect_Data)) # nolint
  expect_identical(DA_PD, DA_PD_check)
  if(._DR_ver_0_8) expect_warning(DA_PD_trafo <- DR_data(Perfect_Data, trafo = TRUE), regexp = "transformation forced", fixed = TRUE) else # nolint
                   expect_no_condition(DA_PD_trafo <- DR_data(Perfect_Data, trafo = TRUE)) # older versions did no issue warning # nolint
  expect_identical(DA_PD_trafo, DA_PD_trafo_check)
  
  # almost perfect data
  expect_warning(DA_aPD_norm <- DR_data(almost_Perfect_Data), regexp = "normalization forced", fixed = TRUE) # nolint
  expect_identical(DA_aPD_norm, DA_aPD_norm_check)
  expect_warning(DA_aPD_norm_trafo <- DR_data(almost_Perfect_Data, trafo = TRUE), regexp = if(._DR_ver_0_8) "normalization forced.+transformation forced" else "normalization forced") # old version produced an inadequate warning (did not mention trnasformation) # nolint
  expect_identical(DA_aPD_norm_trafo, DA_aPD_norm_trafo_check)
})


# TODO check base variations!
