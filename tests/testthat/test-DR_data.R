# some checks require DirichletReg >= 0.8-0
._DR_ver_0_8 <- as.numeric_version(packageDescription("DirichletReg")[["Version"]]) >= as.numeric_version("0.8-0")



################################################################################
### Basic Checks ###############################################################
################################################################################

compodat <- seq(0.1, 0.9, 0.1)

test_that("Checks that should result in errors", {
  if(._DR_ver_0_8) expect_error(DR_data(compodat, trafo = FALSE), regexp = "\"trafo\" must be a small number > 0 or TRUE.", fixed = TRUE)
  if(._DR_ver_0_8) expect_error(DR_data(compodat, trafo =   0.2), regexp = "\"trafo\" must be a small number > 0 or TRUE.", fixed = TRUE)
  if(._DR_ver_0_8) expect_error(DR_data(compodat, trafo =  -1.0), regexp = "\"trafo\" must be a small number > 0 or TRUE.", fixed = TRUE)
  
  if(._DR_ver_0_8) expect_error(DR_data(compodat, norm_tol = TRUE), regexp = "\"norm_tol\" must be a small number > 0.", fixed = TRUE)
  if(._DR_ver_0_8) expect_error(DR_data(compodat, norm_tol =  0.2), regexp = "\"norm_tol\" must be a small number > 0.", fixed = TRUE)
  if(._DR_ver_0_8) expect_error(DR_data(compodat, norm_tol = -1.0), regexp = "\"norm_tol\" must be a small number > 0.", fixed = TRUE)
})



################################################################################
### Test Vector Input ##########################################################
################################################################################

beta_data     <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_1NA <- c(NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_0   <- c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_1   <- c(1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_m1  <- c(-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta_data_p2  <- c(+2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

bd1_check <- structure(c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data", "beta_data")), Y.original = structure(list(Y.original = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = c("1 - beta_data", "beta_data"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = FALSE, base = 1, class = "DirichletRegData") # nolint
bd2_check <- structure(c(NA, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_1NA", "beta_data_1NA")), Y.original = structure(list(Y.original = c(NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = c("1 - beta_data_1NA", "beta_data_1NA"), obs = 9L, valid_obs = 8L, normalized = FALSE, transformed = FALSE, base = 1L, class = "DirichletRegData") # nolint
bd3_check <- structure(c(0.944444444444444, 0.766666666666667, 0.677777777777778, 0.588888888888889, 0.5, 0.411111111111111, 0.322222222222222, 0.233333333333333, 0.144444444444444, 0.0555555555555556, 0.233333333333333, 0.322222222222222, 0.411111111111111, 0.5, 0.588888888888889, 0.677777777777778, 0.766666666666667, 0.855555555555556), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_0", "beta_data_0")), Y.original = structure(list(Y.original = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = c("1 - beta_data_0", "beta_data_0"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = 1, class = "DirichletRegData") # nolint
bd4_check <- structure(c(0.0555555555555556, 0.766666666666667, 0.677777777777778, 0.588888888888889, 0.5, 0.411111111111111, 0.322222222222222, 0.233333333333333, 0.144444444444444, 0.944444444444444, 0.233333333333333, 0.322222222222222, 0.411111111111111, 0.5, 0.588888888888889, 0.677777777777778, 0.766666666666667, 0.855555555555556), dim = c(9L, 2L), dimnames = list(NULL, c("1 - beta_data_1", "beta_data_1")), Y.original = structure(list(Y.original = c(1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = c("1 - beta_data_1", "beta_data_1"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = 1, class = "DirichletRegData") # nolint

test_that("Testing Well-Behaved Beta-Distributed Data", {
  expect_message(bd1 <- DR_data(beta_data), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_equal(bd1, bd1_check)
  expect_s3_class(bd1, "DirichletRegData")
  expect_no_error(print(bd1))
  expect_no_error(summary(bd1))
})

if(._DR_ver_0_8){ # versions before 0.8-0 could not handle this scenario
  test_that("Testing a Beta-Distributed Data with Missings", {
    expect_message(bd2 <- DR_data(beta_data_1NA), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
    expect_equal(bd2, bd2_check)
    expect_s3_class(bd2, "DirichletRegData")
    expect_no_error(print(bd2))
    expect_no_error(summary(bd2))
  })
}

test_that("Testing a Beta-Distributed Data with a Zero", {
  expect_warning(bd3 <- DR_data(beta_data_0), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_equal(bd3, bd3_check)
  expect_s3_class(bd3, "DirichletRegData")
  expect_no_error(print(bd3))
  expect_no_error(summary(bd3))
})

test_that("Testing a Beta-Distributed Data with a One", {
  expect_warning(bd4 <- DR_data(beta_data_1), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_equal(bd4, bd4_check)
  expect_s3_class(bd4, "DirichletRegData")
  expect_no_error(print(bd4))
  expect_no_error(summary(bd4))
})

test_that("Testing a Beta-Distributed Data LESS THAN ZERO or GREATER THAN ONE", {
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
mbd1_check <- structure(c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), dim = c(9L, 2L), dimnames = list( NULL, if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_matrix", "beta_matrix")), Y.original = structure(list( beta_data = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 )), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_matrix", "beta_matrix"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = FALSE, base = 1, class = "DirichletRegData") # nolint # nolint
mbd2_check <- structure(c(NA, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), dim = c(9L, 2L), dimnames = list( NULL, if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_matrix_1NA", "beta_matrix_1NA")), Y.original = structure(list( beta_data_1NA = c(NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_matrix_1NA", "beta_matrix_1NA"), obs = 9L, valid_obs = 8L, normalized = FALSE, transformed = FALSE, base = 1, class = "DirichletRegData") # nolint
mbd3_check <- structure(c(0.944444444444444, 0.766666666666667, 0.677777777777778, 0.588888888888889, 0.5, 0.411111111111111, 0.322222222222222, 0.233333333333333, 0.144444444444444, 0.0555555555555556, 0.233333333333333, 0.322222222222222, 0.411111111111111, 0.5, 0.588888888888889, 0.677777777777778, 0.766666666666667, 0.855555555555556), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_matrix_0", "beta_matrix_0" )), Y.original = structure(list(beta_data_0 = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_matrix_0", "beta_matrix_0" ), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = 1, class = "DirichletRegData") # nolint
mbd4_check <- structure(c(0.0555555555555556, 0.766666666666667, 0.677777777777778, 0.588888888888889, 0.5, 0.411111111111111, 0.322222222222222, 0.233333333333333, 0.144444444444444, 0.944444444444444, 0.233333333333333, 0.322222222222222, 0.411111111111111, 0.5, 0.588888888888889, 0.677777777777778, 0.766666666666667, 0.855555555555556), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_matrix_1", "beta_matrix_1" )), Y.original = structure(list(beta_data_1 = c(1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_matrix_1", "beta_matrix_1" ), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = 1, class = "DirichletRegData") # nolint

test_that("Testing Well-Behaved Beta-Distributed Data (as a Matrix)", {
  expect_message(mbd1 <- DR_data(beta_matrix), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_equal(mbd1, mbd1_check)
  expect_s3_class(mbd1, "DirichletRegData")
  expect_no_error(print(mbd1))
  expect_no_error(summary(mbd1))
})

# works for versions before 0.8-0 too
test_that("Testing a Beta-Distributed Data with Missings (as a Matrix)", {
  expect_message(mbd2 <- DR_data(beta_matrix_1NA), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_equal(mbd2, mbd2_check)
  expect_s3_class(mbd2, "DirichletRegData")
  expect_no_error(print(mbd2))
  expect_no_error(summary(mbd2))
})

test_that("Testing a Beta-Distributed Data with a Zero (as a Matrix)", {
  expect_warning(mbd3 <- DR_data(beta_matrix_0), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_equal(mbd3, mbd3_check)
  expect_s3_class(mbd3, "DirichletRegData")
  expect_no_error(print(mbd3))
  expect_no_error(summary(mbd3))
})

test_that("Testing a Beta-Distributed Data with a One (as a Matrix)", {
  expect_warning(mbd4 <- DR_data(beta_matrix_1), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_equal(mbd4, mbd4_check)
  expect_s3_class(mbd4, "DirichletRegData")
  expect_no_error(print(mbd4))
  expect_no_error(summary(mbd4))
})

test_that("Testing a Beta-Distributed Data LESS THAN ZERO or GREATER THAN ONE (as a Matrix)", {
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
dbd1_check <- structure(c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), dim = c(9L, 2L), dimnames = list( NULL, if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_df", "beta_df")), Y.original = structure(list( beta_data = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 )), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data", "beta_data") else c("1 - beta_df", "beta_df"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = FALSE, base = 1, class = "DirichletRegData") # nolint
dbd2_check <- structure(c(NA, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), dim = c(9L, 2L), dimnames = list( NULL, if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_df_1NA", "beta_df_1NA")), Y.original = structure(list( beta_data_1NA = c(NA, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), row.names = c(NA, -9L), class = "data.frame"), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1NA", "beta_data_1NA") else c("1 - beta_df_1NA", "beta_df_1NA"), obs = 9L, valid_obs = 8L, normalized = FALSE, transformed = FALSE, base = 1, class = "DirichletRegData") # nolint
dbd3_check <- structure(c(0.944444444444444, 0.766666666666667, 0.677777777777778, 0.588888888888889, 0.5, 0.411111111111111, 0.322222222222222, 0.233333333333333, 0.144444444444444, 0.0555555555555556, 0.233333333333333, 0.322222222222222, 0.411111111111111, 0.5, 0.588888888888889, 0.677777777777778, 0.766666666666667, 0.855555555555556), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_df_0", "beta_df_0")), Y.original = structure(list( beta_data_0 = c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 )), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_0", "beta_data_0") else c("1 - beta_df_0", "beta_df_0"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = 1, class = "DirichletRegData") # nolint
dbd4_check <- structure(c(0.0555555555555556, 0.766666666666667, 0.677777777777778, 0.588888888888889, 0.5, 0.411111111111111, 0.322222222222222, 0.233333333333333, 0.144444444444444, 0.944444444444444, 0.233333333333333, 0.322222222222222, 0.411111111111111, 0.5, 0.588888888888889, 0.677777777777778, 0.766666666666667, 0.855555555555556), dim = c(9L, 2L), dimnames = list(NULL, if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_df_1", "beta_df_1")), Y.original = structure(list( beta_data_1 = c(1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 )), class = "data.frame", row.names = c(NA, -9L)), dims = 2L, dim.names = if(._DR_ver_0_8) c("1 - beta_data_1", "beta_data_1") else c("1 - beta_df_1", "beta_df_1"), obs = 9L, valid_obs = 9L, normalized = FALSE, transformed = TRUE, base = 1, class = "DirichletRegData") # nolint

test_that("Testing Well-Behaved Beta-Distributed Data (as a Data Frame)", {
  expect_message(dbd1 <- DR_data(beta_df), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_equal(dbd1, dbd1_check)
  expect_s3_class(dbd1, "DirichletRegData")
  expect_no_error(print(dbd1))
  expect_no_error(summary(dbd1))
})

# works for versions before 0.8-0 too
test_that("Testing a Beta-Distributed Data with Missings (as a Data Frame)", {
  expect_message(dbd2 <- DR_data(beta_df_1NA), regexp = "beta-distribution assumed", fixed = TRUE) # nolint
  expect_equal(dbd2, dbd2_check)
  expect_s3_class(dbd2, "DirichletRegData")
  expect_no_error(print(dbd2))
  expect_no_error(summary(dbd2))
})

test_that("Testing a Beta-Distributed Data with a Zero (as a Data Frame)", {
  expect_warning(dbd3 <- DR_data(beta_df_0), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_equal(dbd3, dbd3_check)
  expect_s3_class(dbd3, "DirichletRegData")
  expect_no_error(print(dbd3))
  expect_no_error(summary(dbd3))
})

test_that("Testing a Beta-Distributed Data with a One (as a Data Frame)", {
  expect_warning(dbd4 <- DR_data(beta_df_1), regexp = "transformation forced", fixed = TRUE) # nolint
  expect_equal(dbd4, dbd4_check)
  expect_s3_class(dbd4, "DirichletRegData")
  expect_no_error(print(dbd4))
  expect_no_error(summary(dbd4))
})

test_that("Testing a Beta-Distributed Data LESS THAN ZERO or GREATER THAN ONE (as a Data Frame)", {
  expect_error(DR_data(beta_df_m1), regexp = if(._DR_ver_0_8) "beta distribution cannot safely be assumed" else '"Y" contains values < 0', fixed = TRUE)
  expect_error(DR_data(beta_df_p2), regexp = "beta distribution cannot safely be assumed", fixed = TRUE)
})

rm(list = c("dbd1_check", "dbd2_check", "dbd3_check", "dbd4_check",
  "beta_df", "beta_df_0", "beta_df_1", "beta_df_1NA", "beta_df_m1", "beta_df_p2"))

rm(list = c("beta_data", "beta_data_1NA", "beta_data_0", "beta_data_1", "beta_data_m1", "beta_data_p2"))










DR_AL_check <- structure(c(0.775, 0.719, 0.507, 0.523570712136409, 0.7, 0.665, 0.431, 0.534, 0.155, 0.317, 0.657, 0.704, 0.174, 0.106, 0.382, 0.108, 0.184, 0.046, 0.156, 0.319, 0.095, 0.171, 0.105, 0.0477611940298507, 0.026, 0.114, 0.067, 0.069, 0.04, 0.0740740740740741, 0.048, 0.045, 0.066, 0.0670670670670671, 0.0740740740740741, 0.06, 0.063, 0.025, 0.02, 0.195, 0.249, 0.361, 0.410230692076229, 0.265, 0.322, 0.553, 0.368, 0.544, 0.415, 0.278, 0.29, 0.536, 0.698, 0.431, 0.527, 0.507, 0.474, 0.504, 0.451, 0.535, 0.48, 0.554, 0.544278606965174, 0.452, 0.527, 0.469, 0.497, 0.449, 0.516516516516517, 0.495, 0.485, 0.521, 0.473473473473474, 0.456456456456456, 0.489, 0.538, 0.48, 0.478, 0.03, 0.032, 0.132, 0.0661985957873621, 0.035, 0.013, 0.016, 0.098, 0.301, 0.268, 0.065, 0.006, 0.29, 0.196, 0.187, 0.365, 0.309, 0.48, 0.34, 0.23, 0.37, 0.349, 0.341, 0.407960199004975, 0.522, 0.359, 0.464, 0.434, 0.511, 0.409409409409409, 0.457, 0.47, 0.413, 0.459459459459459, 0.469469469469469, 0.451, 0.399, 0.495, 0.502), dim = c(39L, 3L), dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"), c("sand", "silt", "clay")), Y.original = structure(list( sand = c(0.775, 0.719, 0.507, 0.522, 0.7, 0.665, 0.431, 0.534, 0.155, 0.317, 0.657, 0.704, 0.174, 0.106, 0.382, 0.108, 0.184, 0.046, 0.156, 0.319, 0.095, 0.171, 0.105, 0.048, 0.026, 0.114, 0.067, 0.069, 0.04, 0.074, 0.048, 0.045, 0.066, 0.067, 0.074, 0.06, 0.063, 0.025, 0.02), silt = c(0.195, 0.249, 0.361, 0.409, 0.265, 0.322, 0.553, 0.368, 0.544, 0.415, 0.278, 0.29, 0.536, 0.698, 0.431, 0.527, 0.507, 0.474, 0.504, 0.451, 0.535, 0.48, 0.554, 0.547, 0.452, 0.527, 0.469, 0.497, 0.449, 0.516, 0.495, 0.485, 0.521, 0.473, 0.456, 0.489, 0.538, 0.48, 0.478 ), clay = c(0.03, 0.032, 0.132, 0.066, 0.035, 0.013, 0.016, 0.098, 0.301, 0.268, 0.065, 0.006, 0.29, 0.196, 0.187, 0.365, 0.309, 0.48, 0.34, 0.23, 0.37, 0.349, 0.341, 0.41, 0.522, 0.359, 0.464, 0.434, 0.511, 0.409, 0.457, 0.47, 0.413, 0.459, 0.469, 0.451, 0.399, 0.495, 0.502)), class = "data.frame", row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39")), dims = 3L, dim.names = c("sand", "silt", "clay"), obs = 39L, valid_obs = 39L, normalized = TRUE, transformed = FALSE, base = 1, class = "DirichletRegData") # nolint
DR_AL_withNA_check <- structure(c(NA, NA, 0.507, 0.523570712136409, 0.7, 0.665, 0.431, 0.534, 0.155, 0.317, 0.657, 0.704, 0.174, 0.106, 0.382, 0.108, 0.184, 0.046, 0.156, 0.319, 0.095, 0.171, 0.105, 0.0477611940298507, 0.026, 0.114, 0.067, 0.069, 0.04, 0.0740740740740741, 0.048, 0.045, 0.066, 0.0670670670670671, 0.0740740740740741, 0.06, 0.063, 0.025, 0.02, NA, NA, 0.361, 0.410230692076229, 0.265, 0.322, 0.553, 0.368, 0.544, 0.415, 0.278, 0.29, 0.536, 0.698, 0.431, 0.527, 0.507, 0.474, 0.504, 0.451, 0.535, 0.48, 0.554, 0.544278606965174, 0.452, 0.527, 0.469, 0.497, 0.449, 0.516516516516517, 0.495, 0.485, 0.521, 0.473473473473474, 0.456456456456456, 0.489, 0.538, 0.48, 0.478, NA, NA, 0.132, 0.0661985957873621, 0.035, 0.013, 0.016, 0.098, 0.301, 0.268, 0.065, 0.006, 0.29, 0.196, 0.187, 0.365, 0.309, 0.48, 0.34, 0.23, 0.37, 0.349, 0.341, 0.407960199004975, 0.522, 0.359, 0.464, 0.434, 0.511, 0.409409409409409, 0.457, 0.47, 0.413, 0.459459459459459, 0.469469469469469, 0.451, 0.399, 0.495, 0.502), dim = c(39L, 3L), dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"), c("sand", "silt", "clay")), Y.original = structure(list( sand = c(NA, NA, 0.507, 0.522, 0.7, 0.665, 0.431, 0.534, 0.155, 0.317, 0.657, 0.704, 0.174, 0.106, 0.382, 0.108, 0.184, 0.046, 0.156, 0.319, 0.095, 0.171, 0.105, 0.048, 0.026, 0.114, 0.067, 0.069, 0.04, 0.074, 0.048, 0.045, 0.066, 0.067, 0.074, 0.06, 0.063, 0.025, 0.02), silt = c(NA, NA, 0.361, 0.409, 0.265, 0.322, 0.553, 0.368, 0.544, 0.415, 0.278, 0.29, 0.536, 0.698, 0.431, 0.527, 0.507, 0.474, 0.504, 0.451, 0.535, 0.48, 0.554, 0.547, 0.452, 0.527, 0.469, 0.497, 0.449, 0.516, 0.495, 0.485, 0.521, 0.473, 0.456, 0.489, 0.538, 0.48, 0.478), clay = c(NA, NA, 0.132, 0.066, 0.035, 0.013, 0.016, 0.098, 0.301, 0.268, 0.065, 0.006, 0.29, 0.196, 0.187, 0.365, 0.309, 0.48, 0.34, 0.23, 0.37, 0.349, 0.341, 0.41, 0.522, 0.359, 0.464, 0.434, 0.511, 0.409, 0.457, 0.47, 0.413, 0.459, 0.469, 0.451, 0.399, 0.495, 0.502)), row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39" ), class = "data.frame"), dims = 3L, dim.names = c("sand", "silt", "clay"), obs = 39L, valid_obs = 37L, normalized = TRUE, transformed = FALSE, base = 1, class = "DirichletRegData") # nolint

DA_PD_check <- structure(c(0.999, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001, 0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2")), Y.original = structure(list(y_2 = c(0.999, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001), y_2 = c(0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = FALSE, transformed = FALSE, base = 1L, class = "DirichletRegData") # nolint
DA_PD_trafo_check <- structure(c(0.953636363636364, 0.863636363636364, 0.772727272727273, 0.681818181818182, 0.590909090909091, 0.5, 0.409090909090909, 0.318181818181818, 0.227272727272727, 0.136363636363636, 0.0463636363636364, 0.0463636363636364, 0.136363636363636, 0.227272727272727, 0.318181818181818, 0.409090909090909, 0.5, 0.590909090909091, 0.681818181818182, 0.772727272727273, 0.863636363636364, 0.953636363636364), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2")), Y.original = structure(list(y_2 = c(0.999, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001), y_2 = c(0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = FALSE, transformed = TRUE, base = 1L, class = "DirichletRegData") # nolint

DA_aPD_norm_check <- structure(c(0.999899909918927, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001, 0.000100090081072966, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2")), Y.original = structure(list(y_2 = c(0.999, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001), y_2 = c(1e-04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = TRUE, transformed = FALSE, base = 1L, class = "DirichletRegData") # nolint
DA_aPD_norm_trafo_check <- structure(c(0.954454463562661, 0.863636363636364, 0.772727272727273, 0.681818181818182, 0.590909090909091, 0.5, 0.409090909090909, 0.318181818181818, 0.227272727272727, 0.136363636363636, 0.0463636363636364, 0.0455455364373391, 0.136363636363636, 0.227272727272727, 0.318181818181818, 0.409090909090909, 0.5, 0.590909090909091, 0.681818181818182, 0.772727272727273, 0.863636363636364, 0.953636363636364), dim = c(11L, 2L), dimnames = list(NULL, c("y_2", "y_2")), Y.original = structure(list(y_2 = c(0.999, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001), y_2 = c(1e-04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999)), class = "data.frame", row.names = c(NA, -11L)), dims = 2L, dim.names = c("y_2", "y_2"), obs = 11L, valid_obs = 11L, normalized = TRUE, transformed = TRUE, base = 1L, class = "DirichletRegData") # nolint

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
  expect_warning(DR_AL <- DR_data(AL), regexp = "normalization forced", fixed = TRUE) # nolint
  expect_equal(DR_AL, DR_AL_check)
  expect_equivalent(DR_AL, suppressWarnings(DR_data(mat_AL)))
  expect_s3_class(DR_AL, "DirichletRegData")
  expect_no_error(print(DR_AL))
  expect_no_error(summary(DR_AL))
  
  expect_warning(DR_AL_withNA <- DR_data(AL_withNA), regexp = "normalization forced", fixed = TRUE) # nolint
  expect_equal(DR_AL_withNA[, ], DR_AL_withNA_check[, ])
  expect_s3_class(DR_AL_withNA, "DirichletRegData")
  expect_no_error(print(DR_AL_withNA))
  expect_no_error(summary(DR_AL_withNA))
  
  expect_warning(DR_AL_withNA_base3 <- DR_data(AL_withNA, base = 3), regexp = "normalization forced", fixed = TRUE) # nolint
  expect_equal(attr(DR_AL_withNA_base3, "base"), 3L)

  # process "perfect" data
  expect_no_condition(DA_PD <- DR_data(Perfect_Data)) # nolint
  expect_equal(DA_PD, DA_PD_check)
  if(._DR_ver_0_8) expect_warning(DA_PD_trafo <- DR_data(Perfect_Data, trafo = TRUE), regexp = "transformation forced", fixed = TRUE) else # nolint
                   expect_no_condition(DA_PD_trafo <- DR_data(Perfect_Data, trafo = TRUE)) # older versions did no issue warning # nolint
  expect_equal(DA_PD_trafo, DA_PD_trafo_check)
  
  # almost perfect data
  if(._DR_ver_0_8) expect_warning(DA_aPD_norm <- DR_data(almost_Perfect_Data), regexp = "normalization forced", fixed = TRUE) # nolint
  if(._DR_ver_0_8) expect_equal(DA_aPD_norm, DA_aPD_norm_check)
  if(._DR_ver_0_8) expect_warning(DA_aPD_norm_trafo <- DR_data(almost_Perfect_Data, trafo = TRUE), regexp = "normalization forced.+transformation forced") # nolint
  if(._DR_ver_0_8) expect_equal(DA_aPD_norm_trafo, DA_aPD_norm_trafo_check)
})
