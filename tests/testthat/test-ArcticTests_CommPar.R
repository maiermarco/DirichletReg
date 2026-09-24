#TO DO:
#anova(resC1)
#residuals(resC1)
#update(resC1)

test_that("Arctic Lake - Check the Original Data", {
  expect_true(exists("ArcticLake"))
  expect_identical(dim(ArcticLake), c(39L, 4L))
  expect_identical(names(ArcticLake), c("sand", "silt", "clay", "depth"))
  expect_true(all(unlist(lapply(ArcticLake, is.numeric))))
})

AL <- ArcticLake[, 4L, drop = FALSE]

test_that("Arctic Lake - Data Transformation", {
  expect_identical(dim(AL), c(39L, 1L))
  expect_s3_class(AL, "data.frame", exact = TRUE)
  expect_warning(AL$Y <<- DR_data(ArcticLake[, 1:3]), regexp = "normalization\\ forced$")
  expect_equal(unname(rowSums(AL$Y)), rep(1.0, 39L))
})

#--- Checks: Common Model ------------------------------------------------------

#                                                                              #
#                                                                             ##
#                                                                              #
#                                                                              #
#                                                                              #

resC1 <- DirichReg(Y ~ 1, AL)

load("testdata/resC1.RData")

test_that("Model Estimation", {
  expect_equal(resC1_mathematica$MLE    , resC1$logLik)
  expect_equal(resC1_mathematica$DEV    , -2.0*resC1$logLik)
  expect_equal(resC1_mathematica$COEFS  , unname(resC1$coefficients))
  expect_equal(resC1_mathematica$SE     , unname(resC1$se))
  expect_equal(resC1_mathematica$Z      , unname(resC1$coefficients / resC1$se))
  expect_equal(resC1_mathematica$P      , 2*pnorm(-abs(unname(resC1$coefficients / resC1$se))))
  expect_equal(resC1_mathematica$HESSIAN, unname(resC1$hessian))
  expect_equal(resC1_mathematica$VCOV   , unname(resC1$vcov))
})

test_that("Methods", {
  expect_equal(resC1_mathematica$NOBS , nobs(resC1))
  expect_equal(resC1_mathematica$MLE  , as.numeric(logLik(resC1)))
  expect_equal(resC1_mathematica$NPAR , attr(logLik(resC1), "df", exact = TRUE))
  expect_equal(resC1_mathematica$AIC  , AIC(resC1))
  expect_equal(resC1_mathematica$BIC  , BIC(resC1))
  expect_equal(resC1_mathematica$COEFS, unname(unlist(coef(resC1))))
  expect_equal(resC1_mathematica$VCOV , unname(vcov(resC1)))
  
  expect_equal(resC1_mathematica$PREDICT$ALPHA, unname(fitted(resC1, alpha=TRUE , phi=FALSE, mu=FALSE)[1L,]))
  expect_equal(resC1_mathematica$PREDICT$PHI  , unname(fitted(resC1, alpha=FALSE, phi=TRUE , mu=FALSE)[1L] ))
  expect_equal(resC1_mathematica$PREDICT$MU   , unname(fitted(resC1, alpha=FALSE, phi=FALSE, mu=TRUE )[1L,]))
  
  expect_equal(resC1_mathematica$PREDICT$ALPHA, unname(predict(resC1, data.frame("depth" = 0), alpha=TRUE , phi=FALSE, mu=FALSE)[1L,]))
  expect_equal(resC1_mathematica$PREDICT$PHI  , unname(predict(resC1, data.frame("depth" = 0), alpha=FALSE, phi=TRUE , mu=FALSE)[1L] ))
  expect_equal(resC1_mathematica$PREDICT$MU   , unname(predict(resC1, data.frame("depth" = 0), alpha=FALSE, phi=FALSE, mu=TRUE )[1L,]))
  
  conf_ints <- confint(resC1, level = c(.99, .95))
  conf_ints <- lapply(1:3, function(listelement){ sort(unlist(lapply(c(list(conf_ints$coefficients), conf_ints$ci), `[[`, listelement))) })
  conf_ints <- unname(t(matrix(unlist(conf_ints), 5)))
  expect_equal(resC1_mathematica$CONFINT, conf_ints)
})

#                                                                           ###
#                                                                          #   #
#                                                                             #
#                                                                           # 
#                                                                          #####

resC2 <- DirichReg(Y ~ depth, AL)

load("testdata/resC2.RData")

test_that("Model Estimation", {
  expect_equal(resC2_mathematica$MLE    , resC2$logLik)
  expect_equal(resC2_mathematica$DEV    , -2.0*resC2$logLik)
  expect_equal(resC2_mathematica$COEFS  , unname(resC2$coefficients))
  expect_equal(resC2_mathematica$SE     , unname(resC2$se))
  expect_equal(resC2_mathematica$Z      , unname(resC2$coefficients / resC2$se))
  expect_equal(resC2_mathematica$P      , 2*pnorm(-abs(unname(resC2$coefficients / resC2$se))))
  expect_equal(resC2_mathematica$HESSIAN, unname(resC2$hessian))
  expect_equal(resC2_mathematica$VCOV   , unname(resC2$vcov))
})

test_that("Methods", {
  expect_equal(resC2_mathematica$NOBS , nobs(resC2))
  expect_equal(resC2_mathematica$MLE  , as.numeric(logLik(resC2)))
  expect_equal(resC2_mathematica$NPAR , attr(logLik(resC2), "df", exact = TRUE))
  expect_equal(resC2_mathematica$AIC  , AIC(resC2))
  expect_equal(resC2_mathematica$BIC  , BIC(resC2))
  expect_equal(resC2_mathematica$COEFS, unname(unlist(coef(resC2))))
  expect_equal(resC2_mathematica$VCOV , unname(vcov(resC2)))
  
 #expect_equal(resC2_mathematica$PREDICT$ALPHA, unname(fitted(resC2, alpha=TRUE, phi=FALSE, mu=FALSE)[1L,]))
 #expect_equal(resC2_mathematica$PREDICT$PHI  , unname(fitted(resC2, alpha=FALSE, phi=TRUE, mu=FALSE)[1L] ))
 #expect_equal(resC2_mathematica$PREDICT$MU   , unname(fitted(resC2, alpha=FALSE, phi=FALSE, mu=TRUE)[1L,]))
  
  expect_equal(resC2_mathematica$PREDICT$ALPHA, unname(predict(resC2, data.frame("depth" = 0:150), alpha = TRUE , phi = FALSE, mu = FALSE)))
  expect_equal(resC2_mathematica$PREDICT$PHI  , unname(predict(resC2, data.frame("depth" = 0:150), alpha = FALSE, phi = TRUE , mu = FALSE)))
  expect_equal(resC2_mathematica$PREDICT$MU   , unname(predict(resC2, data.frame("depth" = 0:150), alpha = FALSE, phi = FALSE, mu = TRUE )))

  # BEWARE! UGLY CODE AHEAD!
  conf_ints <- confint(resC2, level = c(.99, .95))
  conf_mat <- matrix(NA_real_, 6L, 4L); rowz <- rep(1L:2L, 3L); listz <- rep(1L:3L, each = 2L)
  for(i in 1L:6L){
    conf_mat[i,] <- c(conf_ints$ci[[2L]][[listz[i]]][rowz[i], 1L], conf_ints$ci[[1L]][[listz[i]]][rowz[i], ], conf_ints$ci[[2L]][[listz[i]]][rowz[i], 2L])
  }
  conf_ints <- unname(cbind(conf_mat[, 1L:2L], unlist(conf_ints$coefficients), conf_mat[, 3:4]))
  expect_equal(resC2_mathematica$CONFINT, conf_ints)
})

#                                                                          ####
#                                                                              #
#                                                                            ##
#                                                                              #
#                                                                          ####

resC3 <- DirichReg(Y ~ depth + I(depth^2), AL)

load("testdata/resC3.RData")

test_that("Model Estimation", {
  expect_equal(resC3_mathematica$MLE    , resC3$logLik)
  expect_equal(resC3_mathematica$DEV    , -2.0*resC3$logLik)
  expect_equal(resC3_mathematica$COEFS  , unname(resC3$coefficients))
  expect_equal(resC3_mathematica$SE     , unname(resC3$se))
  expect_equal(resC3_mathematica$Z      , unname(resC3$coefficients / resC3$se))
  expect_equal(resC3_mathematica$P      , 2*pnorm(-abs(unname(resC3$coefficients / resC3$se))))
  expect_equal(resC3_mathematica$HESSIAN, unname(resC3$hessian))
  expect_equal(resC3_mathematica$VCOV   , unname(resC3$vcov))
})

test_that("Methods", {
  expect_equal(resC3_mathematica$NOBS , nobs(resC3))
  expect_equal(resC3_mathematica$MLE  , as.numeric(logLik(resC3)))
  expect_equal(resC3_mathematica$NPAR , attr(logLik(resC3), "df", exact = TRUE))
  expect_equal(resC3_mathematica$AIC  , AIC(resC3))
  expect_equal(resC3_mathematica$BIC  , BIC(resC3))
  expect_equal(resC3_mathematica$COEFS, unname(unlist(coef(resC3))))
  expect_equal(resC3_mathematica$VCOV , unname(vcov(resC3)))
  
 #expect_equal(resC3_mathematica$PREDICT$ALPHA, unname(fitted(resC3, alpha=TRUE , phi=FALSE, mu=FALSE)[1L,]))
 #expect_equal(resC3_mathematica$PREDICT$PHI  , unname(fitted(resC3, alpha=FALSE, phi=TRUE , mu=FALSE)[1L] ))
 #expect_equal(resC3_mathematica$PREDICT$MU   , unname(fitted(resC3, alpha=FALSE, phi=FALSE, mu=TRUE )[1L,]))
  
  expect_equal(resC3_mathematica$PREDICT$ALPHA, unname(predict(resC3, data.frame("depth" = 0:150), alpha=TRUE , phi=FALSE, mu=FALSE)))
  expect_equal(resC3_mathematica$PREDICT$PHI,   unname(predict(resC3, data.frame("depth" = 0:150), alpha=FALSE, phi=TRUE , mu=FALSE)))
  expect_equal(resC3_mathematica$PREDICT$MU,    unname(predict(resC3, data.frame("depth" = 0:150), alpha=FALSE, phi=FALSE, mu=TRUE )))
  
  # AGAIN, SORRY FOR THE UGLY CODE...
  conf_ints <- confint(resC3, level = c(.99, .95))
  conf_mat <- matrix(NA_real_, 9L, 4L); rowz <- rep(1L:3L, 3L); listz <- rep(1L:3L, each = 3L)
  for(i in 1L:9L){
    conf_mat[i,] <- c(conf_ints$ci[[2L]][[listz[i]]][rowz[i], 1L], conf_ints$ci[[1L]][[listz[i]]][rowz[i], ], conf_ints$ci[[2L]][[listz[i]]][rowz[i], 2L])
  }
  conf_ints <- unname(cbind(conf_mat[, 1L:2L], unlist(conf_ints$coefficients), conf_mat[, 3L:4L]))
  expect_equal(resC3_mathematica$CONFINT, conf_ints)
})

#                                                                            ##
#                                                                           # # 
#                                                                          #####
#                                                                             # 
#                                                                             #

# Test the update() method

resC4 <- update(resC3, . ~ . - I(depth^2) | . | . - depth - I(depth^2))

load("testdata/resC4.RData")

test_that("Model Estimation", {
  expect_equal(resC4_mathematica$MLE    , resC4$logLik)
  expect_equal(resC4_mathematica$DEV    , -2.0*resC4$logLik)
  expect_equal(resC4_mathematica$COEFS  , unname(resC4$coefficients))
  expect_equal(resC4_mathematica$SE     , unname(resC4$se))
  expect_equal(resC4_mathematica$Z      , unname(resC4$coefficients / resC4$se))
  expect_equal(resC4_mathematica$P      , 2*pnorm(-abs(unname(resC4$coefficients / resC4$se))))
  expect_equal(resC4_mathematica$HESSIAN, unname(resC4$hessian))
  expect_equal(resC4_mathematica$VCOV   , unname(resC4$vcov))
})

test_that("Methods", {
  expect_equal(resC4_mathematica$NOBS , nobs(resC4))
  expect_equal(resC4_mathematica$MLE  , as.numeric(logLik(resC4)))
  expect_equal(resC4_mathematica$NPAR , attr(logLik(resC4), "df", exact = TRUE))
  expect_equal(resC4_mathematica$AIC  , AIC(resC4))
  expect_equal(resC4_mathematica$BIC  , BIC(resC4))
  expect_equal(resC4_mathematica$COEFS, unname(unlist(coef(resC4))))
  expect_equal(resC4_mathematica$VCOV , unname(vcov(resC4)))
  
 #expect_equal(resC4_mathematica$PREDICT$ALPHA, unname(fitted(resC4, alpha=TRUE , phi=FALSE, mu=FALSE)[1L,]))
 #expect_equal(resC4_mathematica$PREDICT$PHI  , unname(fitted(resC4, alpha=FALSE, phi=TRUE , mu=FALSE)[1L] ))
 #expect_equal(resC4_mathematica$PREDICT$MU   , unname(fitted(resC4, alpha=FALSE, phi=FALSE, mu=TRUE )[1L,]))
  
  expect_equal(resC4_mathematica$PREDICT$ALPHA, unname(predict(resC4, data.frame("depth" = 0:150), alpha=TRUE , phi=FALSE, mu=FALSE)))
  expect_equal(resC4_mathematica$PREDICT$PHI  , unname(predict(resC4, data.frame("depth" = 0:150), alpha=FALSE, phi=TRUE , mu=FALSE)))
  expect_equal(resC4_mathematica$PREDICT$MU   , unname(predict(resC4, data.frame("depth" = 0:150), alpha=FALSE, phi=FALSE, mu=TRUE )))
  
  # 🙈
  conf_ints <- confint(resC4, level = c(.99, .95))
  conf_mat <- matrix(NA_real_, 6L, 4L); rowz <- c(1L, 2L, 1L, 2L, 3L, 1L); listz <- rep(1L:3L, c(2L, 3L, 1L))
  for(i in 1L:6L){
    conf_mat[i,] <- c(conf_ints$ci[[2L]][[listz[i]]][rowz[i], 1L], conf_ints$ci[[1L]][[listz[i]]][rowz[i], ], conf_ints$ci[[2L]][[listz[i]]][rowz[i], 2L])
  }
  conf_ints <- unname(cbind(conf_mat[, 1L:2L], unlist(conf_ints$coefficients), conf_mat[, 3L:4L]))
  expect_equal(resC4_mathematica$CONFINT, conf_ints)
})

#                                                                          #####
#                                                                          #    
#                                                                          #### 
#                                                                              #
#                                                                          ####


