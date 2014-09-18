cat("=== Arctic Lake Data Tests =====================================================\n")

context("AL: Original Data")

test_that("Arctic Lake - Original Data Structure", {
	expect_true(exists("ArcticLake"))
	expect_identical(dim(ArcticLake), c(39L, 4L))
	expect_identical(names(ArcticLake), c("sand", "silt", "clay", "depth"))
	expect_true(all(unlist(lapply(ArcticLake, function(colElement){ class(colElement) == "numeric" }))))
})

context("AL: Transformation")

AL <- ArcticLake[, 4, drop=FALSE]

test_that("Arctic Lake - Data Transformation", {
	expect_identical(dim(AL), c(39L, 1L))
	expect_identical(class(AL), "data.frame")
	expect_warning(AL$Y <<- DR_data(ArcticLake[, 1:3]), ".*normalization forced.*")
	expect_equal(unname(rowSums(AL$Y)), rep(1.0, 39L))
})

cat("\n--- Common Model Checks - DirichletReg vs. Mathematica")

context("AL: Common - Null Model ( Y ~ 1 )")

resC1 <- DirichReg(Y ~ 1, AL)

resC1_mathematica <- list(
	MLE     =  39.52929411370061365530,
	DEV     = -79.05858822740122731060,
	COEFS   = c(0.02097861493616783052574, 0.8408687713317540853724, 0.2613372419951986613011),
	SE      = c(0.1662656209409404954810, 0.1726406312466813760947, 0.1682475856662426178722),
	Z       = c(0.1261753020103878320402, 4.870630773645980149831, 1.553289700772411471728),
	P       = c(0.8995931612401563064894, 1.112425458312201743960e-6, 0.1203539396684395661496),
	GRAD    = c(0e-64,0e-62,0e-63),
	HESSIAN = matrix(c(-55.10500066979587876697, 22.20582474400168179970, 12.43883086502252973858, 22.20582474400168179970, -62.21464112109048121996, 28.23926138550186271570, 12.43883086502252973858, 28.23926138550186271570, -58.89205189102017582447), 3L),
	VCOV    = matrix(c(0.02764425670687651061389, 0.01599938109376821138710, 0.01351070157568414676500, 0.01599938109376821138710, 0.02980478755725261791073, 0.01767095470989548347560, 0.01351070157568414676500, 0.01767095470989548347560, 0.02830725008251964840603),3)
)

test_that("Arctic Lake - Data Transformation", {
	expect_equal(resC1_mathematica$MLE, resC1$logLik)
	expect_equal(resC1_mathematica$DEV, -2.0*resC1$logLik)
	expect_equal(resC1_mathematica$COEFS, unname(resC1$coefficients), check.attributes = FALSE)
	expect_equal(resC1_mathematica$SE, unname(resC1$se), check.attributes = FALSE)
	expect_equal(resC1_mathematica$Z, unname(resC1$coefficients / resC1$se), check.attributes = FALSE)
	expect_equal(resC1_mathematica$P, 2*pnorm(-abs(unname(resC1$coefficients / resC1$se))), check.attributes = FALSE)
	expect_equal(resC1_mathematica$HESSIAN, unname(resC1$hessian), check.attributes = FALSE)
	expect_equal(resC1_mathematica$VCOV, unname(resC1$vcov), check.attributes = FALSE)
})
