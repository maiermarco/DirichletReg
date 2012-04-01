pkgname <- "DirichletReg"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('DirichletReg')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ArcticLake")
### * ArcticLake

flush(stderr()); flush(stdout())

### Name: ArcticLake
### Title: Arctic Lake Data (Aitchison)
### Aliases: ArcticLake
### Keywords: datasets

### ** Examples

head(ArcticLake)
AL <- DR_data(ArcticLake[,1:3])
plot(AL)
summary(AL)



cleanEx()
nameEx("BloodSamples")
### * BloodSamples

flush(stderr()); flush(stdout())

### Name: BloodSamples
### Title: Serum Protein Composition In Blood Samples
### Aliases: BloodSamples
### Keywords: datasets

### ** Examples

head(BloodSamples)
Bl <- DR_data(BloodSamples[,1:4])
summary(Bl)



cleanEx()
nameEx("DR_data")
### * DR_data

flush(stderr()); flush(stdout())

### Name: DR_data
### Title: Prepare Compositional Data
### Aliases: DR_data print.DirichletRegData summary.DirichletRegData
### Keywords: manip

### ** Examples

# create a DirichletRegData object from the Arctic Lake data
head(ArcticLake[, 1:3])
AL <- DR_data(ArcticLake[, 1:3])
summary(AL)
head(AL)



cleanEx()
nameEx("Dirichlet")
### * Dirichlet

flush(stderr()); flush(stdout())

### Name: Dirichlet
### Title: The Dirichlet Distribution
### Aliases: rdirichlet ddirichlet

### ** Examples

X1 <- rdirichlet(100, c(5, 5, 10))

a.mat <- cbind(1:10, 5, 10:1)
a.mat
X2 <- rdirichlet(10, a.mat)
# note how the probabilities in the first an last column relate to a.mat
round(X2, 2)

ddirichlet(X1, c(5, 5, 10))
ddirichlet(X2, a.mat)



cleanEx()
nameEx("DirichletReg-package")
### * DirichletReg-package

flush(stderr()); flush(stdout())

### Name: DirichletReg-package
### Title: The 'DirichletReg' Package
### Aliases: DirichletReg DirichletReg
### Keywords: package

### ** Examples

  example(plot.DirichletRegData)
  example(DirichReg)



cleanEx()
nameEx("DirichletRegModel")
### * DirichletRegModel

flush(stderr()); flush(stdout())

### Name: DirichletRegModel
### Title: Methods for the class 'DirichletRegModel'
### Aliases: print.DirichletRegModel summary.DirichletRegModel
###   fitted.DirichletRegModel predict.DirichletRegModel
###   residuals.DirichletRegModel logLik.DirichletRegModel
###   AIC.DirichletRegModel BIC.DirichletRegModel nobs.DirichletRegModel
###   vcov.DirichletRegModel update.DirichletRegModel
###   confint.DirichletRegModel print.DirichletRegConfint

### ** Examples

  ALake <- ArcticLake
  ALake$AL <- DR_data(ArcticLake[, 1:3])
  mod1 <- DirichReg(AL ~ depth + I(depth^2) | depth, data = ALake, model="alternative")
  update(mod1, . ~ . | . + I(depth^2), evaluate = FALSE)
  mod1
  summary(mod1)
  head(fitted(mod1))
  head(predict(mod1))
  head(residuals(mod1))
  confint(mod1)
  logLik(mod1)
  vcov(mod1)



cleanEx()
nameEx("anova.DirichletRegModel")
### * anova.DirichletRegModel

flush(stderr()); flush(stdout())

### Name: anova.DirichletRegModel
### Title: Compare Dirichlet Regression Models using an LRT
### Aliases: anova.DirichletRegModel

### ** Examples

ALake <- ArcticLake
ALake$AL <- DR_data(ArcticLake[,1:3])
mod0 <- DirichReg(AL ~ 1, ALake)
mod1 <- DirichReg(AL ~ depth, ALake)
mod2 <- DirichReg(AL ~ depth + I(depth^2), ALake)
anova(mod1, mod0, mod2, sorted = TRUE)



cleanEx()
nameEx("dirichreg")
### * dirichreg

flush(stderr()); flush(stdout())

### Name: DirichReg
### Title: Fitting a Dirichlet Regression
### Aliases: DirichReg
### Keywords: multivariate models regression

### ** Examples

ALake <- ArcticLake
ALake$Y <- DR_data(ALake[,1:3])

# fit a quadratic Dirichlet regression models ("common")
res1 <- DirichReg(Y ~ depth + I(depth^2), ALake)

# fit a Dirichlet regression with quadratic predictor for the mean and
# a linear predictor for precision ("alternative")
res2 <- DirichReg(Y ~ depth + I(depth^2) | depth, ALake, model="alternative")

# test both models
anova(res1, res2)

res1
summary(res2)



cleanEx()
nameEx("plot.DirichletRegData")
### * plot.DirichletRegData

flush(stderr()); flush(stdout())

### Name: plot.DirichletRegData
### Title: Plot Dirichlet-Distributed Data
### Aliases: plot.DirichletRegData lines.DirichletRegData
### Keywords: hplot

### ** Examples

# plot of "Sand" in the Arctic Lake data set
plot(DR_data(ReadingSkills[, 1]), main="Reading Accuracy")

# ternary plot of Arctic Lake data
plot(DR_data(ArcticLake[, 1:3]), a2d = list(colored = FALSE))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
