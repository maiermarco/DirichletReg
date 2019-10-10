\encoding{UTF-8}
\name{DirichletRegModel}
\alias{print.DirichletRegModel}
\alias{summary.DirichletRegModel}
\alias{fitted.DirichletRegModel}
\alias{predict.DirichletRegModel}
\alias{residuals.DirichletRegModel}
\alias{logLik.DirichletRegModel}
\alias{AIC.DirichletRegModel}
\alias{BIC.DirichletRegModel}
\alias{nobs.DirichletRegModel}
\alias{vcov.DirichletRegModel}
\alias{update.DirichletRegModel}
\alias{confint.DirichletRegModel}
\alias{drop1.DirichletRegModel}
\alias{print.DirichletRegConfint}
\title{Methods for the Class \code{DirichletRegModel}}
\description{
  These are available methods for the results of Dirichlet regression models and objects of class \code{DirichletRegModel}.
  These methods contain functions for \code{print} and \code{summary} of the data, generate \code{fitted} values and predicting new values using \code{predict}.
  Various types of \code{residuals} are implemented and \code{confint} can be used to compute confidence intervals of the parameters.
  Furthermore \code{logLik} extracts the log-likelihood of the model and \code{vcov} extracts the covariance matrix of the parameter estimates.
}
\usage{
\method{print}{DirichletRegModel}(x, digits = max(3, getOption("digits") - 3), \dots)

\method{summary}{DirichletRegModel}(object, \dots)

\method{fitted}{DirichletRegModel}(object, mu = TRUE, alpha = FALSE, phi = FALSE, \dots)

\method{predict}{DirichletRegModel}(object, newdata, mu = TRUE, alpha = FALSE, phi = FALSE, \dots)

\method{residuals}{DirichletRegModel}(object, type = c("standardized", "composite", "raw"), \dots)

\method{confint}{DirichletRegModel}(object, parm, level, \dots, type=c("all", "beta", "gamma"), exp = FALSE)

\method{print}{DirichletRegConfint}(x, digits = 3, \dots)

\method{logLik}{DirichletRegModel}(object, \dots)

\method{AIC}{DirichletRegModel}(object, \dots, k = 2)

\method{BIC}{DirichletRegModel}(object, \dots)

\method{nobs}{DirichletRegModel}(object, \dots)

\method{vcov}{DirichletRegModel}(object, \dots)

\method{update}{DirichletRegModel}(object, formula., \dots, evaluate = TRUE)

\method{drop1}{DirichletRegModel}(object, scope, test = c("LRT", "none"), k = 2, sort = TRUE, \dots)
}
\arguments{
  \item{x}{an object of class \code{DirichletRegModel}}
  \item{object}{an object of class \code{DirichletRegModel} or \code{DirichletRegConfint} for printing an object obtained by \code{confint.DirichletRegModel}}
  \item{alpha}{logical; returns alpha values}
  \item{mu}{logical; returns expected values}
  \item{phi}{logical; returns precision values}
  \item{type}{for \code{residuals}: defines the type of residuals to be computed \code{"standardized"} (i.e., Pearson), \code{"composite"}, or \code{"raw"}

              for \code{confint}: defines the type of parameter (\code{"all"}, \code{"beta"}, or \code{"gamma"}) for which confidence values are returned}
  \item{newdata}{a \code{data.frame} containing new observations}
  \item{k}{number for the weighting of parameters}
  \item{parm}{a vector containing names of the parameters to print}
  \item{level}{(a vector of) confidence level(s), defaults to \code{.95}}
  \item{exp}{logical; returns parameters in exponentiated form}
  \item{digits}{the number of digits in the output}
  \item{formula.}{the new formula to be updated, see \code{\link[=update.formula]{update.formula}} and \code{\link[=update.Formula]{update.Formula}}}
  \item{evaluate}{if \code{FALSE} the updated call will be returned, but not evaluated}
  \item{scope}{defines the scope of variables to be dropped, see \code{\link[=drop1]{drop1}}}
  \item{test}{defines the type of test for \code{drop1}}
  \item{sort}{if \code{TRUE}, p-values will be sorted in decreasing order.}
  \item{\dots}{further arguments}
}
%\details{
%.
%}
\author{Marco J. Maier}
\examples{
ALake <- ArcticLake
ALake$AL <- DR_data(ArcticLake[, 1:3])

mod1 <- DirichReg(AL ~ depth + I(depth^2) | depth, data = ALake, model="alternative")

update(mod1, . ~ . | . + I(depth^2), evaluate = FALSE)
mod1

drop1(mod1)   ### issues a caveat when used for the first time in an R session

summary(mod1)

head(fitted(mod1))

predict(mod1, newdata = data.frame("depth" = seq(10, 100, 10)))

head(residuals(mod1))

confint(mod1)
confint(mod1, exp = TRUE)

logLik(mod1)
round(vcov(mod1), 5)
}
