\name{DirichletRegModel}
\alias{print.DirichletRegModel}
\alias{summary.DirichletRegModel}
\alias{fitted.DirichletRegModel}
\alias{predict.DirichletRegModel}
\alias{residuals.DirichletRegModel}
\alias{logLik.DirichletRegModel}
\alias{AIC.DirichletRegModel}
\alias{BIC.DirichletRegModel}
\alias{vcov.DirichletRegModel}
\alias{confint.DirichletRegModel}
\alias{print.DirichletRegConfint}
\title{Methods for the class \code{DirichletRegModel}}
\description{.
}
\usage{
\method{print}{DirichletRegModel}(x, digits = max(3, getOption("digits") - 3), \ldots)

\method{summary}{DirichletRegModel}(object, digits = max(3, getOption("digits") - 3), \ldots)

\method{fitted}{DirichletRegModel}(object, mu = TRUE, alpha = FALSE, phi = FALSE, \ldots)

\method{predict}{DirichletRegModel}(object, newdata, mu = TRUE, alpha = FALSE, phi = FALSE, \ldots)

\method{residuals}{DirichletRegModel}(object, type = c("standardized", "composite", "raw"), \ldots)

\method{confint}{DirichletRegModel}(object, parm, level, \ldots, type=c("all", "beta", "gamma"), e = FALSE)

\method{print}{DirichletRegConfint}(x, digits = 3, \ldots)

\method{logLik}{DirichletRegModel}(object, \ldots)

%\method{AIC}{DirichletRegModel}(object, \ldots, k = 2)

%\method{BIC}{DirichletRegModel}(object, \ldots)

\method{vcov}{DirichletRegModel}(object, \ldots)
}
\arguments{
  \item{x}{an object of class \code{DirichletRegModel}}
  \item{object}{an object of class \code{DirichletRegModel} or \code{DirichletRegConfint} for printing an object obtained by \code{confint.DirichletRegModel}}
  \item{alpha}{logical; returns alpha values}
  \item{mu}{logical; returns expected values}
  \item{phi}{logical; returns precision values}
  \item{type}{for \code{residuals}: defines the type of residuals to be computed \code{"standardized"} (i.e. Pearson), \code{"composite"}, or \code{"raw"}
  
              for \code{confint}: defines the type of parameter (\code{"all"}, \code{"beta"}, or \code{"gamma"}) for which confidence values are returned}
  \item{newdata}{a \code{data.frame} containing new observations}
%  \item{k}{number for the weighting of parameters}
  \item{parm}{a vector containing names of the parameters to print}
  \item{level}{(a vector of) confidence level(s), defaults to \code{.95}}
  \item{e}{logical; returns parameters in exponentiated form}
  \item{digits}{the number of digits in the output}
  \item{\ldots}{further arguments}
}
\details{
.
}
\author{
  Marco Maier
}
%\note{
%%  ~~further notes~~
%}
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  AL <- DR_data(ArcticLake[, 1:3])
  modl <- DirichReg(AL ~ depth + I(depth^2) | phi ~ depth, data = ArcticLake)
  modl
  summary(modl)
  head(fitted(modl))
  head(predict(modl))
  head(residuals(modl))
  confint(modl)
  logLik(modl)
  vcov(modl)
}