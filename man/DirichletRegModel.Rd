\name{DirichletRegModel}
\alias{print.DirichletRegModel}
\alias{summary.DirichletRegModel}
\alias{fitted.DirichletRegModel}
\alias{predict.DirichletRegModel}
\alias{residuals.DirichletRegModel}
\title{Methods for the class \code{DirichletRegModel}}
\description{.
}
\usage{
\method{print}{DirichletRegModel}(x, digits = max(3, getOption("digits") - 3), \ldots)

\method{summary}{DirichletRegModel}(object, digits = max(3, getOption("digits") - 3), \ldots)

\method{fitted}{DirichletRegModel}(object, mu = TRUE, alpha = FALSE, phi = FALSE, \ldots)

\method{predict}{DirichletRegModel}(object, newdata, alpha = FALSE, \ldots)

\method{residuals}{DirichletRegModel}(object, type = c("score", "raw", "standardized", "composite"), \ldots)
}
\arguments{
  \item{x}{an object}
  \item{object}{an object}
  \item{alpha}{whether alphas or probabilities should be returned}
  \item{mu}{expected values}
  \item{phi}{precision}
  \item{type}{can return different types of residuals}
  \item{newdata}{supply new data for your predictions}
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

}