\name{Dirichlet}
\alias{rdirichlet}
\alias{ddirichlet}
\title{The Dirichlet Distribution}
\description{Density and random generation for the Dirichlet distribution.}
\usage{
rdirichlet(n, alpha)

ddirichlet(x, alpha, log = FALSE, sum.up = FALSE)
}
\arguments{
  \item{n}{number of random observations to draw.}
  \item{x}{a matrix containing observations.}
  \item{alpha}{the Dirichlet distribution's parameters. Can be a vector (one set of parameters for all observations) or a matrix (a different set of parameters for each observation), see details.}
  \item{log}{if TRUE, logarithmic densities are returned.}
  \item{sum.up}{if TRUE, the (log-)likelihood is returned.}
}
\details{
  The Dirichlet distribution is a multidimensional generalization of the Beta distribution where each dimension is governed by an alpha-parameter.
  Usually, \code{alpha} is a vector thus the same parameters will be used for all observations.
  If \code{alpha} is a matrix, a complete set of alpha-parameters must be supplied for each observation.
}
\value{
\item{rdirichlet}{returns a matrix with random numbers according to the supplied alpha vector or matrix.}
\item{ddirichlet}{returns a vector of densities for the given data and alphas.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
X1 <- rdirichlet(100, c(5, 5, 10))

a.mat <- cbind(1:10, 5, 10:1)
a.mat
X2 <- rdirichlet(10, a.mat)

round(X2, 2)

ddirichlet(X1, c(5, 5, 10))
ddirichlet(X2, a.mat)
}
\keyword{distribution}
\keyword{datagen}
