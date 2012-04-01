\name{Dirichlet}
\alias{rdirichlet}
\alias{ddirichlet}
\title{The Dirichlet Distribution}
\description{Density function and random number generation for the Dirichlet distribution}
\usage{
rdirichlet(n, alpha)

ddirichlet(x, alpha, log = FALSE, sum.up = FALSE)
}
\arguments{
  \item{n}{number of random observations to draw}
  \item{x}{a matrix containing observations}
  \item{alpha}{the Dirichlet distribution's parameters. Can be a vector (one set of parameters for all observations) or a matrix (a different set of parameters for each observation), see \dQuote{Details}}
  \item{log}{if \code{TRUE}, logarithmic densities are returned}
  \item{sum.up}{if \code{TRUE}, the (log-)likelihood is returned}
}
\details{
  The Dirichlet distribution is a multidimensional generalization of the Beta distribution where each dimension is governed by an alpha-parameter.
  Formally this is \eqn{\mathcal{D}(\alpha_i)=\left[\left.\Gamma(\sum_{i}\alpha_i)\right/\prod_i\Gamma(\alpha_i)\right]\prod_{i}y_i^{\alpha_i-1}}{D(\alpha)=[\Gamma(\sum\alpha)/\prod\Gamma(\alpha)]\prod y^(\alpha-1)}.
  Usually, \code{alpha} is a vector thus the same parameters will be used for all observations.
  If \code{alpha} is a matrix, a complete set of alpha-parameters must be supplied for each observation.
  
  \code{log} returns the logarithm of the densities (therefore the log-likelihood) and \code{sum.up} returns the product or sum and thereby the likelihood or log-likelihood.
}
\value{
\item{rdirichlet}{returns a matrix with random numbers according to the supplied alpha vector or matrix.}
\item{ddirichlet}{returns a vector of densities (if \code{sum.up = FALSE}) or the (log-)likelihood (if \code{sum.up = TRUE}) for the given data and alphas.}
}
\author{
  Marco Maier
}
\examples{
X1 <- rdirichlet(100, c(5, 5, 10))

a.mat <- cbind(1:10, 5, 10:1)
a.mat
X2 <- rdirichlet(10, a.mat)
# note how the probabilities in the first an last column relate to a.mat
round(X2, 2)

ddirichlet(X1, c(5, 5, 10))
ddirichlet(X2, a.mat)
}
