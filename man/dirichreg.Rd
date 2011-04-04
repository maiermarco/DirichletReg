\name{DirichReg}
\alias{DirichReg}
\title{Fitting a Dirichlet Regression}
\description{See above}
\usage{
DirichReg(formula, data, weights, control, verbosity = 0, \ldots)
}
\arguments{
  \item{formula}{model formula}
  \item{data}{data}
  \item{weights}{vector of weights}
  \item{control}{a list containing control parameter used for the optimization}
  \item{verbosity}{information about the function's progress}
  \item{\ldots}{further arguments}
}
\details{
  Using the \code{control} argument, the settings passed to \code{optim} can be
  altered. To alter starting values for example the code would be \code{control = list(sv=c(1,2,3, ...))}
  where \code{c(1,2,3, ...)} contains the appropriate starting values to be used.
}
\value{
\item{all the values}{description}
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

AL <- DR.data(ArcticLake[,1:3])

res <- DirichReg(AL ~ depth + I(depth^2), ArcticLake)

res
summary(res)
}
\keyword{multivariate}
\keyword{models}
\keyword{regression}