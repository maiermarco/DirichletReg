\name{DirichMixture}
\alias{DirichMixture}
\title{Fitting a Dirichlet Regression}
\description{See above}
\usage{
DirichMixture(Y, classes = 2, EMs = 5, verbosity = 0, ctl =
                 list(em.tol = 1e-08, op.tol = 1e-08))
}
\arguments{
  \item{Y}{a \code{DR_data()} prepared data object}
  \item{classes}{the number of latent classes to extract}
  \item{EMs}{the number of full EM-optimizations}
  \item{verbosity}{detailed information about the current status}
  \item{ctl}{a list of control parameters}
}
\details{
details?
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
RR <- DR_data(Rocks[,1:5])

drm4 <- DirichMixture(RR, classes=4, EMs=1)
drm5 <- DirichMixture(RR, classes=5, EMs=1)

drm4
#summary(drm5)

#anova(drm4, drm5)
}
\keyword{multivariate}
\keyword{models}
\keyword{regression}