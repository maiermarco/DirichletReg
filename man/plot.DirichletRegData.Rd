\encoding{UTF-8}
\name{plot.DirichletRegData}
\alias{plot.DirichletRegData}
\alias{lines.DirichletRegData}
\title{Plot Dirichlet-Distributed Data}
\description{
With this function you can plot Dirichlet-distributed data in 2, 3 and 4 dimensions.
}
\usage{
\method{plot}{DirichletRegData}(x, dims, ticks = TRUE, ref.lines = NULL, dim.labels, a2d = list(colored =
  TRUE, c.grid = TRUE, col.scheme = c("dims", "entropy"), entropy.contours =
  FALSE, entropy.colors = FALSE), a3d = list(rgl = TRUE, ...), rug = TRUE,
  reset_par = TRUE, ...)
%\method{lines}{DirichletRegData}(x, indep, groups, \dots, cumul = FALSE, orig.scale = FALSE)
}
\arguments{
  \item{x}{data prepared with \code{\link{DR_data}}}
  \item{dims}{select two, three, or four Dimensions of your data \code{x} to plot}
  \item{ticks}{display ticks?}
  \item{ref.lines}{.}
  \item{dim.labels}{a character vector giving labels for the dimensions/variables}
  \item{a2d}{a named list of settings for ternary plots (3 variables), see Details}
  \item{a3d}{a named list of settings for quaternary plots (4 variables), see Details}
  \item{rug}{display a rug for a one-dimensional plot (2 variables)}
  \item{reset_par}{reset graphical parameters of \code{\link{DR_data}} after creating a two-dimensional plot (2 variables), see Details}
  \item{\dots}{further graphical arguments as \code{col}, \code{pch}, \code{cex}, \ldots}
%  \item{indep}{.}
%  \item{groups}{.}
%  \item{cumul}{.}
%  \item{orig.scale}{.}
}
%\details{
%No details yet.
%}
\author{Marco J. Maier}
%\note{
%%  ~~further notes~~
%}
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# plot of "Sand" in the Arctic Lake data set
plot(DR_data(ReadingSkills[, 1]), main="Reading Accuracy")

# ternary plot of Arctic Lake data
plot(DR_data(ArcticLake[, 1:3]), a2d = list(colored = FALSE))
}
\keyword{hplot}