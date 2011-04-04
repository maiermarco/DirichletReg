\name{plot.DirichletRegData}
\alias{plot.DirichletRegData}
\alias{lines.DirichletRegData}
\title{Plot Dirichlet-Distributed Data}
\description{
With this function you can plot Dirichlet-distributed data in 2, 3 and 4 dimensions.
}
\usage{
\method{plot}{DirichletRegData}(x, \ldots, dims, explore, model, c.grid = TRUE, ticks =
                 TRUE, colored = TRUE, col.scheme = "classes",
                 entropy.contours = FALSE, entropy.colors = FALSE,
                 pt.size = 0.75, dim.labels, args.3d = list(rgl = TRUE,
                 theta = NULL, phi = NULL, ref.lines = NULL, ...))

\method{lines}{DirichletRegData}(x, indep, groups, \ldots, cumul = FALSE, orig.scale = FALSE)
}
\arguments{
  \item{x}{the data prepared with \code{DR.data}}
  \item{\ldots}{further arguments}
  \item{dims}{.}
  \item{explore}{.}
  \item{model}{.}
  \item{c.grid}{.}
  \item{ticks}{.}
  \item{colored}{.}
  \item{col.scheme}{.}
  \item{entropy.contours}{.}
  \item{entropy.colors}{.}
  \item{pt.size}{.}
  \item{dim.labels}{a character vector giving labels for the dimensions/variables}
  \item{args.3d}{a named list of settings for 3-dimensional plots}

  \item{indep}{.}
  \item{groups}{.}
  \item{cumul}{.}
  \item{orig.scale}{.}
}
%\details{
%No details yet.
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

plot(DR.data(ArcticLake[, 1:3]))
}
\keyword{hplot}