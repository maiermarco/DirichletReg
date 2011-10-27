\name{plot.DirichletRegData}
\alias{plot.DirichletRegData}
\alias{lines.DirichletRegData}
\title{Plot Dirichlet-Distributed Data}
\description{
With this function you can plot Dirichlet-distributed data in 2, 3 and 4 dimensions.
}
\usage{
\method{plot}{DirichletRegData}(x, dims, c.grid = TRUE, ticks = TRUE,
        colored = TRUE, ref.lines = NULL, col.scheme =
        c("dims", "entropy"), entropy.contours = FALSE,
        entropy.colors = FALSE, dim.labels, args.3d = list(rgl
        = TRUE, ...), rug = T, reset_par = TRUE, ...)
%\method{lines}{DirichletRegData}(x, indep, groups, \ldots, cumul = FALSE, orig.scale = FALSE)
}
\arguments{
  \item{x}{data prepared with \code{\link{DR_data}}}
  \item{dims}{select two, three, or four Dimensions of your data \code{x} to plot}
  \item{c.grid}{display a grid?}
  \item{ticks}{display ticks?}
  \item{colored}{colored plot?}
  \item{col.scheme}{coloring scheme can be \code{dims} or \code{entropy}}
  \item{entropy.contours}{plot entropy contours in a ternary plot?}
  \item{entropy.colors}{.}
  \item{rug}{display a rug?}
  \item{reset_par}{reset graphical parameters after creating the plot?}
  \item{ref.lines}{.}
  \item{dim.labels}{a character vector giving labels for the dimensions/variables}
  \item{args.3d}{a named list of settings for 3-dimensional plots}
  \item{\ldots}{further arguments}

%  \item{indep}{.}
%  \item{groups}{.}
%  \item{cumul}{.}
%  \item{orig.scale}{.}
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
# ternary plot of Arctic Lake data
plot(DR_data(ArcticLake[, 1]))
plot(DR_data(ArcticLake[, 1:3]))
plot(DR_data(Rocks[, 1:4]))
}
\keyword{hplot}