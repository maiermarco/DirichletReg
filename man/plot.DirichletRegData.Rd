\encoding{UTF-8}
\name{plot.DirichletRegData}
\alias{plot.DirichletRegData}
\title{Plot Dirichlet-Distributed Data}
\description{
With this function you can plot Dirichlet-distributed data in 2, 3 and 4 dimensions.
}
\usage{
\method{plot}{DirichletRegData}(x, dims, ticks = TRUE, ref.lines = NULL, dim.labels, a2d = list(colored =
  TRUE, c.grid = TRUE, col.scheme = c("dims", "entropy"), entropy.contours =
  FALSE, entropy.colors = FALSE), a3d = list(rgl = TRUE, ...), rug = TRUE,
  reset_par = TRUE, ...)
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
# density plot of "accuracy" in the Reading Skills data set
plot(DR_data(ReadingSkills["accuracy"]), main = "Reading Accuracy")

# ternary plot of Arctic Lake data
plot(DR_data(ArcticLake[, 1:3]), a2d = list(colored = FALSE))

# 10,000 random four-dimensional data points generated with alpha = (2, 2, 2, 2)
set.seed(2222L)
random_4D <- rdirichlet(10000L, c(2, 2, 2, 2))
colnames(random_4D) <- sprintf("Dim. \%s", as.roman(1:4))
# plot the data using the static plot device (not the default "rgl")
plot(DR_data(random_4D), cex = 2/3, a3d = list(rgl = FALSE))
}
\keyword{hplot}
