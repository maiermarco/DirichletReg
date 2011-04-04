\name{ArcticLake}
\alias{ArcticLake}
\docType{data}
\title{Arctic Lake Data (Aitchison)}
\description{Aitchison's `arctic lake' data set. These data show how the
ground's composition of a lake varies depending on the water depth.}
\usage{
data(ArcticLake)
}
\format{
  A data frame with 39 observations on the following 4 variables:
  \describe{
    \item{\code{sand}}{relative frequencies of sand}
    \item{\code{silt}}{relative frequencies of silt}
    \item{\code{clay}}{relative frequencies of clay}
    \item{\code{depth}}{water depth in meters}
  }
}
%\source{
%%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
head(ArcticLake)
AL <- DR.data(ArcticLake[,1:3])
AL
}
\keyword{datasets}