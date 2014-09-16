\encoding{UTF-8}
\name{ArcticLake}
\alias{ArcticLake}
\docType{data}
\title{Arctic Lake Data (Aitchison)}
\description{
These data are taken from Aitchison (2003) and contain information on the relation of sediment composition with depth in an Arctic lake.
}
\usage{ArcticLake}
\format{
  A data frame with 39 observations on the following 4 variables:
  \describe{%
    \item{\code{sand}, \code{silt}, \code{clay}}{relative frequencies of sand, silt, and clay}%
    \item{\code{depth}}{water depth in meters}%
  }%
}%
\source{
Aitchison, J. (2003). \emph{The Statistical Analysis of Compositional Data.} The Blackburn Press, Caldwell, NJ.
}
\examples{
head(ArcticLake)
AL <- DR_data(ArcticLake[,1:3])
plot(AL)
summary(AL)
}
\keyword{datasets}