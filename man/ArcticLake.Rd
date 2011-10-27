\name{ArcticLake}
\alias{ArcticLake}
\docType{data}
\title{Arctic Lake Data (Aitchison)}
\description{
  These data are taken from Aitchison (2003) and contain information on the relation
  of sediment composition with depth in an Arctic lake.
}
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
\source{
  Aitchison, J.\ (2003). \emph{The Statistical Analysis of Compositional Data.}
  The Blackburn Press, Caldwell, NJ.
}
\examples{
head(ArcticLake)
AL <- DR_data(ArcticLake[,1:3])
plot(AL)
AL
}
\keyword{datasets}