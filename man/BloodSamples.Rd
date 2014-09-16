\encoding{UTF-8}
\name{BloodSamples}
\alias{BloodSamples}
\docType{data}
\title{Serum Protein Composition in Blood Samples}
\description{
These data (Aitchison, 2003) list blood samples' compositions of \emph{Albumin}, \emph{Pre-Albumin}, \emph{Globulin A}, and \emph{Globulin B} in relation to two types of diseases. 14 patients suffer from disease A, 16 from disease B and 6 are unclassified.
}
\usage{BloodSamples}
\format{
  A data frame with 36 observations on the following 6 variables.
  \describe{%
    \item{\code{Albumin}, \code{Pre.Albumin}, \code{Globulin.A}, \code{Globulin.B}}{the amounts of Albumin, Pre-Albumin, Globulin A, and Globulin B.}%
    \item{\code{Disease}}{diagnosis of disease \code{A}, \code{B}, or \code{NA} for unclassified observations.}%
    \item{\code{New}}{a factor indicating whether the observations are old and classified (\code{No}) or new and unclassified (\code{Yes}).}%
  }%
}%
\source{%
Aitchison, J. (2003). \emph{The Statistical Analysis of Compositional Data.} The Blackburn Press, Caldwell, NJ.%
}
\examples{
head(BloodSamples)
Bl <- DR_data(BloodSamples[,1:4])
summary(Bl)
}
\keyword{datasets}
