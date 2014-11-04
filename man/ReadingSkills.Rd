\encoding{UTF-8}
\name{Reading Accuracy Data}
\alias{ReadingSkills}
\alias{ReadDisc}
\docType{data}
\title{Pammer and Kevan's Data on Reading Skills}
\description{These data provide transformed reading accuracy scores predicted by \acronym{IQ} and diagnosed dyslexia.}
\usage{ReadingSkills
ReadDisc}
\format{
  A data frame containing 44 observations on 3 variables.
  \describe{
    \item{\code{accuracy}}{reading accuracy score transformed to fit into \eqn{(0,\,1)}{(0, 1)}}
    \item{\code{dyslexia}}{a factor with the diagnosis of dyslexia (\dQuote{yes} or \dQuote{no})}
    \item{\code{iq}}{non-verbal \acronym{IQ} (\eqn{z}{z}-scores; \eqn{\mu=0}{mean = 0}, \eqn{\sigma^2=1}{sigma² = 1})}
  }
}
\source{Example 3 from \url{http://psychology3.anu.edu.au/people/smithson/details/betareg/betareg.html}}
\keyword{datasets}
