\encoding{UTF-8}
\name{Simplex-Transformations}
\alias{toQuaternary}
\alias{toQuaternaryVectors}
\alias{toTernary}
\alias{toTernaryVectors}
\alias{toSimplex}
\title{Transform Compositional Data for a Simplex}
\description{These functions transform a matrix with three or four components to fit into a two- or three-dimensional simplex (triangle or tetrahedron).}
\usage{toSimplex(x)

toTernary(abc)
toTernaryVectors(c1, c2, c3)

toQuaternary(abcd)
toQuaternaryVectors(c1, c2, c3, c4)
}
\arguments{
  \item{x}{a matrix-like object with 3 or 4 columns.}

  \item{abc}{a matrix-like object with 3 columns.}
  \item{abcd}{a matrix-like object with 4 columns.}

  \item{c1}{a numeric vector with values of the first component.}
  \item{c2}{a numeric vector with values of the second component.}
  \item{c3}{a numeric vector with values of the third component.}
  \item{c4}{a numeric vector with values of the fourth component.}
}
\details{Most of these functions are only used internally, but sometimes it might be useful to plot \dQuote{custom} ternary or quaternary graphics.

Note that, apart from \code{toSimplex()}, functions do not have \emph{any} checks, so it is advisable to use this function if elements are added to plots or own graphics are created.
}
\note{In prior versions (up to 0.5-0), an unexported function \code{coord.trafo()} was used internally and could also be accessed via \code{DirichletReg:::coord.trafo()}.

If you have used this in your code, you will get a message that the function is now deprecated and will become defunct in the future.
Use \code{toSimplex()} instead.
}
\value{The function returns a \code{matrix} object with coordinates in two or three dimensions}
\author{Marco J. Maier}
\examples{
# create a DirichletRegData object from the Arctic Lake data
"to be added"
}
\keyword{manip}
