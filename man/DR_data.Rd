\encoding{UTF-8}
\name{DirichletRegData}
\alias{DirichletRegData}
\alias{DR_data}
\alias{print.DirichletRegData}
\alias{summary.DirichletRegData}
%
%
\title{Prepare Compositional Data}
%
\description{This function prepares a matrix with compositional variables for further processing in the \pkg{DirichletReg} package.}
%
%
%
\usage{
DR_data(Y, trafo = sqrt(.Machine$double.eps), base = 1,
    norm_tol = sqrt(.Machine$double.eps), fix_names = TRUE)

\method{print}{DirichletRegData}(x, type = c("processed", "original"), \dots)

\method{summary}{DirichletRegData}(object, \dots)
}
%
%
%
\arguments{%
  \item{Y}{%
    A \code{matrix} or \code{data.frame} with nonnegative values of all compositional variables (in some cases, a vector is also permissible, see \dQuote{Details}).
  }%
  \item{trafo}{%
    Either a small numeric value in (0, .1] or \code{TRUE} (default: \ifelse{latex}{\eqn{\sqrt{\mathtt{.Machine\$double.eps}}\approx1.49\cdot10^8}}{\code{sqrt(.Machine$double.eps) = 1.49e-08}}).
    Transformation of variables causes the values to shrink away from extreme values of 0 and 1, see \dQuote{Details}.
    \cr%
    \code{TRUE} will force transformation, \code{FALSE} results in an error.
    \cr%
    If \code{trafo} is numeric it is used as a \dQuote{threshold}, so transformation will be applied if values in \code{Y} are \eqn{y < \mathtt{trafo}}{y < trafo} or \eqn{y > \left(1-\mathtt{trafo}\right)}{y > (1 - trafo)}.
  }%
  \item{base}{%
    The \dQuote{base} component to use in the reparametrized model
  }%
  \item{norm_tol}{%
    Due to numerical precision, row sums of \eqn{\mathbf{Y}} may not be \emph{exactly} equal to 1.
    Therefore, \code{norm_tol} is a small value in (0, .1] which represents the tolerance when testing for \dQuote{near equality} to 1 (default: \ifelse{latex}{\eqn{\sqrt{\mathtt{.Machine\$double.eps}}\approx1.49\cdot10^8}}{\code{sqrt(.Machine$double.eps) = 1.49e-08}}, see \code{\link[base]{all.equal}}).
  }%
  \item{fix_names}{%
    A logical value that determines whether \code{DR_data()} will try to infer the variable name if \code{Y} is a single column or vector.
    By default, \code{DR_data} will attempt to find a name, if \code{fix_names = FALSE}, whatever is supplied to \code{Y} will be used as is.
  }%
  \item{x}{%
    A \code{DirichletRegData} object
  }%
  \item{type}{%
    Displays either the (possibly normalized or transformed) \code{"processed"} or \code{"original"} data
  }%
  \item{object}{%
    A \code{DirichletRegData} object
  }%
  \item{\dots}{%
    Further arguments
  }%
}
\details{
  \subsection{\code{Y}}{%
    \code{Y} is a \code{matrix} or \code{data.frame} containing compositional variables.
    If they do not sum up to 1 for all observations, normalization is forced where each row entry is divided by the row's sum (a warning will be issued that normalization was applied).\cr
    In case one row-entry (or more) is \code{NA}, the whole row will be returned as \code{NA}.
    Beta-distributed variables can be supplied as a single vector which, however, has to have values in the interval \eqn{[0,\,1]}{[0, 1]}.
    The second variable will be generated (\code{1 - Y}) and a \code{matrix} consisting of the columns \code{1 - Y} and \code{Y} will be returned.
    A message will be issued that a beta-distributed variable was assumed and that this assumption needs to be checked.
  }
  \subsection{\code{trafo}}{%
    The transformation (done if \code{trafo = TRUE}) is a generalization of that proposed by Smithson and Verkuilen (2006) that transforms each component \eqn{y}{y} of \eqn{Y}{Y} by computing \eqn{y^{*}=\frac{y(n-1)+\frac{1}{2}}{n}}{y*=[y(n-1)+1/2]/n} where \eqn{n}{n} is the number of observations in \eqn{Y}{Y} (this approach is also used in the package \pkg{betareg}, see Cribari-Neto & Zeileis, 2010).\cr
    For an arbitrary number of dimensions (or variables) \eqn{d}{d} the transformation is \eqn{y^{*}=\frac{y(n-1)+\frac{1}{d}}{n}}{y*=[y(n-1)+1/d]/n}.
  }
  \subsection{\code{base}}{%
    To set the base (i.e., omitted) component of \code{Y} for the \dQuote{alternative} (mean/precision) model, the argument \code{base} can be used. This is by default set to the first variable in \code{Y} (if a vector is be supplied, the column \code{1 - Y} becomes the base component).\cr
    Note that the definition can be overruled in \code{\link{DirichReg}}.
  }
  \subsection{\code{x} and \code{object}}{%
    Objects created by \code{DR_data}.
  }
  \subsection{\code{type}}{%
    specifies for the print method whether the original or processed data are displayed.
  }
}
\value{The function returns a \code{matrix} object of class \code{DirichletRegData} with the following attributes:
  \item{attr(*, "dimnames")}{a list with two entries, row names (by default \code{NULL}) and column names.}
  \item{attr(*, "Y.original")}{the original data}
  \item{attr(*, "dims")}{number of dimensions of \code{Y} (i.e., number of columns)}
  \item{attr(*, "dim.names")}{the number of components in \code{Y}}
  \item{attr(*, "obs")}{number of observations of \code{Y} (i.e., number of rows)}
  \item{attr(*, "valid_obs")}{number of valid observations}
  \item{attr(*, "normalized")}{a logical value indicating whether the data were normalized}
  \item{attr(*, "transformed")}{a logical value indicating whether the data were transformed}
  \item{attr(*, "base")}{number of the variable used as the base in the reparametrized model}
}
\references{
Smithson, M. & Verkuilen, J. (2006). A Better Lemon Squeezer? Maximum-Likelihood Regression With Beta-Distributed Dependent Variables. \emph{Psychological Methods, 11}(1), 54--71. \doi{10.1037/1082-989X.11.1.54}

Cribari-Neto, F. & Zeileis, A. (2010). Beta Regression in \R. \emph{Journal of Statistical Software, 34}(2), 1--24. \doi{10.18637/jss.v034.i02}
}
\author{Marco J. Maier}
\examples{
# inspect the raw data
head(ArcticLake[, 1:3])

# create a DirichletRegData object from the Arctic Lake data
AL <- DR_data(ArcticLake[, 1:3])

# print a summary
summary(AL)

# here are the first few rows of the transformed data
head(AL)

# you can also supply a vector in [0, 1]
rand_data <- DR_data(runif(10))

# since no names were supplied, DR_data guessed the names
head(rand_data)

# note that here, only one variable is supplied, therefore the variable name
# would be ReadingSkills[["accuracy"]]
read_acc <- DR_data(ReadingSkills[["accuracy"]])

# however, since fix_names = TRUE by default, DR_data() infers that the variable name is
head(read_acc)

# accuracy and names the columns accordingly. To prevent this, set fix_names = FALSE
head(DR_data(ReadingSkills[["accuracy"]], fix_names = FALSE))
}
\keyword{manip}
