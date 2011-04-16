\name{DirichReg}
\alias{DirichReg}
\title{Fitting a Dirichlet Regression}
\description{See above}
\usage{
DirichReg(formula, data, weights, control, verbosity = 0, \ldots)
}
\arguments{
  \item{formula}{model formula, for different specifications see details}
  \item{data}{data}
  \item{weights}{vector of weights}
  \item{control}{a list containing control parameter used for the optimization}
  \item{verbosity}{information about the function's progress}
  \item{\ldots}{further arguments}
}
\details{
  \subsection{Formula Specification}{
    \code{formula} determines the used predictors \emph{and} the parameterization
    that will be employed. The responses \strong{must} be prepared by \link{DR.data}.
    
    There are two different parameterization: the \emph{\sQuote{common}} param. that
    models each alpha by an individual set of predictors and the \emph{\sQuote{alternative}}
    param. that models expected values (as in multinomial logistic regression) and
    precision parameters with two sets of predictors.
    
    \itemize{
    \item To set up a commonly parameterized model, type \code{DV ~ p1 + p2}, where \code{DV}
      is the dependent variable prepared by \link{DR.data} and \code{p1} and \code{p1}
      are predictors. The employed predictors will be the same for every dependent variable.
      
      If you want to use a different set for the dependent variables you can separate
      the sets by \code{|}s like this: \code{DV ~ p1 | p1 + p2 | p1 * p2}.
      This means that the first response in \code{DV} will be modeled using \code{p1},
      the second will be predicted by \code{p1 + p2} and so on.

      Note that if you use the latter approach, the predictors have to be stated
      explicitly for all response variables.
    \item The alternative parameterization works consists of two parts: the
      expected values and their \sQuote{precision}. As in multinomial logistic
      regression, one response variable is omitted (by default the first, but this
      can be changed by the \code{base} argument in \link{DR.data}) and for the
      rest a set of predictors is used with a logit-link. For precisions, a different
      set of predictors can be set up using a log-link.
      
      \code{DV ~ p1 * p2 | phi ~ p1 + p2} will set up a model where the expected
      values are predicted by \code{p1 * p2} and precision are modeled using
      \code{p1 + p2}. Note that the second part of the formula \strong{must} be
      \code{phi ~} to set up the alterantive parameterization.
    }  
  }
  Using the \code{control} argument, the settings passed to \code{optim} can be
  altered. To alter starting values for example the code would be \code{control = list(sv=c(1,2,3, ...))}
  where \code{c(1,2,3, ...)} contains the appropriate starting values to be used.

}
\value{
\item{all the values}{description}
}
%\references{
%% ~put references to the literature/web site here ~
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

AL <- DR.data(ArcticLake[,1:3])

res <- DirichReg(AL ~ depth + I(depth^2), ArcticLake)

res
summary(res)
}
\keyword{multivariate}
\keyword{models}
\keyword{regression}