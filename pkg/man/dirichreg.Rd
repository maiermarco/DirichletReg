\encoding{UTF-8}
\name{DirichReg}
\alias{DirichReg}
\title{Fitting a Dirichlet Regression}
\description{
  This function allows for fitting Dirichlet regression models using two different parametrizations.
}
\usage{
DirichReg(formula, data, model = c("common", "alternative"),
          subset, sub.comp, base, weights, control, verbosity = 0)
}
\arguments{
  \item{formula}{the model formula (for different specifications see \dQuote{Details})}
  \item{data}{a \code{data.frame} containing independent \strong{and} dependent variables}
  \item{model}{specifies whether the \code{"common"} (\eqn{\alpha\mathrm{s}}{alphas}) or \code{"alternative"} (\eqn{\mu/\phi}{mean/precision}) parametrization is employed (see \dQuote{Details})}
  \item{subset}{estimates the model for a subset of the data}
  \item{sub.comp}{analyze a subcomposition by selecting specific components (see \dQuote{Details})}
  \item{base}{redefine the base variable}
  \item{weights}{frequency weights}
  \item{control}{a list containing control parameters used for the optimization}
  \item{verbosity}{prints information about the function's progress, see Details}
}
\details{
  \subsection{Formula Specification and Models}{
    \code{formula} determines the used predictors.
    The responses \strong{must} be prepared by \code{\link{DR_data}} and can be optionally stored in the object containing all covariates which is then specified as the argument \code{data} (\dQuote{on-the-fly} processing is not implemented yet).
    
    There are two different parametrization (controlled by the argument \code{model}, see below):
    \itemize{
      \item the \emph{\dQuote{common}} param. that models each \eqn{\alpha}{alpha} by an (possibly individual) set of predictors, and
      \item the \emph{\dQuote{alternative}} param. that models expected values (\eqn{\mu}{mu}; as in multinomial logistic regression) and precision parameters (\eqn{\phi}{phi}) with two sets of predictors.
    }
    
    As the two models offer different modeling strategies, the specification of their formulae differ:
    
    \subsection{Formulae for the \dQuote{Common} Model}{
      The simplest possible model here is to include only an intercept for all components.
      If \code{DV} is the \emph{\sQuote{dependent variable}} (i.e., compositional data) with three components, we can request this null-model by \code{DV ~ 1}.
      We always have at least two dependent variables, so simple formulae as the one given above will be expanded to \code{DV ~ 1 | 1 | 1}, because \code{DV} hast three components.
      Likewise, it is possible to specify a common set of predictors for all components, as in \code{DV ~ p1 * p2}, where \code{p1} and \code{p2} are predictors.
      
      If the covariates of the components shall differ, one has to set up a complete formula for each subcomposition, using \code{|} as separators between the components, for example, \code{DV ~ p1 | p1 + p2 | p1 * p2} will lead to a model where the first response in \code{DV} will be modeled using \code{p1}, the second will be predicted by \code{p1 + p2} and the third by \code{p1 * p2}.
      Note that if you use the latter approach, the predictors have to be stated
      explicitly for all response variables.
    }
    \subsection{Formulae for the \dQuote{Alternative} Model}{
      The simplest possible model here is to include an intercept for all components (except the base) and an intercept for precision. This can be achieved by \code{DV ~ 1}, which is expanded to \code{DV ~ 1 | 1}. The part modeling the \sQuote{mean} (first element on the right-hand side) is mandatory, if no specification for precision is included, an intercept will be added. Note that you need to set \code{model = "alternative"} to use this parametrization!
      
      The alternative parametrization consists of two parts: modeled expected values (\eqn{\mu}{mu}) and their \sQuote{precision} (\eqn{\phi}{phi}). 
      As in multinomial logistic regression, one response variable is omitted (by default the first, but this can be changed by the \code{base} argument in \code{\link{DR_data}} or \code{DirichReg}) and for the rest a set of predictors is used with a multinomial logit-link.
      For precisions, a different set of predictors can be set up using a log-link.
      
      \code{DV ~ p1 * p2 | p1 + p2} will set up a model where the expected values are predicted by \code{p1 * p2} and precision are modeled using \code{p1 + p2}.
    }
  }
  \subsection{Data Preparation}{
    The \code{data} argument accepts a \code{data.frame} that \strong{must} include the dependent variable as a named element (see examples how to do this).
  }
  \subsection{Changing the Base Component and Analyzing Subcompositions}{
    The base-component (i.e., omitted component) is initially set during the stage of data preparation \code{\link{DR_data}}, but can easily be changed using the argument \code{base} which takes integer values from 1 to the maximum number of components.
    
    If a data set contains a large number of components, of which only a few are relevant, the latter can be \sQuote{sorted out} and the irrelevant (i.e., not selected) components will be aggregated into a single variable (row sums) that automatically becomes the base category for the model, unless specified otherwise by \code{base}. The positioning of variables will necessarily change: the aggregated variable takes the first column and the others are appended in their order of selection.
  }
  \subsection{Subsets and Weights}{
    Using \code{subset}, the model can be fitted only to a part of the data, for more information about this functionality, see \code{\link{subset}}.\cr
    Note that, unlike in \code{\link{glm}}, \code{weights} are \strong{not} treated as prior weights, but as frequency weights!
  }
  \subsection{Optimization and Verbosity}{
    Using the \code{control} argument, the settings passed to the optimizers can be altered.
    This argument takes a named list.
    To supply user-defined starting values, use \code{control = list(sv=c(...))} and supply a vector containing initial values for all parameters.
    Optimizer-specific options include the number of iterations (\code{iterlim = 1000}) and convergence criteria for the BFGS- and NR-optimization ((\code{tol1 = 1e-5}) and (\code{tol2 = 1e-10})).
    
    Verbosity takes integer values from \code{0} to \code{4}.
    \code{0}, no information is printed (default).
    \code{1} prints information about 3 stages (preparation, starting values, estimation).
    \code{2} prints little information about optimization (\code{verbosity} values greater than one are passed to \code{print.default = verbosity - 1} of \code{\link[maxLik]{maxBFGS}} and \code{\link[maxLik]{maxNR}}).
    \code{3} prints more information about optimization.
    \code{4} prints all information about optimization.
  }
}
\value{
\item{call}{[\code{language}] function call}
\item{parametrization}{[\code{character}] used parametrization}
\item{varnames}{[\code{character}] components' names}
\item{n.vars}{[\code{numeric}] vector with the number of parameters per set of predictors}
\item{dims}{[\code{numeric}] number of components}
\item{Y}{[\code{numeric}] used components}
\item{X}{[\code{numeric list}] sets of predictors}
\item{Z}{[\code{numeric list}] sets of predictors (only for the alternative parametrization)}
%\item{subset}{[\code{logical}] a vector determining in/exclusion of observations}
\item{sub.comp}{[\code{numeric}] vector of single components}
\item{base}{[\code{numeric}] base (only for the alternative parametrization)}
\item{weights}{[\code{numeric}] vector of frequency weights}
\item{orig.resp}{[\code{DirichletRegData}] the original response}
\item{data}{[\code{data.frame}] original data}
\item{d}{[\code{data.frame}] used data}
\item{formula}{[\code{Formula}] expanded formula}
\item{mf_formula}{[\code{language}] expression for generating the model frame}
\item{npar}{[\code{numeric}] number of parameters}
\item{coefficients}{[\code{numeric}] named vector of parameters}
\item{coefnames}{[\code{character}] names of the parameters}
\item{fitted.values}{[\code{list of matrices}] list containing alpha's, mu's, phi's for the observations}
\item{logLik}{[\code{numeric}] the log-likelihood}
\item{vcov}{[\code{matrix}] covariance-matrix of parameter estimates}
\item{hessian}{[\code{matrix}] (observed) Hessian}
\item{se}{[\code{numeric}] vector of standard errors}
\item{optimization}{[\code{list}] contains details about the optimization process provided by \code{\link[maxLik]{maxBFGS}} and \code{\link[maxLik]{maxNR}}}
}
\author{Marco J. Maier}
\examples{
ALake <- ArcticLake
ALake$Y <- DR_data(ALake[,1:3])

# fit a quadratic Dirichlet regression models ("common")
res1 <- DirichReg(Y ~ depth + I(depth^2), ALake)

# fit a Dirichlet regression with quadratic predictor for the mean and
# a linear predictor for precision ("alternative")
res2 <- DirichReg(Y ~ depth + I(depth^2) | depth, ALake, model="alternative")

# test both models
anova(res1, res2)

res1
summary(res2)
}
\keyword{multivariate}
\keyword{models}
\keyword{regression}