\encoding{UTF-8}
\name{anova.DirichletRegModel}
\alias{anova.DirichletRegModel}
\title{Compare Dirichlet Regression Models using an LRT}
\description{
This function allows for pairwise tests of Dirichlet regression models using a likelihood ratio test (\acronym{LRT}).
}
\usage{
\method{anova}{DirichletRegModel}(object, \dots, sorted = FALSE)
}
\arguments{
\item{object}{the model to be compared against those listed in \dots}
\item{\dots}{models to be tested against the one specified as \code{object}}
\item{sorted}{should the models be sorted according to their number of parameters?}
}
\details{
The test statistic is computed \eqn{LR=-2\left[\log(L_a)-\log(L_b)\right]}{LR = -2 * [log(L_a) - log(L_b)]} where \eqn{L_i}{L_i} is the likelihood of model \eqn{i}{i} with \eqn{df}{df} equal to the difference of the number of parameters in the models.
}
\value{%
  An object of class \code{"anova_DirichletRegModel"} is a list containing these components:
  
  \item{Deviance}{a vector of deviances (\eqn{-2L_i}{-2 * L_i}) of the models}
  \item{N. par}{a vector with the numbers of parameters}
  \item{Difference}{the likelihood ratios}
  \item{df}{the degrees of freedom}
  \item{Pr(>Chi)}{the \eqn{p}{p}-values}
  \item{sorting}{sorting indices}
  \item{n.mods}{number of models}
  \item{calls}{a list containing the calls of the tested models}
}%
\author{Marco J. Maier}
\examples{
ALake <- ArcticLake
ALake$AL <- DR_data(ArcticLake[,1:3])

# compute three models
mod0 <- DirichReg(AL ~ 1, ALake)
mod1 <- DirichReg(AL ~ depth, ALake)
mod2 <- DirichReg(AL ~ depth + I(depth^2), ALake)

# sorted = TRUE orders the models by number of parameters
anova(mod1, mod0, mod2, sorted = TRUE)
}
