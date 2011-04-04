\name{anova.DirichletRegModel}
\alias{anova.DirichletRegModel}
\title{Compare Dirichlet Regression Models using an LRT}
\description{
This function allows for pairwise tests of Dirichlet regression using a likelihood-ratio-test.
}
\usage{
\method{anova}{DirichletRegModel}(object, \ldots, sorted = FALSE)
}
\arguments{
  \item{object}{a model to be tested}
  \item{\ldots}{more models to be included in the LRT}
  \item{sorted}{should the models be sorted according to their numbers or parameters?}
}
\details{
no details.
}
\author{
  Marco Maier
}
\examples{
AL <- DR.data(ArcticLake[,1:3])
mod1 <- DirichReg(AL ~ depth, ArcticLake)
mod2 <- DirichReg(AL ~ depth + I(depth^2), ArcticLake)
anova(mod1, mod2, sorted = TRUE)
}