\encoding{UTF-8}
\name{anova.DirichletRegModel}
\alias{anova.DirichletRegModel}
\title{Compare Dirichlet Regression Models using an LRT}
\description{
This function allows for pairwise tests of Dirichlet regression models using a likelihood ratio test (LRT).
}
\usage{
\method{anova}{DirichletRegModel}(object, \ldots, sorted = FALSE)
}
\arguments{
\item{object}{the model to be compared against those listed in \ldots}
\item{\ldots}{models to be tested against the one specified as \code{object}}
\item{sorted}{should the models be sorted according to their numbers or parameters?}
}
\details{
The test statistic is computed \eqn{LR=-2\left[\log(L_a)-\log(L_b)\right]}{LR=-2[log(La)-log(Lb)]} where \eqn{L_i}{Li} is the likelihood of model \eqn{i}{i} with \eqn{df}{df} equal to the difference of the number of parameters in the models.
}
\author{
Marco Maier
}
\examples{
ALake <- ArcticLake
ALake$AL <- DR_data(ArcticLake[,1:3])
mod0 <- DirichReg(AL ~ 1, ALake)
mod1 <- DirichReg(AL ~ depth, ALake)
mod2 <- DirichReg(AL ~ depth + I(depth^2), ALake)
anova(mod1, mod0, mod2, sorted = TRUE)
}