\name{sim1}
\alias{sim1}
\docType{data}
\title{Simulated dataset to demonstrate Bayesian Survival Trajectory Analysis (BaSTA)}
\description{
This dataset was created by stochastically simulating a hypothetical population with different mortality patterns between males and females and with proportional decreases in mortality as a function of a hypothetical continuous covariate (e.g. birth weight, average adult weight, etc.). The population was simulated for 100 years, at each one of which 100 individuals were born. The number of females per generation was randomly drawn from a binomial distribution using function \code{rbinom} with probability of 0.5 (i.e. 1:1 sex ratio). The individual continuous covariate was randomly drawn from a random normal distribution (with function \code{\link{rnorm}}) with mean parameter equal to 0 (e.g. anomaly of weights) and standard deviation equal to 1.  The time of death for each individual was inversed sampled from a Gompertz CDF of ages at death. The Gompertz parameters for females were: \eqn{\alpha = -4} and \eqn{\beta = 0.15}; and for males at \eqn{\alpha = -3} and \eqn{\beta = 0.15}. The gamma parameter for the continuous covariate was \eqn{\gamma = 0.2}.

The study was assumed to start at year 51 and to finish at year 70. Recapture probability was set to 0.6 and thus each observation per individual was randomly drawn from a Bernoulli trial with parameter \eqn{p = 0.6}. Captures at birth and recoveries were randomly drawn from a Bernoulli trial with parameters \eqn{p = 0.5} and \eqn{p = 0.2}, respectively. 

Therefore, the resulting dataset includes individuals where the data are left-truncated and/or right-censored. This is typical of capture mark recovery datasets.

}
\usage{sim1}
\format{RData file, data frame.}
\keyword{datasets}