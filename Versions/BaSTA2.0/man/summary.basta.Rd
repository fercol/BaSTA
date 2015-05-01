\name{summary.basta}
\alias{summary.basta}
\alias{print.basta}
\alias{plot.basta}

\title{
Summarizing and plotting Bayesian Survival Trajectory Analysis (BaSTA) model outputs.
}
\description{
These functions are all generic methods for class \code{basta}.
}
\usage{
\method{summary}{basta}(object, \dots)
\method{print}{basta}(x, \dots)
\method{plot}{basta}(x, plot.trace = TRUE, trace.name = "theta", \dots)
}
\arguments{
  \item{object }{An object of class \code{basta}.}
  
  \item{x }{An object of class \code{basta}.}

  \item{plot.trace}{A logical argument. If \code{TRUE} the raw parameter traces are plotted, else, the predictive intervals for the resulting survival probability and mortality rates are plotted.}
  
  \item{trace.name }{Character string indicating the set of parameters or posteriors to be plotted. The options are: "\code{theta}" to plot the survival model parameters; "\code{gamma}" to plot the proportional hazards parameters (if it applies, else \code{plot.basta} returns an error); "\code{pi}" to plot the recapture probabilities; and "post" to plot the conditional posteriors for the parameters and for the latent ages at death, and the full posterior for the model.}
  
  \item{\dots }{Additional arguments passed to functions \code{print}, \code{summary} and \code{plot}.}
}
\details{
For objects of class \code{basta}, the \code{print} function returns three summary elements describing the model and its results, namely: \code{call}, \code{run} and \code{coefficients}. \code{call} describes the basic model used (i.e. exponential, Gompertz, Weibull or logistic), the shape chosen, \dQuote{\code{simple}}, \dQuote{\code{Makeham}} or \dQuote{\code{bathtub}}, the covariate structure chosen, \dQuote{\code{fused}}, \dQuote{\code{prop.haz}} or \dQuote{\code{all.in.mort}} and which covariates (if any) were categorical and which continuous. 

The summary element \code{run} describes whether all of the simulations ran for  all of the iterations (specified by \code{niter}) specified by the user in the \code{\link{basta}} function. If not all of the runs were completed, it outlines which of them have failed.

The summary element \code{coefficients} prints out the estimated coefficients for all parameters in the model, as well as their standard errors and 95\% upper and lower credible intervals. It also includes a measure of serial autocorrelation for each parameter calculated from the thinned parameter chains, an update rate per parameter (commonly the same for all survival and proportional hazards parameters), and the potential scale reduction factor for each parameter as a measure of convergence (Gelman \emph{et al.} 2004). 

Function \code{summary} includes all the previous elements, as well as a summary description of the priors and jump standard deviations for all survival parameters, a calibration of the Kullback-Liebler discrepancy as a measure of parameter similarities for those parameters associated to categorical covariates (McCulloch 1989), and a measure of model fit based on the deviance information criterion (DIC) (Spiegelhalter \emph{et al.} 2002). 

Function \code{plot} takes objects of class \code{basta} to create trace plots or, if the argument for \code{plot.trace} is set to \code{FALSE}, it plots estimated survival probabilities and mortality rates with their 95\% predictive intervals. 
}
\references{
Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B. (2004) \emph{Bayesian data analysis}. 2nd edn. Chapman & Hall/CRC, Boca Raton, Florida, USA.

McCulloch, R.E. (1989) Local model influence. \emph{Journal of the American Statistical Association}, 84, 473-478.

Spiegelhalter, D.J., Best, N.G., Carlin, B.P. and Van Der Linde, A. (2002) Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B} 64, 583-639.
}
\author{
Fernando Colchero \email{colchero@demogr.mpg.de}, Owen R. Jones \email{jones@demogr.mpg.de} and Maren Rebke \email{rebke@demogr.mpg.de}
}

\seealso{
\code{\link{basta}}
}

\examples{
## Load BaSTA output:
data("sim1Out", package = "BaSTA")

## Print summary output:
summary(sim1Out)

## Plot traces for mortality parameters (theta):
plot(sim1Out)

## Plot traces for proportional hazards parameters (gamma):
plot(sim1Out, trace.name = "gamma")

## Plot predicted mortality and survival:
plot(sim1Out, plot.trace = FALSE)
}

