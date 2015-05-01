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
  
  \item{\dots }{Additional arguments passed to functions \code{print}, \code{summary} and \code{plot} (see details).}
}
\details{
For objects of class \code{basta}, the \code{print} function returns three summary elements describing the model and its results, namely: \code{call}, \code{run}, \code{coefficients} and, if convergence was reached, the \code{DIC} values for model fit. \code{call} describes the basic model used (i.e. exponential, Gompertz, Weibull or logistic), the shape chosen, \dQuote{\code{simple}}, \dQuote{\code{Makeham}} or \dQuote{\code{bathtub}}, the covariate structure chosen, \dQuote{\code{fused}}, \dQuote{\code{prop.haz}} or \dQuote{\code{all.in.mort}} and which covariates (if any) were categorical and which continuous. Argument \code{digits} can be used for number formatting (see \code{\link{summary}}() or \code{\link{signif}}() for details).  

The summary element \code{coefficients} prints out the estimated coefficients for all parameters in the model, as well as their standard errors and 95\% upper and lower credible intervals. It also includes a measure of serial autocorrelation for each parameter calculated from the thinned parameter chains, an update rate per parameter, and the potential scale reduction factor for each parameter as a measure of convergence (Gelman \emph{et al.} 2004). 

Function \code{summary} includes all the previous elements, as well as a summary description of the priors and jump standard deviations for all survival parameters, a calibration of the Kullback-Leibler discrepancy as a measure of parameter similarities for those parameters associated to categorical covariates (McCulloch 1989), and a measure of model fit based on the deviance information criterion (DIC) (Spiegelhalter \emph{et al.} 2002). 

Function \code{plot} takes objects of class \code{basta} to create trace plots or, if the argument for \code{plot.trace} is set to \code{FALSE}, it plots estimated survival probabilities and mortality rates with their 95\% predictive intervals. If argument \code{plot.trace} is set to \code{FALSE}, argument \code{xlim} can be used to define a range of ages to visualize survival and mortality trends. Also, if logical argument \code{noCI} is set to \code{TRUE}, credible intervals around survival and mortality are not plotted, leaving only the mean trends. This can be handy when several categorical covariates have been evaluated and the plots get too crowded.  
}

\value{ 
Function \code{summary}() outputs the folowing values:

\item{coefficients }{A matrix with estimated coefficients (i.e. mean values per parameter on the thinned sequences after burnin), which includes standard errors, upper and lower 95\% credible intervals, update rates per parameter (commonly the same for all survival and proportional hazards parameters), serial autocorrelation on the thinned sequences and the potential scale reduction factor for convergence (see \code{Convergence} value below). }
  
\item{DIC }{Basic deviance information criterion (DIC) calculations to be used for model selection (Spiegelhalter \emph{et al.} 2002).}
	
\item{KullbackLeibler }{List with Kullback-Leibler discrepancy matrices between pair of parameters for categorical covariates (McCulloch 1989, Burnham and Anderson 2001) and McCulloch's (1989) calibration measure. If only one simulation was ran or if no convergence was reached, then the returned value is \dQuote{\code{Not calculated}}.}

\item{convergence }{A matrix with convergence coefficients based on potential scale reduction as described by Gelman \emph{et al.} (2004). If only one simulation was ran, then the returned value is \dQuote{\code{Not calculated}}.}
    
\item{modelSpecs }{Model specifications inidicating the \code{model}, the \code{shape} and the covariate structure that were specified by the user.}

\item{settings }{A vector indicating the number of iterations for each MCMC, the burn in sequence, the thinning interval, and the number of simulations that were run.}
	
}

\references{
Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B. (2004) \emph{Bayesian data analysis}. 2nd edn. Chapman & Hall/CRC, Boca Raton, Florida, USA.

McCulloch, R.E. (1989) Local model influence. \emph{Journal of the American Statistical Association}, 84, 473-478.

Spiegelhalter, D.J., Best, N.G., Carlin, B.P. and Van Der Linde, A. (2002) Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B} 64, 583-639.

See also:

Colchero, F. and J.S. Clark (2012) Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. \emph{Journal of Animal Ecology}. 81(1):139-149.

Colchero, F., O.R. Jones and M. Rebke. (2012) BaSTA: an R package for Bayesian estimation of age-specific survival from incomplete mark-recapture/recovery data with covariates. \emph{Method in Ecology and Evolution}. DOI: 10.1111/j.2041-210X.2012.00186.x
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

## Plot traces for recapture probability(ies) (pi):
plot(sim1Out, trace.name = "pi")

## Plot predicted mortality and survival:
plot(sim1Out, plot.trace = FALSE)

## Plot predicted mortality and survival between 2 and 8 years of age:
plot(sim1Out, plot.trace = FALSE, xlim = c(2, 8))
}

