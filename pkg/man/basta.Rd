\name{basta}
\alias{basta}
\alias{basta.default}
\title{
Parametric Bayesian estimation of age-specific survival for left-truncated and right-censored capture-recapture/recovery data.
}

\description{
This function performs multiple Markov Chain Monte Carlo (MCMC) simulations for the Bayesian estimation of age-specific mortality and survival trends when a large proportion of (or all) records have unknown times of birth and/or death. Survival parameters and unknown (i.e. latent) birth and death times are estimated, allowing the user to test a range of mortality patterns, and to test the effect of continuous and/or discrete covariates following Colchero and Clark's (2012) approach.}

\usage{
basta(object, \dots)

\method{basta}{default}(object, studyStart, studyEnd, 
    minAge = 0, model = "GO", shape = "simple", 
    covarsStruct = "fused", niter = 11000, burnin = 1001, 
    thinning = 20, recaptTrans = studyStart, nsim = 1, 
    parallel = FALSE, ncpus = 2, lifeTable = TRUE, 
    updateJumps = TRUE, \dots)
}

\arguments{
  
\item{object }{A \code{data.frame} to be used as an input data file for BaSTA. The first column is a vector of individual unique IDs, the second and third columns are birth and death years respectively. Columns 4-(nt-1) represent the observation window of nt years.  This is followed (optionally) by columns for categorical and continuous covariates (see \code{details}). }
 
\item{studyStart }{The first year of the study.}
  
\item{studyEnd }{The last year of the study.}
  
\item{minAge }{Age at which the analysis should start (see \code{details})}
  
\item{model }{The underlying mortality model to be used. "EX" = exponential, "GO"= Gompertz, "WE" = Weibull and "LO" = logistic (see \code{details}).}
  
\item{shape }{The overall shape of the model. Values are: \code{simple} = no extra parameters added; \code{Makeham} = a constant parameter is added to the mortality rate; and \code{bathtub} = a Gompertz declining mortality for early ages and a constant parameter are added to the mortality model (see \code{details}).}
  
\item{covarsStruct }{Character string that indicates how covariates should be evaluated. The options are: \dQuote{\code{fused}}, which defines all categorical variables as covariates for each mortality parameter and all continuous covariates under a proportional hazards structure; \dQuote{\code{prop.haz}}, which puts all covariates under a proportional hazards structure; and \dQuote{\code{all.in.mort}} puts all covariates as a multilevel function of the mortality parameters (see \code{details}).}
  
\item{niter }{The total number of MCMC steps.}
  
\item{burnin }{The number of iterations for the burn in (see \code{details}).}
  
\item{thinning }{The number of skipped MCMC steps to minimize serial autocorrelation (see \code{details}).}
  
\item{recaptTrans }{A vector (of length npi) defining the recapture probability transition times (RPTP). These are points (years) where the recapture probability is thought to change.  The default setting is for the recapture probability to be constant throughout the study, so the \code{recaptTrans} is simply defined as a single element vector of the first year of the observation period (e.g. c(1985)). If recapture probabilities were known to change at year say, 1990, the RPTP should be defined as c(1985, 1990).}
    
\item{nsim }{A numerical value for the number of simulations to be run.}
  
\item{parallel }{A logical argument indicating whether the multiple simulations should be run in parallel or not. If \code{TRUE}, package \pkg{snowfall} is called and multiple simulations are run in parallel. If \pkg{snowfall} is not installed, the model is ran in series.}
  
\item{ncpus }{a numerical value that indicates the number of cpus to be used if \code{parallel} is \code{TRUE} and package \pkg{snowfall} is installed. The default is 2 cpus. If package pkg{snowfall} is not installed, the simulations are run in series.}
  
\item{lifeTable }{A logical argument indicating whether or not to produce life tables. If \code{TRUE}, a cohort life table is calculated using function \code{\link{MakeLifeTable}}.}
  
\item{updateJumps }{A logical argument indicating wheter to update jump standard deviations (adaptive independent Metropolis) until an update rate of 0.25 is achieved (see \code{details}). }

\item{\dots }{Additional arguments to be passed to function \code{basta} (see \code{details})}

}

\details{

To construct the input data \code{object} the function \code{\link{CensusToCaptHist}} can be used to build the capture-recapture matrix, while the covariate (design) matrix can be constructed with the  \code{\link{MakeCovMat}} function.

\code{basta} uses parametric mortality functions to estimate age-specific mortality (survival) from capture-recapture/recovery data. The mortality rate function describes how the risk of mortality changes with age, and is defined as \eqn{\mu(x | \theta)}, where \eqn{x} corresponds to age and \eqn{\theta} is a vector of parameters to be estimated. 

The \code{model} argument allows the user to choose between four basic mortality rate functions: 

(a) Exponential (\dQuote{\code{EX}}; Cox and Oakes 1974), with constant mortality rate:

\deqn{\mu_b(x | b) = b}
with \eqn{b > 0}.

(b) Gompertz (\dQuote{\code{GO}}; Gompertz 1925, Pletcher 1999): 

\deqn{\mu_b(x | b) = exp(b_0 + b_1 x)}
with \eqn{-\infty < b_0, b_1 < \infty}.


(c) Weibull (\dQuote{\code{WE}}; Pinder III \emph{et al.} 1978):

\deqn{\mu_b(x | b) = b_0 b_1^(b_0) x^(b_0 -1)}
with \eqn{b_0, b_1 > 0}.


(d) logistic (\dQuote{\code{LO}}; Pletcher 1999):

\deqn{\mu_b(x | b) = exp(b_0 + b_1 x) / (1 + b_2 exp(b_0)/b_1 (exp(b_1 x)-1))}
with \eqn{b_0, b_1, b_2 > 0}.

The \code{shape} argument allows the user to extend these models in order to explore more complex mortality shapes. The default value is \dQuote{\code{simple}} which leaves the model as defined above. With value \dQuote{\code{Makeham}}, a constant is added to the mortality rate, making the model equal to \eqn{\mu_0(x | \theta)= \mu_b(x | b) + c}, where \eqn{\theta = [c, b]}. With value \dQuote{\code{bathtub}}, concave shapes in mortality can be explored. This is achieved by adding a declining Gompertz term and a constant parameter to the basic mortality model, namely:

\deqn{\mu_0(x | \theta) = exp(a_0 - a_1 x) + c + \mu_b(x | b)}.
with \eqn{-\infty < a_0 < \infty}, \eqn{a_1 \leq 0} and \eqn{c \leq 0}.

To incorporate covariates into the inference process, the mortality model is further extended by including a proportional hazards structure, of the form: 

\deqn{\mu(x | \theta, \Gamma, Z_a, Z_c) = \mu_0(x | \theta, Z_a) exp(\Gamma Z_c)}

where \eqn{\mu_0(x | \theta, Z_a)} represents the mortality section as defined above, while the second term \eqn{exp(\Gamma Z_c)} corresponds to the proportional hazards function.  \eqn{Z_a} and \eqn{Z_c} are covariate (design) matrices for categorical and continuous covariates, respectively.

When covariates are included in the dataset, the \code{basta} function provides three different ways in which these can be evaluated by using argument \code{covarsStruct}: 

1. \dQuote{\code{fused}} will make the mortality parameters linear functions of all categorical covariates (analogous to a generalised linear model (GLM) structure) and will put all continuous covariates under a proportional hazards structure. Thus, for a simple exponential mortality rate of the form \eqn{\mu_0(x | b) = b}, the parameter is equal to \eqn{b = b_0 + b_1 z_1 + \dots}, where \eqn{[b_0, \dots, b_k]} are paramters that link the mortality parameter \eqn{b} with the categorical covariates \eqn{[z_1,\dots,z_k]}.

2. \dQuote{\code{prop.haz}} will put all covariates under a proportional hazards structure irrespective of the type of variable.

3. \dQuote{\code{all.in.mort}} will put all covariates as linear functions of the survival parameters as explained above. Since most models require the lower bounds for the mortality parameters to be equal to 0, the only model that can be used for this test is Gompertz with \code{shape} set to \dQuote{\code{simple}}. In case these arguments are specified deferently, a warning message is printed noting that \code{model} will be forced to be \dQuote{\code{GO}} and \code{shape} will be set to \dQuote{\code{simple}}.

The \code{burnin} argument represents the number of steps at the begining of the MCMC run that is be discarded. This sequence commonly corresponds to the non-converged section of the MCMC sequence. Convergence and model selection measures are calculated from the remaining thinned parameter chains if multiple simulations are run, and all if all of them run to completion.

The \code{thinning} argument specifies the number of steps to be skipped in order to reduce serial autocorrelation. The thinned sequence, which only includes steps after burn in, is then used to calculate convergence statistics and model for selection.

The \code{updateJumps} argument specifies wether to run a simulation to find appropriate jump standard deviations for theta and gamma parameters. If argument \dQuote{\code{nsim}} is set to 1, then the simulation runs with the update jumps routine active. If \dQuote{\code{nsim}} is larger than 1, then an initial simulation is ran to find apropriate jumps before the main analysis is ran. 

Additional arguments for priors, jumps and start values can be passed on the \dots section. For instance, argument \code{thetaStart} can be specified as a vector defining the initial values for each parameter in the survival model. If this argument is not specified, a set of random parameters is generated for each simulation. Similarly, argument \code{gammaStart} can be specified for all parameters in the proportional hazards section of the model. Jump standard deviations (i.e. the standard error in the Metropolis step) can be specified with arguments \code{thetaJumps} and \code{gammaJumps}. As with \code{thetaStart}, default values are assigned if these arguments are not specified. Arguments \code{thetaPriorMean}, \code{thetaPriorSd}, \code{gammaPriorMean} and \code{gammaPriorSd} can be used to specify prior means and standard errors for each survival and proportional hazards parameters. Again, if not specified, default values are assigned.

The number of parameters in \code{thetaStart}, \code{thetaJumps}, \code{thetaPriorMean} and \code{thetaPriorSd} should be a vector or matrix for the parameters in the mortality function. The number of parameters will depend on the model chosen with \code{model} (see above). If the number of parameters specified does not match the number of parameters inherent to the model and shape selected, the function returns an error. 

As described above, the number of parameters for \code{gammaStart}, \code{gammaJumps}, \code{gammaPriorMean} and \code{gammaPriorSd} arguments (i.e. section b), namely the proportional hazards section, will be a function of the number of continuous covariates if argument \code{covarsStruct} is \dQuote{\code{fused}}, or to the total number of covariates when \code{covarsStruct} is \dQuote{\code{prop.haz}}.

}

\value{  
\item{coefficients }{A matrix with estimated coefficients (i.e. mean values per parameter on the thinned sequences after burnin), which includes standard errors, upper and lower 95\% credible intervals, update rates per parameter, serial autocorrelation on the thinned sequences and the potential scale reduction factor for convergence (see \code{Convergence} value below). }
	
\item{DIC }{Basic deviance information criterion (DIC) calculations to be used for model selection (Spiegelhalter \emph{et al.} 2002). Small differences between values should only be used a reference (see comments in Spiegelhalter \emph{et al.} 2002). If all or some of the simulations failed, then the returned value is \dQuote{\code{Not calculated}}.}
	
\item{convergence }{If requested, a matrix with convergence coefficients based on potential scale reduction as described by Gelman \emph{et al.} (2004). If only one simulation was ran, then the returned value is \dQuote{\code{Not calculated}}.}
  
\item{KullbackLeibler }{If called by \code{\link{summary}}, list with Kullback-Leibler discrepancy matrices between pair of parameters for categorical covariates (McCulloch 1989, Burnham and Anderson 2001) and McCulloch's (1989) calibration measure. If only one simulation was ran or if no convergence is reached, then the returned value is \dQuote{\code{Not calculated}}.}
  
\item{params }{If requested, a matrix with the thinned, converged parameter traces of all runs. This matrix is used to calculate quantiles for parameters, survival probability and mortality (see below).}

\item{settings }{If called by \code{\link{summary}}, this is a vector indicating the number of iterations for each MCMC, the burn in sequence, the thinning interval, and the number of simulations that were run.}
	
\item{modelSpecs }{Model specifications inidicating the \code{model}, the \code{shape}, the covariate structure and the minimum age that were specified by the user.}
	
\item{jumpPriors }{If requested or called by functions \code{summary} or \code{\link{summary.basta}}, a matrix with the jump standard deviations and prior means and standard deviation for the parameters.}
	
\item{birthQuant }{If requested, summary matrix of estimated times of birth including median and upper and lower 95\% predictive intervals.}
	
\item{deathQuant }{If requested, summary matrix of estimated times of birth including median and upper and lower 95\% predictive intervals.}
	
\item{agesQuant }{If requested, summary matrix of estimated ages at death including median and upper and lower 95\% predictive intervals.}
  	
\item{mortQuant }{If requested or called by functions \code{plot} or \code{\link{plot.basta}} median and 95\% predictive intervals for the estimated mortality rates.}

\item{survQuant }{If requested or called by functions \code{plot} or \code{\link{plot.basta}} median and 95\% predictive intervals for the estimated survival probability.}
	
\item{parsForPlot }{If requested or called by functions \code{plot} or \code{\link{plot.basta}} thinned sequences, including burn-in, of parameter traces used to plot outputs.}

\item{lifeTable }{If requested and specified in the argument \code{lifeTable}, a cohort life table calculated from the estimated ages at death.}
}

\references{
Burnham, K.P. and Anderson, D.R. (2001) Kullback-Leibler information as a basis for strong inference in ecological studies. \emph{Widlife Research}, 28, 111-119.

Colchero, F. and J.S. Clark (2012) Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. \emph{Journal of Animal Ecology}. 81, 139-149.

Colchero, F., O.R. Jones and M. Rebke. (2012) BaSTA: an R package for Bayesian estimation of age-specific survival from incomplete mark-recapture/recovery data with covariates. \emph{Method in Ecology and Evolution}. 3, 466-470. 

Cox, D. R., and Oakes D. (1984) \emph{Analysis of Survival Data}. Chapman and Hall, London.

Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B. (2004) \emph{Bayesian data analysis}. 2nd edn. Chapman & Hall/CRC, Boca Raton, Florida, USA.

Gompertz, B. (1825) On the nature of the function expressive of the law of human mortality, and on a new mode of determining the value of life contingencies. \emph{Philosophical Transactions of the Royal Society of London}, 115, 513-583.

King, R. and Brooks, S.P. (2002) Bayesian model discrimination for multiple strata capture-recapture data. \emph{Biometrika}, 89, 785-806. 

McCulloch, R.E. (1989) Local model influence. \emph{Journal of the American Statistical Association}, 84, 473-478.

Pinder III, J.E., Wiener, J.G. and Smith, M.H. (1978) The Weibull distribution: a new method of summarizing survivorship data. \emph{Ecology}, 59, 175-179.

Spiegelhalter, D.J., Best, N.G., Carlin, B.P. and van der Linde, A. (2002) Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B}, 64, 583-639.
}

\author{
Fernando Colchero \email{colchero@imada.sdu.dk}, Owen R. Jones \email{jones@biology.sdu.dk}, and Maren Rebke \email{maren.rebke@avitec-research.de}
}

\seealso{
\code{\link{summary.basta}}, \code{\link{print.basta}}, \code{\link{plot.basta}} to visualise summary outputs for objects of class \dQuote{\code{basta}}.
\code{\link{CensusToCaptHist}} and \code{\link{MakeCovMat}} for raw data formatting.
}
\examples{
## Load data:
data("sim1", package = "BaSTA")

## Check data consistency:
new.dat  <- DataCheck(sim1, studyStart = 51, 
            studyEnd = 70, autofix = rep(1,7))

## Run short version of BaSTA on the data:
out <- basta(sim1, studyStart = 51, studyEnd = 70, 
             niter = 200, burnin = 11, thinning = 10,
             updateJumps = FALSE)

## Print results:
summary(out, digits = 3)

## Plot traces for survival parameters:
plot(out)

## Plot traces for proportional hazards parameter:
plot(out, trace.name = "gamma")

## Plot survival and mortality curves:
plot(out, plot.trace = FALSE)
}

\keyword{ methods }
