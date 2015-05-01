\name{multibasta}
\alias{multibasta}
\alias{summary.multibasta}
\alias{print.multibasta}
\title{
Function to run multiple BaSTA models on the same dataset.
}

\description{
These functions allow users to run models with different functional forms on the same dataset, to perform model comparison and to visualise the results of the multiple runs.}

\usage{
multibasta(object, studyStart, studyEnd, models, shapes, \dots)
\method{summary}{multibasta}(object, \dots)
\method{print}{multibasta}(x, \dots)
}

\arguments{
  
\item{object }{For function \code{multibasta}, a \code{data.frame} to be used as an input data file for BaSTA. The first column is a vector of individual unique IDs, the second and third columns are birth and death years respectively. Columns 4-(nt-1) represent the observation window of nt years.  This is followed (optionally) by columns for categorical and continuous covariates (see \code{details}). For function \code{summary}, a \code{multibasta} output object of class ``multibasta.''}
 
\item{x }{a \code{multibasta} output object of class ``multibasta.''}

\item{studyStart }{The first year of the study.}
  
\item{studyEnd }{The last year of the study.}
  
\item{models }{A character vector specifying the models to be tested (e.g. ``GO'', ``WE'', ``LO'')}
  
\item{shapes }{A character vector specifying the shapes to be tested (e.g. ``simple'', ``Makeham'', ``bathtub'').}
  
\item{\dots }{Additional arguments to be passed to function \code{basta} (see \code{details})}

}

\value{  
\item{runs }{A list with the \code{basta} outputs for the models and shapes tested.}
	
\item{DICs }{A matrix with the DIC values for each model, sorted from the model with lowest DIC to the model with highest DIC.}
	
\item{models }{A summary table showing the models and shaped tested.}
  
}

\author{
Fernando Colchero \email{colchero@demogr.mpg.de}, Owen R. Jones \email{jones@biology.sdu.dk}, and Maren Rebke \email{maren.rebke@avitec-research.de}
}

\seealso{
\code{\link{basta}}.
\code{\link{CensusToCaptHist}} and \code{\link{MakeCovMat}} for raw data formatting.
}

\examples{
## Load data:
data("sim1", package = "BaSTA")

## Run short version of BaSTA on the data:
multiout <- multibasta(sim1, studyStart = 51, studyEnd = 70, 
            models = c("GO", "WE"), shapes = "simple",
             niter = 200, burnin = 11, thinning = 10, nsim = 1,
             updateJumps = FALSE, parallel = FALSE)

## Print results:
summary(multiout, digits = 3)
}

\keyword{ methods }
