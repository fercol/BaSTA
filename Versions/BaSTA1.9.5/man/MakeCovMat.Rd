\name{MakeCovMat}
\alias{MakeCovMat}

\title{
Function to build a matrix of covariates (i.e. design matrix) for a Bayesian Survival Trajectory Analysis (BaSTA) analysis.
}
\description{
This function creates a matrix of covariates or design matrix appropriate for BaSTA from raw individual level covariate data. The function identifies categorical and continuous covariates and organizes them accordingly. 
}
\usage{
MakeCovMat(x, data)
}

\arguments{
 
 \item{x }{A character string vector or a numerical vector indicating the columns to be included or an object of class \code{formula}.}
 
 \item{data }{A data frame of n rows (n = number of individuals in dataset) and nz columns (number of general covariates) including categorical covariates as factors (e.g. sex with individual labels "f","m", location, etc.) and/or continuous individual covariates (e.g. weight at birth, general weather conditions at each location, etc.).}
  }

\details{
The \code{x} argument can be of class \code{character}, \code{numeric} or \code{formula} as long as the elements described correspond to the column names in the \code{data} data frame. 
The data frame specified in argument \code{data} needs to explicitly differentiate between categorical and numerical variables. The elements in the column of a categorical variable must be coerced to be \code{factors}.  
}

\value{
The function returns a new covariate matrix to be collated to a matrix that includes a column for individual ID, a column for time of birth, and a column for times of death, plus the full recapture matrix.
}

\author{
Owen R. Jones \email{jones@demogr.mpg.de} and Fernando Colchero \email{colchero@demogr.mpg.de}
}

\seealso{
\code{\link{basta}}
}
\examples{
## Simulated sex and weight data for 5 individuals:
sex        <- sample(c("f", "m"), 5, replace = TRUE)
weight     <- rnorm(5, mean = 10, sd = 1)
raw.mat    <- data.frame(sex, weight)
new.mat    <- MakeCovMat(~sex + weight, data = raw.mat)
}
