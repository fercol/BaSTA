\name{CensusToCaptHist}
\alias{CensusToCaptHist}

\title{
Constructs capture-history matrix from census matrix to be used in Bayesian Survival Trajectory Analysis (BaSTA).
}
\description{
This function takes a common census table, consisting of a pair of vectors with ID and observation date, and converts it into a capture-history matrix to be used as part of a BaSTA data input.
}
\usage{
CensusToCaptHist(ID, d, dformat = "\%Y", timeInt = "Y")
}

\arguments{
  \item{ID}{
A vector with individual IDs.
}
  \item{d}{
A vector of dates when each individual was observed. 
}
  \item{dformat}{
Defines the date format for \code{d} when \code{d} is of class \code{character}.
}
  \item{timeInt }{
A one character string specifying which time interval should be used between capture occasions. Arguments are \dQuote{\code{Y}} for years, \dQuote{\code{M}} for months, \dQuote{\code{W}} for weeks and \dQuote{\code{D}} for days. 
}
}
\details{
The \code{d} argument can be specified as an object of class \code{POSIXct} or \code{POSIXlt}, as a vector of integer time intervals or as a character string indicating the day, month and year. (e.g. \code{dd/mm/yyyy}, \code{mmddyyyy}, \code{mm-dd-yyyy} etc.). When \code{d} is of class \code{character} then argument \code{dformat} needs to be specified using the same conventions as in function \code{\link{format.POSIXct}} for objects of class \code{POSIXct} or \code{POSIXlt}.
}

\author{
Owen R. Jones \email{jones@demogr.mpg.de} and Maren Rebke \email{rebke@demogr.mpg.de}
}

\seealso{
\code{\link{MakeCovMat}}, which formats a covariate matrix compatible with this output. 
}
\examples{
id.vec <- sort(sample(1:5, size = 15, replace = TRUE))
d.vec <- rep(0, length(id.vec))
for(i in unique(id.vec)){
  svec <- which(id.vec == i)
  d.vec[svec] <- sort(sample(1990:1995, length(svec)))
}
Y <- CensusToCaptHist(ID = id.vec,  d = d.vec)
}
