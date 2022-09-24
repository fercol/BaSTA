\name{DataCheck}
\alias{DataCheck}
\title{A function to check the input data file for a Bayesian Survival Trajectory Analysis (BaSTA) analysis.}

\description{This function performs some basic error checking on the input files for a BaSTA analysis. A number of checks are performed including; (1) individuals dying before the observation window starts; (2) individuals with no observations of any kind (i.e. phantom individuals); (3) individuals with birth date recorded as being after death date; (4) individuals with observations recorded after death; (5) individuals with observations before birth; (6) years of birth must appear as 0 in the observation matrix; (7) years of death must appear as 0 in the observation matrix.}

\usage{DataCheck(object, studyStart, studyEnd, autofix = rep(0, 7), 
       silent = TRUE)}

\arguments{

  \item{object }{A data.frame to be used as an input data file for BaSTA. The first column is the individual's ID, the second and third columns are birth and death years respectively. Columns 4 to \code{nt+3} represent the observation window of nt years. This is followed (optionally) by columns for covariate.}

  \item{studyStart }{The start year of the observation window.}

  \item{studyEnd }{The end year of the observation window.}

  \item{autofix }{A vector argument with a length of 7 indicating whether to automatically fix any errors (see \code{details}). This should be used with extreme caution. We recommend going back to the individual-based data and fixing each error "by hand".}

  \item{silent }{A logical argument indicating whether to print a detailed report to the screen or not.}

}

\details{
Argument \code{autofix} allows the user to fix the potential errors by specifying a code for each fix. Below are the descriptions of the actions that are taken depending on the error type and the fix code:

Type 1: 0 = do nothing; 1 = remove from dataframe.

Type 2: 0 = do nothing; 1 = remove from dataframe.

Type 3: 0 = do nothing; 1 = replace death records with 0; 2 = replace birth records with 0; 3 = replace both birth and death records with 0.

Type 4: 0 = do nothing; 1 = remove spurious post-death observations.

Type 5: 0 = do nothing; 1 = remove observations that pre-date year of birth.

Type 6: 0 = do nothing; 1 = replace birth year element of observation matrix with 0.

Type 7: 0 = do nothing; 1 = replace death year element of observation matrix with 0.
}

\value{
  \item{ok }{A logical indicator that indicates if the data are free of errors or not. i.e. \code{TRUE} = the data have no apparent errors, and \code{FALSE} = there is at leat one error.}
  \item{newData }{A corrected data frame.}
  \item{type1 }{A vector of row numbers in the original data frame where there are deaths occurring before the study starts.}
  \item{type2 }{A vector of row numbers in the original data frame where there are no birth/death AND no obervations.}
  \item{type3 }{A vector of row numbers in the original data frame where there are births recorded after death.}
  \item{type4 }{A vector of row numbers in the original data frame where there are observations (i.e. recaptures) after death.}
  \item{type5 }{A vector of row numbers in the original data frame where there are observations (i.e. recaptures) before birth.}
  \item{type6 }{A vector of row numbers in the original data frame where the year of birth is not a zero in the recapture matrix.}
  \item{type7 }{A vector of row numbers in the original data frame where the year of death is not a zero in the recapture matrix.}
}

\author{Owen R. Jones \email{jones@demogr.mpg.de},  Maren Rebke \email{rebke@demogr.mpg.de}, and Fernando Colchero \email{colchero@demogr.mpg.de}}

\seealso{\code{\link{basta}}}


