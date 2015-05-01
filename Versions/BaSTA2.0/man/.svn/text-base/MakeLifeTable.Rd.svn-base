\name{MakeLifeTable}
\alias{MakeLifeTable}

\title{A function for calculating standard cohort-based life tables from age at death data.
}

\description{
This function calculates a standard cohort-based life table from death ages.
}

\usage{MakeLifeTable(DeathAges, ax = 0.5, n = 1)}

\arguments{
  \item{DeathAges}{A vector of ages at death. These can be precise (e.g. 5.4 yrs) or rounded to the integer year (e.g. 5 yrs).}

  \item{ax}{A vector describing the averaged individual-years lived in the interval by those dying in the interval.   It is only necessary to put the first few values - further values are extrapolated to the correct dimensions. e.g. if \code{ax} is set as \code{c(0.3, 0.5)} and there are 5 rows in the life tables, \code{ax} will become \code{c(0.3, 0.5, 0.5, 0.5, 0.5)}.}

  \item{n}{Defines the interval width. The default is 1 (yr), but it can be set to be any numerical value(s). For example, if you want small intervals for early ages, and larger intervals for older ages, you can define n as a sequence such as \code{c(1,1,5)} to indicate that the first 2 intervals are 1 year long, then the rest of the intervals are 5 years long.}
}

\value{
The function returns a \code{data.frame} in the form of a life table with the following column headings.
Start and end of the interval, where start of the interval is defined as \emph{x}, \emph{lx} (number alive at age \emph{x} (i.e. at the start of the interval)), \emph{dx} (number dying between ages \emph{x} and \emph{x+n}), \emph{mx} (death rate in the cohort between ages \emph{x} and \emph{x+n}), \emph{ax} (average person-years lived in the interval by those dying in the interval), \emph{qx} (probability of dying between ages \emph{x} and \emph{x+n}), \emph{px} (probability of surviving between ages \emph{x} and \emph{x+n}), \emph{Lx} (`person' years lived between ages \emph{x} and \emph{x+n}), \emph{Tx} (`person' years lived above age \emph{x}), \emph{ex} (life expectancy at age \emph{x}).
}

\references{
Preston, S.H., Heuveline, P. and Guillot, M. (2001) \emph{Demography: Measuring and modeling population processes.} Blackwell Publishers, Oxford, UK.
}

\author{Owen Jones \email{jones@mpg.demogr.de}}

\seealso{
\code{\link{basta}}
}

\examples{
MakeLifeTable(rpois(100,3))
}

\keyword{ misc }
\keyword{ datagen }