# TODO: Add comment
# 
# Author: fernando
###############################################################################
.PrepCovDat <- function(object, ...) UseMethod(".PrepCovDat")

.PrepCovDat.noCovs <- function(object, ...) {
  covars <- 0
  class(covars) <- c("noCovs", "noClass")
}

.PrepCovDat.fixedCovs <- function(object, covarsStruct) {
  covars <- .PrepCovsFixed(object, covarsStruct)
  covClass <- class(covars)
  class(covars) <- c("fixedCovs", covClass)
  return(covars)
}

.PrepCovDat.bothCovs <- function(object, covarsStruct) {
  covars <- .PrepCovsFixed(object, covarsStruct)
  covClass <- class(covars)
  covars$prHazTime$mat <- object$timeCovs
  if (is.array(object$timeCovs)) {
    covars$prHazTime$length <- dim(object$timeCovs)[3]
    covars$prHazTime$names <- unlist(dimnames(object$timeCovs)[3])
  } else {
    covars$prHazTime$length <- 1
    covars$prHazTime$names <- "TimeCov"
  }
  covars$timeMat <- matrix(1:ncol(object$timeCovs) + (timeCovStart - 1),
      nrow = nrow(object$timeCovs), ncol = ncol(object$timeCovs), 
      byrow = TRUE)
  covars$cumMat <- matrix(t(apply(diag(1, ncol(object$timeCovs), 
                  ncol(object$timeCovs)), 1, cumsum)), ncol(object$timeCovs), 
      ncol(object$timeCovs), 
      dimnames = list(colnames(object$timeCovs), colnames(object$timeCovs)))
  class(covars) <- c("bothCovs", covClass)
  return(covars)
}

.PrepCovs.datTimeCovs <- function(object, ...) {
  covars <- list(inMort = list(mat = matrix(0, 1, 1, 
              dimnames = list(NULL, "NoCov")),
          length = 1,
          names = "NoCov"),
      prHazFix = list(mat = matrix(0, 1, 1, 
              dimnames = list(NULL, "NoCov")),
          length = 1,
          names = "NoCov"))
  covars$prHazTime$mat <- object$timeCovs
  if (is.array(object$timeCovs)) {
    covars$prHazTime$length <- dim(object$timeCovs)[3]
    covars$prHazTime$names <- unlist(dimnames(object$timeCovs)[3])
  } else {
    covars$prHazTime$length <- 1
    covars$prHazTime$names <- "TimeCov"
  }
  covars$timeMat <- matrix(1:ncol(object$timeCovs) + (timeCovStart - 1),
      nrow = nrow(object$timeCovs), ncol = ncol(object$timeCovs), 
      byrow = TRUE)
  covars$cumMat <- matrix(t(apply(diag(1, ncol(object$timeCovs), 
                  ncol(object$timeCovs)), 1, cumsum)), ncol(object$timeCovs), 
      ncol(object$timeCovs), 
      dimnames = list(colnames(object$timeCovs), colnames(object$timeCovs)))
  class(covars) <- c("timeCovs", "noClass")
  return(covars)
}

.PrepCovsFixed <- function(object, covarsStruct) {
  covMat <- object$fixedCovs
  covarsType <- .FindCovariateType(covMat)
  if (!is.null(covarsType$cont)) {
    rangeCont <- apply(matrix(covMat[, covarsType$cont],
            ncol = length(covarsType$cont)), 2, range)
    idChange <- which(rangeCont[1, ] * rangeCont[2, ] > 0)
    if (length(idChange) > 0) {
      for (i in covarsType$cont[idChange]) {
        covMat[, i] <- covMat[, i] - mean(covMat[, i])
      }
    }
  }
  if (covarsStruct == "fused") {
    if (is.null(covarsType$cat)) {
      covMatInMort <- matrix(0, 1, 1, dimnames = list(NULL, "NoCov"))
    } else {
      covMatInMort <- covMat[, covarsType$cat]
    }
    if (is.null(covarsType$cont)) {
      covMatPrHaz <- matrix(0, 1, 1, dimnames = list(NULL, "NoCov"))
    } else {
      covMatPrHaz <- covMat[, covarsType$cont]
    }
    covars <- list(inMort = list(mat = covMatInMort,
            length = ncol(covMatInMort),
            names = colnames(covMatInMort)),
        prHazFix = list(mat = covMatPrHaz,
            length = ncol(covMatPrHaz),
            names = colnames(covMatPrHaz)))
    if (is.null(covarsType$cat)) {
      class(covars) <- "propHaz"
    } else if (is.null(covarsType$cont)) {
      class(covars) <- "inMort"
    } else {
      class(covars) <- "fused"
    }
  } else if (covarsStruct == "prop.haz") {
    covars <- list(inMort = list(mat = matrix(1, 1, 1, 
                dimnames = list(NULL, "NoCov")),
            length = 1,
            names = "NoCov"), 
        prHazFix = list(mat = covMat,
            length = ncol(covMat),
            names = colnames(covMat)))
    class(covars) = "propHaz"
  } else {
    covars <- list(inMort = list(mat = covMat,
            length = ncol(covMat),
            names = colnames(covMat)),
        prHazFix = list(mat = matrix(0, 1, 1, 
                dimnames = list(NULL, "NoCov")),
            length = 1,
            names = "NoCov"))
    class(covars) = "inMort"
  }
  return(covars)
}

# 3.1. Covariate type (i.e. categorical and continuous):
.FindCovariateType   <- function(covMat) {
  lu  <- apply(covMat, 2, function(x) length(unique(x)))
  ru  <- apply(covMat, 2, range)
  idcat <- which(lu == 2 & apply(ru, 2, sum) == 1)
  if (length(idcat) == 0) {
    idcat <- NULL
  }
  idint <- which(lu == 1)
  if (length(idint) == 0) {
    idint <- NULL
  }
  idcon <- which(lu > 2)
  if (length(idcon) == 0) {
    idcon <- NULL
  }
  return(list(int  = idint, 
          cat  = idcat, 
          cont = idcon))
}

