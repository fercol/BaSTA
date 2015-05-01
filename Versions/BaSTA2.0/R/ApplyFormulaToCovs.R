# TODO: Add comment
# 
# Author: fernando
###############################################################################

.CreateFinalCovList <- function(object, form = NULL) {
  if (class(object)[2] == "noCovs") {
    covars <- 0
    class(covars) <- c("covars", "noCovs")
  } else {
    if (is.null(form)) {
      if (class(object)[2] == "fixedCovs") {
        form <- paste("~", paste(colnames(object$fixedCovs), collapse = "+"))
      } else if (class(object)[2] == "timeCovs") {
        form <- paste("~", paste(colnames(object$timeCovs), collapse = "+"))
      } else {
        form <- paste("~", paste(c(colnames(object$fixedCovs), 
                    names(object$timeCovs)), collapse = "+"))
      }
      class(covars) <- c("covars", class(object)[2])
      form <- as.formula(form)
      
    } else {
      .VerifyCovsInForm(object, form)
      covars <- .CreateModelList(object, x = form)
      class(covars) <- c("covars", class(covars))
    }
  }
  return(covars)
} 


.VerifyCovsInForm <- function(object, form = NULL) {
  allCovNames <- c(object$nameFCovs, object$nameTCovs)
  nCovs <- length(allCovNames)
  if (is.null(form)) {
    covNames <- allCovNames
  } else {
    if (is.numeric(form)) {
      if (!all(form <= nCovs)) {
        stop("Some arguments in 'form' do not match the number of covariates",
            "in 'object'.\n", call. = FALSE)
      } else {
        covNames <- allCovNames[form]
      }
    } else if (is.character(form)) {
      if (!all(form %in% allCovNames)) {
        stop("Some arguments in 'form' do not match the covariate names in ",
            "'object'.\n", call. = FALSE)
      }
      covNames <- form
    } else if (class(form) == 'formula') {
      covNames <- labels(terms(form))
      covNames <- unique(unlist(strsplit(covNames, split = ":")))
      if (!all(covNames %in% allCovNames)) {
        stop("Some arguments in 'form' do not match the covariate names in ",
            "'object'.\n", call. = FALSE)
      } else {
        covNames <- allCovNames[which(allCovNames %in% covNames)]
      }
    }
  }
#  return(covNames)
}

.CreateModelList <- function(object, ...) UseMethod(".CreateModelList")

.CreateModelList.fixedCovs <- function(object, x) {
  fixedCovsMat <- model.matrix(object = x, data = object$fixedCovs)
  allCovList <- list(fixedCovs = cbind(fixedCovsMat))
  class(allCovList) <- "fixedCovs"
  return(allCovList)
}

.CreateModelList.timeCovs <- function(object, x) {
  n <- nrow(object$fixedCovs)
  modFrame <- .ExtractModelFrame(object, x)
  modelNames <- colnames(model.matrix(x, 
          data = data.frame(modFrame$basicCovList)))
  modelNameList <- strsplit(modelNames, split = ":")
  
  timeCovList <- list()
  for (i in 1:length(modelNames)) {
    localCovs <- modelNameList[[i]]
    nLocCovs <- length(localCovs)
    if (nLocCovs == 1) {
      if (localCovs == "(Intercept)") {
        timeCovList[["(Intercept)"]] <- object$timeCovs[[1]] * 0 + 1
      } else {
        idCov <- which(modFrame$covNames == localCovs)
        timeCovList[[localCovs]] <- modFrame$covList[[idCov]]
      }
    } else {
      timeCovList[[modelNames[i]]] <- modFrame$covList[[1]] * 
          0 + 1
      for(jj in 1:nLocCovs) {
        timeCovList[[modelNames[i]]] <- timeCovList[[modelNames[i]]] *
            modFrame$covList[[localCovs[jj]]]
      }
    }
  }
  allCovList <- list(timeCovs = timeCovList)
  class(allCovList) <- "timeCovs"
  return(allCovList)
}

.CreateModelList.bothCovs <- function(object, x) {
  modFrame <- .ExtractModelFrame(object, x)
  if (!("time" %in% modFrame$covType)) {
    fixedCovsMat <- model.matrix(object = x, data = object$fixedCovs)
    allCovList <- list(fixedCovs = cbind(fixedCovsMat))
    class(allCovList) <- "fixedCovs"
  } else {
    n <- nrow(object$fixedCovs)
    modelNames <- colnames(model.matrix(object = x, 
            data = data.frame(modFrame$basicCovList)))
    modelNameList <- strsplit(modelNames, split = ":")
    
    fixedCovMat <- matrix(0, n, 0)
    timeCovList <- list()
    for (i in 1:length(modelNames)) {
      localCovs <- modelNameList[[i]]
      nLocCovs <- length(localCovs)
      if (nLocCovs == 1) {
        if (localCovs == "(Intercept)") {
          if (("fixed" %in% modFrame$covType)) {
            fixedCovMat <- matrix(1, n, ncol(fixedCovMat) + 1)
            colnames(fixedCovMat) <- "(Intercept)"
          } else {
            timeCovList[["(Intercept)"]] <- object$timeCovs[[1]] * 0 + 1 
          }
        } else {
          idCov <- which(modFrame$covNames == localCovs)
          if (modFrame$covType[idCov] == "fixed") {
            namesFixed <- colnames(fixedCovMat)
            fixedCovMat <- cbind(fixedCovMat, modFrame$covList[[idCov]])
            colnames(fixedCovMat) <- c(namesFixed, localCovs)
          } else {
            timeCovList[[localCovs]] <- modFrame$covList[[idCov]]
          }
        }
      } else {
        logicTcov <- which(modFrame$covType[localCovs] == "time")
        if(length(logicTcov) > 0) {
          timeCovList[[modelNames[i]]] <- modFrame$covList[[logicTcov[1]]] * 
              0 + 1
          for(jj in 1:nLocCovs) {
            timeCovList[[modelNames[i]]] <- timeCovList[[modelNames[i]]] *
                modFrame$covList[[localCovs[jj]]]
          }
        } else {
          fixedInt <- rep(1, n)
          for(jj in 1:nLocCovs) {
            fixedInt <- fixedInt * modFrame$covList[[localCovs[jj]]]
          }
          fixedCovMat <- cbind(fixedCovMat, fixedInt)
          colnames(fixedCovMat) <- c(colnames(fixedCovMat)[-ncol(fixedCovMat)], 
              modelNames[i])
        }
      }
    }
    if (length(timeCovList) == 0) {
      allCovList <- list(fixedCovs = fixedCovMat)
      class(allCovList) <- "fixedCovs"
    } else {
      allCovList <- list(fixedCovs = fixedCovMat, timeCovs = timeCovList)
      class(allCovList) <- "bothCovs"
    }
  }
  return(allCovList)
}


.ExtractModelFrame <- function(object, x) {
  labelsCov <- labels(terms(x))
  xcov <- strsplit(labels(terms(x)), split = ":")
  uniqueCov <- unique(unlist(strsplit(labels(terms(x)), split = ":")))
  covList <- list()
  basicCovList <- list()
  covType <- c()
  for (i in 1:length(uniqueCov)) {
    idCovType <- which(c(is.element(uniqueCov[i], object$nameFCovs), 
            is.element(uniqueCov[i], object$nameTCovs)))
    if (idCovType == 1) {
      tempCov <- object$fixedCov[, uniqueCov[i]]
      basicCovList[[uniqueCov[i]]] <- tempCov
      if (is.factor(tempCov) | is.character(tempCov)) {
        tempCov <- tapply(rep(1, length(tempCov)), 
            data.frame(id = 1:length(tempCov), tempCov), sum)
        tempCov[is.na(tempCov)] <- 0
        colnames(tempCov) <- paste(uniqueCov[i], colnames(tempCov), sep = "")
        for (jj in 1:ncol(tempCov)) {
          covList[[colnames(tempCov)[jj]]] <- tempCov[, jj]
          covType <- c(covType, "fixed")
        }
      } else {
        covList[[uniqueCov[i]]] <- tempCov
        covType <- c(covType, "fixed")
      }
    } else {
      covType <- c(covType, "time")
      covList[[uniqueCov[i]]] <- object$timeCovs[[uniqueCov[i]]]
      basicCovList[[uniqueCov[i]]] <- object$timeCovs[[uniqueCov[i]]][, 1]
    }
  }
  covNames <- names(covList)
  names(covType) <- covNames
  return(list(covList = covList, basicCovList = basicCovList, 
          covType = covType, covNames = covNames))
}
