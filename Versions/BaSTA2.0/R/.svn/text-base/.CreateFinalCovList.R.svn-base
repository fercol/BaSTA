# TODO: Add comment
# 
# Author: fernando
###############################################################################

.CreateCovList <- function(object, form = NULL) {
  if (class(object)[2] == "noCovs") {
    covars <- list(names = NULL, length = 0, class = "intercept")
    class(covars) <- c("covars", "noCovs")
  } else {
    if (is.null(form)) {
      if (class(object)[2] == "fixedCovs") {
        form <- paste("~", paste(colnames(object$fixedCovs), collapse = "+"))
      } else if (class(object)[2] == "timeCovs") {
        form <- paste("~", paste(names(object$timeCovs), collapse = "+"))
      } else {
        form <- as.formula(paste("~", paste(c(colnames(object$fixedCovs), 
                    names(object$timeCovs)), collapse = "+")))
      }
    }
    .VerifyCovsInForm(object, form)
    covars <- .CreateModelList(object, x = form)
    covars$names <- c(colnames(covars$fixedCovs), names(covars$timeCovs))
    covars$length <- length(covars$names)
    covars$class <- .FindCovarsClass(covars)
    class(covars) <- c("covars", class(covars))
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
}

.CreateModelList <- function(object, ...) UseMethod(".CreateModelList")

.CreateModelList.fixedCovs <- function(object, x) {
  fixedCovsMat <- model.matrix(object = x, data = object$fixedCovs)
  allCovList <- list(fixedCovs = cbind(fixedCovsMat))
  class(allCovList) <- "fixedCovs"
  return(allCovList)
}

.CreateModelList.timeCovs <- function(object, x) {
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
  ncolTime <- ncol(timeCovList[[1]])
  colnamesTime <- colnames(timeCovList[[1]])
  cumMat <- matrix(t(apply(diag(1, ncolTime, ncolTime), 1, cumsum)), 
      ncolTime, ncolTime, dimnames = list(colnamesTime, colnamesTime))
  timeMat <- t(t(timeCovList[[1]] * 0) + 
          0:(ncolTime - 1) + object$firstCohort)
  allCovList <- list(timeCovs = timeCovList, cumMat = cumMat, 
      timeMat = timeMat)
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
    modelNames <- colnames(model.matrix(object = x, 
            data = data.frame(modFrame$basicCovList)))
    modelNameList <- strsplit(modelNames, split = ":")
    
    fixedCovMat <- matrix(0, object$N, 0)
    timeCovList <- list()
    for (i in 1:length(modelNames)) {
      localCovs <- modelNameList[[i]]
      nLocCovs <- length(localCovs)
      if (nLocCovs == 1) {
        if (localCovs == "(Intercept)") {
          if (("fixed" %in% modFrame$covType)) {
            fixedCovMat <- matrix(1, object$N, ncol(fixedCovMat) + 1)
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
          fixedInt <- rep(1, object$N)
          for(jj in 1:nLocCovs) {
            fixedInt <- fixedInt * modFrame$covList[[localCovs[jj]]]
          }
          fixedCovMat <- cbind(fixedCovMat, fixedInt)
          colnames(fixedCovMat) <- c(colnames(fixedCovMat)[-ncol(fixedCovMat)], 
              modelNames[i])
        }
      }
    }
    if (length(timeCovList) > 0) {
      ncolTime <- ncol(timeCovList[[1]])
      colnamesTime <- colnames(timeCovList[[1]])
      cumMat <- matrix(t(apply(diag(1, ncolTime, ncolTime), 1, cumsum)), 
          ncolTime, ncolTime, dimnames = list(colnamesTime, colnamesTime))
      timeMat <- t(t(timeCovList[[1]] * 0) + 
              0:(ncolTime - 1) + object$firstCohort)
    }
    if (ncol(fixedCovMat) == 0) {
      allCovList <- list(timeCovs = timeCovList, cumMat = cumMat, 
          timeMat = timeMat)
      class(allCovList) <- "fixedCovs"
    } else {
      allCovList <- list(fixedCovs = fixedCovMat, timeCovs = timeCovList, 
          cumMat = cumMat, timeMat = timeMat)
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

.FindCovarsClass <- function(covars) {
  covClass <- rep("categorical", covars$length)
  names(covClass) <- covars$names
  if (class(covars) %in% c("fixedCovs", "bothCovs")) {
    nUniCovs  <- apply(covars$fixedCovs, 2, function(x) length(unique(x)))
    idint <- which(nUniCovs == 1)
    idcon <- which(nUniCovs > 2)
    if (length(idint) > 0) covClass[names(idint)] <- "intercept"
    if (length(idcon) > 0) covClass[names(idcon)] <- "numeric"
  }
  if (class(covars) %in% c("timeCovs", "bothCovs")) {
    nTimeCovs <- length(covars$timeCovs)
    nUniCovs  <- sapply(1:nTimeCovs, function(x) 
          length(unique(covars$timeCovs[[x]])))
    names(nUniCovs) <- names(covars$timeCovs)
    idint <- which(nUniCovs == 1)
    idcon <- which(nUniCovs > 2)
    if (length(idint) > 0) covClass[names(idint)] <- "intercept"
    if (length(idcon) > 0) covClass[names(idcon)] <- "numeric"
  }
  return(covClass)
}
