.CreateCovList <- function(object, form = NULL, cohortCov) {
  if (class(object)[2] == "noCovs") {
    covars <- list(names = 'nc', length = 1, class = "intercept")
    class(covars) <- c("covars", "noCovs")
  } else {
    if (is.null(form)) {
      if (class(object)[2] == "fixedCovs") {
        if (is.null(cohortCov)) {
          form <- as.formula(paste("~", paste(colnames(object$fixedCovs), 
                      collapse = "+")))
        } else {
          noCoh <- which(colnames(object$fixed) != cohortCov)
          reordCovs <- colnames(object$fixedCovs)[noCoh]
          form <- as.formula(paste("~", paste(reordCovs, collapse = "+")))
        }
      } else if (class(object)[2] == "timeCovs") {
        form <- as.formula(paste("~", paste(names(object$timeCovs), 
                    collapse = "+")))
      } else {
        if (is.null(cohortCov)) {
          form <- as.formula(paste("~", paste(c(colnames(object$fixedCovs), 
                          names(object$timeCovs)), collapse = "+")))
        } else {
          noCoh <- which(colnames(object$fixed) != cohortCov)
          reordCovs <- c(colnames(object$fixedCovs)[noCoh], 
              names(object$timeCovs))
          form <- as.formula(paste("~", paste(reordCovs, collapse = "+")))
        }
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
      timeMat = timeMat, start = object$firstCohort)
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
              1:ncolTime - 1 + object$firstCohort)
    }
    if (ncol(fixedCovMat) == 0) {
      allCovList <- list(timeCovs = timeCovList, cumMat = cumMat, 
          timeMat = timeMat)
      class(allCovList) <- "fixedCovs"
    } else {
      allCovList <- list(fixedCovs = fixedCovMat, timeCovs = timeCovList, 
          cumMat = cumMat, timeMat = timeMat, start = object$firstCohort)
      class(allCovList) <- "bothCovs"
    }
  }
  return(allCovList)
}

.ExtractModelFrame <- function(object, x) {
  labelsCov <- labels(terms(x))
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




.DefineParams <- function(covars, times, model, shape, covarsStruct, 
    default = TRUE, cohortCov) {
  defaultTheta <- .SetDefaultTheta(model, shape)
  if (!default) {
    cat("\n----------------------------------------------------------\n",
        "Here you will be asked to define starting values, jumps\n",
        "and priors for the parameters in your BaSTA analysis.\n",
        "(If you wish to stop this function simply click 'esc').",
        "\n----------------------------------------------------------\n")
    # a) theta:
    cat("\n1. THETA PARAMETERS (IN MORTALITY):", 
        "\n-----------------------------------")
  }
  thetaPars <- .AssignPars(covars, defaultTheta, "theta", 
      covarsStruct, default)
  parList <- list(theta = thetaPars)
  if (covarsStruct == "prop.haz" & class(covars)[2] != "noCovs") {
    if (!default) {
      cat("\n2. GAMMA PARAMETERS (PROPORTIONAL HAZARDS):", 
          "\n-------------------------------------------")
    }
    defaultGamma <- list(start = 0.1, jump = 0.01, priorMean = 0, 
        priorSd = 0.001, length = 1)
    gammaPars <- .AssignPars(covars, defaultGamma, "gamma", 
        covarsStruct, default)
    parList$gamma <- gammaPars
  }
  parList$theta$indExp <- defaultTheta$indExp
  class(parList) <- ifelse(class(covars)[2] == "noCovs", 
      "noCovs", ifelse(covarsStruct == "prop.haz", "propHaz", "inMort"))
  return(parList)
}

.AssignPars <- function(covars, defaultPar, genParName, 
    covarsStruct, default) {
  parCats <- c("start", "jump", "priorMean", "priorSd")
  parCatsText <- c("starting", "jump", "prior mean", "prior sd")
  if (genParName == "theta") {
    if (covarsStruct == "all.in.mort") {
      parMat <- data.frame(matrix(0, covars$length, defaultPar$length, 
              dimnames = list(covars$names, defaultPar$name)))
      nCovRows <- covars$length
    } else {
      parMat <- data.frame(matrix(0, 1, defaultPar$length, 
              dimnames = list(NULL, defaultPar$name)))
      nCovRows <- 1
    }
  } else if (covarsStruct == "prop.haz" | class(covars)[2] == "noCovs"){
    parMat <- matrix(0, covars$length, 1, 
        dimnames = list(covars$names, "gamma"))
    nCovRows <- covars$length
  }
  parList <- list(start = parMat, jump = parMat, priorMean = parMat, 
      priorSd = parMat)  
  if (!default) {
    specAllPar <- .AskYesNo(paste("Use default ",
            genParName, " parameters", sep = ""))
  } else {
    specAllPar <- "y"
  }
  for (nCovs in 1:nCovRows) {
    if (specAllPar == "n") {
      if (nCovRows == 1) {
        specCovPar <- "n"
      } else {
        specCovPar <- .AskYesNo(paste("Use default values for", 
                covars$names[nCovs]))
      }
    } else {
      specCovPar <- "y"
    }
    for (par in 1:4) {
      if (specCovPar == "n") {
        specPar <- .AskYesNo(paste("Use default", parCatsText[par], 
                "values"))
      } else {
        specPar <- "y"
      }
      if (specPar == "y") {
        if (covars$class[nCovs] == "intercept") {
          indPar <- 1
        } else if (covars$class[nCovs] == "categorical") {
          indPar <- ifelse(par %in% c(1, 3), 
              ifelse("intercept" %in% covars$class, 0, 1), 1)
        } else {
          indPar <- ifelse(par %in% c(1, 3), 0, 1)
        }
        parList[[parCats[par]]][nCovs, ] <- defaultPar[[parCats[par]]] * 
            indPar
      } else {
        parUser <- .SpecifyParNum(defaultPar$length, parCatsText[par],
            ":")
        parList[[par]][nCovs, ] <- parUser
      }
    }
  }
  return(parList)
}

.SpecifyParNum <- function(length, cat, name) {
  parNum <- readline(paste("Type ", 
          length, " ", cat, 
          " parameter(s) (use comas between parameters)", name, sep = ""))
  op <- options()
  options(warn = -1)
  parNum <- as.numeric(unlist(strsplit(parNum, split = ",")))
  options(op)
  while (length(parNum) != length | any(is.na(parNum))) {
    if (any(is.na(parNum))) {
      cat("\nSpecify numerical values.\n")
    } else {
      cat("\nWrong number of parameters.\n")
    }
    parNum <- readline(paste("Type ", 
            length, " ", cat, 
            " parameter(s) (use comas between parameters)", name, sep = ""))
    op <- options()
    options(warn = -1)
    parNum <- as.numeric(unlist(strsplit(parNum, split = ",")))
  }
  return(parNum)
}

.AskYesNo <- function(question) {
  answer <- readline(paste("\n", question, "? (y or n): ", sep = ""))
  while(!is.element(answer, c("y", "n"))) {
    answer <- readline(paste("Please, type 'y' or 'n' or use the esc",
            "key to terminate: "))
  }
  return(answer)
}

.DefineGamma <- function(type, length, parmat, parCatsText) {
  default <- .AskYesNo(paste("Use default values for", 
          type, "covariates"))
  if (default == 'n') {
    for (par in 1:3) {
      defaultPar <- .AskYesNo(paste("Use default", parCatsText[par], 
              "values"))
      if (defaultPar == 'n') {
        parGa <- .SpecifyParNum(length, 
            parCatsText[par], ":")
        parmat[par, ] <- parGa
      }
    }
  }
  return(parmat)
}
