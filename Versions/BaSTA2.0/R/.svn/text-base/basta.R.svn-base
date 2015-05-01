# TODO: Add comment
# 
# Author: fernando
###############################################################################
basta <-
    function(object, ... ) UseMethod("basta")

# 0.- FIND ERRORS (unfinished):
.FindErrors <- function(model, shape, covarsStruct, 
    niter, burnin, thinning) {
  # a) Data errors:
#	data.check  <- DataCheck(object, studyStart, autofix = rep(0, 7), 
#      silent = TRUE)
#	if (!data.check[[1]]) {
#		stop("\nYou have an error in Dataframe 'object',\nplease use function ",
#				"'DataCheck'\n", call. = FALSE)
#	}
  
  # b) Check that niter, burnin, and thinning are compatible.
  if (burnin > niter) {
    stop("Object 'burnin' larger than 'niter'.", call. = FALSE)
  }
  if (thinning > niter) {
    stop("Object 'thinning' larger than 'niter'.", call. = FALSE)
  }
  
  # c) Model type, shape and covariate structure:
  if (!is.element(model, c("EX", "GO", "WE", "LO"))) {
    stop("Model misspecification: specify available models", 
        " (i.e. 'EX', 'GO', 'WE' or 'LO')\n", call. = FALSE)
  }
  if (!is.element(shape, c("simple", "Makeham", "bathtub"))) {
    stop("shape misspecification. Appropriate arguments are:", 
        " 'simple', 'Makeham' or 'bathtub'.\n", call. = FALSE)
  }
  if (!is.element(covarsStruct, c("fused", "prop.haz", "all.in.mort"))) {
    stop("Covariate structure misspecification. Appropriate arguments are:", 
        " 'fused', 'prop.haz' or 'all.in.mort'.\n", call. = FALSE)
  }
  if (model == "EX" & shape != "simple") {
    stop("Model misspecification: EX model can only be fitted with a", 
        " simple shape", call. = FALSE)
  }
  if (model == "EX" & covarsStruct != "fused") {
    stop("Model misspecification: EX model can only be fitted with a", 
        " fused covariate structure", call. = FALSE)
  }
  if (covarsStruct == "all.in.mort" & sum(model == "GO", shape == 
          "simple") < 2) {
    stop("Model misspecification: all.in.mort is only available with", 
        " Gompertz (GO) models and simple shape.", call. = FALSE)
  }
}

# 1.- PREPARE MORTALITY DATA AND COVARIATES:
# 1.1.- Mortality Data:
.PrepMortDat <- 
    function(object) UseMethod(".PrepMortDat")

.PrepMortDat.fullCens <- function(object) {
  # 1. Extract times of birth and death:
  bi <- object$birthDeath[, 1]
  di <- object$birthDeath[, 2]
  classAge <- ifelse(any(di - bi != round(di - bi)), 
      "contAge", "discAge")
  
  # 2. Construct data list:
  mortDat <- list(birth = bi, death = di, N = object$N)
  class(mortDat) <- c(classAge, "fullCens")
  return(mortDat)
}

.PrepMortDat.census <- function(object) {
  # 1. Extract times of birth and death:
  bi <- object$birthDeath[, 1]
  di <- object$birthDeath[, 2]
  classAge <- ifelse(any(di - bi != round(di - bi)), 
      "contAge", "discAge")
  
  # 2. Calculate first and last time observed 
  #    and total number of times observed:
  firstObs <- object$capHistMat[, 1]
  lastObs <- object$capHistMat[, 2]
  
  # 3. Construct data list and define class:
  idNoB <- which(bi == 0)
  idNoD <- which(di == 0)
  mortDat <- list(birth = bi, death = di, N = object$N, 
      firstObs = firstObs, lastObs = lastObs, 
      start = object$studyStart, end = object$studyEnd, 
      idNoB = idNoB, idNoD = idNoD)
  class(mortDat) <- c(classAge, "census")
  return(mortDat)
}

.PrepMortDat.capRecov <- function(object) {
  # 1. Extract times of birth and death:
  bi <- object$birthDeath[, 1]
  di <- object$birthDeath[, 2]
  classAge <- ifelse(any(di - bi != round(di - bi)), 
      "contAge", "discAge")
  
  # 2. Calculate first and last time observed 
  #    and total number of times observed:
  firstObs <- object$capHistMat[, 1]
  lastObs <- object$capHistMat[, 2]
  idNoB <- which(bi == 0)
  idNoD <- which(di == 0)
  mortDat <- list(birth = bi, death = di, N = object$N, 
      firstObs = firstObs, lastObs = lastObs, 
      start = object$studyStart, end = object$studyEnd, 
      idNoB = idNoB, idNoD = idNoD)
  class(mortDat) <- c(classAge, "capRecov") 
  return(mortDat)
}

.PrepMortDat.capRecap <- function(object) {
  # 1. Extract times of birth and death:
  bi <- object$birthDeath[, 1]
  di <- object$birthDeath[, 2]
  classAge <- ifelse(any(di - bi != round(di - bi)), 
      "contAge", "discAge")
  
  # 2. Calculate first and last time observed 
  #    and total number of times observed:
  Y <- as.matrix(object$capHistMat)
  nTimes <- ncol(Y)
  times <- object$studyStart + 1:nTimes - 1
  ytemp <- t(t(Y) * times)
  lastObs <- c(apply(ytemp, 1, max))
  ytemp[ytemp == 0] <- max(times) * 2
  firstObs <- c(apply(ytemp, 1, min))
  firstObs[firstObs == max(times) * 2]  <- 0
  nCaps <- Y %*% rep(1, nTimes)
  
  # 3. Define study duration:
  Dx <- (times[2] - times[1])
  Tm <- matrix(times, object$N, nTimes, byrow=TRUE)
  knownFirstObs <- firstObs
  id <- which(bi > 0 & bi >= object$studyStart)
  knownFirstObs[id] <- bi[id] + 1
  id <- which(bi > 0 & bi < object$studyStart)
  knownFirstObs[id]  <- object$studyStart
  knownLastObs <- lastObs
  id <- which(di > 0 & di <= object$studyEnd)
  knownLastObs[id] <- di[id] - 1
  id <- which(di > 0 & di > object$studyEnd)
  knownLastObs[id] <- object$studyEnd
  detectIdentMat <- .MakeIndicatorMat(firstObs, lastObs, 
      Tm)
  idNoB <- which(bi == 0)
  idNoD <- which(di == 0)
  mortDat <- list(birth = bi, death = di, N = object$N, 
      firstObs = firstObs, lastObs = lastObs, Y = Y, 
      nCaps = nCaps, Dx = Dx, studyTimeMat = Tm, 
      detecMat = detectIdentMat, start = object$studyStart, 
      end = object$studyEnd, idNoB = idNoB, idNoD = idNoD)
  class(mortDat) <- c(classAge, "capRecap") 
  return(mortDat)
}

# 1.2.- List of Covariates:
# Note: rename to ".CreateCovObj"
.CreateCovList <- function(object, form = NULL) {
  if (class(object)[2] == "noCovs") {
    covars <- list(names = 'nc', length = 1, class = "intercept")
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

# 2.- DEFINE PARAMETERS:
# 2.1.- User Defined or Default Parameter Values:
.DefineParams <- function(covars, model, shape, covarsStruct, 
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

.SetDefaultTheta  <- function(model, shape) {
  if (model == "EX") {
    nTh <- 1
    startTh <- -2 
    jumpTh <- 0.01 
    priorMean <- -4
    priorSd <- 0.1
    nameTh <- "b"
    indExp <- "b"
  } else if (model == "GO") {
    nTh <- 2 
    startTh <- c(-2, 0.01) 
    jumpTh <- c(0.05, 0.025)
    priorMean <- c(-3, 0.01)
    priorSd <- c(0.1, 0.001)
    nameTh <- c("b0", "b1")
    indExp <- c()
  } else if (model == "WE") {
    nTh <- 2
    startTh <- c(0.4, -2.5) 
    jumpTh <- c(0.05, 0.01)
    priorMean <- c(0.2, -4)
    priorSd <- c(0.01, 0.1)
    nameTh <- c("b0", "b1")
    indExp = c("b0", "b1")
  } else if (model == "LO") {
    nTh <- 3 
    startTh <- c(-3, -4, -10) 
    jumpTh <- c(0.001, 0.001, 0.001) 
    priorMean <- c(-3, -4, -20)
    priorSd <- c(0.1, 0.1, 0.5)
    nameTh <- c("b0", "b1", "b2")
    indExp <- c("b1", "b2")
  }
  if (shape == "Makeham") {
    nTh <- nTh + 1 
    startTh <- c(0, startTh) 
    jumpTh <- c(0.01, jumpTh) 
    priorMean <- c(0, priorMean)
    priorSd <- c(0.1, priorSd)
    nameTh <- c("c", nameTh)
  } else if (shape == "bathtub") {
    nTh <- nTh + 3 
    startTh <- c(-0.1, -0.6, 0, startTh)
    jumpTh <- c(0.001, 0.001, 0.01, jumpTh) 
    priorMean <- c(-2, -2, 0, priorMean)
    priorSd <- c(0.1, 0.1, 0.1, priorSd)
    nameTh <- c("a0", "a1", "c", nameTh)
    indExp <- c("a1", indExp)
  }
  defaultTheta  <- list(length = nTh, start = startTh, jump = jumpTh, 
      priorMean = priorMean, priorSd = priorSd, name = nameTh, 
      indExp = indExp)
  attr(defaultTheta, "model") = model
  attr(defaultTheta, "shape") = shape
  return(defaultTheta)
}

# 2.2.- Extract Parameters for Analysis:
.DefineParsG <- function(params) {
  parsG <- list(theta = params$theta$start, indExp = params$theta$indExp)
  if (class(params) == "propHaz") {
    parsG$gamma <- params$gamma$start
  }
  class(parsG) <- class(params)
  return(parsG)
}

# 3.- DEFINE MORTALITY AND SURVIVAL FUNCTIONAL FORMS:
.DefineSimpleMort <- function(model) {
  if (model == "EX") {
    function(x, theta) theta$b0
  } else if (model == "GO") {
    CalcSimpleMort <- function(x, theta) {
      exp(theta$b0 + theta$b1 * x)
    }
  } else if (model == "WE") {
    function(x, theta) {
      theta$b0 * theta$b1^theta$b0 * 
          x^(theta$b0 - 1)
    }
  } else if (model == "LO") {
    function(x, theta) {
      exp(theta$b0 + theta$b1 * x) / 
          (1 + theta$b2 * exp(theta$b0) / 
            theta$b1 * (exp(theta$b1 * x) - 1))
    }
  }
}

.DefineMort <- function(shape, CalcSimpleMort) {
  if (shape == "simple") {
    function(x, theta) {
      CalcSimpleMort(x, theta)
    }
  } else if (shape == "Makeham") {
    function(x, theta) {
      theta$c + CalcSimpleMort(x, theta)
    }
  } else if (shape == "bathtub") {
    function(x, theta) {
      exp(theta$a0 - theta$a1 * x) + theta$c + 
          CalcSimplecMort(x, theta)
    }
  }
}

.DefineSimpleSurv <- function(model, shape) {
  if (model == "EX") {
    function(x, theta) exp(- theta[, 1] * x)
  } else if (model == "GO") {
    CalcBasicSurv <- function(x, theta) {
      exp(exp(theta$b0) / theta$b1 * 
              (1 - exp(theta$b1 * x)))
    }
  } else if (model == "WE") {
    function(x, theta) {
      exp(-(theta$b1 * x)^theta$b0)
    }
  } else if (model == "LO") {
    function(x, theta) {
      (1 + theta$b2 * exp(theta$b0) / theta$b1 * 
            (exp(theta$b1 * x) - 1))^(-1 / theta$b2)
    }
  }
}

.DefineSurv <- function(shape, CalcSimpleSurv) {
  if (shape == "simple") {
    function(x, theta) {
      CalcSimpleSurv(x, theta)
    }
  } else if (shape == "Makeham") {
    function(x, theta) {
      exp(-theta$c * x) * 
          CalcSimpleSurv(x, theta)
    }
  } else if (shape == "bathtub") {
    function(x, theta) {
      exp(exp(theta$a0) / theta$a1 * (exp(-theta$a0 * x) - 1) - 
              theta$c * x) * CalcSimpleSurv(x, theta)
    }
  }
}

# 4.- FUNCTION TO CALCULATE PARAMETER VALUES WITH OR WITHOUT COVARIATES:
.CalcParams <- function(params, ...) UseMethod(".CalcParams")

.CalcParams.propHaz <- function(params, covars) {
  thetaPars <- params$theta
  if (!is.null(params$indExp)) {
    thetaPars[params$indExp] <- exp(thetaPars[params$indExp])
  }
  gammaPars <- .CalcGamma(covars, params$gamma)
  parList <- list(theta = data.frame(thetaPars), gamma = gammaPars)
  class(parList) <- class(params)
  return(parList)
}

.CalcParams.inMort <- function(params, covars) {
  thetaPars <- .CalcThetaInMort(covars, params$theta, params$indExp)
  parList <- list(theta = .CalcThetaInMort(covars, params$theta, 
          params$indExp))
  class(parList) <- class(params)
  return(parList)
}

# 4.1.- Calculate Thetas:
.CalcThetaInMort <- function(covars, ...) UseMethod(".CalcThetaInMort")

.CalcThetaInMort.noCovs <- function(covars, theta, indExp) {
  if (!is.null(indExp)) {
    theta[, indExp] <- exp(theta[, indExp])
  }
  return(theta)
}

.CalcThetaInMort.fixedCovs <- function(covars, theta, indExp) {
  theta <- as.matrix(theta)
  thetaPars <- covars$fixedCovs %*% theta
  if ("c" %in% colnames(theta)) {
    cVec <- theta[1, "c"] * 
        exp(covars$fixedCovs[, -1] %*% theta[-1, "c"])
    thetaPars[, "c"] <- cVec
  }
  if (!is.null(indExp)) {
    thetaPars[, indExp] <- exp(thetaPars[, indExp])
  }
  thetaPars <- data.frame(thetaPars)
  return(thetaPars)
}

.CalcThetaInMort.timeCovs <- function(covars, theta, indExp) {
  thetaPars <- list()
  thetaNames <- colnames(theta)
  for (iTh in thetaNames[which(thetaNames != "c")]) {
    thetaPars[[iTh]] <- covars$timeCovs[[1]] * 0
    for (jCov in 1:length(covars$timeCovs)) {
      thetaPars[[iTh]] <- thetaPars[[iTh]] + 
          covars$timeCovs[[jCov]] * theta[jCov, iTh]
    }
    if (iTh %in% indExp) {
      thetaPars[[thName]] <- exp(thetaPars[[thName]])
    }
  }
  if ("c" %in% thetaNames) {
    thetaPars[["c"]] <- covars$timeCovs[[1]] * 0
    for (jCov in 2:length(covars$timeCovs)) {
      thetaPars[["c"]] <- thetaPars[["c"]] + 
          covars$timeCovs[[jCov]] * theta[jCov, "c"]
    }
    thetaPars[["c"]] <- theta[1, "c"] * exp(thetaPars[["c"]])
  }
  return(thetaPars)
}

.CalcThetaInMort.bothCovs <- function(covars, theta, indExp) {
  thetaPars <- list()
  thetaNames <- colnames(theta)
  fThetaPars <- covars$fixedCovs %*% 
      as.matrix(theta[colnames(covars$fixedCovs), ])
  if ("c" %in% thetaNames) {
    fThetaPars[, "c"] <- covars$fixedCovs[, -1] %*% 
        as.matrix(theta[colnames$fixedCovs[-1], "c"])
  }
  for (iTh in thetaNames[which(thetaNames != "c")]) {
    thName <- colnames(theta)[iTh]
    thetaPars[[iTh]] <- covars$timeCovs[[1]] * 0 + fThetaPars[, iTh]
    tCovNames <- names(covars$timeCovs)
    for (jCov in tCovNames) {
      thetaPars[[iTh]] <- thetaPars[[iTh]] + 
          covars$timeCovs[[jCov]] * theta[jCov, iTh]
    }
    if (iTh %in% indExp) {
      thetaPars[[iTh]] <- exp(thetaPars[[iTh]])
    }
  }
  if ("c" %in% thetaNames) {
    thetaPars[["c"]] <- covars$timeCovs[[1]] * 0 + fThetaPars[, "c"]
    for (jCov in 1:length(covars$timeCovs)) {
      thetaPars[["c"]] <- thetaPars[["c"]] + 
          covars$timeCovs[[jCov]] * theta[jCov, "c"]
    }
    thetaPars[["c"]] <- theta[1, "c"] * exp(thetaPars[["c"]])
  }
  return(thetaPars)
}

# 4.2.- Calculate Gammas:
.CalcGamma <- function(covars, ...) UseMethod(".CalcGamma")

.CalcGamma.fixedCovs <- function(covars, gamma) {
  gammaPars <- covars$fixedCovs %*% c(gamma)
  return(exp(gammaPars))
}

.CalcGamma.timeCovs <- function(covars, gamma) {
  gammaPars <- covars$timeCovs[[1]] * 0
  for (iGa in names(covars$timeCovs)) {
    gammaPars <- gammaPars + covars$timeCovs[[iGa]] * gamma[iGa, ]
  }
  return(exp(gammaPars))
}

.CalcGamma.bothCovs <- function(covars, gamma) {
  fGammaPars <- covars$fixedCovs %*% c(gamma[colnames(covars$fixedCovs), ])
  gammaPars <- covars$timeCovs[[1]] * 0 + c(fGammaPars)
  for (iGa in names(covars$timeCovs)) {
    gammaPars <- gammaPars + covars$timeCovs[[iGa]] * gamma[iGa, ]
  }
  return(exp(gammaPars))
}

# 5.- FUNCTIONS TO CALCULATE PRIOR AGE DISTRIBUTION:
.SetPriorAgeDist <- function(params) {
  priorRange <- list()
  priorMean <- params$theta$priorMean[1, ]
  priorSd <- params$theta$priorSd[1, ]
  nPriors <- ncol(params$theta$priorMean)
  Dprior <- rep(0, nPriors)
  names(Dprior) <- colnames(priorMean)
  priorRange <- lapply(1:nPriors, function(th) 
        seq(qnorm(0.001, mean = priorMean[1, th], sd = priorSd[1, th]),
            qnorm(0.999, mean = priorMean[1, th], sd = priorSd[1, th]), 
            length = 100))
  Dprior <- sapply(1:nPriors, function(th) priorRange[[th]][2] -
            priorRange[[1]][1])
  
  priorCombs <- expand.grid(priorRange)
  if (!is.null(params$indExp)) {
    priorCombs[, params$indExp] <- exp(priorCombs[, params$indExp])
  }
  colnames(priorCombs) <- colnames(priorMean)
  dx   <- 0.5
  xv   <- seq(0, 1000, dx)
  meanLifeExpec <- sapply(1:nrow(priorCombs), 
      function(ii) 
        sum(CalcSurv(xv, priorCombs[ii, ]) * dx))
  
  priorDens <- priorCombs * 0
  for(i in 1:nPriors) {
    priorDens[, i] <- dnorm(priorCombs[, i], mean = priorMean[1, i], 
        sd = priorSd[1, i])
  }
  
  xvdisc <- 0:1000
  xMat <- matrix(xvdisc, nrow(priorCombs), length(xvdisc), byrow = TRUE)
  vxMat <- CalcSurv(xMat, priorCombs) / meanLifeExpec * 
      apply(priorDens, 1, prod) 
  priorAgeDist <- apply(vxMat, 2, sum) * prod(Dprior)
  priorAgeDist <- priorAgeDist / sum(priorAgeDist)
  return(priorAgeDist)
}

.CalcPriorAgeDist <- function(x, priorAgeDist) priorAgeDistr[x + 1]

# 6.- SET INITIAL VALUES FOR TIMES OF BIRTH AND DEATH, AGES, ETC.:
# 6.1.- Main method:
.SetIniTimeDat <- function(covars, ...) UseMethod(".SetIniTimeDat") 

.SetIniTimeDat.noCovs <- function(covars, mortdat, minAge) {
  iniTimesDat <- .SetIniTimes(mortdat, minAge)
  class(iniTimesDat) <- c(class(iniTimesDat)[2:3], class(covars)[2])
  return(iniTimesDat)
}

.SetIniTimeDat.fixedCovs <- .SetIniTimeDat.noCovs

.SetIniTimeDat.timeCovs <- function(covars, mortdat, minAge) {
  iniTimesDat <- .SetIniTimes(mortdat, minAge)
  idb <- (iniTimesDat$times$birth - covars$start) * 
      mortdat$N + 1:mortdat$N
  idd <- (iniTimesDat$times$death - covars$start) * 
      mortdat$N + 1:mortdat$N
  iniTimesDat$times <- data.frame(cbind(iniTimesDat$times, idb, idd))
  timeMat <- .MakeIndicatorMat(iniTimesDat$time$birth, iniTimesDat$time$death, 
      covars$timeMat)
  iniTimesDat$timeMat <- timeMat
  class(iniTimesDat) <- c(class(iniTimesDat)[2:3], class(covars)[2])
  return(iniTimesDat)
}

.SetIniTimeDat.bothCovs <- .SetIniTimeDat.timeCovs

# 6.2.- Method to set initial times of birth and death, 
#       and detection indicator matrix for capRecap data.
.SetIniTimes <- function(mortdat, ...) UseMethod(".SetIniTimes")

.SetIniTimes.fullCens <- function(mortdat, minAge) {
  idTr <- mortdat$birth * 0
  idTr[mortdat$birth < mortdat$start] <- 1
  iniTimes <- data.frame(birth = mortdat$birth, death = mortdat$death, 
      age = mortdat$death - mortdat$birth, idTr = idTr)
  if (minAge > 0) {
    iniTimes <- .SplitByMinAgeIni(iniTimes, mortdat, minAge)
    iniTimesClass <- c(class(mortdat), "minAge")
  } else {
    iniTimesClass <- c(class(mortdat), "noMinAge")
  }
  iniTimesList <- list(times = iniTimes)
  class(iniTimesList) <- iniTimesClass
  return(iniTimesList)
}

.SetIniTimes.census <- function(mortdat, minAge) {
  iniTimes <- .SetIniBirthDeath(mortdat)
  if (minAge > 0) {
    iniTimes <- .SplitByMinAgeIni(iniTimes, mortdat, minAge)
    iniTimesClass <- c(class(mortdat), "minAge")
  } else {
    iniTimesClass <- c(class(mortdat), "noMinAge")
  }
  iniTimes$idNew <- c()
  class(iniTimes) <- iniTimesClass
  return(iniTimes)
}

.SetIniTimes.capRecov <- .SetIniTimes.census

.SetIniTimes.capRecap <- function(mortdat, minAge) {
  iniTimes <- .SetIniBirthDeath(mortdat)
  if (minAge > 0) {
    iniTimes <- .SplitByMinAgeIni(iniTimes, mortdat, minAge)
    iniTimesClass <- c(class(mortdat), "minAge")
  } else {
    iniTimesClass <- c(class(mortdat), "noMinAge")
  }
  firstTimeInStudy <- c(apply(cbind(mortdat$start, iniTimes$times$birth + 1), 
          1, max))
  lastTimeInStudy <- c(apply(cbind(mortdat$end, iniTimes$times$death - 1), 
          1, min))
  aliveIdentMat <- .MakeIndicatorMat(firstTimeInStudy, lastTimeInStudy, 
      mortdat$studyTimeMat)
  iniTimes$idNew <- c()
  iniTimes$aliveMat <- aliveIdentMat
  class(iniTimes) <- iniTimesClass
  return(iniTimes)
}

# 6.3.- Function to set initial times of birth and death:
.SetIniBirthDeath <- function(mortdat) {
  birth <- mortdat$birth
  if (length(mortdat$idNoB) > 0) {
    idFirstObs <- mortdat$idNoB[mortdat$firstObs[mortdat$idNoB] > 0]
    if (length(idFirstObs) > 0) {
      birth[idFirstObs] <- mortdat$firstObs[idFirstObs] - 1
    }
    idDeath <- mortdat$idNoB[mortdat$firstObs[mortdat$idNoB] == 0 & 
            mortdat$death[mortdat$idNoB] > 0]
    if (length(idDeath) > 0) {
      birth[idDeath] <- mortdat$death[idDeath] - 1
    }
  }
  death <- mortdat$death
  if (length(mortdat$idNoD) > 0) {
    idLastObs <- mortdat$idNoD[mortdat$lastObs[mortdat$idNoD] > 0]
    if (length(idLastObs) > 0) {
      death[idLastObs] <- mortdat$lastObs[idLastObs] + 1
    }
    idBirth <- mortdat$idNoD[mortdat$lastObs[mortdat$idNoD] == 0]
    if (length(idBirth) > 0) {
      death[idBirth] <- birth[idBirth] + 1
    }
  }
  idTr <- which(birth < mortdat$start)
  iniBD <- list(times = data.frame(birth = birth, death = death, 
          ages = death - birth), idTr = idTr)
  return(iniBD)
}

# 6.4.- Function to separate juveniles from adults:
.SplitByMinAgeIni <- function(times, mortdat, minAge) {
  baseTimes <- times$times
  indAd <- baseTimes$birth * 0
  indJu <- indAd
  ageAdTrunc <- indAd
  ageJuTrunc <- indAd
  indAd[baseTimes$age >= minAge] <- 1
  indJu[baseTimes$age < minAge] <- 1
  ageAd <- baseTimes$age - minAge
  ageAd[baseTimes$age < minAge] <- 0
  ageJu <- baseTimes$age
  ageJu[baseTimes$age > minAge] <- minAge
  idTrAd <- which(baseTimes$birth + minAge < mortdat$start)
  if (length(idTrAd) > 0) {
    ageAdTrunc[idTrAd] <- mortdat$start - 
        (baseTimes$birth[idTrAd] + minAge)
  }
  idTrJu <- which(baseTimes$birth < mortdat$start &
              mortdat$start - baseTimes$birth < minAge)
  if (length(idTrJu) > 0) {
    ageJuTrunc[idTrJu] <- mortdat$start - baseTimes$birth[idTrJu]
  }
  timesMinAge <- data.frame(baseTimes, ageAd = ageAd, ageJu = ageJu, 
      ageAdTr = ageAdTrunc, ageJuTr = ageJuTrunc, indAd = indAd, indJu = indJu)
  times$times <- timesMinAge
  times$idTrAd <- idTrAd
  times$idTrJu <- idTrJu
  return(times)
}

# 7.- UPDATE TIMES OF BIRTH & DEATH:
# 7.1.- Main function:
.SampleTimeDat <- function(times, ...) UseMethod(".SampleTimeDat") 

.SampleTimeDat.noCovs <- function(times, mortdat, minAge, covars) {
  times <- .UpdateTimes(times, mortdat, minAge)
  return(times)
}

.SampleTimeDat.fixedCovs <- .SampleTimeDat.noCovs

.SampleTimeDat.timeCovs <- function(times, mortdat, minAge, covars) {
  times <- .UpdateTimes(times, mortdat, minAge)
  if (length(times$idNew) > 0) {
    times$times$idb[times$idNew] <- (times$times$birth[times$idNew] - 
          covars$start) * mortdat$N + c(1:mortdat$N)[times$idNew]
    times$times$idd[times$idNew] <- (times$times$death[times$idNew] - 
          covars$start) * mortdat$N + c(1:mortdat$N)[times$idNew]
    times$timeMat[times$idNew, ] <- 
        .MakeIndicatorMat(times$times$birth[times$idNew], 
            times$times$death[times$idNew], covars$timeMat[times$idNew, ])
  }
  return(times)
}

.SampleTimeDat.bothCovs <- .SampleTimeDat.timeCovs

# 7.2.- Functions to update additional data based on proposed ages:
.UpdateTimes <- function(times, ...) UseMethod(".UpdateTimes")

.UpdateTimes.fullCens <- function(times, ...) {
  return(times)
}

.UpdateTimes.census <- function(times, mortdat, minAge) {
  times <- .ProposeTimes(times, mortdat, minAge)
  return(times)
}

.UpdateTimes.capRecov <- .UpdateTimes.census

.UpdateTimes.capRecap <- function(times, mortdat, minAge) {
  times <- .ProposeTimes(times, mortdat, minAge)
  if (length(times$idNew) > 0) {
    firstTimeInStudy <- c(apply(cbind(mortdat$start, 
                times$times$birth[times$idNew] + 1), 1, max))
    lastTimeInStudy <- c(apply(cbind(mortdat$end, 
                times$times$death[times$idNew] - 1), 1, min))
    times$aliveMat[times$idNew, ] <- .MakeIndicatorMat(firstTimeInStudy, 
        lastTimeInStudy, mortdat$studyTimeMat[times$idNew, ])
  }
  return(times)
}

# 7.3.- Functions to sample ages and recalculate min ages:
.ProposeTimes <- function(times, ...) UseMethod(".ProposeTimes")

.ProposeTimes.noMinAge <- function(times, mortdat, ...) {
  newTimes <- .ProposeAges(times, mortdat)
  return(newTimes)
}

.ProposeTimes.minAge <- function(times, mortdat, minAge) {
  newTimes <- .ProposeAges(times, mortdat)
  newTimes <- .SplitByMinAge(newTimes, mortdat, minAge)
  return(newTimes)
}

# 7.4.- Function to propose ages:
.ProposeAges <- function(times, mortdat) {
  if (length(mortdat$idNoB) > 0) {
    birthNew <- times$times$birth
    birthNew[mortdat$idNoB] <- times$times$birth[mortdat$idNoB] + sample(-1:1, 
        length(mortdat$idNoB), replace = TRUE)
    idCap <- mortdat$idNoB[which(mortdat$nCaps[mortdat$idNoB] > 0)]
    birthNew[idCap] <- apply(cbind(birthNew[idCap], 
            mortDat$firstObs[idCap] - 1), 1, min)
    idNCap <- mortdat$idNoB[which(mortdat$nCaps[mortdat$idNoB] == 0)]
    birthNew[idNCap] <- apply(cbind(birthNew[idNCap], 
            mortDat$death[idNCap] - 1), 1, min)
    idNew <- mortdat$idNoB[which(times$times$birth[mortdat$idNoB] != 
                birthNew[mortdat$idNoB])]
    times$times$birth[mortdat$idNoB] <- birthNew[mortdat$idNoB]
  }
  if (length(mortdat$idNoD) > 0) {
    deathNew <- times$times$death
    deathNew[mortdat$idNoD] <- times$times$death[mortdat$idNoD] + sample(-1:1,
        length(mortdat$idNoD), replace = TRUE)
    idCh <- mortdat$idNoD[which(deathNew[mortdat$idNoD] != 
                times$death[mortdat$idNoD])]
    deathNew[idCh] <- apply(cbind(deathNew[idCh], birthNew[idCh], 
            mortdat$lastObs[idCh] + 1), 1, max)
    idNew <- sort(unique(c(idNew, 
                mortdat$idNoD[which(times$times$death[mortdat$idNoD] != 
                            deathNew[mortdat$idNoD])])))
    times$times$death[mortdat$idNoD] <- deathNew[mortdat$idNoD]
  }
  if (length(idNew) > 0) {
    times$idNew <- idNew
    times$times$ages[idNew] <- deathNew[idNew] - birthNew[idNew]
  }
  return(times)
}

# 7.5.- Function to update adults and juveniles:
.SplitByMinAge <- function(times, mortdat, minAge) {
  if (length(times$idNew) > 0) {
    splitTimes <- times$times[times$idNew, ]
    splitTimes[, c("ageAd", "ageJu", "ageAdTr", "ageJuTr", 
            "indAd", "indJu")] <- 0
    splitTimes$indAd[splitTimes$ages >= minAge] <- 1
    splitTimes$indJu[splitTimes$ages < minAge] <- 1
    splitTimes$ageAd <- splitTimes$ages - minAge
    splitTimes$ageAd[splitTimes$ages < minAge] <- 0
    splitTimes$ageJu <- splitTimes$ages
    splitTimes$ageJu[splitTimes$ages > minAge] <- minAge
    idTrAdNew <- which(splitTimes$birth + minAge < mortdat$start)
    splitTimes$ageAdTr[idTrAdNew] <- mortdat$start - 
        (splitTimes$birth[idTrAdNew] + minAge)
    idTrJuNew <- which(splitTimes$birth < mortdat$start &
            mortdat$start - splitTimes$birth < minAge)
    splitTimes$ageJuTr[idTrJuNew] <- mortdat$start - 
        splitTimes$birth[idTrJuNew]
    times$times[times$idNew, ] <- splitTimes
    times$idTrAd <- 
        sort(unique(c(times$idTrAd[!(times$idTrAd %in% times$idNew)], 
                    idTrAdNew)))
    times$idTrJu <- 
        sort(unique(c(times$idTrJu[!(times$idTrJu %in% times$idNew)], 
                    idTrJuNew)))
  }
  return(times)
}

# CALCULATE INDIVIDUAL SURVIVAL AND MORTALITY (unfinished):
.CalcMortFuns <- function(pars, ...) UseMethod(".CalcMortFuns")

.CalcMortFuns.parnoCovs <- function(parMats, ages) {
#  hazRate <- CalcMort(ages$x, parMats$theta)
  intCensMort <- CalcSurv(ages$x, parMats$theta) - 
      CalcSurv(ages$x + Dx, parMats$theta)
  survProbTr <- CalcSurv(ages$xTr, parMats$theta)
  survProbTr[-ages$idTr] <- 1
  return(cbind(intCensMort, survProbTr))
}

.CalcMortFuns.parfixedCovs <- function(parMats, ages) {
#  hazRate <- CalcMort(ages$x, parMats$theta) * exp(parMats$gamma)
  intCensMort <- (CalcuSurv(ages$x, parMats$theta) - 
        CalcuSurv(ages$x + Dx, parMats$theta))^exp(parMats$gamma)
  survProbTr <- CalcSurv(ages$xTr, parMats$theta)^exp(parMats$gamma)
  survProbTr[-ages$idTr] <- 1
  return(cbind(intCensMort, survProbTr))
}

.CalcMortFuns.parbothCovs <- function(parMats, ages) {
  hazRate <- CalcMort(ages$x, parMats$theta) * exp(parMats$gamma)
  hazRate <- (hazRate - hazRate[ages$idb]) * ages$aliveMat
  survProb <- exp(-hazRate %*% rep(1, ncol(ages$aliveMat)))
  survProbTr <- exp(-(hazRate * ages$xTr) %*% rep(1, ncol(ages$aliveMat)))
  return(cbind(intCensMort, survProbTr))
}

.CalcMortFuns.partimeCovs <- function(parMats, ages) {
  hazRate <- CalcMort(ages$x, parMats$theta) * exp(parMats$gamma)
  hazRate <- (hazRate - hazRate[ages$idb]) * ages$aliveMat
  survProb <- exp(-hazRate %*% rep(1, ncol(ages$aliveMat)))
  survProbTr <- exp(-(hazRate * ages$xTr) %*% rep(1, ncol(ages$aliveMat)))
  return(cbind(hazRate[ages$idd], survProb, survProbTr))
}

# MAKE ALIVE MATRIX (unfinished):
.MakeIndicatorMat <- function(first, last, timeMat) {
  firstMat <- timeMat - first
  firstMat[firstMat >= 0] <- 1
  firstMat[firstMat < 0] <- 0
  lastMat <- timeMat - last
  lastMat[lastMat <= 0] <- -1
  lastMat[lastMat > 0]  <- 0
  return(firstMat * (-lastMat))	
}

# CALCULATE LOWER BOUND FOR c PARAMETER 
# (unfinished and possibly deprecated):
.CalculateLowC              <- function(theta) {
  if (shape == "Makeham") {
    if (model == "GO") {
      c.low <- ifelse(theta[3] > 0, -exp(theta[2]), 0)
    } else if (model == "WE") {
      c.low <- 0
    } else if (model == "LO") {
      c.low <- ifelse(theta[2] > theta[3] * exp(theta[1]), 
          -exp(theta[2]), 0)
    } 
  }
  if (shape == "bathtub") {
    if (model == "GO") {
      x.min <- (theta[1] + log(theta[2]) - theta[4] - 
            log(theta[5])) / (theta[2] + theta[5])
    } else if (model == "LO" | model == "WE") {
      x.vec <- seq(0, 100, 0.1)
      mort  <- .CalculateMort(x.vec, matrix(theta, 
              length(x.vec), lengthTheta, 
              byrow = TRUE), 0)
      x.min <- x.vec[which(mort == min(mort))[1]]
    }
    c.low <- -exp(theta[1] - theta[2] * x.min) - 
        .CalculateBasicMort(x.min, 
            matrix(theta[-c(1:3)], 
                1, lengthBasicTheta))
  }
  return(c.low)
}

# DYNAMIC UPDATE OF JUMP SD'S (verify):
UpdateJumps <- function(jObject, updateVec, targetUpdate, g, 
    updateInt, nPar, updateLen) {
  gUpdate <- which(updateInt == g)
  if (nPar > 1) {
    parCount <- seq(nPar, nPar * 100, nPar) - gUpdate
    parCount <- nPar - parCount[which(parCount >= 0)][1]
  } else {
    parCount <- 1
  }
  updateRate <- sum(updateVec[g + c(-(updateLen - 1):0)]) / updateLen
  if (updateRate == 0) updateRate <- 1e-3
  if (gUpdate == 1) {
    jObject$jumpsMat <- jObject$startJump
    jObject$updateRateVec <- updateRate
  } else {
    jObject$jumpsMat <- rbind(jObject$jumpsMat, jObject$jump)
    jObject$updateRateVec <- c(jObject$updateRateVec, updateRate)
  }
  if (gUpdate > nPar) {
    idTarget <- which(jObject$updateRateVec[gUpdate - (nPar - 1):0] > 
            targetUpdate * 0.9 &
            jObject$updateRateVec[gUpdate - (nPar - 1):0] < 
            targetUpdate * 1.1)
  } else {
    idTarget <- 0
  }
  if (length(idTarget) < 3) {
    if (parCount == 1) {
      updateDiff <- abs(targetUpdate - jObject$shortUpdVec)
      jObject$updateOrder <- sort.int(updateDiff, index.return = TRUE)$ix
    }
    if (gUpdate <= nPar) {
      jObject$jump <- jObject$startJump
    } 
    if (gUpdate > 1) {
      jObject$shortUpdVec[jObject$parNum] <- updateRate
    }
    jObject$parNum <- jObject$updateOrder[parCount]
    jObject$jump[jObject$parNum] <- jObject$jump[jObject$parNum] * 
        updateRate / targetUpdate
    jObject$gUpdate <- gUpdate
  } else {
    jObject$update <- FALSE
  }
  return(jObject)
}




