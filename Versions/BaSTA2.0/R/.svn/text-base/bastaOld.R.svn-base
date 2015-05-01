basta <-
function(object, ... ) UseMethod("basta")

# 1. Survival analysis functions:
# 1.1. Mortality and survival functions:
.CalculateMort <- function(model, shape) {
  nIniTh <- if (shape == "simple") 0 else if (shape == "Makeham") 1 else 3
    if (model == "EX") {
      CalcBasicMort <- function(x, theta) theta[, 1]
    } else if (model == "GO") {
      CalcBasicMort <- function(x, theta) {
        exp(theta[, nIniTh + 1] + theta[, nIniTh + 2] * x)
      }
    } else if (model == "WE") {
      CalcBasicMort <- function(x, theta) {
        theta[, nIniTh + 1] * theta[, nIniTh + 2]^theta[, nIniTh + 1] * 
        x^(theta[, nIniTh + 1] - 1)
      }
	  } else if (model == "LO") {
      CalcBasicMort <- function(x, theta) {
        exp(theta[, nIniTh + 1] + theta[, nIniTh + 2] * x) / 
  				(1 + theta[, nIniTh + 3] * exp(theta[, nIniTh + 1]) / 
          theta[, nIniTh + 2] * (exp(theta[, nIniTh + 2] * x) - 1))
    }
	}
  if (shape == "simple") {
    function(x, theta) {
      CalcBasicMort(x, theta)
    }
	} else if (shape == "Makeham") {
    function(x, theta) {
      theta[, 1] + CalcBasicMort(x, theta)
    }
	} else if (shape == "bathtub") {
    function(x, theta) {
      exp(theta[, 1] - theta[, 2] * x) + theta[, 3] + 
      CalcBasicMort(x, theta)
    }
	}
}

.CalculateSurv <- function(model, shape) {
  nIniTh <- if (shape == "simple") 0 else if (shape == "Makeham") 1 else 3
  if (model == "EX") {
  	CalcBasicSurv <- function(x, theta) exp(- theta[, 1] * x)
	} else if (model == "GO") {
		CalcBasicSurv <- function(x, theta) {
			exp(exp(theta[, nIniTh + 1]) / theta[, nIniTh + 2] * 
        (1 - exp(theta[, nIniTh + 2] * x)))
		}
	} else if (model == "WE") {
		CalcBasicSurv <- function(x, theta) {
			theta[, nIniTh + 1] * theta[, nIniTh + 2]^theta[, nIniTh + 1] * 
        x^(theta[, nIniTh + 1] - 1)
		}
	} else if (model == "LO") {
		CalcBasicSurv <- function(x, theta) {
			(1 + theta[, nIniTh + 3] * exp(theta[, nIniTh + 1]) / 
        theta[, nIniTh + 2] * 
        (exp(theta[, nIniTh + 2] * x) - 1))^(-1 / theta[, nIniTh + 3])
		}
	}
	if (shape == "simple") {
		function(x, theta) {
			CalcBasicSurv(x, theta)
		}
	} else if (shape == "Makeham") {
		function(x, theta) {
			exp(-theta[, 1] * x) * 
						CalcBasicSurv(x, theta)
		}
	} else if (shape == "bathtub") {
		function(x, theta) {
			exp(exp(theta[, 1]) / theta[, 2] * (exp(-theta[, 2] * x) - 1) - 
									theta[, 3] * x) * 
						CalcBasicSurv(x, theta)
		}
	}
}

# 1.2. Main survival analysis function:
.CalcAges <- function(covars, ...) UseMethod(".CalcAges")

.CalcAges.noCovs <- function(covars, b, d) {
  x <- d - b
  idTr <- which(b < studyStart)
  xTr <- x * 0
  xTr[idTr] <- studyStart - b[idTr]
  return(list(x = x, xTr = xTr, idTr = idTr))
}

.CalcAges.fixedCovs <- function(covars, b, d) {
  x <- d - b
  idTr <- which(b < studyStart)
  xTr <- x * 0
  xTr[idTr] <- studyStart - b[idTr]
  return(list(x = x, xTr = xTr, idTr = idTr))
}

.CalcAges.timeCovs <- function(covars, b, d) {
  aliveMat <- .MakeAliveMatrix(b, d + 1, covars$timeMat)
  x <- (aliveMat %*% covars$cumMat - 1) * aliveMat
  x[x<0] = 0
  idTr <- which(b < studyStart)
  xTr <- x * 0
  xTr[idTr] <- studyStart - b[idTr]
  idb <- (b - timeCovStart) * n + 1:n
  idd <- (d - timeCovStart) * n + 1:n
  return(list(x = x, xTr = xTr, idTr = idTr, aliveMat = aliveMat, 
          idb = idb, idd = idd))
}

.CalcAges.bothCovs <- function(covars, b, d) {
  aliveMat <- .MakeAliveMatrix(b, d + 1, covars$timeMat)
  x <- (aliveMat %*% covars$cumMat - 1) * aliveMat
  x[x<0] = 0
  idTr <- which(b < studyStart)
  xTr <- x * 0
  xTr[idTr] <- studyStart - b[idTr]
  idb <- (b - timeCovStart) * n + 1:n
  idd <- (d - timeCovStart) * n + 1:n
  return(list(x = x, xTr = xTr, idTr = idTr, aliveMat = aliveMat, 
          idb = idb, idd = idd))
}



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

# 1.3. Function to calculate lower bounds for 'c' parameter 
# (only relevant for 'Makeham' and 'bathtub' shapes):
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

# 1.4. Matrix of yearly indicators for time alive:
.MakeAliveMatrix <- function(f, l, Tm) {
  Fm <- Tm - f
	Fm[Fm >= 0] <- 1
	Fm[Fm < 0] <- 0
	Lm <- Tm - l
	Lm[Lm <= 0] <- -1
	Lm[Lm > 0]  <- 0
	return(Fm * (-Lm))	
}

# 2. Parameter functions:
# 2.1. Default values for theta parameters:
.SetDefaultTheta  <- function(model, shape) {
	if (model == "EX") {
		basicTheta <- list(length = 1, 
				low = -Inf, 
				start = 0.01, 
				jump = 0.005, 
				prior = 0.01,
				name = "b")
	} else if (model == "GO") {
		basicTheta <- list(length = 2, 
				low = c(-Inf, -Inf), 
				start = c(-3, 0.01), 
				jump = c(0.05, 0.025), 
				prior = c(-3, 0.01),
				name = paste("b", (1:2) - 1, sep=""))
	} else if (model == "WE") {
		basicTheta <- list(length = 2, 
				low = c(0, 0), 
				start = c(1.5, 0.1), 
				jump = c(0.01, 0.001), 
				prior = c(1, 0.01),
				name = paste("b", (1:2) - 1, sep=""))
	} else if (model == "LO") {
		basicTheta <- list(length = 3, 
				low = c(-Inf, 0, 0), 
				start = c(-3, 0.01, 0.0001), 
				jump = c(0.001, 0.001, 0.001), 
				prior = c(-3, 0.01, 1e-10),
				name = paste("b", (1:3) - 1, sep=""))
	}
	
	if (shape == "simple") {
		defaultTheta  <- basicTheta
	} else if (shape == "Makeham") {
		defaultTheta  <- list(length = basicTheta$length + 1, 
				low = c(-Inf, basicTheta$low), 
				start = c(0, basicTheta$start), 
				jump = c(0.01, basicTheta$jump), 
				prior = c(0, basicTheta$prior),
				name = c("c", basicTheta$name))
		if (model == "GO") {
			defaultTheta$low <- c(-Inf, -Inf, 0)
		}
	} else if (shape == "bathtub") {
		defaultTheta  <- list(length = basicTheta$length + 3, 
				low = c(-Inf, 0, -Inf, basicTheta$low), 
				start = c(-0.1, 0.5, 0, basicTheta$start), 
				jump = c(0.001, 0.001, 0.01, basicTheta$jump), 
				prior = c(-2, 0.01, 0, basicTheta$prior),
				name = c("a0", "a1", "c", basicTheta$name))
		if (model == "GO") {
			defaultTheta$low <- c(-Inf, 0, -Inf, -Inf, 0)
		}
	}
  attr(defaultTheta, "model") = model
  attr(defaultTheta, "shape") = shape
	return(defaultTheta)
}

# 2.2. Check if input values for parameters, jumps and priors are consistent:
# a) Check and returns mortality parameters:
.CheckParsMort <- function(userPars, covars, 
    model, shape) {
  if (!("defaultTheta" %in% ls())) {
    defaultTheta <- .SetDefaultTheta(model, shape)
  }
  parNames <- c("start", "jump", "prior")
  parId <- c("starting", "jumps", "priors")
  parValAllInMort <- c(0, 0.001, 0)
  parMat <- list()
  par <- 0
  Error <- FALSE
  while(!Error & par <= 2) {
    par <- par + 1
    if (is.null(userPars[[par]])) {
      par.mat <- matrix(defaultTheta[[parNames[par]]], covars$inMort$length, 
          defaultTheta$le, byrow = TRUE, 
          dimnames = list(covars$inMort$name, 
              defaultTheta$name))
      if (class(covars)[2] %in% c("inMort", "fused")) {
        covType <- .FindCovariateType(covars$inMort$mat)
        if (!is.null(covType$cont)) {
          par.mat[covType$cont, ] = parValAllInMort[par]
        }
      }
    } else {
      lengthPar <- length(userPars[[par]])
      if (!is.element(lengthPar, 
          c(defaultTheta$length, defaultTheta$length * 
                  covars$inMort$length))) {
        par.mat <- NULL
        if (class(covars)[2] != "inMort") {
          stop(paste("\nIncorrect length or dimensions for ", parId[par],
                  " parameters\nfor the mortality model. Provide a single\n",
                  " vector of length ",
                  defaultTheta$length, " (number of parameters\n",
                  " for model = '", model,"' and shape = '",shape, "').",
                  sep = ""), call. = FALSE)
        } else {
          if (!is.null(dim(userPars[[par]]))) {
            stop(paste("\nDimensions of ", parId[par], 
                    " matrix for the mortality model\nparameters ",
                    "are incorrect. Provide a single vector\nof",
                    " length ", defaultTheta$length, " or a matrix of",
                    " dimensions ", covars$inMort$length ," times ", 
                    defaultTheta$length, "\n(i.e. number of fixed covariates ", 
                    "in mortality section times number\nof parameters",
                    " for model = '", model,"' and shape = '",shape, "').",
                    sep = ""), call. = FALSE)
          } else {
            stop(paste("\nLength of ", parId[par], 
                    " vector for the mortality model\nparameters is",
                    " incorrect. Provide a single vector\nof",
                    " length ", defaultTheta$length, " or a matrix of",
                    " dimensions ", covars$inMort$length ," times ", 
                    defaultTheta$length, "\n(i.e. number of fixed covariates ", 
                    "in mortality section times number\nof parameters",
                    " for model = '", model,"' and shape = '",shape, "').",
                    sep=""), call. = FALSE)
          }
       }
       Error <- TRUE
     } else {
        if (!is.null(dim(userPars))) {
          par.mat <- userPars[[par]]
        } else {
          par.mat <- matrix(userPars[[par]], covars$inMort$length, 
              defaultTheta$length, byrow = TRUE)
        }
        dimnames(par.mat) <- list(covars$inMort$names, 
            defaultTheta$name)
      }
    }
    parMat[[parNames[par]]] <- par.mat
  }
  return(parMat)
}

# b) Check and return proportional hazards parameters:
.CheckParsPH <- function(userPars, covars){
  if(class(covars)[1] == "noCovs" | 
      (class(covars)[1] == "fixedCovs" & class(covars)[2] == "inMort")) {
    parList <- NULL
  } else {
    parNames <- c("start", "jumps", "priors")
    parId <- c("starting", "jump", "prior")
    parList <- list()
    if (class(covars)[1] == "fixedCovs") {
      covTypeVec <- 1
    } else if (class(covars)[1] == "timeCovs") {
      covTypeVec <- 2
    } else {
      if (class(covars)[2] %in% c("propHaz", "fused")) {
        covTypeVec <- 1:2
      } else {
        covTypeVec <- 2
      }
    }
    for (covt in covTypeVec) {
      nGam <- covars[[covt + 1]]$length
      namesCovs <- covars[[covt + 1]]$name
      if (is.null(namesCovs)) {
        namesCovs <- paste("Cov", 1:nGam, sep = "")
      }
      parMat <- matrix(0, 3, nGam, dimnames = list(parNames, namesCovs))
      for (par in 1:3) {
        if (is.null(userPars[[covt]][[par]])) {
          parVec <- rep(c(0, 0.001, 0)[par], nGam)
        } else {
          if (length(userPars[[covt]][[par]]) != nGam) {
            stop(paste("\nLength of ", parId[par], 
                    " parameters for ", c("fixed", "time")[covt], 
                    " covariates\nin prop. hazards section is not ",
                    "equal to number\nof covariates (n = ", nGam, ").", 
                    sep = ""), call. = FALSE)
          } else {
            parVec <- userPars[[covt]][[par]]
          }
        }
        parMat[par, ] = parVec
      }
      parList[[c("fixed", "time")[covt]]] = parMat
    }
  }
  return(parList)
}

# Construct final parameter matrices:
# a) For fixed covariates:
.CalcFixPars <- function(object, ...) UseMethod(".CalcFixPars")

.CalcFixPars.fused <- function(object, params) {
  thetaMat <- object$inMort$mat %*% params$theta
  gaMatFix <- c(object$prHazFix$mat %*% params$gamma$fixed)
  return(list(theta = thetaMat, gamma = gaMatFix))
}

.CalcFixPars.inMort <- function(object, params) {
  thetaMat <- object$inMort$mat %*% params$theta
  return(list(theta = thetaMat, gamma = 0))
}

.CalcFixPars.propHaz <- function(object, params) {
  thetaMat <- matrix(params$theta, 1, length(params$theta), 
      dimnames = list(NULL, names(params$theta)))
  gaMatFix <- c(object$prHazFix$mat %*% params$gamma$fixed)
  return(list(theta = thetaMat, gamma = gaMatFix))
}


# b) all parameters:
.CalcParMats <- function(object, ...) UseMethod(".CalcParMats")

.CalcParMats.noCovs <- function(object, params) {
  parmat <- list(theta = matrix(params$theta, 1, length(params$theta), 
      dimnames = list(NULL, names(params$theta))))
  class(parmat) <- paste("par", class(object)[1], sep = "")
  return(parmat)
}

.CalcParMats.fixedCovs <- function(object, params) {
  parmat<- .CalcFixPars(object, params)
  class(parmat) <- paste("par", class(object)[1], sep = "")
  return(parmat)
}

.CalcParMats.timeCovs <- function(object, params) {
  thetaMat <- object$inMort$mat %*% params$theta
  gamMat <- object$prHazTime$mat[, , 1] * 0
  for(ga in 1:object$prHazTime$length){
    gamMat <- gamMat + object$prHazTime$mat[, , ga] * params$gamma$time[ga]
  }
  parmat <- list(theta = thetaMat, gamma = gamMat)
  class(parmat) <- paste("par", class(object)[1], sep = "")
  return(parmat)
}

.CalcParMats.bothCovs <- function(object, params) {
  parList <- .CalcFixPars(object, params)
  gamMat <- object$prHazTime$mat[, , 1] * 0
  for(ga in 1:object$prHazTime$length){
    gamMat <- gamMat + object$prHazTime$mat[, , ga] * params$gamma$time[ga]
  }
  gamMat <- gamMat + parList$gamma
  parmat <- list(theta = parList$theta, gamma = gamMat)
  class(parmat) <- paste("par", class(object)[1], sep = "")
  return(parmat)
}

# 2.3. Final parameter names:
.DefineParNames <- function(defaultTheta, covars, recaptTrans) {
  nameTheta  <- paste(rep(defaultTheta$name, 
                              each = covars$inMort$length), 
                          "[", rep(colnames(covars$inMortMat),
                                   defaultTheta$length), 
                          "]", sep = "")
  if (covars$inMort$length == 1) {
    nameTheta <- defaultTheta$name
	}
	nThetaFull <- length(nameTheta)
	nameFixedGamma <- paste("gamma[", colnames(covars$prHazFix$mat), "]", sep="")
	if (covars$nGamFixed == 1) {
		nameFixedGamma  <- "gamma"
	}
  if (class(covars) == "fixedCovs") {
    nameTimeGamma <- NULL
  } else {
    nameTimeGamma <- colnames(covars$prHazTime$mat)
  }
	namePi <- ifelse(length(recaptTrans) == 1, "pi", 
			paste("pi[", recaptTrans, "]", sep = ""))
	namePost <- c("post[theta,gamma]", "post[X0]", "post[full]")
  return(list(theta = nameTheta,
              gamFixed = nameFixedGamma,
              gamTime = nameTimeGamma,
              pi = namePi,
              post = namePost))
}


# 3. Data handling functions:
.CreateBastaData <- function(object, studyStart, fixedCovs = NULL, 
    timeCovs = NULL, ...) {
  bd <- as.matrix(object[, 2:3])
  colnames(bd) <- c("birth", "death")
  Y <- as.matrix(object[, 4:ncol(object)])
  nt <- ncol(Y)
  colnames(Y) <- 1:nt + studyStart - 1
  bastadat <- list (bd = bd,
      Y = Y,
      fixedCovs = fixedCovs,
      timeCovs = timeCovs,
      n = nrow(object), 
      studyEnd = nt + studyStart - 1,
      years = 1:nt + studyStart - 1,
      nYears = nt)
  if (is.null(fixedCovs) & is.null(timeCovs)) {
    class(bastadat) <- "datNoCovs"
  } else {
    if (is.null(timeCovs)) {
      class(bastadat) <- "datFixedCovs"
    } else if (is.null(fixedCovs)) {
      class(bastadat) <- "datTimeCovs"
    } else {
      class(bastadat) <- "datBothCovs"
    }
  }
  return(bastadat)
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

# 3.2. Build final list of covariates and their attributes:
.PrepCovs <- function(object, ...) UseMethod(".PrepCovs")

.PrepCovs.datNoCovs <- function(object, ...) {
  covars <- 0
  class(covars) <- c("noCovs", "noClass")
}

.PrepCovs.datFixedCovs <- function(object, covarsStruct) {
  covars <- .PrepCovsFixed(object, covarsStruct)
  covClass <- class(covars)
  class(covars) <- c("fixedCovs", covClass)
  return(covars)
}

.PrepCovs.datBothCovs <- function(object, covarsStruct) {
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


# 4. Error checking:
.FindErrors <- function(model, shape, covarsStruct, niter, burnin, thinning) {
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


