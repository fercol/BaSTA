DefineParams <- function(object, fixedCovs = NULL, timeCovs = NULL, 
    covars = NULL, model = "GO", shape = "simple", 
    covarsStruct = "fused", default = TRUE) {
  if (substr(class(object), 1, 3) != "dat") {
    object <- .CreateBastaData(object, fixedCovs, timeCovs)
  }
  parCats <- c("start", "jump", "prior")
  parCatsText <- c("starting", "jump", "prior")
  if (is.null(covars)) {
    covars <- .PrepCovs(object, covarsStruct)
  }
  defaultTheta <- .SetDefaultTheta(model, shape)
  thetaPars <- .CheckParsMort(list(start = NULL, jump = NULL, prior = NULL), 
      covars, model, shape)
  gammaPars <- .CheckParsPH(list(fixed = list(start = NULL, jump = NULL, 
              prior = NULL), time = list(start = NULL, jump = NULL, 
              prior = NULL)), covars)
  # Use default values for all parameters:
  if (!default) {
    cat("\n----------------------------------------------------------\n",
        "With this function you can define starting values, jumps\n",
        "and priors for the parameters in your BaSTA analysis.\n",
        "(If you wish to stop this function simply click 'esc').",
        "\n----------------------------------------------------------\n")
    # a) theta:
    cat("\n1. THETA PARAMETERS (IN MORTALITY):", 
        "\n-----------------------------------")
    specTheta <- .AskYesNo("Use default theta (mortality) parameters")
    if(specTheta == 'n') {
      if(class(covars)[2] %in% c("fused", "inMort")) {
        diffTheta <- .AskYesNo(paste("Use the same parameters for covariates",
                paste(paste(covars$inMort$names[-covars$inMort$length], 
                        collapse = ", "), 
                    covars$inMort$names[covars$inMort$length], 
                    sep = " and ")))
        if (diffTheta == 'y') {
          dTh <- 1
        } else {
          dTh <- covars$inMort$length
        }
      } else {
        dTh <- 1
      }
      for(par in 1:3) {
        specThPar <- .AskYesNo(paste("Use default", parCatsText[par], 
                "values"))
        if (specThPar == 'n') {
          if (dTh == 1) {
            parTh <- .SpecifyParNum(defaultTheta$length, parCatsText[par],
                ":")
            parMat <- matrix(parTh, covars$inMort$length,
                defaultTheta$length,
                dimnames = dimnames(thetaPars[[parCats[par]]]),
                byrow = TRUE)
          } else {
            parMat <- matrix(0, covars$inMort$length,
                defaultTheta$length,
                dimnames = dimnames(thetaPars[[parCats[par]]]),
                byrow = TRUE)
            for (i in 1:dTh) {
              covName <- paste(" for ", covars$inMort$names[i], ":",
                  sep = "")
              parTh <- .SpecifyParNum(defaultTheta$length, parCatsText[par],
                  covName)
              parMat[i, ] <- parTh[1:defaultTheta$length]
            }
          }
          thetaPars[[parCats[par]]] <- parMat
        }
      }
    }
    # b) gamma:
    if (class(covars)[1] == "noCovs" | 
        (class(covars)[1] == "fixedCovs" & class(covars)[2] == "inMort")) {
      cat("\nNo proportional hazards parameters are ",
          "required for this model.")
    } else {
      cat("\n2. GAMMA PARAMETERS (PROPORTIONAL HAZARDS):", 
          "\n-------------------------------------------")
      if (class(covars)[2] %in% c("fused", "propHaz")) {
        if (class(covars)[1] %in% c("bothCovs", "timeCovs")) 
          cat("\n2.1 GAMMA PARAMETERS FOR FIXED COVARIATES:")
        gammaPars$fixed <- .DefineGamma("fixed", 
            covars$prHazFix$length, gammaPars$fixed, parCatsText)
      }
      if (class(covars)[1] %in% c("bothCovs", "timeCovs")) {
        if (class(covars)[2] %in% c("fused", "propHaz")) 
          cat("\n2.2 GAMMA PARAMETERS FOR TIME CHANGING COVARIATES:")
        gammaPars$time <- .DefineGamma("time changing", 
            covars$prHazTime$length, gammaPars$time, parCatsText)
      }
    }
    cat("\nList of parameters built.\n")
  }
  parList <- list(theta = thetaPars, gamma = gammaPars)
  class(parList) <- "bastaParams"
  return(parList)
}

.SpecifyParNum <- function(length, cat, name) {
  parNum <- readline(paste("Type the ", 
          length, " ", cat, 
          " paramters separated by comas", name, sep = ""))
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
    parNum <- readline(paste("Type the ", 
            length, " ", cat, 
            " paramters separated by comas", name, sep = ""))
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
