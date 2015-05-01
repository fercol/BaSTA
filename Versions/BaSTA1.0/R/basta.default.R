basta.default <- 
function(object, studyStart, studyEnd, minAge = 0, model = "GO", 
         shape = "simple", covarsStruct = "fused", niter = 50000, 
         burnin = 5001, thinning = 50, recaptTrans = studyStart, 
         thetaStart = NULL, thetaJumps = NULL, thetaPriors = NULL, 
         gammaStart = NULL, gammaJumps = NULL, gammaPriors = NULL, 
         nsim = 1, parallel = FALSE, ncpus = 2, lifeTable = TRUE, 
         progrPlots = FALSE, ...) {
         	
  # This function estimates age-specific mortality from capture-recapture/
  # recovery (CRR) data when a large proportion of (or all) the records have
  # unknown times of birth and death. It uses the framework described by
  # Colchero & Clark (2012) Journal of Animal Ecology. The code is organized
  # as follows:
  # - Section 1 (line   24): load packages
  # - Section 2 (line   27): Define functions
  # - Section 3 (line  290): Error checking
  # - Section 4 (line  341): Data formatting
  # - Section 5 (line  449): MCMC prep
  # - Section 6 (line  594): Define multi-MCMC function
  # - Section 7 (line  856): Run multiple BaSTA MCMCs
  # - Section 8 (line  928): Calculate diagnostics
  # - Section 9 (line 1202): Create output object
  
  # 1. Load package msm:
  require(msm)

  # 2. Initial error checking:
  # 2.1 Data errors:
  data.check                 <- DataCheck(object, studyStart, studyEnd, 
                                          silent = TRUE)
  if (!data.check[[1]]) {
    stop("\nYou have an error in Dataframe 'object',\nplease use function ",
         "'DataCheck'\n", call. = FALSE)
  }

  # 2.2 Check that niter, burnin, and thinning are compatible.
  if (burnin > niter) {
    stop("\nObject 'burnin' larger than 'niter'.", call. = FALSE)
  }
  if (thinning > niter) {
    stop("\nObject 'thinning' larger than 'niter'.", call. = FALSE)
  }
    
  # 2.3 Model type, shape and covariate structure:
  if (!is.element(model, c("EX","GO","WE","LO"))) {
    stop("\nModel misspecification: specify available models", 
         "(i.e. 'EX', 'GO', 'WE' or 'LO')\n", call. = FALSE)
  }
  if (!is.element(shape, c("simple","Makeham","bathtub"))) {
    stop("\nshape misspecification. Appropriate arguments are:",
         " 'simple', 'Makeham' or 'bathtub'.\n", call. = FALSE)
  }
  if (!is.element(covarsStruct, c("fused", "prop.haz", "all.in.mort"))) {
    stop("\nCovariate structure misspecification. Appropriate arguments are:",
         " 'fused', 'prop.haz' or 'all.in.mort'.\n", call. = FALSE)
  }
    

  # 3. Functions:
  # 3.1 Survival, mort, pdf:
  # a) Basic model:
  if (model == "EX") {
    CalculateBasicMx         <- function(x, theta) theta
    CalculateBasicSx         <- function(x, theta) exp(- theta * x)
    length.theta0            <- 1
    low.theta0               <- -Inf
    ini.theta0               <- 0.01
    jump.theta0              <- 0.005
    prior.theta0             <- 0.01
  } else if (model == "GO") {
    CalculateBasicMx         <- function(x, theta) {
      exp(theta[, 1] + theta[, 2] * x)
    }
    CalculateBasicSx         <- function(x, theta) {
      exp(exp(theta[, 1]) / theta[, 2] * (1 - exp(theta[, 2] * x)))
    }
    length.theta0            <- 2
    low.theta0               <- c(-Inf, -Inf)
    ini.theta0               <- c(-3, 0.01)
    jump.theta0              <- c(0.05, 0.025)
    prior.theta0             <- c(-3, 0.01)
  } else if (model == "WE") {
    CalculateBasicMx         <- function(x, theta) {
      theta[, 1] * theta[, 2]^theta[, 1] * x^(theta[, 1] - 1)
    }
    CalculateBasicSx         <- function(x, theta) {
      exp(-(theta[, 2] * x)^theta[, 1])
    }
    length.theta0            <- 2
    low.theta0               <- c(0, 0)
    ini.theta0               <- c(1.5, 0.1)
    jump.theta0              <- c(0.01, 0.001)
    prior.theta0             <- c(1, 0.01)
  } else if (model == "LO") {
    CalculateBasicMx         <- function(x, theta) {
      exp(theta[, 1] + theta[, 2] * x) / 
          (1 + theta[, 3] * exp(theta[, 1]) / theta[, 2] * 
          (exp(theta[, 2] * x) - 1))
    }
    CalculateBasicSx         <- function(x, theta) {
      (1 + theta[, 3] * exp(theta[, 1]) / theta[, 2] * 
      (exp(theta[, 2] * x) - 1))^(-1 / theta[, 3])
    }
    length.theta0            <- 3
    low.theta0               <- c(-Inf, 0, 0)
    ini.theta0               <- c(-3, 0.01, 0.0001)
    jump.theta0              <- c(0.001, 0.001, 0.001)
    prior.theta0             <- c(-3, 0.01, 1e-10)
  }
  name.theta0                <- paste("b", (1:length.theta0) - 1, sep="")

  # b) Extended model for different mortality shapes:
  if (model == "EX") {
    shape                    <- "simple"
  }
  if (shape == "simple") {
    CalculateFullMx          <- function(x, theta, gamma) {
      CalculateBasicMx(x, theta) * exp(gamma)
    }
    CalculateFullSx          <- function(x, theta, gamma) {
      CalculateBasicSx(x, theta)^exp(gamma)
    }
    length.theta             <- length.theta0
    low.theta                <- low.theta0
    ini.theta                <- ini.theta0
    jump.theta               <- jump.theta0
    prior.theta              <- prior.theta0
    name.theta               <- name.theta0
  } else if (shape == "Makeham") {
    CalculateFullMx          <- function(x, theta, gamma) {
      (theta[, 1] + CalculateBasicMx(x, matrix(theta[, -1], 
       ncol = length.theta0))) * exp(gamma)
    }
    CalculateFullSx          <- function(x, theta, gamma) {
      (exp(-theta[, 1] * x) * CalculateBasicSx(x, matrix(theta[, -1], 
      ncol = length.theta0)))^exp(gamma)
    }
    length.theta             <- length.theta0 + 1
    ini.theta                <- c(0, ini.theta0)
    jump.theta               <- c(0.01, jump.theta0)
    prior.theta              <- c(0, prior.theta0)
    low.theta                <- c(-Inf, low.theta0)
    if (model == "GO") {
      low.theta              <- c(-Inf, -Inf, 0)
    }
    name.theta               <- c("c", name.theta0)
  } else if (shape == "bathtub") {
    CalculateFullMx          <- function(x, theta, gamma) {
      (exp(theta[, 1] - theta[, 2] * x) + theta[, 3] + 
       CalculateBasicMx(x, matrix(theta[, -c(1:3)], 
       ncol = length.theta0))) * exp(gamma)
    }
    CalculateFullSx          <- function(x, theta, gamma) {
      (exp(exp(theta[, 1]) / theta[, 2] * (exp(-theta[, 2] * x) - 1) - 
       theta[, 3] * x) * CalculateBasicSx(x, matrix(theta[, -c(1:3)], 
       ncol = length.theta0)))^exp(gamma)
    } 
    length.theta             <- length.theta0 + 3
    ini.theta                <- c(-0.1, 0.5, 0, ini.theta0)
    jump.theta               <- c(0.001, 0.001, 0.01, jump.theta0)
    prior.theta              <- c(-2, 0.01, 0, prior.theta0)
    low.theta                <- c(-Inf, 0, -Inf, low.theta0)
    if (model == "GO") {
      low.theta              <- c(-Inf, 0, -Inf, -Inf, 0)
    }
    name.theta               <- c("a0", "a1", "c", name.theta0)
  }

  CalculateFullFx            <- function(x, theta, gamma) {
    CalculateFullMx(x, theta, gamma) * CalculateFullSx(x, theta, gamma)
  }

  CalculateMultiSx           <- function(theta) {
    CalculateFullSx(xv, matrix(theta, ncol = length.theta), gaa)
  }
                        
  CalculateMultiMx           <- function(theta) {
    CalculateFullMx(xv, matrix(theta, ncol = length.theta), gaa)
  }

  # 3.2 Covariate type (i.e. categorical and continuous):
  FindCovariateType          <- function(Z) {
    # This functions finds and returns if an intercecpt was included 
    # and which covariates are categorical or continuous.
    if (!is.null(Z)) {
      lu                     <- apply(Z, 2, function(x) length(unique(x)))
      ru                     <- apply(Z, 2, range)
      idcat                  <- which(lu == 2 & apply(ru, 2, sum) == 1)
      if (length(idcat) == 0) {
        idcat                <- NULL
      }
      idint                  <- which(lu == 1)
      if (length(idint) == 0) {
        idint                <- NULL
      }
      idcon                  <- which(lu > 2)
      if (length(idcon) == 0) {
        idcon                <- NULL
      }
    } else {
      idcat                  <- NULL
      idint                  <- NULL
      idcon                  <- NULL
    }
		
    return(list(int  = idint, 
                cat  = idcat, 
                cont = idcon))
  }

  # 3.3 Matrix of indicators for time alive:
  BuildAliveMatrix           <- function(f, l, Tm) {
    Fm             <- Tm - f
    Fm[Fm >= 0]    <- 1
    Fm[Fm < 0]     <- 0
    Lm             <- Tm - l
    Lm[Lm <= 0]    <- -1
    Lm[Lm > 0]     <- 0
    return(Fm * (-Lm))	
  }

  # 3.4 Calculate lower bound for 'c' parameter:
  CalculateLowC              <- function(theta) {
    if (shape == "Makeham") {
      if (model == "GO") {
        c.low                <- ifelse(theta[3] > 0, -exp(theta[2]), 0)
      } else if (model == "WE") {
        c.low                <- 0
      } else if (model == "LO") {
        c.low                <- ifelse(theta[2] > theta[3] * exp(theta[1]), 
                                       -exp(theta[2]), 0)
      } 
    }
    if (shape == "bathtub") {
      if (model == "GO") {
        x.min                <- (theta[1] + log(theta[2]) - theta[4] - 
                                log(theta[5])) / (theta[2] + theta[5])
      } else if (model == "LO" | model == "WE") {
        x.vec                <- seq(0, 100, 0.1)
        mort                 <- CalculateFullMx(x.vec, matrix(theta, 
                                                length(x.vec), length.theta, 
                                                byrow = TRUE), 0)
        x.min                <- x.vec[which(mort == min(mort))[1]]
      }
      c.low                  <- -exp(theta[1] - theta[2] * x.min) - 
                                     CalculateBasicMx(x.min, 
                                    matrix(theta[-c(1:3)], 
                                    1, length.theta0))
    }
    return(c.low)
  }

  # 3.5 Check if input values for parameters, jumps and priors are consistent:
  # a) Mortality:
  CheckParsMort              <- function(par, user.par, par.name) {
    if (is.null(user.par)) {
      par.mat           <- matrix(par, length.cat, length.theta, byrow = TRUE, 
                              dimnames = dimnames(low.full.theta))
    } else {
      length.par        <- length(user.par)
      if (!is.element(length.par, c(length.theta, length.theta * length.cat))) {
        par.mat         <- NULL
        if (!is.null(dim(user.par))) {
          stop(paste("\nDimensions of ", par.name, 
                     " matrix for the mortality model parameters are incorrect.",
                      "\n\nProvide a single vector of length ", length.theta,
                      " or a matrix of dimensions ", length.cat ," times ", 
                      length.theta, ".(i.e. number of categorical ", 
                      "covariates times number of parameters for model ", 
                      model," with ", shape, " shape).", sep = ""), 
                      call. = FALSE)
        } else if (is.null(dim(user.par))) {
          stop(paste("\nLength of ", par.name, 
                     " vector for the mortality model parameters is incorrect.",
                     "\n\nProvide a single vector of length ", length.theta, 
                     " or a matrix of dimensions ", length.cat ," times ", 
                     length.theta, ".(i.e. number of categorical ",
                     "covariates times number of parameters for model ", 
                     model," with ",shape, " shape).", sep=""), call. = FALSE)
        }
      } else {
        if (!is.null(dim(user.par))) {
          par.mat         <- user.par
        } else {
          par.mat         <- matrix(user.par, length.cat, length.theta, 
                                    byrow = TRUE)
        }
        dimnames(par.mat) <- dimnames(low.full.theta)
      }
    }
    return(par.mat)
  }

  # b) Proportional hazards:
  CheckParsPH                <- function(user.par, par.name){
    if (is.null(user.par)) {
      par               <- rep(0, length.cont)
      names(par)        <- colnames(Zcont)
    } else {
      length.par        <- length(user.par)
      if (length.par != length.cont) {
        stop(paste("\nLength of ", par.name, 
                   " for prop. hazards section is not equal to ",
                   "number of covariates (n = ", length.cont, ").", 
                   sep = ""), call. = FALSE)
      } else {
        par             <- user.par
        names(par)      <- colnames(Zcont)
      }
    }
    return(par)
  }

  parallelVars               <- c("CalculateBasicMx", "CalculateBasicSx", 
                                  "CalculateFullFx", "CalculateFullMx", 
                                  "CalculateFullSx", "BuildAliveMatrix", 
                                  "CalculateLowC", "name.theta", 
                                  "length.theta0", "length.theta", 
                                  "studyStart", "studyEnd", "recaptTrans", 
                                  "progrPlots")
	
  # 3.6 Check when all covariates in mortality:
  if (covarsStruct == "all.in.mort") {
    if (!is.null(covariate.type$cont)) {
      if (model != "GO") {
        warning("For effects of all covariate types on mortality parameters ",
                "only simple Gompertz (GO) model can be used. ",
                "Model and shape arguments were changed to 'GO' and 'simple',",
                " respectively.\n", call. = FALSE)
      }
      model                  <- "GO"
      shape                  <- "simple"
    } else {
      warning("No continuous covariates were included in the data. Argument",
              " 'covarsStruct' will be set to 'fused'.\n", call. = FALSE)
      covarsStruct           <- 'fused'
    }
  }
  
  parallelVars               <- c(parallelVars, "model", "shape", "covarsStruct")
  
  # 4 Data formatting:
  # 4.1 Extract raw data and create BaSTA data tables:
  # 4.1.1 Study duration, birth & death matrix and recapture matrix:
  study.years                <- studyStart:studyEnd
  study.length               <- length(study.years)
  n                          <- nrow(object)
  Y                          <- as.matrix(object[, 1:study.length + 3])
  colnames(Y)                <- study.years
  bd                         <- as.matrix(object[, 2:3])
  bi                         <- bd[, 1]
  di                         <- bd[, 2]

  # 4.1.2 Calculate first and last time observed 
  #       and total number of times observed:
  ytemp                      <- t(t(Y) * study.years)
  last.obs                   <- c(apply(ytemp,1,max))
  ytemp[ytemp == 0]          <- 1e4
  first.obs                  <- c(apply(ytemp,1,min))
  first.obs[first.obs == 1e4]<- 0
  oi                         <- Y %*% rep(1, study.length)
  rm("ytemp")

  # 4.1.3 Define study duration:
  Dx                         <- (study.years[2] - study.years[1])
  Tm                         <- matrix(study.years, n, study.length, byrow=TRUE)

  parallelVars               <- c(parallelVars, "niter", "burnin", "thinning", 
                                  "bi", "di", "Dx", "Tm", "last.obs", 
                                  "first.obs", "study.years", "study.length", 
                                  "n", "bd", "Y", "oi") 

  # 4.2 Extract covariates:
  # Find if there are covariates:
  if (ncol(object) > study.length + 3) {
    Z                        <- as.matrix(object[, 
                                          (study.length + 4):ncol(object)])
    covariate.type           <- FindCovariateType(Z)
    if (covarsStruct == "prop.haz") {
      Zcont                  <- Z
      Zcat                   <- matrix(1, n, 1)
      colnames(Zcat)         <- "NoCat"
      Cat                    <- FALSE
      Cont                   <- TRUE
    } else {
      # If all covariates should be included in the mortality section:
      if (covarsStruct == "all.in.mort" & !is.null(covariate.type$cont)) {
        Zcat                 <- Z
        Zac                  <- matrix(Z[, covariate.type$cont], nrow(Z), 
                                       length(covariate.type$cont))
        colnames(Zac)        <- names(covariate.type$cont)
        Zcat[, covariate.type$cont] <- t(t(Zac) - apply(Zac, 2, mean))
        Zcat                 <- Zcat[, c(covariate.type$int, 
                                       covariate.type$cat, 
                                       covariate.type$cont)]
        if (is.null(covariate.type$cat) & is.null(covariate.type$int)) {
          Zcat               <- cbind(1,Zcat)
          colnames(Zcat)     <- c("Intercept", names(covariate.type$cont))
        }
        covariate.type       <- FindCovariateType(Zcat)
        Zcont                <-  matrix(0,n,1)
        colnames(Zcont)      <- "NoCont"
        Cont                 <- FALSE
      } else {

        # Find if there are continuous covariates:
        if (!is.null(covariate.type$cont) & covarsStruct != "all.in.mort") {
          Zcont              <- matrix(Z[, covariate.type$cont], nrow(Z), 
                                       length(covariate.type$cont))
          rangeZc            <- apply(Zcont, 2, range)
          idNoChange         <- which(rangeZc[1, ] * rangeZc[2, ] < 0)
          if (length(idNoChange) > 0) {
            for (i in idNoChange) {
            	 Zcont[, i]      <- Zcont[, i] - mean(Zcont[, i])
            }
          }
#          Zcont              <- t(t(Zcont) - apply(Zcont, 2, mean))
          colnames(Zcont)    <- names(covariate.type$cont)
          Cont               <- TRUE
        } else {
          Zcont              <- matrix(0,n,1)
          colnames(Zcont)    <- "NoCont"
          Cont               <- FALSE
        }
		
        # Find if there are categorical covariates:
        if (!is.null(covariate.type$cat)) {
          Zcat               <- Z[, covariate.type$cat]
          if (!is.null(covariate.type$int)) {
            Zcat             <- cbind(1, Zcat)
            colnames(Zcat)   <- c("Intercept", colnames(Zcat)[-1])
          }
          Cat                <- TRUE
        } else {
          Zcat               <- matrix(1, n, 1); colnames(Zcat) <- "NoCat"
          Cat                <- FALSE
        }
      }
    }
  } else {
    Z                        <- NULL
    Zcat                     <- matrix(1, n, 1)
    colnames(Zcat)           <- "NoCat"
    Zcont                    <- Zcat
    colnames(Zcont)          <- "NoCont"
    Cat                      <- FALSE
    Cont                     <- FALSE
    covariate.type           <- FindCovariateType(NULL)
  }
  length.cat                 <- ncol(Zcat)
  length.cont                <- ncol(Zcont)
  
  parallelVars               <- c(parallelVars, "Zcat", "Zcont", "Cont", 
                                  "length.cat", "length.cont", "covariate.type") 


  # 5. MCMC prep:
  # 5.1 Parameter names and lower bounds for theta parameters:
  name.full.theta            <- paste(rep(name.theta, each = length.cat), 
                                      "[", rep(colnames(Zcat), length.theta), 
                                      "]", sep = "")
  if (length.cat == 1) {
    name.full.theta          <- name.theta	
  }
  length.full.theta          <- length(name.full.theta)
  low.full.theta             <- matrix(low.theta, length.cat, length.theta, 
                                       byrow = TRUE, 
                                       dimnames = list(colnames(Zcat), 
                                       name.theta))
  name.gamma                 <- paste("gamma[", colnames(Zcont), "]", sep="")
  if (length.cont == 1) {
    name.gamma               <- "gamma"
  }
  name.pi                    <- ifelse(length(recaptTrans) == 1, "pi", 
                                       paste("pi[", recaptTrans, "]", sep = ""))
  name.post                  <- c("post[theta,gamma]", "post[X0]", "post[full]")

  parallelVars               <- c(parallelVars, "length.theta", 
                                  "length.full.theta", "name.full.theta", 
                                  "name.gamma", "name.pi", "name.post", 
                                  "low.full.theta") 

  # 5.2 Verify thetaJumps, initial parameters and thetaPriors: 
  # 5.2.1 Survival model:
  # a) Initial parameters:
  theta.g                    <- CheckParsMort(ini.theta, thetaStart, "theta")
  if (covarsStruct == "all.in.mort" & is.null(thetaStart)) {
    theta.g[names(covariate.type$cont), ] <- 0
  }
  thetaStart                 <- theta.g

  # b) Jumps:
  theta.jump                 <- CheckParsMort(jump.theta, thetaJumps, "jumps")
  if (covarsStruct == "all.in.mort" & is.null(thetaJumps)) {
    theta.jump[names(covariate.type$cont), ] <- 0.001
  }

  # c) Priors:
  theta.prior                <- CheckParsMort(prior.theta, thetaPriors, "priors")
  if (covarsStruct == "all.in.mort" & is.null(thetaPriors)) {
    theta.prior[names(covariate.type$cont), ] <- 0
  }
	
	
  # 5.2.2 Proportional hazards section:
  if (Cont) {
    # a) Initial parametes:
    gamma.g                  <- CheckParsPH(gammaStart, "gamma parameters")
    gammaStart               <- gamma.g
		
    # b) Jumps:
    gamma.jump               <- CheckParsPH(gammaJumps, "gamma jumps")
    		
    # c) Priors:
    gamma.prior              <- CheckParsPH(gammaPriors, "gamma priors")
  } else {
    gamma.g                  <- 0
    gammaStart               <- gamma.g
    gamma.jump               <- 0
    gamma.prior              <- 0
  }


  # 5.3 Define final priors:
  # a) Survival parameters:
  Ztheta.p                   <- Zcat %*% theta.prior
  theta.sd                   <- 0.5
		
  # b) Prop. hazards section:
  Zgamma.p                   <- Zcont %*% gamma.prior
  gamma.sd                   <- 1

  # c) Age distribution:
  dxx                        <- 0.001
  xx                         <- seq(0,100,dxx)
  zza                        <- cbind(1, matrix(0, length(xx), length.cat - 1))
  zzc                        <- sum(apply(Zcont, 2, mean) * gamma.prior) 
  Ex                         <- sum(xx * CalculateFullFx(xx, zza %*% theta.prior, 
                                   zzc) * dxx)
  v.x                        <- function(x) {
    CalculateFullSx(x, Ztheta.p, Zgamma.p) / Ex
  }
  rm(list=c("dxx", "xx", "zza", "zzc"))

  # d) Detection probability:
  idpi                       <- findInterval(study.years, recaptTrans)
  names(idpi)                <- study.years
  npi                        <- length(unique(idpi))
  rho1                       <- 0.1
  rho2                       <- 0.1
  
  parallelVars               <- c(parallelVars, "Ztheta.p", "theta.sd", 
                                  "Zgamma.p", "gamma.sd", "Ex", "v.x", 
                                  "idpi", "npi", "rho1", "rho2") 

  # 5.4 Starting values:
  # a) Matrix of survival parameters:
  Ztheta.g                   <- Zcat %*% theta.g
		
  # b) Matrix of prop. hazards parameters:
  Zgamma.g                   <- Zcont %*% gamma.g

  # c) Recapture probability:
  pi.g                       <- rep(0.5, npi)
  Pig                        <- pi.g[idpi]

  # d) Times of birth and death:
  bi0                        <- which(bi==0)
  bg                         <- bi
  bg[bi == 0 & first.obs > 0]<- first.obs[bi == 0 & first.obs > 0] - 1
  bg[bi == 0 & first.obs == 0 & di > 0] <- 
    di[bi == 0 & first.obs == 0 & di > 0] - 1

  di0                        <- which(di==0)
  dg                         <- di
  dg[di == 0 & last.obs > 0] <- last.obs[di == 0 & last.obs > 0] + 1
  dg[di == 0 & last.obs == 0]<- bi[di == 0 & last.obs == 0] + 1
  dg[dg < studyStart]<- studyStart + 1

  xg                         <- dg - bg

  parallelVars               <- c(parallelVars, "theta.g", "theta.jump", 
                                  "theta.prior", "gamma.g", "gamma.jump", 
                                  "gamma.prior", "pi.g", "Pig", "bg", "dg", 
                                  "xg", "bi0", "di0") 

  # 5.5 Full observation matrix:
  Fg                         <- c(apply(cbind(studyStart, bg+1), 1, max))
  Lg                         <- c(apply(cbind(studyEnd, dg-1), 1, min))
  Og                         <- BuildAliveMatrix(Fg, Lg, Tm)
  fii                        <- first.obs
  fii[bi > 0 & bi >= studyStart] <- bi[bi > 0 & bi >= studyStart] + 1
  fii[bi > 0 & bi < studyStart]  <- studyStart
  lii                        <- last.obs
  lii[di > 0 & di <= studyEnd]   <- di[di > 0 & di <= studyEnd] - 1
  lii[di > 0 & di > studyEnd]    <- studyEnd
  lfi                        <- BuildAliveMatrix(fii, lii, Tm)

  parallelVars               <- c(parallelVars, "Fg", "Lg", "Og", "lfi") 


  # 6.  Multiple MCMC function:
  multiMCMC  <- function(sim) {
    if (parallel){
      for(ii in 1:(sim * 2)) {}
    } 
    nlow                     <- low.full.theta
    if (nsim > 1) {
      thetaJitter            <- theta.g * 0 + 0.5
      thetaJitter[theta.jump==0] <- 0
      if (covarsStruct=="all.in.mort") {
        thetaJitter[covariate.type$cont, ] <- 0.05
      }
      theta.n                <- matrix(rtnorm(length.full.theta, theta.g, 
                                       thetaJitter, lower=nlow), length.cat, 
                                       length.theta, dimnames=dimnames(theta.g))
      if (shape!="simple") {
        nlow[, 'c']          <- apply(theta.n, 1, CalculateLowC)
        idc.low              <- which(theta.n[, 'c'] < nlow[, 'c'])
        if (length(idc.low)>0) {
          for(cc in idc.low) {
            theta.n[cc,'c']  <- c(rtnorm(1, theta.g[cc, 'c'], 0.5, 
                                  lower = nlow[cc, 'c']))
          }
        }
      }
      if (Cont) {
        gamma.g              <- rnorm(length.cont, gamma.g, 0.5)      
      }
      theta.g                <- theta.n
    }
		
	
    # Output tables:
    thin.seq                 <- seq(burnin, niter, by = thinning)
    theta.mat                <- matrix(NA,niter,length.full.theta)
    colnames(theta.mat)      <- name.full.theta
    gamma.mat                <- matrix(0, niter, length.cont)
    colnames(gamma.mat)      <- name.gamma
    pi.mat                   <- matrix(NA, niter, npi)
    colnames(pi.mat)         <- name.pi
    bi.mat                   <- matrix(NA,length(thin.seq),n)
    di.mat                   <- bi.mat
    la.vec                   <- rep(NA, niter)
    posterior.mat            <- matrix(NA, niter, 3)
    colnames(posterior.mat)  <- name.post
    theta.mat[1, ]           <- theta.g
    pi.mat[1, ]              <- pi.g
    if (Cont) {
      gamma.mat[1, ]         <- gamma.g
    }
    Ztheta.g                 <- Zcat %*% theta.g
    Zgamma.g                 <- Zcont %*% gamma.g
    lag                      <- 0.01
    
    # Juvenile and adult ages:
    IminAge                  <- ifelse(minAge > 0, 1, 0)
    Iag                      <- rep(0, n)
    Ijg                      <- Iag
    Iag[xg > minAge]         <- 1
    Ijg[xg <= minAge]        <- 1
    xjg                      <- xg
    xjg[xg > minAge]         <- minAge
    xag                      <- xg - minAge
    xag[xg < minAge]         <- 0
    xjtg                     <- xg * 0
    idtr                     <- which(bg < studyStart & 
                                      studyStart - bg < minAge)
    xjtg[idtr]               <- studyStart - bg[idtr]
    xatg                     <- xg * 0
    idtr                     <- which(bg + minAge < studyStart)
    xatg[idtr]               <- studyStart - (bg[idtr] + minAge)

    # Run Gibbs sampler:
    naflag                   <- FALSE
    g                        <- 2
    gg                       <- 1
    if (progrPlots) {
      if (.Platform$OS.type=="unix") {
        devtype              <- quartz
      } else {
        devtype              <- windows
      }
      devtype(width = 2, height = 0.5)
      progrpl                <- dev.cur()
      par(mar = rep(0,4))
    }
    while(g <= niter & !naflag) {
		
      # 1.- SAMPLING:
      # a) Sample survival parameters:
      # i) Metropolis draw for params:
      theta.n                <- matrix(rtnorm(length.full.theta, theta.g, 
                                       theta.jump, lower = low.full.theta), 
                                       length.cat, length.theta, 
                                       dimnames = dimnames(theta.g))
      if (shape!="simple") {
        nlow[, 'c']          <- apply(theta.n, 1, CalculateLowC)
        idc.low              <- which(theta.n[, 'c'] < nlow[, 'c'])
        if (length(idc.low)>0) {
          for(cc in idc.low) {
            theta.n[cc,'c']  <- c(rtnorm(1, theta.g[cc,'c'], 0.5, 
                                  lower=nlow[cc,'c']))
          }
        }
      }
      
      if (Cont) {
        gamma.n              <- rnorm(length.cont, gamma.g, gamma.jump) 
      } else {
        gamma.n              <- gamma.g
      }
      
      # ii) Build individual parameter matrices:
      Ztheta.n               <- Zcat %*% theta.n
      Zgamma.n               <- Zcont %*% gamma.n
            
      # iii) Calculate conditional posteriors:
      # - Current parameters:
      p.thg                  <- (log(CalculateFullFx(xag + 
                                0.5 * Dx, Ztheta.g, Zgamma.g)) -
                                log(CalculateFullSx(xatg + 0.5 * Dx, 
                                Ztheta.g, Zgamma.g))) * Iag 
      # - Correct for precision limit:
      p.thg[p.thg==-Inf]     <- -1e200
      # - Priors:
      p.thg                  <- sum(p.thg) + sum(dtnorm(c(theta.g), 
                                c(theta.prior), theta.sd, 
                                lower = low.full.theta, log = TRUE)) + 
                                sum(dnorm(gamma.g, gamma.prior, gamma.sd, 
                                log = TRUE))

      # - Proposed parameters:
      p.thn                  <- (log(CalculateFullFx(xag + 
                                0.5 * Dx, Ztheta.n, Zgamma.n)) -
                                log(CalculateFullSx(xatg + 0.5 * Dx, 
                                Ztheta.n, Zgamma.n))) * Iag
      # - Correct for precision limit:
      p.thn[p.thn==-Inf]     <- -1e200
      # - Priors:
      p.thn                  <- sum(p.thn) + sum(dtnorm(c(theta.n), 
                                c(theta.prior), theta.sd, 
                                lower = low.full.theta, log=TRUE)) + 
                                sum(dnorm(gamma.n, gamma.prior, gamma.sd,
                                log=TRUE))

      r                      <- exp(p.thn-p.thg)
      z                      <- runif (1,0,1)
      
      if (is.na(r) & g==1) {
        if (progrPlots) dev.off(progrpl)
        naflag               <- TRUE 
      } else if (is.na(r) & g > 1) {
      	 if (progrPlots) dev.off(progrpl)
      	 naflag               <- TRUE 
      } else {
      	 if (r>z) {
      	   theta.g            <- theta.n
      	   Ztheta.g           <- Ztheta.n
      	   p.thg              <- p.thn
      	   gamma.g            <- gamma.n
      	   Zgamma.g           <- Zgamma.n
      	 }
      }
      
      # b) Sample times of birth and death:
      # i) New times of birth and death:
      bn                     <- bg 
      bn[bi0]                <- bg[bi0] + sample(-1:1, length(bi0), 
                                replace = TRUE) 
      bn[bi0][oi[bi0] > 0]   <- apply(cbind(bn[bi0][oi[bi0] > 0],
                                      first.obs[bi0][oi[bi0] > 0] - 1), 1, min)
      bn[bi0][oi[bi0] == 0]  <- apply(cbind(bn[bi0][oi[bi0] == 0],
                                      dg[bi0][oi[bi0] == 0] - 1), 1, min)
      dn                     <- dg 
      dn[di0]                <- dg[di0] + sample(-1:1, length(di0), 
                                replace = TRUE) 
      dn[di0]                <- apply(cbind(dn[di0],bn[di0],last.obs[di0] + 1), 
                                      1, max) 
      xn                     <- dn - bn
      
      # ii) New full alive matrices:
      Fn                     <- c(apply(cbind(studyStart, bn + 1), 1, max))
      Ln                     <- c(apply(cbind(studyEnd, dn - 1), 1, min))
      On                     <- BuildAliveMatrix(Fn, Ln, Tm)

      # iii) New juvenile and adult ages:
      Ian                    <- rep(0, n)
      Ijn                    <- Ian
      Ian[xn > minAge]       <- 1
      Ijn[xn <= minAge]      <- 1
      xjn                    <- xn
      xjn[xn > minAge]       <- minAge
      xan                    <- xn - minAge
      xan[xn < minAge]       <- 0
      xjtn                   <- xn * 0
      idtr                   <- which(bn < studyStart & 
                                      studyStart - bn < minAge)
      xjtn[idtr]             <- studyStart - bn[idtr]
      xatn                   <- xn * 0
      idtr                   <- which(bn + minAge < studyStart)
      xatn[idtr]             <- studyStart - (bn[idtr] + minAge)
      
      # iiii) Calculate conditional posteriors:
      # - Current ages:
      p.bdg                  <- (log(lag) * Ijg - lag * xjg) * IminAge + 
                                log(CalculateFullFx(xag + 
                                0.5 * Dx, Ztheta.g, Zgamma.g)) * Iag
      p.bdg[p.bdg==-Inf]     <- -1e200
      p.bdg                  <- p.bdg + (Og - lfi) %*% log(1 - Pig) + 
                                log(v.x(xag + 0.5 * Dx)) * Iag
      # - New ages:
      p.bdn                  <- (log(lag) * Ijn - lag * xjn) * IminAge  + 
                                log(CalculateFullFx(xan + 
                                0.5 * Dx, Ztheta.g, Zgamma.g)) * Ian
      p.bdn[p.bdn==-Inf]     <- -1e200
      p.bdn                  <- p.bdn + (On - lfi) %*% log(1 - Pig) + 
                                log(v.x(xan + 0.5 * Dx)) * Ian
      
      r                      <- exp(p.bdn-p.bdg)
      if (length(which(is.na(r)))>0) {
        if (progrPlots) dev.off(progrpl) 
        naflag               <- TRUE 
      } else {
        z                    <- runif (n, 0, 1)
        idrz                 <- which(r > z)
        bg[idrz]             <- bn[idrz]
        dg[idrz]             <- dn[idrz]
        xg[idrz]             <- xn[idrz]
        p.bdg[idrz]          <- p.bdn[idrz]
        Og[idrz, ]           <- On[idrz, ]
        xjg[idrz]            <- xjn[idrz]
        xjtg[idrz]           <- xjtn[idrz]
        xag[idrz]            <- xan[idrz]
        xatg[idrz]           <- xatn[idrz]
        Iag[idrz]            <- Ian[idrz]
        Ijg[idrz]            <- Ijn[idrz]
      }
      
      # c) Sample recapture probability(ies):
      rho1g                  <- rho1 + t(t(Y) %*% rep(1, n))
      rho2g                  <- rho2 + t(t(Og - Y) %*% rep(1, n))
      Rho1                   <- tapply(rho1g, idpi, sum)
      Rho2                   <- tapply(rho2g, idpi, sum)
      pi.g                   <- rbeta(npi, Rho1, Rho2)
      if (1 %in% pi.g) {
        pi.g[pi.g==1]        <- 1-1e-5
        warning("Some recapture probabilities are equal to 1.",
                "\nThey have been constraint to be fractionally less than 1 ",
                "for computational reasons\n", call. = FALSE)
      }
      Pig                    <- pi.g[idpi]
      
      # d) if minAge > 0, sample lambda:
      if (minAge > 0) {
        lan                   <- rtnorm(n = 1, mean = lag, 
                                       sd = 0.001, lower = 0)
        p.lag                <- sum(log(lag) * Ijg - lag * xjg + 
                                    lag * xjtg) +
                                dtnorm(lag, mean = 0.01, sd = 1, lower = 0)
        p.lan                <- sum(log(lan) * Ijg - lan * xjg + 
                                    lan * xjtg) +
                                dtnorm(lan, mean = 0.01, sd = 1, lower = 0)
        
        r                    <- exp(p.lan - p.lag)
        z                    <- runif(1, 0, 1)
        if (r > z) {
          lag                 <- lan
        }
      }
      
      # 2.- STORE RESULTS:
      # Parameters and latent states:
      theta.mat[g, ]         <- theta.g
      pi.mat[g, ]            <- pi.g
      if (Cont) {
      	 gamma.mat[g, ]       <- gamma.g
      }
      if (minAge > 0) {
        la.vec[g]            <- lag
      }
      if (g %in% thin.seq) {
      	 bi.mat[gg, ]         <- bg
      	 di.mat[gg, ]         <- dg
      	 gg                   <- gg + 1
      }
      
      # Conditional posteriors: 
      posterior.mat[g, ]     <- c(p.thg, sum(p.bdg), p.thg + 
                                 sum((Og - lfi) %*% log(1 - Pig)))
  
      # Progress plot:
      if (g %in% round(seq(1, niter, length = 100)) & progrPlots) {
      	 par(mar        = rep(0, 4))
      	 plot(x         = c(0, niter * 1.1), 
      	      y         = c(0, 1), 
      	      axes      = FALSE, 
      	      col       = NA, 
      	      xlab      = "", 
      	      ylab      = "")
      	 polygon(x      = c(0, niter, niter, 0), 
      	         y      = c(0.35, 0.35, 0.65, 0.65), 
      	         col    = NA, 
      	         border = 'dark red')
      	 polygon(x      = c(0, g, g, 0), 
      	         y      = c(0.35, 0.35, 0.65, 0.65), 
      	         col    = 'dark red', 
      	         border = 'dark red')
      	 text(x         = niter / 2, 
      	      y         = 0.85, 
      	      labels    = paste("MCMC progress (Sim. ", sim, ")", sep = ""), 
      	      cex       = 0.9)
      	 text(x         = g, 
      	      y         = 0.15, 
      	      labels    = paste(round(g / niter * 100), "%", sep = ""), 
      	      cex       = 0.8)
      }
      g                      <- g + 1
    }
    if (g == niter + 1) {
      g                      <- niter
    }
    if (progrPlots) {
      dev.off(progrpl)
    }
    # Return results:
    return(list(theta  = theta.mat, 
                gamma  = gamma.mat, 
                pi     = pi.mat,
                la     = la.vec, 
                bi     = bi.mat, 
                di     = di.mat, 
                post   = posterior.mat, 
                g      = g, 
                naflag = naflag))
  }
  
  # 7. Run (multi) MCMC:
  # 7.1 Run simulations either in series or in parallel:
  if (nsim == 1) {
    parallel <- FALSE
  }
  if (nsim > 1) {
    cat("Multiple simulations started...\n\n") 
  } else {
    cat("Simulation started...\n\n")
  }
  Start                      <- Sys.time()
  if (parallel) {
    avail.pkgs               <- available.packages()
    if (!is.element("snowfall", avail.pkgs)) {
      warning("\nPackage 'snowfall' is not installed.\nSimulations ",
              "will not be ran in parallel (computing time will ",
              "be longer...)\n")
      basta.out              <- lapply(1:nsim, multiMCMC)
    } else {
      require(snowfall)
      sfInit(parallel = TRUE, cpus = ncpus);
      sfExport(list = c(parallelVars, "parallel", "nsim", "minAge"))
      sfLibrary(msm)
      basta.out              <- sfClusterApplyLB(1:nsim, multiMCMC)
      sfStop()
    }
  } else {
    basta.out                <- lapply(1:nsim, multiMCMC)
  }
  End                        <- Sys.time()
  
  # 7.2 Report if all simulations ran through:
  simNames                   <- paste("Sim.", (1:nsim), sep="")
  full.runs                  <- rep(0,nsim)
  names(full.runs)           <- simNames
  last.steps                 <- full.runs
  for(i in 1:nsim) {
    last.steps[i]            <- basta.out[[i]]$g	
    full.runs[i]             <- ifelse(last.steps[i] == niter, 1, 0)
  } 
  id.failed                  <- which(full.runs == 0)
  all.ran                    <- FALSE
  if (nsim==1) {
    if (full.runs==1) {
      cat("MCMC finished running\n")
      cat(paste("Total MCMC computing time: ", round(as.numeric(julian(End) - 
          julian(Start)) * 24 * 60, 2), " minutes\n\n", sep=""))
      all.ran                <- TRUE
    } else {
      cat(paste("MCMC stopped at step ", basta.out[[1]]$g,
          "\nPdf of ages at death equal to 0 for some individuals.",
          "\nReduce jumps to avoid pdf of ages at death equal to 0.\n", 
          sep = ""))
    }
  } else {
    if (length(id.failed)>0 & length(id.failed)<nsim) {
      cat("\nOne or more simulations failed\nConvergence diagnostics ",
          "and model selection will not be calculated.\n",
          "Reduce jumps to avoid pdf of ages at death equal to 0.\n")
    } else if (length(id.failed)==nsim) {
      cat("\nAll simulations failed\nConvergence diagnostics and model ",
          "selection will not be calculated.\n",
          "Reduce jumps to avoid pdf of ages at death equal to 0.\n")
    } else {
      all.ran                <- TRUE
      cat("\nMultiple simulations finished.\n")
      cat(paste("Total MCMC computing time: ", 
          round(as.numeric(julian(End) - julian(Start)) * 24 * 60, 2), 
          " minutes\n\n", sep = ""))
    }
  }	
  
  # 8. Diagnostics:
  thin.seq                   <- seq(burnin, niter, thinning)
  nthin                      <- length(thin.seq)
  
  # 8.1 Thinned result matrices:
  if (Cont) {
    name.out.mat             <- c(name.full.theta, name.gamma, name.pi, 
                                  name.post) 
  } else {
    name.out.mat             <- c(name.full.theta, name.pi, name.post)
  }
  out.mat                    <- matrix(NA, niter * nsim, length(name.out.mat))
  dimnames(out.mat)          <- list(rep(simNames, each = niter), name.out.mat)
  Bimat                      <- matrix(NA, nthin * nsim, n)
  rownames(Bimat)            <- rep(simNames, each = nthin)
  Dimat                      <- Bimat
  idthin                     <- rep(0, niter*nsim)
  if (minAge > 0) {
    la.mat                   <- matrix(NA, nrow = niter, ncol = nsim, 
                                       dimnames = list(NULL, simNames))
  } else {
    la.mat                   <- matrix(NA, nrow = 1, ncol = nsim,
                                       dimnames = list(NULL, simNames))
  }
  
  for(i in 1:nsim) {
    Idsim                    <- which(rownames(out.mat) == simNames[i])
    if (Cont) {
      out.mat[Idsim, ]       <- cbind(basta.out[[i]]$theta, 
                                      basta.out[[i]]$gamma, 
                                      basta.out[[i]]$pi, 
                                      basta.out[[i]]$post)
    } else {
      out.mat[Idsim, ]       <- cbind(basta.out[[i]]$theta, 
                                      basta.out[[i]]$pi, 
                                      basta.out[[i]]$post)
    }
    idthin[Idsim[thin.seq]]  <- 1
    Idsim                    <- which(rownames(Bimat) == simNames[i])
    Bimat[Idsim, ]           <- basta.out[[i]]$bi
    Dimat[Idsim, ]           <- basta.out[[i]]$di
    if (minAge > 0){
      la.mat[, i]            <- basta.out[[i]]$la
    } 
  }
  
  # 8.2 Basic summary statistics for parameters:
  par.mat                    <- out.mat[idthin == 1, 
                                        -(c(ncol(out.mat) - c(2:0)))]
  coef                       <- cbind(apply(par.mat, 2, mean, na.rm=TRUE), 
                                      apply(par.mat, 2, sd, na.rm=TRUE), 
                                      t(apply(par.mat, 2, quantile, 
                                      c(0.025, 0.975), na.rm = TRUE)), 
                                      NA, NA, NA)
  colnames(coef)             <- c("Estimate", "StdErr", "Lower95%CI", 
                                  "Upper95%CI", "SerAutocor", "UpdateRate", 
                                  "PotScaleReduc")
  if (length(id.failed) < nsim) {
    idfix                    <- which(theta.jump==0)
    if (length(idfix) > 0) {
      coef[-idfix,"SerAutocor"]<- apply(par.mat[, -c(idfix)], 2, 
                                        function(x) cor(x[-1], x[-length(x)], 
                                        use = "complete.obs"))
      coef[idfix,"SerAutocor"] <- 1
    } else {
      coef[, "SerAutocor"]   <- apply(par.mat, 2, function(x) cor(x[-1], 
                                     x[-length(x)], use = "complete.obs"))
    }
    coef[, "UpdateRate"]     <- apply(out.mat[, -c(ncol(out.mat) - c(2:0))], 
                                      2, function(x) 
                                      length(which(diff(x[!is.na(x)]) != 0)) / 
                                      length(x[!is.na(x)]))
  }
  out.mat                    <- cbind(idthin, out.mat)
  
  # 8.3 Convergence and model selection:
  if (all.ran) {
    if (nsim>1) {
      # 8.3.1 Convergence diagnostics (potential scale reduction):
      Means                  <- apply(par.mat, 2, function(x) 
                                      tapply(x, rownames(par.mat), mean))
      Vars                   <- apply(par.mat, 2, function(x) 
                                      tapply(x, rownames(par.mat), var))
      meanall                <- apply(Means, 2, mean)
      B                      <- nthin / (nsim - 1) * 
                                apply(t((t(Means) - meanall)^2), 2, sum)
      W                      <- 1 / nsim * apply(Vars, 2, sum)
      Varpl                  <- (nthin - 1) / nthin * W + 1 / nthin * B
      Rhat                   <- sqrt(Varpl / W)
      Rhat[Varpl==0]         <- 1
      conv                   <- cbind(B, W, Varpl, Rhat)
      rownames(conv)         <- colnames(par.mat)
      coef[, ncol(coef)]     <- conv[, 'Rhat']
      
      # Report if convergence was reached:
      idnconv                <- which(conv[, 'Rhat']< 0.95 | conv[, 'Rhat']>1.1)
      if (length(idnconv) > 0) {
        modSel               <- NULL
        kl.list              <- NULL
        warning("Convergence not reached for some survival parameters.",
                "\nDIC could not be calculated.\n", call. = FALSE)
      } else {
        # 8.3.2 Model selection (DIC, if convergence was reached):
        posterior            <- out.mat[idthin == 1, ncol(out.mat)]
        L                    <- length(posterior)
        Dm                   <- -2 * posterior
        Dmode                <- -2 * posterior[which(posterior == 
                                max(posterior))[1]]
        Dave                 <- mean(Dm)
        pD                   <- Dave - Dmode
        k                    <- npi + length.full.theta
        if (Cont) {
          k                  <- k + length.cont
        }
        DIC                  <- 2 * Dave - Dmode
        modSel               <- c(Dave, Dmode, pD, k, DIC)
        names(modSel)        <- c("D.ave", "D.mode", "pD", "k", "DIC")
        cat("Survival parameters converged appropriately.",
            "\nDIC was calculated\n.")
        # 8.3.3 Inference on parameter estimates:
        # Kullback-Liebler distances for categorical covariates:
        if (is.null(covariate.type$cat)) {
          kl.list            <- NULL
        } else {
          name.cat           <- names(covariate.type$cat)
          n.cat              <- length(name.cat)
          n.comb             <- (n.cat - 1)^2 - 
                                ((n.cat - 1)^2 - (n.cat - 1)) / 2
          covar.comb         <- matrix(0, n.comb, 2, 
                                       dimnames = list(NULL, c("cov1", "cov2")))
          ij                 <- 1
          for (i in 1:(n.cat - 1)) {
            for (j in 2:n.cat) {
            	if (i < j) {
            	  covar.comb[ij, ]<- c(name.cat[i], name.cat[j])
            	  ij           <- ij + 1
            	}
            }
          }
          if (covarsStruct == "fused") {
            kl.12            <- matrix(NA, nrow = n.comb, ncol = length.theta0, 
                                       dimnames = list(paste(covar.comb[, 1], 
                                       "-", covar.comb[, 2], sep = ""), name.theta0))
            kl.21            <- matrix(NA, nrow = n.comb, ncol = length.theta0, 
                                       dimnames = list(paste(covar.comb[, 2], 
                                       "-", covar.comb[, 1], sep = ""), name.theta0))
            p.low            <- low.full.theta[1,]
          } else {
            kl.12            <- matrix(NA, nrow = n.comb, ncol = 1, 
                                       dimnames = list(paste(covar.comb[,1], 
                                       "-", covar.comb[,2], sep = ""), "gamma"))
            kl.21            <- matrix(NA, nrow = n.comb, ncol = 1, 
                                       dimnames = list(paste(covar.comb[, 2], 
                                       "-", covar.comb[, 1], sep = ""), "gamma"))
            p.low            <- -Inf
          }
          q.12               <- kl.12
          q.21               <- kl.21
          kl.pars            <- colnames(kl.12)
          for(i in 1:ncol(kl.12)){
            for(j in 1:n.comb){
              p1             <- par.mat[, paste(kl.pars[i], "[", 
                                        covar.comb[j, 1], "]", sep = "")]
              mean.p1        <- mean(p1)
              sd.p1          <- sd(p1)
              p2             <- par.mat[, paste(kl.pars[i], "[", 
                                        covar.comb[j, 2], "]", sep = "")]
              mean.p2        <- mean(p2)
              sd.p2          <- sd(p2)
              full.range     <- range(c(mean.p1 + c(-4, 4) * sd.p1, 
                                        mean.p2 + c(-4, 4) * sd.p2))
              full.range[1]  <- max(full.range[1], p.low[i])
              p.vec          <- seq(full.range[1], full.range[2], length = 500)
              dp             <- p.vec[2] - p.vec[1]
              dens.p1        <- dtnorm(p.vec, mean = mean(p1), sd = sd(p1), 
                                       low = p.low[i])  
              dens.p2        <- dtnorm(p.vec, mean = mean(p2), sd = sd(p2), 
                                       low = p.low[i])  
              kl.12[j, i]    <- sum(dens.p1*log(dens.p1/dens.p2) * 
                                     dp)
              kl.21[j, i]    <- sum(dens.p2*log(dens.p2/dens.p1) * 
                                     dp)
              q.12[j, i]     <- (1 + (1 - exp(-2 * kl.12[j, i])^(1 / 2))) / 2
              q.21[j, i]     <- (1 + (1 - exp(-2 * kl.21[j, i])^(1 / 2))) / 2
            }
          }
          kl.list            <- list(kl12 = kl.12, 
                                     kl21 = kl.21, 
                                     q12  = q.12, 
                                     q21  = q.21)
        }
      }
    } else {
      conv                   <- NULL
      modSel                 <- NULL
      kl.list                <- NULL
    }
    
    # 8.3.3 Summary times of birth and ages at death:
    xq                       <- apply(Dimat - Bimat, 2, quantile, 
                                      c(0.5, 0.025, 0.975))
    bq                       <- apply(Bimat, 2, quantile, 
                                      c(0.5, 0.025, 0.975))
    
    # 8.3.4 Summary Survival and mortality:
    thmat                    <- par.mat[, name.full.theta]
    if (Cont) {
      rzc                    <- apply(Zcont, 2, quantile, c(0.5, 0.025, 0.975))
      rownames(rzc)          <- c("Med.", "Lower", "Upper")
      gave                   <- apply(as.matrix(par.mat[, 
                                      which(substr(colnames(par.mat), 
                                      1, 2) == "gamma")]), 2, mean)
    } else if (covarsStruct == "all.in.mort") {
      rzc                    <- apply(as.matrix(Zcat[, covariate.type$cont]), 
                                      2, quantile, c(0.5, 0.025, 0.975))
      dimnames(rzc)          <- list(c("Med.", "Lower", "Upper"), 
                                     names(covariate.type$cont))
      gave                   <- 0
      zcname                 <- colnames(rzc)
    } else {
      rzc                    <- matrix(0, 1, 1, dimnames = list("nc", "nc"))
      gave                   <- 0
      zcname                 <- c("")
    }
    Sxq                      <- list()
    mxq                      <- list()
    xvec                     <- list()
    zaname                   <- c(names(covariate.type$int), 
                                  names(covariate.type$cat))
    if (is.null(zaname)) {
      zaname                 <- "NoCov"
    } 
    for(i in 1:length(zaname)) {
      if (zaname[i] == "NoCov") {
        idza                 <- 1:n
      } else {
        idza                 <- which(Z[, zaname[i]] == 1)
      }
      xv                     <- seq(0, ceiling(max(xq[1, idza]) * 1.1), 0.1)
      xvec[[zaname[i]]]      <- xv + minAge
      for(j in 1:ncol(rzc)) {
        Sxq[[zaname[i]]][[colnames(rzc)[j]]] <- 
                array(0, dim = c(length(xv), 3, nrow(rzc)), 
                dimnames = list(NULL, c("50%", "2.5%", "97.5%"), 
                rownames(rzc)))
        
        mxq[[zaname[i]]][[colnames(rzc)[j]]] <- 
          Sxq[[zaname[i]]][[colnames(rzc)[j]]]
        for(k in 1:nrow(rzc)) {
          gaa                <- sum(gave * rzc[k, j])
          Cols               <- paste(name.theta, "[",zaname[i], "]", sep = "")
          if (length.cat == 1) Cols <- name.theta
          Thm                <- thmat[, Cols]
          if (covarsStruct == "all.in.mort") {
            Thm              <- Thm + thmat[, paste(name.theta, 
                                "[",names(covariate.type$cont)[j],"]", 
                                sep = "")] * rzc[k, j]
          }
          Sxq[[zaname[i]]][[colnames(rzc)[j]]][, , k] <- 
                t(apply(apply(Thm, 1 ,CalculateMultiSx), 1, quantile, 
                c(0.5, 0.025, 0.975)))
          mxq[[zaname[i]]][[colnames(rzc)[j]]][, , k] <- 
                t(apply(apply(Thm, 1 ,CalculateMultiMx), 1, quantile, 
                c(0.5, 0.025, 0.975)))
        }
      }
    }
    
    # 8.3.5 Calculate life table from estimated ages at death:
    if (lifeTable) {
      LT                     <- list()
      for(i in 1:length(zaname)) {
        if (zaname[i] == "NoCov") {
          idza               <- which(bq[1, ] >= studyStart)
        } else {
          idza               <- which(Z[, zaname[i]] == 1 &  
                                      bq[1, ] >= studyStart)
        }
        x                    <- xq[1, idza]
        LT[[zaname[i]]]      <- MakeLifeTable(x, ax = 0.5, n = 1)
      }
    } else {
      LT                     <- NULL
    }
  } else {
    conv                     <- NULL
    modSel                   <- NULL
    kl.list                  <- NULL
    xq                       <- NULL
    bq                       <- NULL
    Sxq                      <- NULL
    mxq                      <- NULL
    xvec                     <- NULL
    if (lifeTable) {
      LT                     <- NULL
    }
  }
  
  # 9. Return a list object of class 'basta':
  Settings                   <- c(niter, burnin, thinning, nsim)
  names(Settings)            <- c("niter", "burnin", "thinning", "nsim") 
  ModelSpecs                 <- c(model, shape, covarsStruct, 
                                  paste(names(covariate.type$cat), 
                                  collapse = ", "), 
                                  paste(names(covariate.type$cont), 
                                  collapse = ", "))
  names(ModelSpecs)          <- c("model", "shape", "Covar. structure", 
                                  "Categorical", "Continuous")
  JumpPriors                 <- cbind(c(theta.jump, gamma.jump), 
                                      c(theta.prior, gamma.prior))
  dimnames(JumpPriors)       <- list(c(name.full.theta, name.gamma), 
                                     c("Jump.sd", "Mean.priors"))
  if (!Cont) {
    JumpPriors               <- JumpPriors[-nrow(JumpPriors), ]
  }
  output                     <- list()
  output$coefficients        <- coef
  output$modSel              <- modSel
  output$Convergence         <- conv
  output$KullbackLiebler     <- kl.list
  output$settings            <- Settings
  output$ModelSpecs          <- ModelSpecs
  output$JumpPriors          <- JumpPriors
  output$Params              <- out.mat
  output$Bis                 <- Bimat
  output$Dis                 <- Dimat
  output$Bq                  <- bq
  output$Xq                  <- xq
  output$Sx                  <- Sxq
  output$mx                  <- mxq
  output$xv                  <- xvec
  output$bd                  <- bd
  output$Y                   <- Y
  output$Zcat                <- Zcat
  output$Zcont               <- Zcont
  output$studyStart          <- studyStart
  output$studyEnd            <- studyEnd
  output$finished            <- full.runs
  if (lifeTable) {
    output$lifeTable         <- LT
  }
  if (minAge > 0) {
    output$Lambda            <- la.mat
  }
  class(output)              <- "basta"
  return(output)
}
