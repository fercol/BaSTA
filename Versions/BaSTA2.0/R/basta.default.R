#  File src/library/BaSTA/R/basta.R
#  Part of the R package, http://www.R-project.org
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#  This function estimates age-specific mortality from capture-recapture/
#  recovery (CRR) data when a large proportion of (or all) the records have
#  unknown times of birth and death. It uses the framework described by
#  Colchero & Clark (2011) Journal of Animal Ecology. 
basta.default <- 
		function(object, form = NULL, minAge = 0, model = "GO", 
        shape = "simple", covarsStruct = "prop.haz", 
        cohortCov = NULL, defaultPars = TRUE, niter = 50000, 
        burnin = 5001, thinning = 50, nsim = 1, parallel = FALSE, 
        ncpus = 2, lifeTable = TRUE, progrPlots = FALSE, ...) {
	
	# 1. Load packages:
	require(msm)
	
	# 2. Run error checks:
	errors <- .FindErrors(model, shape, covarsStruct, niter, burnin, thinning)
	
	# 3. Default values:
  mortDat <- .PrepMortDat(object)
  covList <- .CreateCovList(object, form = form)
  
  # Initial values for times of birth and death:
  timesg <- .SetIniTimeDat(covList, mortDat, minAge)
  parsMort <- .DefineParams(covList, model, shape, covarsStruct, default)
  parsg <- .DefineParsG(parsMort)
  parCovsg <- .CalcParams(parsg, covList)
  
  # Set survival and mortality functions:
  CalcSimpleSurv <- .DefineSimpleSurv(model)
  CalcSurv <- .DefineSurv(shape, CalcSimpleSurv)
  CalcSimpleMort <- .DefineSimpleMort(model)
  CalcMort <- .DefineMort(shape, CalcSimpleMort)

  # Define prior for age distribution:
  priorAgeDistr <- .SetPriorAgeDist(parsMort)
  
	# d) Detection probability:
	idpi                       <- findInterval(bastaDat$years, recaptTrans)
	names(idpi)                <- bastaDat$years
	npi                        <- length(unique(idpi))
	rho1                       <- 0.1
	rho2                       <- 0.1
	
  # c) Recapture probability:
  pi.g <- rep(0.5, npi)
  Pig <- pi.g[idpi]
  
  # d) Times of birth and death:
  iniTimes <- .SetIniTimes(mortDat, minAge)
  
  parallelVars <- c("mortDat", "covList", "parsMort", "parsg", 
      "parCovsg", "CalcSurv", "CalcMort", "CalcSimpleMort", 
      "CalcSimpleSurv", "priorAgeDistr", "iniTimes") 
  
  # 5.5 Full observation matrix:
  Fg <- c(apply(cbind(studyStart, bg+1), 1, max))
  Lg <- c(apply(cbind(studyEnd, dg-1), 1, min))
  Og <- BuildAliveMatrix(Fg, Lg, Tm)
  fii <- first.obs
  fii[bi > 0 & bi >= studyStart] <- bi[bi > 0 & bi >= studyStart] + 1
  fii[bi > 0 & bi < studyStart]  <- studyStart
  lii <- last.obs
  lii[di > 0 & di <= studyEnd] <- di[di > 0 & di <= studyEnd] - 1
  lii[di > 0 & di > studyEnd] <- studyEnd
  lfi <- BuildAliveMatrix(fii, lii, Tm)
  
  parallelVars <- c(parallelVars, "Fg", "Lg", "Og", "lfi") 
  
  # 6.  Multiple MCMC function:
  multiMCMC  <- function(sim) {
    if (parallel){
      for(ii in 1:(sim * 2)) {}
    }     
    
    # Output tables:
    thin.seq <- seq(burnin, niter, by = thinning)
    theta.mat <- matrix(NA,niter,length.full.theta)
    colnames(theta.mat) <- name.full.theta
    gamma.mat <- matrix(0, niter, length.cont)
    colnames(gamma.mat) <- name.gamma
    pi.mat <- matrix(NA, niter, npi)
    colnames(pi.mat) <- name.pi
    bi.mat <- matrix(NA,length(thin.seq),n)
    di.mat <- bi.mat
    la.vec <- rep(NA, niter)
    posterior.mat <- matrix(NA, niter, 3)
    colnames(posterior.mat) <- name.post
    theta.mat[1, ] <- theta.g
    pi.mat[1, ] <- pi.g
    if (Cont) {
      gamma.mat[1, ] <- gamma.g
    }
    updVec <- rep(0, niter)
    updateLen <- 500
    updateInt <- 1000 + 0:(nPar * 15 - 1) * updateLen
    
    # Juvenile and adult ages:
    IminAge <- ifelse(minAge > 0, 1, 0)
    Iag <- rep(0, n)
    Ijg <- Iag
    Iag[xg >= minAge] <- 1
    Ijg[xg < minAge] <- 1
    xjg <- xg
    xjg[xg > minAge] <- minAge
    xag <- xg - minAge
    xag[xg < minAge] <- 0
    xjtg <- xg * 0
    idtr <- which(bg < studyStart & 
            studyStart - bg < minAge)
    xjtg[idtr] <- studyStart - bg[idtr]
    xatg <- xg * 0
    idtr <- which(bg + minAge < studyStart)
    xatg[idtr] <- studyStart - (bg[idtr] + minAge)
    
    # Start parameters from different values:
    nlow <- low.full.theta
    theta.n <- theta.g
    if (Cont) gamma.n <- gamma.g
    if (nsim > 1) {
      
      thetaJitter <- theta.g * 0 + 0.5
      thetaJitter[theta.jump == 0] <- 0
      InfPost <- TRUE
      while (InfPost) {
        theta.n[idUpdJump] <- rtnorm(length(idUpdJump), theta.g[idUpdJump], 
            thetaJitter[idUpdJump], lower=nlow[idUpdJump])
        if (shape != "simple") {
          nlow[, 'c'] <- apply(theta.n, 1, CalculateLowC)
          idc.low <- which(theta.n[, 'c'] < nlow[, 'c'])
          if (length(idc.low) > 0) {
            for(cc in idc.low) {
              theta.n[cc,'c'] <- c(rtnorm(1, theta.g[cc, 'c'], thetaJitter[cc, 'c'], 
                      lower = nlow[cc, 'c']))
            }
          }
        }
        if (Cont) {
          gamma.n <- rnorm(length.cont, gamma.g, 0.5)      
        } else {
          gamma.n <- gamma.g
        }
        Ztheta.n <- Zcat %*% theta.n
        Zgamma.n <- Zcont %*% gamma.n
        
        post <- sum((log(CalculateFullFx(xag, 
                          Ztheta.n, Zgamma.n)) -
                  log(CalculateFullSx(xatg, 
                          Ztheta.n, Zgamma.n))) * Iag)
        if(post == -Inf) {
          InfPost <- TRUE
        } else {
          InfPost <- FALSE
        }
        
      }
    }
    theta.g <- theta.n
    if(Cont) gamma.g <- gamma.n
    Ztheta.g <- Zcat %*% theta.g
    Zgamma.g <- Zcont %*% gamma.g
    lag <- 0.01
    
    # Run Gibbs sampler:
    naflag <- FALSE
    gg <- 1
    if (progrPlots) {
      if (.Platform$OS.type=="unix") {
        devtype <- quartz
      } else {
        devtype <- windows
      }
      devtype(width = 2, height = 0.5)
      progrpl <- dev.cur()
      par(mar = rep(0,4))
    }
    for(g in 1:niter) {
      
      # 1.- SAMPLING:
      # a) Sample survival parameters:
      # i) Metropolis draw for params:
      theta.n <- matrix(rtnorm(length.full.theta, theta.g, 
              theta.jump, lower = low.full.theta), 
          length.cat, length.theta, 
          dimnames = dimnames(theta.g))
      if (shape!="simple") {
        nlow[, 'c'] <- apply(theta.n, 1, CalculateLowC)
        idc.low <- which(theta.n[, 'c'] < nlow[, 'c'])
        if (length(idc.low)>0) {
          for(cc in idc.low) {
            theta.n[cc,'c'] <- c(rtnorm(1, theta.g[cc,'c'], 0.5, 
                    lower=nlow[cc,'c']))
          }
        }
      }
      
      if (Cont) {
        gamma.n <- rnorm(length.cont, gamma.g, gamma.jump) 
      } else {
        gamma.n <- gamma.g
      }
      
      # ii) Build individual parameter matrices:
      Ztheta.n <- Zcat %*% theta.n
      Zgamma.n <- Zcont %*% gamma.n
      
      # iii) Calculate conditional posteriors:
      # - Current parameters:
      p.thg <- (log(CalculateFullFx(xag, 
                    Ztheta.g, Zgamma.g)) -
            log(CalculateFullSx(xatg, 
                    Ztheta.g, Zgamma.g))) * Iag 
      
      # - Priors:
      p.thg <- sum(p.thg) + sum(dtnorm(c(theta.g[idUpdJump]), 
                  c(theta.prior[idUpdJump]), theta.sd, 
                  lower = low.full.theta[idUpdJump], log = TRUE)) + 
          sum(dnorm(gamma.g, gamma.prior, gamma.sd, 
                  log = TRUE))
      
      # - Proposed parameters:
      p.thn <- (log(CalculateFullFx(xag, 
                    Ztheta.n, Zgamma.n)) -
            log(CalculateFullSx(xatg, 
                    Ztheta.n, Zgamma.n))) * Iag
      
      # - Priors:
      p.thn <- sum(p.thn) + sum(dtnorm(c(theta.n[idUpdJump]), 
                  c(theta.prior[idUpdJump]), theta.sd, 
                  lower = low.full.theta[idUpdJump], log=TRUE)) + 
          sum(dnorm(gamma.n, gamma.prior, gamma.sd,
                  log=TRUE))
      
      r <- exp(p.thn-p.thg)
      z <- runif (1,0,1)
      
      if (is.na(r)) {
        naflag <- TRUE 
      } else {
        if (r>z) {
          theta.g <- theta.n
          Ztheta.g <- Ztheta.n
          p.thg <- p.thn
          gamma.g <- gamma.n
          Zgamma.g <- Zgamma.n
          updVec[g] <- 1
        }
      }
      
      # b) Sample times of birth and death:
      # i) New times of birth and death:
      bn <- bg 
      bn[bi0] <- bg[bi0] + sample(-1:1, length(bi0), replace = TRUE) 
      bn[bi0][oi[bi0] > 0] <- apply(cbind(bn[bi0][oi[bi0] > 0],
              first.obs[bi0][oi[bi0] > 0] - 1), 1, min)
      bn[bi0][oi[bi0] == 0] <- apply(cbind(bn[bi0][oi[bi0] == 0],
              dg[bi0][oi[bi0] == 0] - 1), 1, min)
      dn <- dg 
      dn[di0] <- dg[di0] + sample(-1:1, length(di0), 
          replace = TRUE) 
      dn[di0] <- apply(cbind(dn[di0],bn[di0],last.obs[di0] + 1), 
          1, max) 
      xn <- dn - bn
      
      # ii) New full alive matrices:
      Fn <- c(apply(cbind(studyStart, bn + 1), 1, max))
      Ln <- c(apply(cbind(studyEnd, dn - 1), 1, min))
      On <- BuildAliveMatrix(Fn, Ln, Tm)
      
      # iii) New juvenile and adult ages:
      Ian <- rep(0, n)
      Ijn <- Ian
      Ian[xn >= minAge] <- 1
      Ijn[xn < minAge] <- 1
      xjn <- xn
      xjn[xn > minAge] <- minAge
      xan <- xn - minAge
      xan[xn < minAge] <- 0
      xjtn <- xn * 0
      idtr <- which(bn < studyStart & 
              studyStart - bn < minAge)
      xjtn[idtr] <- studyStart - bn[idtr]
      xatn <- xn * 0
      idtr <- which(bn + minAge < studyStart)
      xatn[idtr] <- studyStart - (bn[idtr] + minAge)
      
      # iiii) Calculate conditional posteriors:
      # - Current ages:
      p.bdg <- (log(lag) * Ijg - lag * xjg) * IminAge + 
          log(CalculateFullFx(xag, 
                  Ztheta.g, Zgamma.g)) * Iag
      p.bdg <- p.bdg + (Og - lfi) %*% log(1 - Pig) + 
          log(v.x(xag + 0.5 * Dx)) * Iag
      
      # - New ages:
      p.bdn <- (log(lag) * Ijn - lag * xjn) * IminAge  + 
          log(CalculateFullFx(xan, 
                  Ztheta.g, Zgamma.g)) * Ian
      p.bdn <- p.bdn + (On - lfi) %*% log(1 - Pig) + 
          log(v.x(xan + 0.5 * Dx)) * Ian
      
      r <- exp(p.bdn-p.bdg)
      idNa <- which(!is.na(r))
      z <- runif (length(idNa), 0, 1)
      idrz <- which(r[idNa] > z)
      bg[idNa[idrz]] <- bn[idNa[idrz]]
      dg[idNa[idrz]] <- dn[idNa[idrz]]
      xg[idNa[idrz]] <- xn[idNa[idrz]]
      p.bdg[idNa[idrz]] <- p.bdn[idNa[idrz]]
      Og[idNa[idrz], ] <- On[idNa[idrz], ]
      xjg[idNa[idrz]] <- xjn[idNa[idrz]]
      xjtg[idNa[idrz]] <- xjtn[idNa[idrz]]
      xag[idNa[idrz]] <- xan[idNa[idrz]]
      xatg[idNa[idrz]] <- xatn[idNa[idrz]]
      Iag[idNa[idrz]] <- Ian[idNa[idrz]]
      Ijg[idNa[idrz]] <- Ijn[idNa[idrz]]
      
      # c) Sample recapture probability(ies):
      rho1g <- rho1 + t(t(Y) %*% rep(1, n))
      rho2g <- rho2 + t(t(Og - Y) %*% rep(1, n))
      Rho1 <- tapply(rho1g, idpi, sum)
      Rho2 <- tapply(rho2g, idpi, sum)
      pi.g <- rbeta(npi, Rho1, Rho2)
      if (1 %in% pi.g) {
        pi.g[pi.g==1] <- 1-1e-5
        warning("Some recapture probabilities are equal to 1.",
            "\nThey have been constraint to be fractionally less than 1 ",
            "for computational reasons\n", call. = FALSE)
      }
      Pig <- pi.g[idpi]
      
      # d) if minAge > 0, sample lambda:
      if (minAge > 0) {
        lan <- rtnorm(n = 1, mean = lag, 
            sd = 0.001, lower = 0)
        p.lag <- sum(log(lag) * Ijg - lag * xjg + 
                    lag * xjtg) +
            dtnorm(lag, mean = 0.01, sd = 1, lower = 0)
        p.lan <- sum(log(lan) * Ijg - lan * xjg + 
                    lan * xjtg) +
            dtnorm(lan, mean = 0.01, sd = 1, lower = 0)
        
        r <- exp(p.lan - p.lag)
        z <- runif(1, 0, 1)
        if (r > z) {
          lag <- lan
        }
      }
      
      # 2.- STORE RESULTS:
      # Parameters and latent states:
      theta.mat[g, ] <- theta.g
      pi.mat[g, ] <- pi.g
      if (Cont) {
        gamma.mat[g, ] <- gamma.g
      }
      if (minAge > 0) {
        la.vec[g] <- lag
      }
      if (g %in% thin.seq) {
        bi.mat[gg, ] <- bg
        di.mat[gg, ] <- dg
        gg <- gg + 1
      }
      
      # Conditional posteriors: 
      posterior.mat[g, ] <- c(p.thg, sum(p.bdg), p.thg + 
              sum((Og - lfi) %*% log(1 - Pig)))
      
      # Update Jumps:
      if (g %in% updateInt & jumpObject$update & updateJumps) {
        jumpObject <- UpdateJumps(jObject = jumpObject, updateVec = updVec, 
            targetUpdate = 0.2, g, 
            updateInt, nPar, updateLen)
        theta.jump[idUpdJump] <- jumpObject$jump[1:length(idUpdJump)]
        if (Cont) {
          gamma.jump <- jumpObject$jump[-c(1:length(idUpdJump))]
        }
      }
      # Progress plot:
      if (g %in% round(seq(1, niter, length = 100)) & progrPlots) {
        par(mar = rep(0, 4))
        plot(x  = c(0, niter * 1.1), y = c(0, 1), axes = FALSE, col = NA, 
            xlab = "", ylab = "")
        polygon(x = c(0, niter, niter, 0), y = c(0.35, 0.35, 0.65, 0.65), 
            col = NA, border = 'dark red')
        polygon(x = c(0, g, g, 0), y = c(0.35, 0.35, 0.65, 0.65), 
            col = 'dark red', border = 'dark red')
        text(x  = niter / 2, y = 0.85, labels = paste("MCMC progress (Sim. ", 
                sim, ")", sep = ""), cex = 0.9)
        text(x = g, y = 0.15, labels = paste(round(g / niter * 100), 
                "%", sep = ""), cex = 0.8)
      }
    }
    if (progrPlots) {
      dev.off(progrpl)
    }
    # Return results:
    return(list(theta = theta.mat, gamma = gamma.mat, pi = pi.mat, 
            la = la.vec, bi = bi.mat, di = di.mat, post = posterior.mat, 
            g = g, naflag = naflag, jObject = jumpObject))
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
  Start <- Sys.time()
  if (parallel) {
    avail.pkgs <- available.packages()
    if (!is.element("snowfall", avail.pkgs)) {
      warning("\nPackage 'snowfall' is not installed.\nSimulations ",
          "will not be ran in parallel (computing time will ",
          "be longer...)\n")
      basta.out <- lapply(1:nsim, multiMCMC)
    } else {
      require(snowfall)
      sfInit(parallel = TRUE, cpus = ncpus);
      sfExport(list = c(parallelVars, "parallel", "nsim", "minAge"))
      sfLibrary(msm)
      basta.out <- sfClusterApplyLB(1:nsim, multiMCMC)
      sfStop()
    }
  } else {
    basta.out <- lapply(1:nsim, multiMCMC)
  }
  End <- Sys.time()
  
  # 7.2 Report if all simulations ran through:
  simNames <- paste("Sim.", (1:nsim), sep="")
  full.runs <- rep(0,nsim)
  names(full.runs) <- simNames
  last.steps <- full.runs
  for(i in 1:nsim) {
    last.steps[i] <- basta.out[[i]]$g	
    full.runs[i] <- ifelse(last.steps[i] == niter, 1, 0)
  } 
  id.failed <- which(full.runs == 0)
  all.ran <- FALSE
  if (nsim==1) {
    if (full.runs==1) {
      cat("MCMC finished running\n")
      cat(paste("Total MCMC computing time: ", 
              round(as.numeric(julian(End) - julian(Start)) * 24 * 60, 2), 
              " minutes\n\n", sep=""))
      all.ran <- TRUE
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
      all.ran <- TRUE
      cat("\nMultiple simulations finished.\n")
      cat(paste("Total MCMC computing time: ", 
              round(as.numeric(julian(End) - julian(Start)) * 24 * 60, 2), 
              " minutes\n\n", sep = ""))
    }
  }	
  
  # 8. Diagnostics:
  thin.seq <- seq(burnin, niter, thinning)
  nthin <- length(thin.seq)
  
  # 8.1 Thinned result matrices:
  if (Cont) {
    name.out.mat <- c(name.full.theta, name.gamma, name.pi, 
        name.post) 
  } else {
    name.out.mat <- c(name.full.theta, name.pi, name.post)
  }
  out.mat <- matrix(NA, niter * nsim, length(name.out.mat))
  dimnames(out.mat) <- list(rep(simNames, each = niter), name.out.mat)
  Bimat <- matrix(NA, nthin * nsim, n)
  rownames(Bimat) <- rep(simNames, each = nthin)
  Dimat <- Bimat
  idthin <- rep(0, niter*nsim)
  if (minAge > 0) {
    la.mat <- matrix(NA, nrow = niter, ncol = nsim, 
        dimnames = list(NULL, simNames))
  } else {
    la.mat <- matrix(NA, nrow = 1, ncol = nsim,
        dimnames = list(NULL, simNames))
  }
  
  for(i in 1:nsim) {
    Idsim <- which(rownames(out.mat) == simNames[i])
    if (Cont) {
      out.mat[Idsim, ] <- cbind(basta.out[[i]]$theta, basta.out[[i]]$gamma, 
          basta.out[[i]]$pi, basta.out[[i]]$post)
    } else {
      out.mat[Idsim, ] <- cbind(basta.out[[i]]$theta, basta.out[[i]]$pi, 
          basta.out[[i]]$post)
    }
    idthin[Idsim[thin.seq]] <- 1
    Idsim <- which(rownames(Bimat) == simNames[i])
    Bimat[Idsim, ] <- basta.out[[i]]$bi
    Dimat[Idsim, ] <- basta.out[[i]]$di
    if (minAge > 0){
      la.mat[, i] <- basta.out[[i]]$la
    } 
  }
  
  # 8.2 Basic summary statistics for parameters:
  par.mat <- out.mat[idthin == 1, 
      -(c(ncol(out.mat) - c(2:0)))]
  coef <- cbind(apply(par.mat, 2, mean, na.rm=TRUE), apply(par.mat, 2, 
          sd, na.rm=TRUE), t(apply(par.mat, 2, quantile, c(0.025, 0.975), 
              na.rm = TRUE)), NA, NA, NA)
  colnames(coef) <- c("Estimate", "StdErr", "Lower95%CI", "Upper95%CI", 
      "SerAutocor", "UpdateRate", "PotScaleReduc")
  if (length(id.failed) < nsim) {
    idfix <- which(theta.jump==0)
    if (length(idfix) > 0) {
      coef[-idfix,"SerAutocor"] <- apply(par.mat[, -c(idfix)], 2, 
          function(x) cor(x[-1], x[-length(x)], use = "complete.obs"))
      coef[idfix,"SerAutocor"] <- 1
    } else {
      coef[, "SerAutocor"] <- apply(par.mat, 2, function(x) cor(x[-1], 
                x[-length(x)], use = "complete.obs"))
    }
    coef[, "UpdateRate"] <- apply(out.mat[, -c(ncol(out.mat) - c(2:0))], 
        2, function(x) 
          length(which(diff(x[!is.na(x)]) != 0)) / 
              length(x[!is.na(x)]))
  }
  out.mat <- cbind(idthin, out.mat)
  
  # 8.3 Convergence and model selection:
  if (all.ran) {
    if (nsim>1) {
      # 8.3.1 Convergence diagnostics (potential scale reduction):
      Means <- apply(par.mat, 2, function(x) 
            tapply(x, rownames(par.mat), mean))
      Vars <- apply(par.mat, 2, function(x) 
            tapply(x, rownames(par.mat), var))
      meanall <- apply(Means, 2, mean)
      B <- nthin / (nsim - 1) * apply(t((t(Means) - meanall)^2), 2, sum)
      W <- 1 / nsim * apply(Vars, 2, sum)
      Varpl <- (nthin - 1) / nthin * W + 1 / nthin * B
      Rhat <- sqrt(Varpl / W)
      Rhat[Varpl==0] <- 1
      conv <- cbind(B, W, Varpl, Rhat)
      rownames(conv) <- colnames(par.mat)
      coef[, ncol(coef)] <- conv[, 'Rhat']
      
      # Report if convergence was reached:
      idnconv <- which(conv[, 'Rhat']< 0.95 | conv[, 'Rhat']>1.1)
      if (length(idnconv) > 0) {
        modSel <- NULL
        kl.list <- NULL
        warning("Convergence not reached for some survival parameters.",
            "\nDIC could not be calculated.\n", call. = FALSE)
      } else {
        # 8.3.2 Model selection (DIC, if convergence was reached):
        posterior <- out.mat[idthin == 1, ncol(out.mat)]
        L <- length(posterior)
        Dm <- -2 * posterior
        Dmode <- -2 * posterior[which(posterior ==  max(posterior))[1]]
        Dave <- mean(Dm)
        pD <- Dave - Dmode
        k <- npi + length.full.theta
        if (Cont) {
          k <- k + length.cont
        }
        DIC <- 2 * Dave - Dmode
        modSel <- c(Dave, Dmode, pD, k, DIC)
        names(modSel) <- c("D.ave", "D.mode", "pD", "k", "DIC")
        cat("Survival parameters converged appropriately.",
            "\nDIC was calculated.\n")
        
        # 8.3.3 Inference on parameter estimates:
        # Kullback-Leibler distances for categorical covariates:
        if (is.null(covariate.type$cat)) {
          kl.list <- NULL
        } else {
          name.cat <- names(covariate.type$cat)
          n.cat <- length(name.cat)
          if (n.cat > 1) {
            n.comb <- (n.cat - 1)^2 - 
                ((n.cat - 1)^2 - (n.cat - 1)) / 2
            covar.comb <- matrix(0, n.comb, 2, 
                dimnames = list(NULL, c("cov1", "cov2")))
            ij <- 1
            for (i in 1:(n.cat - 1)) {
              for (j in 2:n.cat) {
                if (i < j) {
                  covar.comb[ij, ]<- c(name.cat[i], name.cat[j])
                  ij <- ij + 1
                }
              }
            }
            if (covarsStruct == "fused") {
              kl.12 <- matrix(NA, nrow = n.comb, 
                  ncol = length.theta0, 
                  dimnames = list(paste(covar.comb[, 1], 
                          "-", covar.comb[, 2], sep = ""),
                      name.theta0))
              kl.21 <- matrix(NA, nrow = n.comb, 
                  ncol = length.theta0, 
                  dimnames = list(paste(covar.comb[, 2], 
                          "-", covar.comb[, 1], sep = ""),
                      name.theta0))
              p.low <- low.full.theta[1,]
            } else {
              kl.12 <- matrix(NA, nrow = n.comb, ncol = 1, 
                  dimnames = list(paste(covar.comb[,1], 
                          "-", covar.comb[,2], sep = ""), "gamma"))
              kl.21 <- matrix(NA, nrow = n.comb, ncol = 1, 
                  dimnames = list(paste(covar.comb[, 2], 
                          "-", covar.comb[, 1], 
                          sep = ""), "gamma"))
              p.low <- -Inf
            }
            q.12 <- kl.12
            q.21 <- kl.21
            kl.pars <- colnames(kl.12)
            for(i in 1:ncol(kl.12)){
              for(j in 1:n.comb){
                p1 <- par.mat[, paste(kl.pars[i], "[", 
                        covar.comb[j, 1], "]", sep = "")]
                mean.p1 <- mean(p1)
                sd.p1 <- sd(p1)
                p2 <- par.mat[, paste(kl.pars[i], "[", 
                        covar.comb[j, 2], "]", sep = "")]
                mean.p2 <- mean(p2)
                sd.p2 <- sd(p2)
                full.range <- range(c(mean.p1 + c(-4, 4) * sd.p1, 
                        mean.p2 + c(-4, 4) * sd.p2))
                full.range[1] <- max(full.range[1], p.low[i])
                p.vec <- seq(full.range[1], full.range[2], 
                    length = 500)
                dp <- p.vec[2] - p.vec[1]
                dens.p1 <- dtnorm(p.vec, mean = mean(p1), sd = sd(p1), 
                    lower = p.low[i])  
                dens.p2 <- dtnorm(p.vec, mean = mean(p2), sd = sd(p2), 
                    lower = p.low[i])  
                kl.12[j, i] <- sum(dens.p1*log(dens.p1/dens.p2) * 
                        dp)
                kl.21[j, i] <- sum(dens.p2*log(dens.p2/dens.p1) * 
                        dp)
                q.12[j, i] <- (1 + (1 - exp(-2 * kl.12[j, i])^(1 / 2))) / 2
                q.21[j, i] <- (1 + (1 - exp(-2 * kl.21[j, i])^(1 / 2))) / 2
              }
            }
            kl.list <- list(kl12 = kl.12, kl21 = kl.21, 
                q12  = q.12, q21  = q.21)
          } else {
            kl.list <- NULL
          }
        }
      }
    } else {
      conv <- NULL
      modSel <- NULL
      kl.list <- NULL
    }
    
    # 8.3.3 Summary times of birth and ages at death:
    xq <- apply(Dimat - Bimat, 2, quantile, c(0.5, 0.025, 0.975))
    bq <- apply(Bimat, 2, quantile, c(0.5, 0.025, 0.975))
    
    # 8.3.4 Summary Survival and mortality:
    thmat <- matrix(par.mat[, name.full.theta], ncol = length.full.theta, 
        dimnames = list(NULL, name.full.theta))
    if (Cont) {
      rzc <- apply(Zcont, 2, quantile, c(0.5, 0.025, 0.975))
      rownames(rzc) <- c("Med.", "Lower", "Upper")
      gave <- apply(as.matrix(par.mat[, which(substr(colnames(par.mat), 
                          1, 2) == "gamma")]), 2, mean)
    } else if (covarsStruct == "all.in.mort") {
      rzc <- apply(as.matrix(Zcat[, covariate.type$cont]), 
          2, quantile, c(0.5, 0.025, 0.975))
      dimnames(rzc) <- list(c("Med.", "Lower", "Upper"), 
          names(covariate.type$cont))
      gave <- 0
      zcname <- colnames(rzc)
    } else {
      rzc <- matrix(0, 1, 1, dimnames = list("nc", "nc"))
      gave <- 0
      zcname <- c("")
    }
    Sxq <- list()
    mxq <- list()
    xvec <- list()
    zaname <- c(names(covariate.type$int), names(covariate.type$cat))
    if (is.null(zaname)) {
      zaname <- "NoCov"
    } 
    for(i in 1:length(zaname)) {
      if (zaname[i] == "NoCov") {
        idza <- 1:n
      } else {
        idza <- which(Z[, zaname[i]] == 1)
      }
      xv <- seq(0, ceiling(max(xq[1, idza]) * 1.1), 0.1)
      xvec[[zaname[i]]] <- xv + minAge
      for(j in 1:ncol(rzc)) {
        Sxq[[zaname[i]]][[colnames(rzc)[j]]] <- 
            array(0, dim = c(length(xv), 3, nrow(rzc)), 
                dimnames = list(NULL, c("50%", "2.5%", "97.5%"), 
                    rownames(rzc)))
        
        mxq[[zaname[i]]][[colnames(rzc)[j]]] <- 
            Sxq[[zaname[i]]][[colnames(rzc)[j]]]
        for(k in 1:nrow(rzc)) {
          gaa <- sum(gave * rzc[k, j])
          Cols <- paste(name.theta, "[",zaname[i], "]", sep = "")
          if (length.cat == 1) Cols <- name.theta
          Thm <- matrix(thmat[, Cols], ncol = length(Cols))
          if (covarsStruct == "all.in.mort") {
            Thm <- Thm + thmat[, paste(name.theta, 
                    "[",names(covariate.type$cont)[j],"]", 
                    sep = "")] * rzc[k, j]
          }
          Sxq[[zaname[i]]][[colnames(rzc)[j]]][, , k] <- 
              t(apply(apply(Thm, 1 ,CalculateMultiSx), 1, quantile, 
                      c(0.5, 0.025, 0.975)))
          if (model == "EX") {
            Thmm <- matrix(Thm, nrow = nrow(Thm), 
                    ncol = length(xv), byrow = TRUE) *
                exp(gaa)
            mxq[[zaname[i]]][[colnames(rzc)[j]]][, , k] <- 
                t(apply(Thmm,2, quantile, c(0.5, 0.025, 0.975)))
          } else {
            mxq[[zaname[i]]][[colnames(rzc)[j]]][, , k] <- 
                t(apply(apply(Thm, 1 ,CalculateMultiMx),
                        1, quantile, c(0.5, 0.025, 0.975)))
          }
        }
      }
    }
    
    # 8.3.5 Calculate life table from estimated ages at death:
    if (lifeTable) {
      LT  <- list()
      for(i in 1:length(zaname)) {
        if (zaname[i] == "NoCov") {
          idza <- which(bq[1, ] >= studyStart)
        } else {
          idza <- which(Z[, zaname[i]] == 1 &  bq[1, ] >= studyStart)
        }
        x <- xq[1, idza]
        tempLT <- MakeLifeTable(x, ax = 0.5, n = 1)
        tempLT <- subset(tempLT, tempLT$StartAge >= minAge)
        rownames(tempLT) <- NULL
        LT[[zaname[i]]] <- tempLT
      }
    } else {
      LT <- NULL
    }
  } else {
    conv <- NULL
    modSel <- NULL
    kl.list <- NULL
    xq <- NULL
    bq <- NULL
    Sxq <- NULL
    mxq <- NULL
    xvec <- NULL
    if (lifeTable) {
      LT <- NULL
    }
  }
  
  # 9. Return a list object of class 'basta':
  Settings <- c(niter, burnin, thinning, nsim)
  names(Settings) <- c("niter", "burnin", "thinning", "nsim") 
  ModelSpecs <- c(model, shape, covarsStruct, 
      paste(names(covariate.type$cat), 
          collapse = ", "), 
      paste(names(covariate.type$cont), 
          collapse = ", "))
  names(ModelSpecs) <- c("model", "shape", "Covar. structure", 
      "Categorical", "Continuous")
  Priors <- c(theta.prior)
  jumpPriorName <- name.full.theta
  if (Cont) {
    Priors <- c(Priors, gamma.prior)
    jumpPriorName <- c(jumpPriorName, name.gamma)
  }
  Jumps <- matrix(0, length(jumpPriorName), nsim)
  
  for (JJ in 1:nsim) {
    Jumps[idUpdJump, JJ] <- basta.out[[JJ]]$jObject$jump[1:length(idUpdJump)]
    if (Cont) {
      Jumps[-c(1:length.full.theta)] <- 
          basta.out[[JJ]]$jObject$jump[-c(1:length(idUpdJump))]
    }
  }
  JumpPriors <- cbind(Priors, Jumps)
  dimnames(JumpPriors) <- list(jumpPriorName, 
      c("Mean.priors", paste("Jump.sd.sim.", 1:nsim, sep = "")))
  output <- list()
  output$coefficients <- coef
  output$DIC <- modSel
  output$Convergence <- conv
  output$KullbackLeibler <- kl.list
  output$settings <- Settings
  output$ModelSpecs <- ModelSpecs
  output$JumpPriors <- JumpPriors
  output$Params <- out.mat
  output$Bis <- Bimat
  output$Dis <- Dimat
  output$Bq <- bq
  output$Xq <- xq
  output$Sx <- Sxq
  output$mx <- mxq
  output$xv <- xvec
  output$bd <- bd
  output$Y <- Y
  output$Zcat <- Zcat
  output$Zcont <- Zcont
  output$studyStart <- studyStart
  output$studyEnd <- studyEnd
  output$finished <- full.runs
  if (lifeTable) {
    output$lifeTable <- LT
  }
  if (minAge > 0) {
    output$Lambda <- la.mat
  }
  class(output) <- "basta"
  
  ## Check Makeham terms, return warning if the Makeham (c) terms overlap 0.
#  if(shape=="Makeham"){
#    MakehamCoefRow <- substr(rownames(output$coef), 1, 1) == "c"
#    MakehamLower95pcCI <- output$coef[MakehamCoefRow, 3]
#    MakehamUpper95pcCI <- output$coef[MakehamCoefRow, 4]
  # Number of coefficients where Lower95%CI < 0 & Upper95%CI > 0?
#    NLower <- length(which(MakehamLower95pccCI <= 0 & 
#                MakehamUpper95pccCI >= 0))
#    NLower <- sum(MakehamLower95pcCI <= 0)
#    if (NLower > 0) {
#    warning(paste(NLower, " of the Makeham coefficients have 95% CI\n",
#              "overlapping 0. Makeham models may not be appropriate", sep = ""), 
#          call. = FALSE)
#    }
#  }
  return(output)
  }
  