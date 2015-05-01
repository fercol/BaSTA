plot.basta <-
    function(x, plot.trace = TRUE, trace.name = "theta", ...){
  
  if (x$settings['nsim'] <= 8) {
    Palette <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
        '#A65628', '#F781BF', '#999999')
  } else {
    Palette <- grey(1:x$settings['nsim'] / x$settings['nsim'])
  }
  # 1. Trace plots for parameters:
  if(plot.trace){
    # OLD VERSION (To be deprecated in the future):
    if ("mx" %in% names(x)) {
      length.cat <- ncol(x$Zcat)
      # 1. Trace plots for parameters:
      traces          <- c("theta","gamma","pi", "post")
      if(!is.element(trace.name, traces)) {
        stop(paste("Wrong 'trace.name' argument. Valid arguments are:", 
                paste(paste("'", traces, "'", sep = ""), collapse = ", "), 
                ".\n"), call. = FALSE)
      }
      varNames <- substr(colnames(x$Par)[-1], 1, 2)
      id.th  <- is.element(substr(varNames, 1, 1), c("a", "b", "c"))
      varNames[id.th] <- "th"
      nTraces <- table(varNames)[c("th", "ga", "pi", "po")]
      names(nTraces) <- c("th", "ga", "pi", "po")
      nTraces[is.na(nTraces)] <- 0
      if (trace.name == "gamma" & nTraces["ga"] == 0) {
        stop("\nTrace plots cannot be drawn for 'gamma' parameters.",
            "\nNo proportional hazards arguments were evaluated.", 
            call. = FALSE)
      }
      niter <- x$settings['niter']
      nsim  <- x$settings['nsim' ]
      simNames <- unique(rownames(x$Par))
      idpl <- which(varNames == substr(trace.name, 1, 2))
      X <- as.matrix(x$Par[, -1][, idpl])
      colnames(X) <- colnames(x$Par)[-1][idpl]
      if (niter > 1e3) {
        xThin <- round(seq(1, niter, length = 1e3))
      } else {
        xThin <- 1:niter
      }
      p <- which(traces == trace.name)
      Cols <- Palette[round(seq(1, 12, length = nsim))]
      model <- as.character(x$ModelSpecs['model'])
      shape <- as.character(x$ModelSpecs['shape'])
      if (model == "EX") {
        nthm <- 1
      } else if (model == "GO" | model == "WE") {
        nthm <- 2
      }  else {
        nthm  <- 3
      }
      if (shape == "Makeham") {
        nthm <- nthm + 1
      } else if(shape == "bathtub") {
        nthm <- nthm + 3
      }
      ydim <- c(nthm, ceiling(nTraces['ga'] / 2), 
          ceiling(nTraces['pi'] / 2), 2)
      xdim <- c(length.cat, 2, 2, 2)
      for(ii in 2:3) if (ydim[ii] == 1) xdim[ii] <- 1
      op <- par(mfrow = c(ydim[p], xdim[p]), 
          mar   = c(3, 3, 3, 1))
      for(i in 1:nTraces[p]){
        x <- X[, i]
        yl <- range(x, na.rm = TRUE)
        plot(x = c(1, niter), y = yl, col = NA, xlab = "Iteration", 
            ylab = "", main = colnames(X)[i], frame.plot = FALSE)
        for(j in 1:nsim) {
          lines(x = xThin, y = x[which(names(x) == simNames[j])][xThin],
              type = 'l', col = Cols[j], lwd = 1.5)
        }
      }
      
    } else {
      # NEW VERSION:
      traces <- c("theta","gamma","pi")
      parNames <- substr(colnames(x$params), 1, 2)
      if (trace.name == "theta") {
        idPars <- which(!(parNames %in% c("ga", "pi")))
      } else if (trace.name == "gamma") {
        if (trace.name == "gamma" & x$modelSpecs[3] %in% c("inMort", "noCov")) {
          stop("\nTrace plots cannot be drawn for 'gamma' parameters.",
              "\nNo proportional hazards arguments were evaluated.", 
              call. = FALSE)
        } else {
          idPars <- which(parNames == "ga")
        }
      } else if (trace.name == "pi") {
        idPars <- which(parNames == "pi")
      } else {
        stop("Wrong 'trace.name' argument. Valid arguments are:\n", 
            "'theta', 'gamma' or 'pi'.\n", call. = FALSE)
      }
      
      nRows <- ifelse(trace.name == "theta", length(unique(parNames[idPars])),
          ceiling(length(idPars) / 2))
      nCols <- ifelse(trace.name == "theta", table(parNames[idPars])[1], 
          ceiling(length(idPars) / nRows))
      yLim <- sapply(colnames(x$params)[idPars], function(par) 
            range(sapply(1:x$settings["nsim"], function(sim)
                      range(x$parsForPlot[[sim]][, par]))))
      op <- par(mfrow = c(nRows, nCols), mar = c(3, 3, 2.5, 1))
      for (i in 1:length(idPars)) {
        plot(c(1, x$settings['niter']), col = NA, yLim[, i], xlab = "", 
            ylab = "", frame.plot = FALSE, main = colnames(x$params)[idPars[i]])
        for (j in 1:x$settings['nsim']) {
          lines(seq(1, x$settings['niter'], x$settings['thinning']), 
              x$parsForPlot[[j]][, idPars[i]], col = Palette[j])
        }
      }
    }
    par(op)
    
    # 2. Plot survival and mortality:
  } else {
    # OLD VERSION:
    if ("mx" %in% names(x)) {
      length.cat <- ncol(x$Zcat)
      if(is.null(x$Sx)){
        stop("MCMC runs on BaSTA did not finish.\n Survival and mortality ",
            "plots cannot be constructed, verify model input and run again.\n",
            call. = FALSE)
      }
      zname <- names(x$Sx)
      length.cat <- length(zname)
      Bord <- Palette[round(seq(1, 12, length = length.cat))]
      Cols <- adjustcolor(Bord, alpha.f = 0.5)
      op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
      
      # Plot survival probability:
      mxv <- ceiling(max(unlist(x$xv)) / 5) * 5
      plot(x = c(0, mxv), y = range(0, 1), col = NA, xlab = "", 
          ylab  = expression(S(x)), main = "Survival probability", 
          frame.plot  = FALSE, ...)
      for(i in 1:length.cat){
        xv <- x$xv[[i]]
        polygon(x = c(xv, rev(xv)), 
            y = c(x$Sx[[i]][[1]][, 2, 1], rev(x$Sx[[i]][[1]][, 3, 1])), 
            col = Cols[i], border = Bord[i])
        lines(x = xv, y = x$Sx[[i]][[1]][,1,1], col = Bord[i], lty = 3)
      }
      if (length.cat > 1) {
        legend(x = 'topright', legend = zname, pch = 15, pt.cex = 3, 
            cex = 1.25, col = Cols, bty = 'n')
      }
      
      # Plot mortality rates:
      ylmx <- c(0, round(max(unlist(x$mx))))
      plot(x = c(0,mxv), y = ylmx, col = NA, ylim = ylmx, xlab  = "Age (x)", 
          ylab = expression(mu(x)), main = "Mortality", frame.plot  = FALSE, 
          ...)
      for(i in 1:length.cat){
        xv <- x$xv[[i]]
        polygon(x = c(xv, rev(xv)), 
            y  = c(x$mx[[i]][[1]][,2,1], rev(x$mx[[i]][[1]][,3,1])), 
            col = Cols[i], 
            border = Bord[i])
        lines(x = xv, y = x$mx[[i]][[1]][, 1, 1], col = Bord[i], lty = 3)
      }
      
      # NEW VERSION:  
    } else {
      ncovs <- length(x$survQuant)
      if (ncovs <= 8) {
        Bord <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
            '#A65628', '#F781BF', '#999999')[round(seq(1, 8, 
                    length = length(x$mortQuant)))]
      } else {
        Bord <- rainbow(ncovs)
      }
      Cols <- Bord
      for(cl in 1:length(Cols)) {
        Cols[cl] <- adjustcolor(Bord[cl], alpha.f = seq(0.1, 0.4, 
                length = ncovs)[cl])
      }
      op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
      xv <- lapply(1:ncovs, function(idcovs) as.numeric(colnames(x$mortQuant[[idcovs]])))
      if (!("xlim" %in% names(list(...)))) {
        xlim <- c(0, max(unlist(xv)))
        idRange <- lapply(1:ncovs, function(idcovs) length(xv[[idcovs]]))
      } else {
        xlim <- list(...)$xlim
        idRange <- lapply(1:ncovs, function(idcovs) which(abs(xv[[idcovs]] - xlim[2]) == min(abs(xv[[idcovs]] - xlim[2])))[1])
        xv <- lapply(1:ncovs, function(idcovs) xv[[idcovs]][1:idRange[[idcovs]]])
      }
      if ("noCI" %in% names(list(...))) {
        noCI <- list(...)$noCI
        if (class(noCI) != "logical") {
          noCI <- FALSE
        }
      } else {
        noCI <- FALSE
      }
      for (dem in c("survQuant", "mortQuant")) {
        if (noCI) {
          idQuants <- 1
        } else {
          idQuants <- 1:3
        }
        ylim <- c(0, max(sapply(1:length(x[[dem]]), function(idem) 
                      max(x[[dem]][[idem]][idQuants, 2:idRange[[idem]]]))))
        xlab <- ifelse(dem == "mortQuant", "Age", "")
        ylab <- ifelse(dem == "mortQuant", expression(mu(x)), expression(S(x)))
        main <- ifelse(dem == "mortQuant", "Mortality", "Survival")
        plot(xlim, ylim, col = NA, xlab = xlab, ylab = ylab, 
            frame.plot = FALSE, main = main, ylim = c(0, ylim[2]))
        if (!noCI) {
          for (cov in 1:length(x$mortQuant)) {
            polygon(c(xv[[cov]][1:idRange[[cov]]], rev(xv[[cov]][1:idRange[[cov]]])), 
                c(x[[dem]][[cov]][2, 1:idRange[[cov]]], 
                    rev(x[[dem]][[cov]][3, 1:idRange[[cov]]])),
                col = Cols[cov], border = Bord[cov], lwd = 0.5, lty = 1)
          }
        }
        for (cov in 1:length(x$mortQuant)) {
          lines(xv[[cov]][1:idRange[[cov]]], x[[dem]][[cov]][1, 1:idRange[[cov]]], col = Bord[cov],
              lwd = 3)
        }
        if (dem == "survQuant" & length(x$mortQuant) > 1) {
          legend('topright', substr(names(x$mortQuant), 2, 
                  nchar(names(x$mortQuant))), 
              col = Bord, lwd = 3, bty = 'n')
        }
      }
    }
    par(op)
  }
}
