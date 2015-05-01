summary.basta <-
function(object,...){
  extraArgs       <- list(...)
  cat("\nCall:\n")
  cat(paste("Model             \t\t: ", object$ModelSpecs[1], "\n", sep = ""))
  cat(paste("Shape             \t\t: ", object$ModelSpecs[2], "\n", sep = ""))
  cat(paste("Covars. structure \t\t: ", object$ModelSpecs[3], "\n", sep = ""))
  cat(paste("Cat. covars.      \t\t: ", object$ModelSpecs[4], "\n", sep = ""))
  cat(paste("Cont. covars.     \t\t: ", object$ModelSpecs[5], "\n", 
            collapse = ""))

  cat("\nModel settings:\n")
  print(object$set)

  cat("\nRuns:\n")
  id.failed      <- which(object$finished == 0)
  if (object$set['nsim'] == 1){
    if (length(id.failed) == 0) {
      cat("The simulation finished.\n") 
    } else {
      cat("The simulation failed.\n")
    }
  } else {
    if (sum(object$finished) == length(object$finished)) {
      cat("All simulations finished.\n")
    } else if (length(id.failed) == 1) {
      cat(paste("Simulation number ", id.failed, " failed.\n", sep = ""))
    } else {
      cat(paste("Simulations number ", 
          paste(id.failed[-length(id.failed)], collapse = ", "),
          " and ", id.failed[length(id.failed)], " failed.\n", sep = ""))
    }
  }
	
  cat("\nJumps and priors:\n")
  print(object$JumpP)
  
  if (length(extraArgs) > 0) {
    if (!is.element('digits', names(extraArgs))){
      digits      <- 4
    } else {
      digits      <- extraArgs$digits
    }
  } else {
    digits        <- 4
  }
  cat("\nMean Kullback-Liebler\ndiscrepancy calibration (KLDC):\n")
  if (!is.null(object$K)){
    mean.q        <- (object$K$q12 + object$K$q21) / 2
    print.default(mean.q, digits = digits)
  } else {
    if (object$set['nsim'] == 1) {
      cat("KLDC was not calculated due to insufficient number",
          " of simulations to estimate convergence.\n")
    } else {
      cat("KLDC was not calculated due to lack of convergence.\n")
    }
  }
  

  cat("\nCoefficients:\n")
  print.default(object$coefficients, digits = digits)

  cat("\nConvergence:\n")
  if (is.null(object$modSel)){
    if (object$set['nsim'] == 1) {
      cat("\nConvergence calculations require more than one run.",
          "\nTo estimate potential scale reduction run at least two simulations.\n")
    } else {
      cat("\nWarning: Convergence not reached for some parameters",
          " (i.e. 'PotScaleReduc' values larger than 1.1).",
          "\nThese estimates should not be used for inference.\n")
    }
  } else {
    cat("Appropriate convergence reached for all parameters.\n")
  } 
  cat("\nDIC:\n")
  if (!is.null(object$modSel)){
    cat(object$modSel["DIC"])
  } else {
    if (object$set['nsim'] == 1) {
      cat("DIC was not calculated due to insufficient number",
          " of simulations to estimate convergence.\n")
    } else {
      cat("DIC was not calculated due to lack of convergence.\n")
    }
  }
}

