print.basta <-
function(x, ...){
  cat("\nCall:\n")
  cat(paste("Model             \t\t: ", x$ModelSpecs[1], "\n", sep = ""))
  cat(paste("Shape             \t\t: ", x$ModelSpecs[2], "\n", sep = ""))
  cat(paste("Covars. structure \t\t: ", x$ModelSpecs[3], "\n", sep = ""))
  cat(paste("Cat. covars.      \t\t: ", x$ModelSpecs[4], "\n", sep = ""))
  cat(paste("Cont. covars.     \t\t: ", x$ModelSpecs[5], "\n", collapse = ""))

  cat("\nRuns:\n")
  id.failed      <- which(x$finished == 0)
  if (x$set['nsim'] == 1){
    if (length(id.failed) == 0) {
      cat("The simulation finished.\n") 
    } else {
      cat("The simulation failed.\n")
    }
  } else {
    if (sum(x$finished) == length(x$finished)) {
      cat("All simulations finished.\n")
    } else if (length(id.failed) == 1) {
      cat(paste("Simulation number ", id.failed, " failed.\n", sep = ""))
    } else {
      cat(paste("Simulations number ", 
          paste(id.failed[-length(id.failed)], collapse = ", "),
          " and ", id.failed[length(id.failed)], " failed.\n", sep = ""))
    }
  }
  cat("\nCoefficients:\n")
  print.default(x$coefficients, ...)
  if (is.null(x$modSel)){
    if (x$set['nsim'] == 1) {
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
}

