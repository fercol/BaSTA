# TODO: Add comment
# 
# Author: fernando
###############################################################################


.PrepMortDat <- 
    function(object) UseMethod(".PrepMortDat")

.PrepMortDat.fullCens <- function(object) {
  # 1. Extract times of birth and death:
  bi <- object$birthDeath[, 1]
  di <- object$birthDeath[, 2]
  classAge <- ifelse(any(di - bi != round(di - bi)), 
      "contAge", "discAge")
  
  # 2. Construct data list:
  mortDat <- list(birth = bi, death = di)
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
  mortDat <- list(birth = bi, death = di, firstObs = firstObs, 
      lastObs = lastObs)
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
  mortDat <- list(birth = bi, death = di, firstObs = firstObs, 
      lastObs = lastObs)
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
  nTimes <- ncol(object$capHistMat)
  times <- object$studyStart + 1:nTimes - 1
  ytemp <- t(t(object$capHistMat) * times)
  lastObs <- c(apply(ytemp, 1, max))
  ytemp[ytemp == 0] <- max(times) * 2
  firstObs <- c(apply(ytemp, 1, min))
  firstObs[firstObs == max(times) * 2]  <- 0
  oi <- object$capHistMat %*% rep(1, nTimes)
  
  # 3. Define study duration:
  Dx <- (times[2] - times[1])
  Tm <- matrix(times, object$N, nTimes, byrow=TRUE)
  mortDat <- list(birth = bi, death = di, firstObs = firstObs, 
      lastObs = lastObs, Y = object$capHistMat, oi = oi, Dx = Dx, 
      aliveMat = Tm)
  class(mortDat) <- c(classAge, "capRecap") 
  return(mortDat)
}

