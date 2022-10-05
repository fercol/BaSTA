CensusToCaptHist <-
    function(ID, d, dformat = "%Y", timeInt = "Y") {
  # Check data
  if(!inherits(ID, "character")) {
    ID <- as.character(ID)
  } 
  if (is.numeric(d)) {
    if (length(which(round(d) != d)) > 0) {
      stop("Please provide integer values or Date class ", 
          "values for argument 'd'.", call. = FALSE)
    } else {
      int <- d
    }
  } else if (is.character(d) | inherits(d, "Date")) {
    if (is.character(d)) {
      d <- as.Date(d, format = dformat)
      if (length(which(is.na(d)))) {
        stop("Wrong 'dformat' argument or wrong 'd' values.", 
            call. = FALSE)
      }
    }
    if (timeInt == "Y"){
      int <- as.numeric(format(d, format = "%Y"))
    } else if (timeInt == "M") {
      int <- as.numeric(format(d, format = "%m")) + 
          12 * (as.numeric(format(d, format = "%Y")) - 
            min(as.numeric(format(d, format = "%Y"))))
    } else if (timeInt == "D" | timeInt == "W") {
      jul <- julian(d, origin = min(d)) + 1
      if (timeInt == "W"){
        int <- ceiling(jul / 7)
      } else {
        int <- jul
      }
    }
  } else {
    stop("Wrong class for argument 'd'. Values\n", 
        "need to be of class 'integer', 'character' or 'Date'.", 
        call. = FALSE)
  }  
  
  # Construct capture-recapture matrix:
  dint <- min(int):max(int)
  ndint <- length(dint)
  uniID <- sort(unique(ID))
  n <- length(uniID)
  mat <- matrix(0, n, ndint, dimnames = list(NULL, dint))
  for (i in 1:ndint) {
    idt <- which(d == dint[i])
    idd <- which(uniID %in% ID[idt])
    mat[idd, i] <- 1
  }

  # Create data.frame with ID and Y:
  dmat <- data.frame(ID = uniID, mat)
  
  return(dmat)
}