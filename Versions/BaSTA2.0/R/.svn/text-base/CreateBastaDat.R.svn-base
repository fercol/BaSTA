# TODO: Add comment
# 
# Author: fernando
###############################################################################

CreateBastaDat <- function(object, studyStart, studyEnd, 
    dataType = "capRecap", capHistMat = NULL, fixedCovs = NULL, 
    timeCovs = NULL, firstCohort = studyStart, Id = "col") {
  if (!(dataType %in% c("capRecap", "capRecov", "census", "fullCens"))) {
    stop("Wrong 'dataType' argument. Options are: 'capRecap', 'capRecov',\n", 
        "'census' or 'fullCens'.", call. = FALSE)
  }
  if (!(Id %in% c("col", "name", "nrow"))) {
    stop("Id argument needs to be either 'name' (i.e. individual IDs\n",
        "in row names), 'col' (i.e. IDs in first column), or\n",
        "'nrow' (i.e. number of rows).", call. = FALSE)
  } else if (Id == 'nrow') {
    idObj <- 1:nrow(object)
    birthDeath <- object[, 1:2]
  } else if (Id == "col") {
    idObj <- object[, 1]
    if (length(unique(idObj)) < length(idObj)) {
      stop(paste("First column in 'object' does not seem to contain ", 
              "unique IDs.\nVerify or change Id argument to 'nrow' or 'name'.", 
              sep = ""), call. = FALSE)
    }
    birthDeath <- object[, 2:3]
  } else {
    idObj <- rownames(object)
    if (is.null(idObj)) {
      stop("No row names in object. Use Id = 'nrow'.", call. = FALSE)
    }
    birthDeath <- object[, 1:2]
  }
  if (!is.null(capHistMat)) {
    capHist <- .MatchDataMats(object, idObj, Id, capHistMat, "capHistMat")
  } else {
    idCol <- 1:ifelse(Id %in% c("nrow", "names"), 2, 3)
    if (ncol(object) > max(idCol)) {
      capHist <- object[, -idCol]
    } else {
      capHist <- NULL
    }
  }
  if (!is.null(fixedCovs)) {
    fixedCovs <- .MatchDataMats(object, idObj, Id, fixedCovs, "fixedCovs")
    nameFcovs <- colnames(fixedCovs)
  } else {
    nameFcovs <- NULL
  }
  if (!is.null(timeCovs)) {
    newTimeCovs <- list()
    if (is.array(timeCovs)) {
#      timeCovs <- .MatchDataMats(object, idObj, Id, timeCovs, "timeCovs")
      if (length(dim(timeCovs)) == 3) {
        for (tC in 1:dim(timeCovs)[3]) {
          newTimeCovs[[dimnames(timeCovs)[[3]][tC]]] <- .MatchDataMats(object, 
              idObj, Id, timeCovs[,, tC], "timeCovs")
        }
      } else {
        newTimeCovs[["tempCov"]] <- .MatchDataMats(object, idObj, Id, 
            timeCovs, "timeCovs")
      }
      timeCovs <- newTimeCovs
    } else {
      for (tC in 1:length(timeCovs)) {
        dtMat <- as.matrix(timeCovs[[tC]])
        newTimeCovs[[names(timeCovs)[tC]]] <- .MatchDataMats(object, idObj, 
            Id, dtMat, paste(names(timeCovs)[tC], "in timeCovs"))
      }
      timeCovs <- newTimeCovs
    }
    nameTcovs <- names(timeCovs)
  } else {
    nameTcovs <- NULL
  }
  if (is.null(nameFcovs) & is.null(nameTcovs)) {
    covClass <- "noCovs"
  } else if (is.null(nameFcovs) & !is.null(nameTcovs)) {
    covClass <- "timeCovs"
  } else if (!is.null(nameFcovs) & is.null(nameTcovs)) {
    covClass <- "fixedCovs"
  } else {
    covClass <- "bothCovs"
  }
  dataClass <- dataType
  if (is.null(capHist)) {
    if (dataType == "capRecap") {
      stop("Argument 'dataType' is set to 'capRecap' but no capture history\n",
          "matrix was found in the dataset.", call. = FALSE)
    } else if (dataType == "capRecov") {
      stop("Argument 'dataType' is set to 'capRecov' but no tag-recovery\n",
          "matrix was found in the dataset.", call. = FALSE)
    } else if (dataType == "census") {
      stop("Argument 'dataType' is set to 'census' but no first-last\n",
          "detection times matrix was found in the dataset.", call. = FALSE)
    }
  }
  bastaData <- list(birthDeath = birthDeath, capHistMat = capHist, 
      fixedCovs = fixedCovs, timeCovs = timeCovs, nameFCovs = nameFcovs, 
      nameTCovs = nameTcovs, studyStart = studyStart, studyEnd = studyEnd,
      firstCohort = firstCohort, N = nrow(birthDeath))
  class(bastaData) <- c(dataClass, covClass)
  return(bastaData)
}


.MatchDataMats <- function(object, idObj, Id, dataMat, nameDataMat) {
  if (Id == "nrow") {
    idDat <- 1:nrow(dataMat)
  } else if (Id == 'col') {
    idDat <- dataMat[, 1]
    dataMat <- dataMat[, -1]
    if (length(unique(idDat)) < length(idDat)) {
      stop(paste("First column in '", nameDataMat, 
              "' does not seem to contain ", 
              "unique IDs.\nVerify or change Id argument to 'nrow' or 'name'.", 
              sep = ""), call. = FALSE)
    }
  } else {
    idDat <- rownames(dataMat)
    if (is.null(idDat)) {
      stop(paste("No row names in '", nameDataMat, 
              "'. Use Id = 'nrow'.", sep = ""), call. = FALSE)
    }
  }
  idMatch <- match(idObj, idDat)
  if (length(which(is.na(idMatch))) == 0) {
    dataMat <- dataMat[idMatch, ]
  } else {
    if (Id == 'nrow') {
      stop(paste("Number of rows in 'object' and '", nameDataMat, 
              "' do not match.", sep = ""), call. = FALSE)      
    } else {
      stop(paste("Missmatch in records between 'object' and '", nameDataMat, 
              "'.\nPlease check that individual IDs in both tables match.", 
              sep = ""), call. = FALSE)
    }
  }
  return(dataMat)
}



