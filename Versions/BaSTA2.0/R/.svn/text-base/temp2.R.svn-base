source("Branches/BaSTA2.0/R/basta.R")
source("Branches/BaSTA2.0/R/CreateBastaDat.R")
#source("Branches/BaSTA2.0/R/DataCheck.R")

niter <- 1000
burnin <- 100
thinning <- 10
tCovStart <- 1
stStart <- 1
model <- "GO"
shape <- "simple"
covarsStruct <- "all.in.mort"
cohortCov <- "cohort"
default <- TRUE
minAge <- 0
myobj <- cbind(1:5, c(1,2,4,2,3), c(5,3,4,5,4), 
    matrix(rbinom(25, 1, 0.4), 5, 5))
n <- nrow(myobj)
colnames(myobj) <- c("id", "b", "d", 
    paste("Y.", 1:5, sep = ""))
tCovs <- array(rnorm(5*10*2, 0, 1), dim = c(5, 10, 2), 
    dimnames = list(NULL, paste("Y.", 1:10, sep = ""), 
        c('w1', 'w2')))
fCovs <- data.frame(sex = c("f", "f", "f", "m", "m"), 
    loc = c("1", "2", "3", "1", "3"), 
    r1 = rnorm(5, 5, 1), r2 = rnorm(5, -1, 1), cohort = rep(1, 5))
tCovs <- list(w1 = matrix(c(1:5, rnorm(5*10, 0, 1)), 5, 11, 
        dimnames = list(NULL, c("ID", paste("Y.", 1:10, 
                    sep = "")))), 
    w2 = matrix(c(1:5, rnorm(5*10, -2, 1)), 5, 11, 
        dimnames = list(NULL, c("ID", paste("Y.", 1:10, sep = "")))))

myform <- ~ sex + r1 + w1 + sex:w1 + cohort

# Test CreateBastaDat:
bastaDat <- CreateBastaDat(object = myobj[, 2:3], capHistMat = NULL, 
    fixedCovs = fCovs, timeCovs = tCovs, dataType = "fullCens",
    studyStart = 1, studyEnd = 5, Id = "nrow")

bastaDat <- CreateBastaDat(object = myobj[, 2:3], capHistMat = NULL, 
    fixedCovs = NULL, timeCovs = NULL, dataType = "fullCens",
    studyStart = 1, studyEnd = 5, Id = "nrow")

mortDat <- .PrepMortDat(bastaDat)
covList <- .CreateCovList(bastaDat, form = myform)
timesg <- .SetIniTimeDat(covList, mortDat, minAge)
timesn <- .SampleTimeDat(timesg, mortDat, minAge, covList)
parsMort <- .DefineParams(covList, model, shape, covarsStruct, default)
parsg <- .DefineParsG(parsMort)
parCovsg <- .CalcParams(parsg, covList)
CalcSimpleSurv <- .DefineSimpleSurv(model)
CalcSurv <- .DefineSurv(shape, CalcSimpleSurv)
CalcSimpleMort <- .DefineSimpleMort(model)
CalcMort <- .DefineMort(shape, CalcSimpleMort)

priorAgeDistr <- .SetPriorAgeDist(parsMort)


newDat <- CreateBastaDat(object = myobj[, -1], capHistMat = NULL, 
    fixedCovs = fCovs2, studyStart = 1, Id = "nrow")

rownames(myobj) <- 1:nrow(myobj)
rownames(fCovs) <- 1:nrow(fCovs)
rownames(tCovs) <- 1:nrow(tCovs)

newDat <- CreateBastaDat(object = myobj[, -1], capHistMat = NULL, 
    fixedCovs = fCovs, timeCovs = tCovs, studyStart = 1, studyEnd = 5, 
    Id = "name")

nfCovs <- fCovs[sample(1:nrow(fCovs), nrow(fCovs), replace = FALSE), ]

newDat <- CreateBastaDat(object = myobj[, -1], capHistMat = NULL, 
    fixedCovs = nfCovs, timeCovs = tCovs, Id = "name", 
    studyStart = 1, studyEnd = 5)

nfCovs <- fCovs[-1, ]
newDat <- CreateBastaDat(object = myobj[, -1], capHistMat = NULL, 
    fixedCovs = nfCovs, timeCovs = tCovs, Id = "name", 
    studyStart = 1, studyEnd = 5)

ntCovs <- tCovs[-1,, ]
newDat <- CreateBastaDat(object = myobj[, -1], capHistMat = NULL, 
    fixedCovs = fCovs, timeCovs = ntCovs, Id = "name", 
    studyStart = 1, studyEnd = 5)


capHis <- myobj[, -c(1:3)]
newDat <- CreateBastaDat(object = myobj[, 2:3], capHistMat = capHis, 
    fixedCovs = fCovs, timeCovs = tCovs, Id = "name",
    studyStart = 1, studyEnd = 5)

rownames(capHis) <- NULL
newDat <- CreateBastaDat(object = myobj[, 2:3], capHistMat = capHis, 
    fixedCovs = fCovs, timeCovs = tCovs, Id = "name",
    studyStart = 1, studyEnd = 5)

newDat <- CreateBastaDat(object = myobj[, 2:3], capHistMat = capHis, 
    fixedCovs = fCovs, timeCovs = tCovs, Id = "nrow",
    studyStart = 1, studyEnd = 5)

newDat <- CreateBastaDat(object = myobj[, 2:3], capHistMat = capHis, 
    fixedCovs = fCovs, timeCovs = tCovs[,, 1], Id = "nrow",
    studyStart = 1, studyEnd = 5)


fCovs2 <- cbind(1:nrow(fCovs), fCovs)
tCovs2 <- array(c(1:5, rnorm(5*10, 0, 1), 1:5, rnorm(5*10, 0, 1)), 
    dim = c(5, 11, 2), dimnames = list(NULL, c("ID", paste("Y.", 1:10, 
                sep = "")), 
        c('w1', 'w2')))
tCovs2 <- list(w1 = matrix(c(1:5, rnorm(5*10, 0, 1)), 5, 11, 
        dimnames = list(NULL, c("ID", paste("Y.", 1:10, 
                    sep = "")))), 
    w2 = matrix(c(1:5, rnorm(5*10, -2, 1)), 5, 11, 
    dimnames = list(NULL, c("ID", paste("Y.", 1:10, sep = "")))))

newDat2 <- CreateBastaDat(object = myobj, capHistMat = NULL, 
    fixedCovs = fCovs2, timeCovs = tCovs2, Id = "col")

newDat2 <- CreateBastaDat(object = myobj[, 1:3], Id = "col")
mortDat2 <- PrepMortDat(newDat2)
newDat2 <- CreateBastaDat(object = myobj[, c(1:3, 2:3)], Id = "col")
mortDat2 <- PrepMortDat(newDat2)
newDat <- CreateBastaDat(object = myobj[, 1:5], capHistMat = NULL, 
    fixedCovs = fCovs2, timeCovs = tCovs2, Id = "col")

# Test .prepMortDat():
bastaDat <- CreateBastaDat(object = myobj[, 1:3], Id = "col", 
    studyStart = 1, studyEnd = 5, dataType = "fullCens")
mortDat <- .PrepMortDat(bastaDat)
iniTimes <- .SetIniTimes(mortDat, minAge = 0)
newTimes <- .SampleAges(iniTimes, mortDat)

bastaDat <- CreateBastaDat(object = myobj[, c(1:3, 2:3)], Id = "col",
    studyStart = 1, studyEnd = 5, dataType = "census")
mortDat <- .PrepMortDat(bastaDat)
iniTimes <- .SetIniTimes(mortDat, minAge = 1)
newTimes <- .SampleAges(iniTimes, mortDat)


bastaDat <- CreateBastaDat(object = myobj[, 1:5], capHistMat = NULL, 
    fixedCovs = fCovs2, timeCovs = tCovs2, Id = "col")
mortDat <- .PrepMortDat(bastaDat)
iniTimes <- .SetIniTimes(mortDat, minAge = 0)

bastaDat <- CreateBastaDat(object = myobj, capHistMat = NULL, 
    fixedCovs = fCovs2, timeCovs = tCovs2, studyStart = 1, studyEnd = 5,
    Id = "col")
mortDat <- .PrepMortDat(bastaDat)
iniTimes <- .SetIniTimes(mortDat, minAge = 0)
newTimes <- .SampleAges(iniTimes, mortDat)
newFullTimes <- .SplitByMinAge(newTimes, iniTimes, mortDat, minAge = 0)
covList <- .CreateCovList(bastaDat)

covList <- .CreateCovList(bastaDat, ~ sex + loc + w1 + w2)


bastaDat <- CreateBastaDat(object = myobj, capHistMat = NULL, 
    fixedCovs = NULL, timeCovs = NULL, studyStart = 1, studyEnd = 5, 
    Id = "col")
mortDat <- .PrepMortDat(bastaDat)

covList <- .CreateCovList(bastaDat)
covList <- .CreateCovList(bastaDat, form = ~ sex + loc)

params <- .DefineParams(covList, covarsStruct = "prop.haz", 
    default = TRUE)
parsg <- .DefineParsG(params)
CalcSimpleSurv <- .DefineSimpleSurv(model = "G0")
CalcSurv <- .DefineSurv(shape = "simple", CalcSimpleSurv)
CalcSimpleMort <- .DefineSimpleMort(model = "G0")
CalcMort <- .DefineMort(shape = "simple", CalcSimpleMort)

parCovs <- .CalcParams(parsg, covList)
priorAgeDistr <- .SetPriorAgeDist(params)

covList <- .CreateCovList(bastaDat, ~ sex + loc + w1 + w2)
parCovs <- .CalcParams(parsg, covList)


