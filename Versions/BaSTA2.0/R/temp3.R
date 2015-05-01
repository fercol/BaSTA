load("/Users/fernando/FERNANDO/PROJECTS/4.PACKAGES/workspace/basta/pkg/data/sim1.RData")
model <- "GO"
n <- nrow(sim1)
totYears <- 1:100
nT <- length(totYears)
shape <- "simple"
covarsStruct <- "all.in.mort"
default <- TRUE
minAge <- 1
tCovsSim <- matrix(rnorm(nT * n, 0, 1), n, nT)
colnames(tCovsSim) <- totYears
bastaDat <- CreateBastaDat(object = sim1[, -c(1, 24:26)], capHistMat = NULL, 
    fixedCovs = sim1[, 24:26], timeCovs = tCovsSim, 
    studyStart = 51, studyEnd = 70, firstCohort = 1, Id = "nrow")

mortDat <- .PrepMortDat(bastaDat)
covList <- .CreateCovList(bastaDat, form = ~ f + m + weight + tempCov - 1)
parsMort <- .DefineParams(covList, model, shape, covarsStruct, default)
parsg <- .DefineParsG(parsMort)
parCovsg <- .CalcParams(parsg, covList)
CalcSimpleSurv <- .DefineSimpleSurv(model)
CalcSurv <- .DefineSurv(shape, CalcSimpleSurv)
CalcSimpleMort <- .DefineSimpleMort(model)
CalcMort <- .DefineMort(shape, CalcSimpleMort)

priorAgeDistr <- .SetPriorAgeDist(parsMort)
timesg <- .SetIniTimeDat(covList, mortDat, minAge)
timesn <- .SampleTimeDat(timesg, mortDat, minAge, covList)
newTs <- .ProposeAges(timesg, mortDat)
splitTs <- .SplitByMinAge(newTs, mortDat, minAge)
updtTs <- .ProposeTimes(timesg, mortDat, minAge)



# Rook example:
dat <- read.csv("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/BASTA/ANALYSIS/DATA/ROOKS/finalRooks.csv")
model <- "GO"
shape <- "simple"
covarsStruct <- "all.in.mort"
default <- TRUE
minAge <- 1
bastaDat <- CreateBastaDat(object = dat[, -c(1, 14:15)], capHistMat = NULL,
    fixedCovs = dat[, 14:15], timeCovs = NULL, studyStart = 1964, 
    studyEnd = 1973, firstCohort = 1940, Id = "nrow")

mortDat <- .PrepMortDat(bastaDat)
covList <- .CreateCovList(bastaDat, form = ~ X1 + X2 - 1)
#covList <- .CreateCovList(bastaDat, form = NULL)
parsMort <- .DefineParams(covList, model, shape, covarsStruct, default)
parsg <- .DefineParsG(parsMort)
parCovsg <- .CalcParams(parsg, covList)
CalcSimpleSurv <- .DefineSimpleSurv(model)
CalcSurv <- .DefineSurv(shape, CalcSimpleSurv)
CalcSimpleMort <- .DefineSimpleMort(model)
CalcMort <- .DefineMort(shape, CalcSimpleMort)

priorAgeDistr <- .SetPriorAgeDist(parsMort)
timesg <- .SetIniTimeDat(covList, mortDat, minAge)
timesn <- .SampleTimeDat(timesg, mortDat, minAge, covList)

