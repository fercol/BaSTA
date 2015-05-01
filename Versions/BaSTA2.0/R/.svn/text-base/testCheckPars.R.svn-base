myobj <- cbind(1:5, c(1,2,4,2,3), c(5,3,4,6,4), 
    matrix(rbinom(25, 1, 0.4), 5, 5))
colnames(myobj) <- c("id", "b", "d", 
    paste("Y.", 1:5, sep = ""))

# TEST CheckParsMort:
# fixed inMort
fCovs <- matrix(c(1,1,1,0,0,0,0,0,1,1), 5, 2, 
    dimnames = list(NULL, c('f', 'm')))
tCovs <- array(rnorm(25, 0, 1), dim = c(5, 5, 2), 
    dimnames = list(NULL, colnames(myobj)[-c(1:3)], 
        c('w1', 'w2')))

datBC <- .CreateBastaData(myobj, fCovs, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
imUserPars <- list(start = NULL, jump = NULL, prior = NULL)

imPars <- .CheckParsMort(imUserPars, covarsBC, model = "GO", shape = "simple")
imPars

# fixed propHaz:
rm(imPars)
fCovs <- matrix(c(rnorm(5, 0, 1), rnorm(5, 10, 1)), 5, 2, 
    dimnames = list(NULL, c('r1', 'r2')))

datBC <- .CreateBastaData(myobj, fCovs, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
imUserPars <- list(start = NULL, jump = NULL, prior = NULL)

imPars <- .CheckParsMort(imUserPars, covarsBC, model = "GO", shape = "simple")
imPars

# fixed inMort wrong user pars:
rm(imPars)
fCovs <- matrix(c(1,1,1,0,0,0,0,0,1,1), 5, 2, 
    dimnames = list(NULL, c('f', 'm')))

datBC <- .CreateBastaData(myobj, fCovs, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
imUserPars <- list(start = c(1,1,1), jump = NULL, prior = NULL)

imPars <- .CheckParsMort(imUserPars, covarsBC, model = "GO", shape = "simple")
imPars

# fixed propHaz right user pars:
rm(imPars)
fCovs <- matrix(c(rnorm(5, 0, 1), rnorm(5, 10, 1)), 5, 2, 
    dimnames = list(NULL, c('r1', 'r2')))

datBC <- .CreateBastaData(myobj, fCovs, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
imUserPars <- list(start = c(1,1), jump = NULL, prior = NULL)

imPars <- .CheckParsMort(imUserPars, covarsBC, model = "GO", shape = "simple")
imPars



# =====================================================#
# TEST CheckParsPH:
# a) Both covariates (time and fixed) and fixed in mort:
fCovs <- matrix(c(1,1,1,0,0,0,0,0,1,1), 5, 2, 
    dimnames = list(NULL, c('f', 'm')))
tCovs <- array(rnorm(25, 0, 1), dim = c(5, 5, 2), 
    dimnames = list(NULL, colnames(myobj)[-c(1:3)], 
        c('w1', 'w2')))

datBC <- .CreateBastaData(myobj, fCovs, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
phUserPars <- list(fixed = list(start = NULL, jump = NULL, prior = NULL), 
    time = list(start = NULL, jump = NULL, prior = NULL))

phPars <- .CheckParsPH(covarsBC, phUserPars)
phPars

# b) both with fix fused:
fCovs <- matrix(c(1,1,1,0,0,0,0,0,1,1, rnorm(5, 0, 1), rnorm(5, 10, 1)), 5, 4, 
    dimnames = list(NULL, c('f', 'm', 'r1', 'r2')))

datBC <- .CreateBastaData(myobj, fCovs, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
phUserPars <- list(fixed = list(start = NULL, jump = NULL, prior = NULL), 
    time = list(start = NULL, jump = c(0.1, 0.1, 0.1), prior = NULL))

phPars <- .CheckParsPH(covarsBC, phUserPars)
phPars

# c) both with fix propHaz:
fCovs <- matrix(c(rnorm(5, 0, 1), rnorm(5, 10, 1)), 5, 2, 
    dimnames = list(NULL, c('r1', 'r2')))

datBC <- .CreateBastaData(myobj, fCovs, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
phUserPars <- list(fixed = list(start = NULL, jump = NULL, prior = NULL), 
    time = list(start = NULL, jump = NULL, prior = NULL))

phPars <- .CheckParsPH(covarsBC, phUserPars)
phPars

# d) Only fixed:
datBC <- .CreateBastaData(myobj, fCovs, NULL)
covarsBC <- .PrepCovs(datBC, 'fused')
phUserPars <- list(fixed = list(start = NULL, jump = NULL, prior = NULL), 
    time = list(start = NULL, jump = NULL, prior = NULL))

phPars <- .CheckParsPH(covarsBC, phUserPars)
phPars

# e) Only time:
datBC <- .CreateBastaData(myobj, NULL, tCovs)
covarsBC <- .PrepCovs(datBC, 'fused')
phUserPars <- list(fixed = list(start = NULL, jump = NULL, prior = NULL), 
    time = list(start = NULL, jump = NULL, prior = NULL))

phPars <- .CheckParsPH(covarsBC, phUserPars)
phPars


# f) ony fixed in Mort:
fCovs <- matrix(c(1,1,1,0,0,0,0,0,1,1), 5, 2, 
    dimnames = list(NULL, c('f', 'm')))

datBC <- .CreateBastaData(myobj, fCovs, NULL)
covarsBC <- .PrepCovs(datBC, 'fused')
phUserPars <- list(fixed = list(start = NULL, jump = NULL, prior = NULL), 
    time = list(start = NULL, jump = NULL, prior = NULL))

phPars <- .CheckParsPH(covarsBC, phUserPars)
phPars

