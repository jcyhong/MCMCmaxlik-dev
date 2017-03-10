# Supplementary materials for "Sampling-Based Approaches to Maximum 
# Likelihood Estimation for Latent Variable Models":
# Code for the logistic regression example (3 parameters)
#
# Description:
# The following code contains the numerical experiments for the logistic 
# regression model: fixed step size, small fixed step size, adadelta, adam, 
# Newton-Raphson, and 1-D sampling. MCEM (from the R package NIMBLE) is used
# as a benchmark.
#
##########################################################################

# Fix the seed for reproducibility.
set.seed(21745)

# Load the packages ---------------------------------------
library("MCMCmaxlik")

# Model specification -------------------------------------
code <- nimbleCode({
  beta0 ~ dnorm(0, sd=10000)
  beta1 ~ dnorm(0, sd=10000)
  sigma_RE ~ dunif(0, 1000)
  for(i in 1:N) {
    beta2[i] ~ dnorm(0, sd=sigma_RE)
    logit(p[i]) <- beta0 + beta1 * x[i] + beta2[i]
    r[i] ~ dbin(p[i], n[i])
  }
})

## constants, data, and initial values
constants <- list(N=10)

data <- list(
  r=c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46),
  n=c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79),
  x=c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
)

inits <- list(beta0=0, beta1=0, sigma_RE=1)

logreg <- nimbleModel(code=code, constants=constants, 
                      data=data, inits=inits, check=FALSE)

paramNodesLogreg <- logreg$getNodeNames(topOnly=T)
Clogreg <- compileNimble(logreg)


# Testing Algorithms ------------------------------------------------------
# Compile the necessary functions.
compiledFunsLogreg <- buildMCMCmaxlik(logreg, paramNodesLogreg)

init <- c(0, 0, 1)
boundary <- list(c(-40, 40), c(-40, 40), c(0.05, 10))
numMCMCSamples <- 300

# 1. Fixed step size ----------------------------------------
ptm <- proc.time()
resultsLogregFixed <- computeMLE(logreg, paramNodesLogreg,
                                 method="fixed", paramInit=init,
                                 stepsize=0.05,
                                 compiledFuns=compiledFunsLogreg,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=300,
                                 boundary=boundary)
timeLogRegFixed <- proc.time() - ptm ##  2.508 
mean(resultsLogregFixed$param[201:300,1], trim=.2)  ## -0.7522742
mean(resultsLogregFixed$param[201:300,2], trim=.2)  ## 1.25159
mean(resultsLogregFixed$param[201:300,3], trim=.2)  ## 0.8827718
save(resultsLogregFixed, file="logRegFixed.RData")

ptm <- proc.time()

# 2. Small fixed step size ----------------------------------------
resultsLogregSmallFixed <- computeMLE(logreg, paramNodesLogreg,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsLogreg,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary)
timeLogRegSmallFixed <- proc.time() - ptm ## 2.530 
mean(resultsLogregSmallFixed$param[201:300,1], trim=.2)  ## -0.5485272
mean(resultsLogregSmallFixed$param[201:300,2], trim=.2)  ##  1.310605
mean(resultsLogregSmallFixed$param[201:300,3], trim=.2)  ## 0.251749
save(resultsLogregSmallFixed, file="logRegSmallFixed.RData")

# 3. Adadelta ----------------------------------------
ptm <- proc.time()
resultsLogregAdadelta <- computeMLE(logreg, paramNodesLogreg,
                                    method="adadelta", paramInit=init,
                                    compiledFuns=compiledFunsLogreg,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary)
timeLogRegAdadelta <- proc.time() - ptm ## 2.758 
mean(resultsLogregAdadelta$param[201:300,1], trim=.2)  ## -0.555296
mean(resultsLogregAdadelta$param[201:300,2], trim=.2)  ## 1.319334
mean(resultsLogregAdadelta$param[201:300,3], trim=.2)  ## 0.3160231
save(resultsLogregAdadelta, file="logRegAdadelta.RData")

# 4. Adam ----------------------------------------
ptm <- proc.time()
resultsLogregAdam <- computeMLE(logreg, paramNodesLogreg,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                stepsize=0.2,
                                eps=1e-4,
                                maxIter=300,
                                boundary=boundary)
timeLogRegAdam <- proc.time() - ptm ## 2.789 
mean(resultsLogregAdam$param[201:300,1], trim=.2) ## -0.5482783
mean(resultsLogregAdam$param[201:300,2], trim=.2) ##  1.310051
mean(resultsLogregAdam$param[201:300,3], trim=.2) ## 0.2494285
save(resultsLogregAdam, file="logRegAdam.RData")

# 5. Newton-Raphson ----------------------------------------
ptm <- proc.time()
resultsLogregNR <- computeMLE(logreg, paramNodesLogreg,
                              method="NR", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary)
timeLogRegNR <- proc.time() - ptm ##  5.790 
mean(resultsLogregNR$param[201:300,1], trim=.2) ## 0.6831159
mean(resultsLogregNR$param[201:300,2], trim=.2) ## -1.625356
mean(resultsLogregNR$param[201:300,3], trim=.2) ## 10 (boundary)
## stuck at these values from iteration 74 on
save(resultsLogregNR, file="logRegNR.RData")

# 6. 1-D sampling ----------------------------------------
ptm <- proc.time()
resultsLogreg1D <- computeMLE(logreg, paramNodesLogreg,
                              method="ga1D", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=300, numMCMCSamples1D=300, 
                              maxIter=300)
timeLogReg1D <- proc.time() - ptm ## 10.704 
mean(resultsLogreg1D$param[201:300,1], trim=.2)  ## -0.5677065
mean(resultsLogreg1D$param[201:300,2], trim=.2)  ## 1.327946
mean(resultsLogreg1D$param[201:300,3], trim=.2)  ## 0.2804238
save(resultsLogreg1D, file="logReg1D.RData")

# MCEM -----------------------------
source("MCEM_with_output.R")
logreg2 <- logreg$newModel()
latentNodes <- logreg2$getNodeNames(latentOnly=TRUE, stochOnly=TRUE)

modMCEM <- buildMCEM(model=logreg2, latentNodes=latentNodes, theta0=c(0,0,1))

ptm <- proc.time()
modMLE <- modMCEM()
timeLogRegMCEM <- proc.time() - ptm 
timeLogRegMCEM ## 252.181 
modMLE  ## -0.54690186   1.31137460       0.2576080 
save(modMLE, file="MCEM_logReg.RData")

## 30 iter

iterTimeLogRegMCEM <- c(1.034, 3.581, 5.141, 9.036, 11.374,
                        13.631, 15.906, 18.187, 20.537, 22.762,
                        25.072, 27.397, 29.691, 32.084, 34.349,
                        36.607, 38.930, 41.257, 43.555, 45.975,
                        48.269, 58.947, 64.053, 69.666, 74.706,
                        79.703, 104.517, 115.557, 213.341,
                        252.181)

timeInfo=cbind(c(timeLogRegFixed, timeLogRegSmallFixed, timeLogRegAdadelta, 
                 timeLogRegAdam, timeLogRegNR, timeLogReg1D, timeLogRegMCEM),
               c("fixed", "smallFixed", "adadelta", "adam", "NR", "1D", 
                 "mcem"))


write.csv(timeInfo,"timeLogReg.csv",row.names=F)


