# Supplementary materials for "Sampling-Based Approaches to Maximum 
# Likelihood Estimation for Latent Variable Models":
# Code for the pump example (2 parameters)
#
# Description:
# The following code contains the numerical experiments for the pump model:
# fixed step size, small fixed step size, adadelta, adam, newton-raphson,
# and 1-D sampling. MCEM (from the R package NIMBLE) is used as a benchmark.
#
##########################################################################

# Fix the seed for reproducibility.
set.seed(21745)

# Load the packages ---------------------------------------
library("MCMCmaxlik")

# Model specification -------------------------------------

pumpCode <- nimbleCode({
  for (i in 1:N) {
    theta[i] ~ dgamma(alpha, beta)
    lambda[i] <- theta[i] * t[i]
    x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1, 1.0)
})
pumpConsts <- list(N=10,
                   t=c(94.3, 15.7, 62.9, 126, 5.24,
                         31.4, 1.05, 1.05, 2.1, 10.5))
pumpData <- list(x=c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
pumpInits <- list(alpha=1, beta=1,
                  theta=rep(0.1, pumpConsts$N))

# Build the model.
pump <- nimbleModel(code=pumpCode, name='pump', constants=pumpConsts,
                    data=pumpData, inits=pumpInits)

paramNodesPump <- pump$getNodeNames(topOnly=T)
Cpump <- compileNimble(pump)


# Testing Algorithms ------------------------------------------------------

# 1. Fixed step size ----------------------------------------
# Compile the necessary functions.
compiledFunsPump <- buildMCMCmaxlik(pump, paramNodesPump)

init <- c(10, 10)
boundary <- list(c(0.05, 10000), c(0.05, 10000))
numMCMCSamples <- 300

setwd("~/Desktop/MCMCmaxlik_results_paper")
ptm <- proc.time()
resultsPumpFixed <- computeMLE(pump, paramNodesPump,
                               method="fixed", paramInit=init,
                               compiledFuns=compiledFunsPump,
                               stepsize=0.05,
                               numMCMCSamples=numMCMCSamples,
                               maxIter=500,
                               boundary=boundary, tol=0)
timePumpFixed <- proc.time() - ptm
timePumpFixed ##  4.475 
mean(resultsPumpFixed$param[401:500,1], trim=.2) ## 0.8221844
mean(resultsPumpFixed$param[401:500,2], trim=.2) ## 1.260489
save(resultsPumpFixed, file="pumpFixed.RData")

# 2. Small fixed step size ----------------------------------------
ptm <- proc.time()
resultsPumpSmallFixed <- computeMLE(pump, paramNodesPump,
                                    method="fixed", paramInit=init,
                                    compiledFuns=compiledFunsPump,
                                    stepsize=0.005,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=5000,
                                    boundary=boundary, tol=0)
timePumpSmallFixed <- proc.time() - ptm
timePumpSmallFixed ## 46.134 
mean(resultsPumpSmallFixed$param[4901:5000,1], trim=.2) ## 0.8230727
mean(resultsPumpSmallFixed$param[4901:5000,2], trim=.2) ## 1.262879
save(resultsPumpSmallFixed, file="pumpSmallFixed.RData")

# 3. Adadelta ----------------------------------------
ptm <- proc.time()
resultsPumpAdadelta <- computeMLE(pump, paramNodesPump,
                                  method="adadelta", paramInit=init,
                                  compiledFuns=compiledFunsPump,
                                  numMCMCSamples=numMCMCSamples,
                                  maxIter=300,
                                  boundary=boundary)
timePumpAdadelta <- proc.time() - ptm
timePumpAdadelta ## 3.023 
mean(resultsPumpAdadelta$param[201:300,1], trim=.2) ## 1.054878
mean(resultsPumpAdadelta$param[201:300,2], trim=.2) ## 1.811111
save(resultsPumpAdadelta, file="pumpAdadelta.RData")

# 4. Adam ----------------------------------------
ptm <- proc.time()
resultsPumpAdam <- computeMLE(pump, paramNodesPump,
                              method="adam", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              stepsize=0.3,
                              eps=1e-4,
                              maxIter=300,
                              boundary=boundary)
timePumpAdam <- proc.time() - ptm
timePumpAdam ## 2.975 
mean(resultsPumpAdam$param[201:300,1], trim=.2) ## 0.8222032
mean(resultsPumpAdam$param[201:300,2], trim=.2) ## 1.260498
save(resultsPumpAdam, file="pumpAdam.RData")

# 5. Newton-Raphson ----------------------------------------
ptm <- proc.time()
resultsPumpNR <- computeMLE(pump, paramNodes=paramNodesPump,
                            method="NR", paramInit=c(10,10),
                            compiledFuns=compiledFunsPump,
                            numMCMCSamples=numMCMCSamples,
                            tol=1e-20,
                            maxIter=300,
                            boundary=boundary)
timePumpNR <- proc.time() - ptm
timePumpNR ## 3.201 
mean(resultsPumpNR$param[201:300,1], trim=.2) ## 7256.679
mean(resultsPumpNR$param[201:300,2], trim=.2) ## 8739.564
save(resultsPumpNR, file="pumpNR.RData")

# 6. 1-D sampling ----------------------------------------
ptm <- proc.time()
resultsPump1D <- computeMLE(pump, paramNodesPump,
                            method="ga1D", paramInit=init,
                            compiledFuns=compiledFunsPump,
                            numMCMCSamples=300, numMCMCSamples1D=300, 
                            maxIter=300)
timePump1D <- proc.time() - ptm
timePump1D ## 10.838 
mean(resultsPump1D$param[201:300,1], trim=.2) ## 0.8666035
mean(resultsPump1D$param[201:300,2], trim=.2) ## 1.330581
save(resultsPump1D, file="pump1D.RData")

# MCEM -----------------------------
source("MCEM_with_output.R")
pump2 <- pump$newModel()
box <- list(list(c('alpha', 'beta'), c(0, Inf)))
pumpMCEM <- buildMCEM(model=pump2, latentNodes='theta[1:10]',
                      boxConstraints=box, theta0=c(10, 10)) ## new version that tracks

ptm <- proc.time()
resultsPumpMCEM <- pumpMCEM() ## 0.8275096    1.276000
timePumpMCEM <- proc.time() - ptm
timePumpMCEM ## 37.616 

iterTimePumpMCEM <- c(0.872, 3.440, 5.129,
                      6.862, 8.522, 12.686,
                      18.877, 37.616)

## 8 iterations

# Save outputs ------------------------
save(resultsPumpMCEM, file="MCEM_pump.RData")


timeInfo <- cbind(c(timePumpFixed, timePumpSmallFixed,
                    timePumpAdadelta, timePumpAdam,
                    timePumpNR, timePump1D, timePumpMCEM),
                  c("fixed", "smallFixed", "adadelta", "adam",
                    "NR", "1D", "mcem"))

write.csv(timeInfo, "timePump.csv", row.names=F)
