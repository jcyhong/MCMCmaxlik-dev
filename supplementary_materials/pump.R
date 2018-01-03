# Supplementary materials for "Sampling-Based Approaches to Maximum 
# Likelihood Estimation for Latent Variable Models":
# Code for the pump example (2 parameters)
#
# Description:
# The following code contains the numerical experiments for the pump model:
# fixed step size, small fixed step size, adadelta, adam, Newton-Raphson,
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

ptm <- proc.time()
resultsPumpFixed <- computeMLE(pump, paramNodesPump,
                               method="fixed", paramInit=init,
                               compiledFuns=compiledFunsPump,
                               stepsize=0.05,
                               numMCMCSamples=numMCMCSamples,
                               maxIter=300,
<<<<<<< HEAD
                               blockSize=20,
                               trackEffSizeGrad=F,
                               boundary=boundary)
=======
                               boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
>>>>>>> 07e88e0e8ef50a723b60ca9ff1d27acfa39cc7f4
timePumpFixed <- proc.time() - ptm
timePumpFixed ##  2.946 
resultsPumpFixed$iter ## 300
mean(tail(resultsPumpFixed$param[, 1], 20), trim=.2) ## 0.8635993
mean(tail(resultsPumpFixed$param[, 2], 20), trim=.2) ## 1.369763
save(resultsPumpFixed, file="pumpFixed.RData")

ptm <- proc.time()
resultsPumpFixed2 <- computeMLE(pump, paramNodesPump,
                               method="fixed", paramInit=init,
                               compiledFuns=compiledFunsPump,
                               stepsize=0.05,
                               numMCMCSamples=numMCMCSamples,
                               maxIter=300,
                               boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timePumpFixed2 <- proc.time() - ptm
timePumpFixed2 ##  3.037 
resultsPumpFixed2$iter ## 300
mean(tail(resultsPumpFixed2$param[, 1], 20), trim=.2) ## 0.866495
mean(tail(resultsPumpFixed2$param[, 2], 20), trim=.2) ## 1.38381
save(resultsPumpFixed2, file="pumpFixedCC.RData")

# 2. Small fixed step size ----------------------------------------
ptm <- proc.time()
resultsPumpSmallFixed <- computeMLE(pump, paramNodesPump,
                                    method="fixed", paramInit=init,
                                    compiledFuns=compiledFunsPump,
                                    stepsize=0.005,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=3000,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timePumpSmallFixed <- proc.time() - ptm
timePumpSmallFixed ## 29.164 
resultsPumpSmallFixed$iter ## 3000
mean(tail(resultsPumpSmallFixed$param[, 1], 20), trim=.2) ## 0.8460056
mean(tail(resultsPumpSmallFixed$param[, 2], 20), trim=.2) ## 1.325574
save(resultsPumpSmallFixed, file="pumpSmallFixed.RData")

ptm <- proc.time()
resultsPumpSmallFixed2 <- computeMLE(pump, paramNodesPump,
                                    method="fixed", paramInit=init,
                                    compiledFuns=compiledFunsPump,
                                    stepsize=0.005,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=3000,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timePumpSmallFixed2 <- proc.time() - ptm
timePumpSmallFixed2 ## 30.356 
resultsPumpSmallFixed2$iter ## 3000
mean(tail(resultsPumpSmallFixed2$param[, 1], 20), trim=.2) ## 0.8457023
mean(tail(resultsPumpSmallFixed2$param[, 2], 20), trim=.2) ## 1.321205
save(resultsPumpSmallFixed2, file="pumpSmallFixedCC.RData")

ptm <- proc.time()
resultsPumpSmallFixed3 <- computeMLE(pump, paramNodesPump,
                                    method="fixed", paramInit=init,
                                    compiledFuns=compiledFunsPump,
                                    stepsize=0.005,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timePumpSmallFixed3 <- proc.time() - ptm
timePumpSmallFixed3 ## 2.752 
resultsPumpSmallFixed3$iter ## 300
mean(tail(resultsPumpSmallFixed3$param[, 1], 20), trim=.2) ##  4.676562
mean(tail(resultsPumpSmallFixed3$param[, 2], 20), trim=.2) ##  11.2526
save(resultsPumpSmallFixed3, file="pumpSmallFixed300.RData")

ptm <- proc.time()
resultsPumpSmallFixed4 <- computeMLE(pump, paramNodesPump,
                                    method="fixed", paramInit=init,
                                    compiledFuns=compiledFunsPump,
                                    stepsize=0.005,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timePumpSmallFixed4 <- proc.time() - ptm
timePumpSmallFixed4 ## 2.873 
resultsPumpSmallFixed4$iter ## 300
mean(tail(resultsPumpSmallFixed4$param[, 1], 20), trim=.2) ##  4.672564
mean(tail(resultsPumpSmallFixed4$param[, 2], 20), trim=.2) ##11.25466
save(resultsPumpSmallFixed4, file="pumpSmallFixed300CC.RData")

# 3. Adadelta ----------------------------------------
ptm <- proc.time()
resultsPumpAdadelta <- computeMLE(pump, paramNodesPump,
                                  method="adadelta", paramInit=init,
                                  compiledFuns=compiledFunsPump,
                                  numMCMCSamples=numMCMCSamples,
                                  maxIter=300,
                                  boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timePumpAdadelta <- proc.time() - ptm
timePumpAdadelta ## 3.129 
mean(tail(resultsPumpAdadelta$param[, 1], 20), trim=.2) ## 1.075411
mean(tail(resultsPumpAdadelta$param[, 2], 20), trim=.2) ## 1.890135
save(resultsPumpAdadelta, file="pumpAdadelta.RData")

ptm <- proc.time()
resultsPumpAdadelta2 <- computeMLE(pump, paramNodesPump,
                                  method="adadelta", paramInit=init,
                                  compiledFuns=compiledFunsPump,
                                  numMCMCSamples=numMCMCSamples,
                                  maxIter=300,
                                  boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timePumpAdadelta2 <- proc.time() - ptm
timePumpAdadelta2 ## 1.110 
resultsPumpAdadelta2$iter ## 108
mean(tail(resultsPumpAdadelta2$param[, 1], 20), trim=.2) ## 1.089937
mean(tail(resultsPumpAdadelta2$param[, 2], 20), trim=.2) ## 1.904102
save(resultsPumpAdadelta2, file="pumpAdadeltaCC.RData")

# 4. Adam ----------------------------------------
ptm <- proc.time()
resultsPumpAdam <- computeMLE(pump, paramNodesPump,
                              method="adam", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timePumpAdam <- proc.time() - ptm
timePumpAdam ##  3.100 
mean(tail(resultsPumpAdam$param[, 1], 20), trim=.2) ## 0.8193667
mean(tail(resultsPumpAdam$param[, 2], 20), trim=.2) ## 1.257927
save(resultsPumpAdam, file="pumpAdam.RData")

ptm <- proc.time()
resultsPumpAdam2 <- computeMLE(pump, paramNodesPump,
                              method="adam", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timePumpAdam2 <- proc.time() - ptm
timePumpAdam2 ##  3.100  
resultsPumpAdam2$iter ## 176
mean(tail(resultsPumpAdam2$param[, 1], 20), trim=.2) ##  0.8298185
mean(tail(resultsPumpAdam2$param[, 2], 20), trim=.2) ## 1.278673
save(resultsPumpAdam2, file="pumpAdamCC.RData")


# 5. Newton-Raphson ----------------------------------------
ptm <- proc.time()
resultsPumpNR <- computeMLE(pump, paramNodes=paramNodesPump,
                            method="NR", paramInit=init,
                            compiledFuns=compiledFunsPump,
                            numMCMCSamples=numMCMCSamples,
                            tol=1e-20,
                            maxIter=300,
                            boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timePumpNR <- proc.time() - ptm
timePumpNR ## 4.480 
mean(resultsPumpNR$param[201:300,1], trim=.2) ## 8379.158
mean(resultsPumpNR$param[201:300,2], trim=.2) ## 10000
save(resultsPumpNR, file="pumpNR.RData")

ptm <- proc.time()
resultsPumpNR2 <- computeMLE(pump, paramNodes=paramNodesPump,
                            method="NR", paramInit=init,
                            compiledFuns=compiledFunsPump,
                            numMCMCSamples=numMCMCSamples,
                            tol=1e-20,
                            maxIter=300,
                            boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timePumpNR2 <- proc.time() - ptm
timePumpNR2 ## 4.377 
resultsPumpNR2$iter ## 300
mean(resultsPumpNR2$param[201:300,1], trim=.2) ## 9911.425
mean(resultsPumpNR2$param[201:300,2], trim=.2) ## 10000
save(resultsPumpNR2, file="pumpNRCC.RData")

# 6. 1-D sampling ----------------------------------------
ptm <- proc.time()
resultsPump1D <- computeMLE(pump, paramNodesPump,
                            method="ga1D", paramInit=init,
                            compiledFuns=compiledFunsPump,
                            numMCMCSamples=300, numMCMCSamples1D=300, 
                            maxIter=300,skipConvCheck=T)
timePump1D <- proc.time() - ptm
#Error in ga1DMLE(model = model, paramNodes = paramNodes, paramInit = paramInit,  : 
#                   unused argument (trackEffSizeGrad = FALSE)
timePump1D ## 12.128 
mean(tail(resultsPump1D$param[, 1], 20), trim=.2) ## 0.8376302
mean(tail(resultsPump1D$param[, 2], 20), trim=.2) ## 1.292479
save(resultsPump1D, file="pump1D.RData")

ptm <- proc.time()
resultsPump1D2 <- computeMLE(pump, paramNodesPump,
                            method="ga1D", paramInit=init,
                            compiledFuns=compiledFunsPump,
                            numMCMCSamples=300, numMCMCSamples1D=300, 
                            maxIter=300,skipConvCheck=F)
timePump1D2 <- proc.time() - ptm
#Error in ga1DMLE(model = model, paramNodes = paramNodes, paramInit = paramInit,  : 
#                   unused argument (trackEffSizeGrad = FALSE)
timePump1D2 ## 3.188 
resultsPump1D2$iter ## 75
mean(tail(resultsPump1D2$param[, 1], 20), trim=.2) ## 0.8955095
mean(tail(resultsPump1D2$param[, 2], 20), trim=.2) ## 1.494974
save(resultsPump1D2, file="pump1DCC.RData")


# 7. Hybrid
# Run 10 iterations of 1D, and then run Adam.
ptm <- proc.time()
resultsPump1DFirst <- computeMLE(pump, paramNodesPump,
                                 method="ga1D", paramInit=init,
                                 compiledFuns=compiledFunsPump,
                                 numMCMCSamples=numMCMCSamples, 
                                 numMCMCSamples1D=300, 
                                 maxIter=10)
resultsPumpAdamSecond <- computeMLE(pump, paramNodesPump,
                                    method="adam", 
                                    paramInit=tail(resultsPump1DFirst$param, 1),
                                    compiledFuns=compiledFunsPump,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary)
mean(tail(resultsPumpAdamSecond$param[, 1], 20), trim=.2) ## 0.8230988
mean(tail(resultsPumpAdamSecond$param[, 2], 20), trim=.2) ## 1.260072
timePumpHybrid <- proc.time() - ptm
timePumpHybrid

# MCEM -----------------------------
#source("MCEM_with_output.R")
pump2 <- pump$newModel()
box <- list(list(c('alpha', 'beta'), c(0, Inf)))
#pumpMCEM <- buildMCEM(model=pump2, latentNodes='theta[1:10]',
#                      boxConstraints=box, theta0=c(10, 10))
## new version that tracks
#Error: In sizeNimbleFunction: A nimbleFunction method should not be processed here.
#This occurred for: calc_E_llk(ARG2_theta_,ARG3_oldTheta_,1)
#This was part of the call:  svals[r] <- calc_E_llk(ARG2_theta_,ARG3_oldTheta_,1)

## regular
pumpMCEM <- buildMCEM(model=pump2, latentNodes='theta[1:10]',
                      boxConstraints=box)

ptm <- proc.time()
resultsPumpMCEM <- pumpMCEM$run(initM = 1000)
timePumpMCEM <- proc.time() - ptm
timePumpMCEM ## 48.970 

resultsPumpMCEM
#alpha     beta 
#0.820380 1.249634 

# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 1.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   alpha      beta 
# 0.8264483 1.1526601 
# Convergence Criterion: 1.001.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 2.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   alpha      beta 
# 0.8146395 1.2225627 
# Convergence Criterion: 0.01901154.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 3.
# Current number of MCMC iterations: 2188.
# Parameter Estimates: 
#   alpha      beta 
# 0.8163346 1.2401123 
# Convergence Criterion: 0.001217226.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 4.
# Current number of MCMC iterations: 2188.
# Parameter Estimates: 
#   alpha      beta 
# 0.8200225 1.2616166 
# Convergence Criterion: 0.001256843.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 5.
# Current number of MCMC iterations: 2188.
# Parameter Estimates: 
#   alpha      beta 
# 0.8160177 1.2431745 
# Convergence Criterion: 0.001192686.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 6.
# Current number of MCMC iterations: 6197.
# Parameter Estimates: 
#   alpha     beta 
# 0.820380 1.249634 
# Convergence Criterion: 0.0002174223.


#iterTimePumpMCEM <- c(0.872, 3.440, 5.129,
#                      6.862, 8.522, 12.686,
#                      18.877, 37.616)

## 8 iterations

# Save outputs ------------------------
save(resultsPumpMCEM, file="MCEM_pump.RData")


timeInfo <- cbind(c(timePumpFixed, timePumpSmallFixed,
                    timePumpAdadelta, timePumpAdam,
                    timePumpNR, timePump1D, timePumpMCEM),
                  c("fixed", "smallFixed", "adadelta", "adam",
                    "NR", "1D", "mcem"))

write.csv(timeInfo, "timePump.csv", row.names=F)
