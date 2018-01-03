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
                                 boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeLogRegFixed <- proc.time() - ptm 
timeLogRegFixed ## 3.042 
mean(tail(resultsLogregFixed$param[,1],20), trim=.2)  ## -0.856524
mean(tail(resultsLogregFixed$param[,2],20), trim=.2)  ## 1.183854
mean(tail(resultsLogregFixed$param[,3],20), trim=.2)  ## 1.175493
save(resultsLogregFixed, file="logRegFixed.RData")

ptm <- proc.time()
resultsLogregFixed2 <- computeMLE(logreg, paramNodesLogreg,
                                 method="fixed", paramInit=init,
                                 stepsize=0.05,
                                 compiledFuns=compiledFunsLogreg,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=300,
                                 boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeLogRegFixed2 <- proc.time() - ptm 
timeLogRegFixed2 ##  0.778 
mean(tail(resultsLogregFixed2$param[,1],20), trim=.2)  ## -0.4839189
mean(tail(resultsLogregFixed2$param[,2],20), trim=.2)  ## 1.323371
mean(tail(resultsLogregFixed2$param[,3],20), trim=.2)  ## 0.647126
resultsLogregFixed2$iter ## 78
save(resultsLogregFixed2, file="logRegFixedCC.RData")


# 2. Small fixed step size ----------------------------------------
ptm <- proc.time()
resultsLogregSmallFixed <- computeMLE(logreg, paramNodesLogreg,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsLogreg,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeLogRegSmallFixed <- proc.time() - ptm 
timeLogRegSmallFixed ## 2.934 
mean(tail(resultsLogregSmallFixed$param[,1],20), trim=.2)  ## -0.5468611
mean(tail(resultsLogregSmallFixed$param[,2],20), trim=.2)  ##  1.304298
mean(tail(resultsLogregSmallFixed$param[,3],20), trim=.2)  ## 0.2528262
save(resultsLogregSmallFixed, file="logRegSmallFixed.RData")

ptm <- proc.time()
resultsLogregSmallFixed2 <- computeMLE(logreg, paramNodesLogreg,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsLogreg,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeLogRegSmallFixed2 <- proc.time() - ptm 
timeLogRegSmallFixed2 ## 3.141 
mean(tail(resultsLogregSmallFixed2$param[,1],20), trim=.2)  ## -0.5449848
mean(tail(resultsLogregSmallFixed2$param[,2],20), trim=.2)  ##  1.305241
mean(tail(resultsLogregSmallFixed2$param[,3],20), trim=.2)  ## 0.2437223
resultsLogregSmallFixed2$iter ## 300
save(resultsLogregSmallFixed2, file="logRegSmallFixedCC.RData")

ptm <- proc.time()
resultsLogregSmallFixed3 <- computeMLE(logreg, paramNodesLogreg,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsLogreg,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=3000,
                                      boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeLogRegSmallFixed3 <- proc.time() - ptm 
timeLogRegSmallFixed3 ## 33.723 
mean(tail(resultsLogregSmallFixed3$param[,1],20), trim=.2)  ## -0.5395547
mean(tail(resultsLogregSmallFixed3$param[,2],20), trim=.2)  ##  1.311037
mean(tail(resultsLogregSmallFixed3$param[,3],20), trim=.2)  ##0.244578
save(resultsLogregSmallFixed3, file="logRegSmallFixed3000.RData")

ptm <- proc.time()
resultsLogregSmallFixed4 <- computeMLE(logreg, paramNodesLogreg,
                                       method="fixed", paramInit=init,
                                       stepsize=0.005,
                                       compiledFuns=compiledFunsLogreg,
                                       numMCMCSamples=numMCMCSamples,
                                       maxIter=3000,
                                       boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeLogRegSmallFixed4 <- proc.time() - ptm 
timeLogRegSmallFixed4 ##  3.296 
mean(tail(resultsLogregSmallFixed4$param[,1],20), trim=.2)  ##-0.5477317
mean(tail(resultsLogregSmallFixed4$param[,2],20), trim=.2)  ##   1.306503
mean(tail(resultsLogregSmallFixed4$param[,3],20), trim=.2)  ##  0.24909
resultsLogregSmallFixed4$iter ## 150
save(resultsLogregSmallFixed4, file="logRegSmallFixedCC3000.RData")

# 3. Adadelta ----------------------------------------
ptm <- proc.time()
resultsLogregAdadelta <- computeMLE(logreg, paramNodesLogreg,
                                    method="adadelta", paramInit=init,
                                    compiledFuns=compiledFunsLogreg,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeLogRegAdadelta <- proc.time() - ptm 
timeLogRegAdadelta ##  3.190 
mean(tail(resultsLogregAdadelta$param[, 1], 20), trim=.2)  ## -0.575974
mean(tail(resultsLogregAdadelta$param[, 2], 20), trim=.2)  ## 1.3102
mean(tail(resultsLogregAdadelta$param[, 3], 20), trim=.2)  ##  0.3415709
save(resultsLogregAdadelta, file="logRegAdadelta.RData")

ptm <- proc.time()
resultsLogregAdadelta2 <- computeMLE(logreg, paramNodesLogreg,
                                    method="adadelta", paramInit=init,
                                    compiledFuns=compiledFunsLogreg,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeLogRegAdadelta2 <- proc.time() - ptm 
timeLogRegAdadelta2 ## 0.543 
mean(tail(resultsLogregAdadelta2$param[, 1], 20), trim=.2)  ##-0.5785057
mean(tail(resultsLogregAdadelta2$param[, 2], 20), trim=.2)  ##1.308722
mean(tail(resultsLogregAdadelta2$param[, 3], 20), trim=.2)  ##0.4058867
resultsLogregAdadelta2$iter ## 44
save(resultsLogregAdadelta2, file="logRegAdadeltaCC.RData")

# 4. Adam ----------------------------------------
ptm <- proc.time()
resultsLogregAdam <- computeMLE(logreg, paramNodesLogreg,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                stepsize=0.2,
                                eps=1e-4,
                                maxIter=300,
                                boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeLogRegAdam <- proc.time() - ptm 
timeLogRegAdam ## 3.005 
mean(tail(resultsLogregAdam$param[, 1], 20), trim=.2) ## -0.5575068
mean(tail(resultsLogregAdam$param[, 2], 20), trim=.2) ##  1.313344
mean(tail(resultsLogregAdam$param[, 3], 20), trim=.2) ##  0.2433034
save(resultsLogregAdam, file="logRegAdam.RData")


ptm <- proc.time()
resultsLogregAdam2 <- computeMLE(logreg, paramNodesLogreg,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                stepsize=0.2,
                                eps=1e-4,
                                maxIter=300,
                                boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeLogRegAdam2 <- proc.time() - ptm 
timeLogRegAdam2 ##  0.752 
mean(tail(resultsLogregAdam2$param[, 1], 20), trim=.2) ## -0.543162
mean(tail(resultsLogregAdam2$param[, 2], 20), trim=.2) ##  1.328982
mean(tail(resultsLogregAdam2$param[, 3], 20), trim=.2) ## 0.2617962
resultsLogregAdam2$iter ## 71
save(resultsLogregAdam2, file="logRegAdamCC.RData")

# 5. Newton-Raphson ----------------------------------------
ptm <- proc.time()
resultsLogregNR <- computeMLE(logreg, paramNodesLogreg,
                              method="NR", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeLogRegNR <- proc.time() - ptm 
timeLogRegNR ## 8.768 
mean(tail(resultsLogregNR$param[,1],20), trim=.2) ##  -6.389591
mean(tail(resultsLogregNR$param[,2],20), trim=.2) ## -10.43822
mean(tail(resultsLogregNR$param[,3],20), trim=.2) ## 10 (boundary)

save(resultsLogregNR, file="logRegNR.RData")

ptm <- proc.time()
resultsLogregNR2 <- computeMLE(logreg, paramNodesLogreg,
                              method="NR", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeLogRegNR2 <- proc.time() - ptm 
timeLogRegNR2 ##  9.029 
mean(tail(resultsLogregNR2$param[,1],20), trim=.2) ##5.316958
mean(tail(resultsLogregNR2$param[,2],20), trim=.2) ## 12.29115
mean(tail(resultsLogregNR2$param[,3],20), trim=.2) ## 5.34089
resultsLogregNR2$iter ## 300
save(resultsLogregNR2, file="logRegNRCC.RData")

# 6. 1-D sampling ----------------------------------------
ptm <- proc.time()
resultsLogreg1D <- computeMLE(logreg, paramNodesLogreg,
                              method="ga1D", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=300, numMCMCSamples1D=300, 
                              maxIter=300,skipConvCheck=T)
timeLogReg1D <- proc.time() - ptm 
timeLogReg1D ## 11.635 
mean(tail(resultsLogreg1D$param[, 1]), trim=.2)  ## -0.5133427
mean(tail(resultsLogreg1D$param[, 2]), trim=.2)  ## 1.278134
mean(tail(resultsLogreg1D$param[, 3]), trim=.2)  ## 0.2491426
save(resultsLogreg1D, file="logReg1D.RData")

ptm <- proc.time()
resultsLogreg1D2 <- computeMLE(logreg, paramNodesLogreg,
                              method="ga1D", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=300, numMCMCSamples1D=300, 
                              maxIter=300,skipConvCheck=F)
timeLogReg1D2 <- proc.time() - ptm 
timeLogReg1D2 ##  4.36
mean(tail(resultsLogreg1D2$param[, 1]), trim=.2)  ##  -0.5453086
mean(tail(resultsLogreg1D2$param[, 2]), trim=.2)  ## 1.312937
mean(tail(resultsLogreg1D2$param[, 3]), trim=.2)  ##  0.232226
resultsLogreg1D2$iter ## 106
save(resultsLogreg1D2, file="logReg1DCC.RData")

# MCEM -----------------------------
#source("MCEM_with_output.R")
logreg2 <- logreg$newModel()
latentNodes <- logreg2$getNodeNames(latentOnly=TRUE, stochOnly=TRUE)

#modMCEM <- buildMCEM(model=logreg2, latentNodes=latentNodes, theta0=c(0,0,1))
modMCEM <- buildMCEM(model=logreg2, latentNodes=latentNodes)


ptm <- proc.time()
modMLE <-  modMCEM$run(initM = 1000)
timeLogRegMCEM <- proc.time() - ptm 
timeLogRegMCEM ## 119.665 
modMLE  
#beta0      beta1   sigma_RE 
#-0.5512313  1.3105122  0.2521761 
save(modMLE, file="MCEM_logReg.RData")

# 
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 1.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0       beta1    sigma_RE 
# -0.05281449  0.09887732  0.80329093 
# Convergence Criterion: 1.001.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 2.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.0877019  0.2308019  0.6997026 
# Convergence Criterion: 0.4068105.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 3.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.1408130  0.3625702  0.6432272 
# Convergence Criterion: 0.3320016.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 4.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.2015608  0.5172245  0.5807458 
# Convergence Criterion: 0.3511456.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 5.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.2805479  0.6661176  0.5180076 
# Convergence Criterion: 0.3702636.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 6.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.3480750  0.8321879  0.4530430 
# Convergence Criterion: 0.4395037.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 7.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.3985227  0.9589871  0.4019443 
# Convergence Criterion: 0.3082527.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 8.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.4518312  1.0729117  0.3646155 
# Convergence Criterion: 0.2218657.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 9.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.4921702  1.1617262  0.3336366 
# Convergence Criterion: 0.182773.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 10.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5020252  1.1946919  0.3130408 
# Convergence Criterion: 0.09175928.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 11.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5265361  1.2461977  0.3061983 
# Convergence Criterion: 0.05902465.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 12.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5358003  1.2833379  0.2857678 
# Convergence Criterion: 0.08164142.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 13.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5385508  1.3029074  0.2715428 
# Convergence Criterion: 0.04706229.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 14.
# Current number of MCMC iterations: 1250.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5455307  1.3047378  0.2643208 
# Convergence Criterion: 0.01814265.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 15.
# Current number of MCMC iterations: 1625.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5504221  1.3097884  0.2539406 
# Convergence Criterion: 0.0282197.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 16.
# Current number of MCMC iterations: 1625.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5470160  1.3180745  0.2578547 
# Convergence Criterion: 0.01526664.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 17.
# Current number of MCMC iterations: 2188.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5464096  1.3097549  0.2558507 
# Convergence Criterion: 0.006060111.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 18.
# Current number of MCMC iterations: 3032.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5427631  1.3061285  0.2533162 
# Convergence Criterion: 0.003766591.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 19.
# Current number of MCMC iterations: 4298.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5495354  1.3157580  0.2518805 
# Convergence Criterion: 0.003697544.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 20.
# Current number of MCMC iterations: 4298.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5528191  1.3194567  0.2536477 
# Convergence Criterion: 0.002433969.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 21.
# Current number of MCMC iterations: 4298.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5489866  1.3083054  0.2531576 
# Convergence Criterion: 0.004153772.
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
#   Iteration Number: 22.
# Current number of MCMC iterations: 13319.
# Parameter Estimates: 
#   beta0      beta1   sigma_RE 
# -0.5512313  1.3105122  0.2521761 
# Convergence Criterion: 0.0008661106.



# iterTimeLogRegMCEM <- c(1.034, 3.581, 5.141, 9.036, 11.374,
#                         13.631, 15.906, 18.187, 20.537, 22.762,
#                         25.072, 27.397, 29.691, 32.084, 34.349,
#                         36.607, 38.930, 41.257, 43.555, 45.975,
#                         48.269, 58.947, 64.053, 69.666, 74.706,
#                         79.703, 104.517, 115.557, 213.341,
#                         252.181)

timeInfo=cbind(c(timeLogRegFixed, timeLogRegSmallFixed, timeLogRegAdadelta, 
                 timeLogRegAdam, timeLogRegNR, timeLogReg1D, timeLogRegMCEM),
               c("fixed", "smallFixed", "adadelta", "adam", "NR", "1D", 
                 "mcem"))


write.csv(timeInfo,"timeLogReg.csv",row.names=F)

## 
r=c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46)
n=c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79)
x=c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)


xlong<-c()
rlong<-c()
indiv<-c()

for(i in 1:length(r)){
  xlong=c(xlong,rep(x[i],n[i]))
  rlong<-c(rlong,rep(1,r[i]),rep(0,n[i]-r[i]))
  indiv<-c(indiv,rep(i,n[i]))
}

data=cbind.data.frame(xlong,rlong,indiv)

require(lme4)
ptm <- proc.time()
gmre=glmer(rlong~xlong+(1|indiv),family=binomial(),data=data)
glmmTime <- proc.time() - ptm ##  0.093 

summary(gmre)

# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: rlong ~ xlong + (1 | indiv)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 715.2    728.2   -354.6    709.2      551 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.6091 -0.7639  0.6215  0.7669  1.4491 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# indiv  (Intercept) 0.06184  0.2487  
# Number of obs: 554, groups:  indiv, 10
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.5483     0.1703  -3.220  0.00128 ** 
#   xlong         1.3104     0.2460   5.326    1e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# xlong -0.691


#beta0: -0.5483 
#beta1:  1.3104   
#sigma_RE:  0.2460 

