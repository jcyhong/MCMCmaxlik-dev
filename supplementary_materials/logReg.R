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
mean(tail(resultsLogregFixed$param[,1],20), trim=.2)  ## -0.7522742
mean(tail(resultsLogregFixed$param[,2],20), trim=.2)  ## 1.25159
mean(tail(resultsLogregFixed$param[,3],20), trim=.2)  ## 0.8827718
save(resultsLogregFixed, file="logRegFixed.RData")

# 2. Small fixed step size ----------------------------------------
ptm <- proc.time()
resultsLogregSmallFixed <- computeMLE(logreg, paramNodesLogreg,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsLogreg,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary)
timeLogRegSmallFixed <- proc.time() - ptm ## 2.530 
mean(tail(resultsLogregSmallFixed$param[,1],20), trim=.2)  ## -0.5485272
mean(tail(resultsLogregSmallFixed$param[,2],20), trim=.2)  ##  1.310605
mean(tail(resultsLogregSmallFixed$param[,3],20), trim=.2)  ## 0.251749
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
mean(tail(resultsLogregAdadelta$param[, 1], 20), trim=.2)  ## -0.555296
mean(tail(resultsLogregAdadelta$param[, 2], 20), trim=.2)  ## 1.319334
mean(tail(resultsLogregAdadelta$param[, 3], 20), trim=.2)  ## 0.3160231
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
mean(tail(resultsLogregAdam$param[, 1], 20), trim=.2) ## -0.5482783
mean(tail(resultsLogregAdam$param[, 2], 20), trim=.2) ##  1.310051
mean(tail(resultsLogregAdam$param[, 3], 20), trim=.2) ## 0.2494285
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
mean(tail(resultsLogregNR$param[,1],20), trim=.2) ## 0.6831159
mean(tail(resultsLogregNR$param[,2],20), trim=.2) ## -1.625356
mean(tail(resultsLogregNR$param[,3],20), trim=.2) ## 10 (boundary)
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
mean(tail(resultsLogreg1D$param[, 1]), trim=.2)  ## -0.5677065
mean(tail(resultsLogreg1D$param[, 2]), trim=.2)  ## 1.327946
mean(tail(resultsLogreg1D$param[, 3]), trim=.2)  ## 0.2804238
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
# 
# 


#beta0: -0.5483
#beta1: 1.3104     
#sigma_RE: 0.2487  

