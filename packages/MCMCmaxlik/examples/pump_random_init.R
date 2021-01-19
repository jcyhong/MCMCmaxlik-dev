# Random initialization to calibrate MSE ----

# Pump example (2 parameters)

# Load the packages ---------------------------------------
library(nimble)
nimbleOptions(experimentalEnableDerivs = TRUE)
library("MCMCmaxlik")
library('parallel')
library('data.table')
library('dplyr')

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
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                         31.4, 1.05, 1.05, 2.1, 10.5))
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))

# Build the model.
pump <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

paramNodesPump <- pump$getNodeNames(topOnly = T)
compileNimble(pump)

# Contour plot prep work ---------------

fr2 <- function(alpha, beta) {
  x <- pumpData$x
  t <- pumpConsts$t
  first_term <- (-alpha - beta) + (0.1 - 1) * log(beta)
  second_term <- 0
  for (i in 1:length(x)) {
    second_term <- second_term + lgamma(alpha + x[i]) - 
      (alpha + x[i]) * log(beta + t[i])
  }
  third_term <- length(x) * (alpha * log(beta) - lgamma(alpha))
  #-(first_term + second_term + third_term)
  -(second_term + third_term)
}

fr2withprior <- function(alpha, beta) {
  x <- pumpData$x
  t <- pumpConsts$t
  first_term <- (-alpha - beta) + (0.1 - 1) * log(beta)
  second_term <- 0
  for (i in 1:length(x)) {
    second_term <- second_term + lgamma(alpha + x[i]) - 
      (alpha + x[i]) * log(beta + t[i])
  }
  third_term <- length(x) * (alpha * log(beta) - lgamma(alpha))
  -(first_term + second_term + third_term)
}

# Generate matrix for contour plot.
alphaSeq <- seq(.01, 10, length=1000)
betaSeq <- seq(.01, 16, length=1000)
mat1 <- outer(alphaSeq, betaSeq, fr2)
contourdf <- data.frame(cbind(expand.grid(alphaSeq, betaSeq), as.vector(mat1)))
names(contourdf) <- c("alpha", "beta", "density")

# Testing Algorithms ------------------------------------------------------
B <- 30

# Compile the necessary functions.
compiledFunsPump <- buildMCMCmaxlik(pump, paramNodesPump)
boundary <- list(c(0.05, 10000), c(0.05, 10000))
numMCMCSamples <- 300

set.seed(1000)
df_list <- mclapply(1:30, function(b) {
  init <- c(runif(1, min=0.05, max=15), runif(1, min=0.05, max=15))
  resultsPumpFixed <- computeMLE(pump, paramNodesPump,
                                 method="fixed", paramInit=init,
                                 compiledFuns=compiledFunsPump,
                                 stepsize=0.05,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=300,
                                 boundary=boundary)
  
  resultsPumpSmallFixed <- computeMLE(pump, paramNodesPump,
                                      method="fixed", paramInit=init,
                                      compiledFuns=compiledFunsPump,
                                      stepsize=0.005,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary)
  
  resultsPumpNR <- computeMLE(pump, paramNodes=paramNodesPump,
                              method="NR", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              tol=1e-20,
                              maxIter=300,
                              boundary=boundary)
  
  resultsPumpAdam <- computeMLE(pump, paramNodesPump,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsPump,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300)
  
  resultsPump1D <- computeMLE(pump, paramNodesPump,
                              method="ga1D", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              numMCMCSamples1D=numMCMCSamples,
                              maxIter=300)
  
  df <- data.frame(do.call(rbind, list(resultsPumpFixed$MLE,
                      resultsPumpSmallFixed$MLE,
                      resultsPumpNR$MLE,
                      resultsPumpAdam$MLE,
                      resultsPump1D$MLE)))
  df$method <- c('Fixed (0.05)', 'Fixed (0.005)', 
                 'Newton-Raphson', 'Adam', '1D sampling')
  
  names(df) <- c("alpha", "beta", "method")
  df
}, mc.cores=4)

df_combined <- rbindlist(df_list)
df_rmse_300 <- df_combined %>% group_by(method) %>%
  summarize(RMSE_alpha=sqrt(mean((alpha - 0.823)^2)),
            RMSE_beta=sqrt(mean((beta - 1.262)^2)))

# Repeat with MCMC sample sizes = 3000.

numMCMCSamples <- 3000

set.seed(2000)
df_list_3000 <- mclapply(1:30, function(b) {
  init <- c(runif(1, min=0.05, max=15), runif(1, min=0.05, max=15))
  resultsPumpFixed <- computeMLE(pump, paramNodesPump,
                                 method="fixed", paramInit=init,
                                 compiledFuns=compiledFunsPump,
                                 stepsize=0.05,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=300,
                                 boundary=boundary)
  
  resultsPumpSmallFixed <- computeMLE(pump, paramNodesPump,
                                      method="fixed", paramInit=init,
                                      compiledFuns=compiledFunsPump,
                                      stepsize=0.005,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary)
  
  resultsPumpNR <- computeMLE(pump, paramNodes=paramNodesPump,
                              method="NR", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              tol=1e-20,
                              maxIter=300,
                              boundary=boundary)
  
  resultsPumpAdam <- computeMLE(pump, paramNodesPump,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsPump,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300)
  
  resultsPump1D <- computeMLE(pump, paramNodesPump,
                              method="ga1D", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              numMCMCSamples1D=numMCMCSamples,
                              maxIter=300)
  
  df <- data.frame(do.call(rbind, list(resultsPumpFixed$MLE,
                                       resultsPumpSmallFixed$MLE,
                                       resultsPumpNR$MLE,
                                       resultsPumpAdam$MLE,
                                       resultsPump1D$MLE)))
  df$method <- c('Fixed (0.05)', 'Fixed (0.005)', 
                 'Newton-Raphson', 'Adam', '1D sampling')
  
  names(df) <- c("alpha", "beta", "method")
  df
}, mc.cores=4)

df_combined_3000 <- rbindlist(df_list_3000)
df_rmse_3000 <- df_combined_3000 %>% group_by(method) %>%
  summarize(RMSE_alpha=sqrt(mean((alpha - 0.823)^2)),
            RMSE_beta=sqrt(mean((beta - 1.262)^2)))

# MCEM result ----
pumpMCEM <- buildMCEM(model = pump, latentNodes = 'theta',
                      boxConstraints = list( list( c('alpha', 'beta'),
                                                   limits = c(0, Inf) ) ))
MCEM_time <- proc.time()
out <- pumpMCEM$run()
MCEM_time <- proc.time() - MCEM_time

abs(0.825 - 0.823)
abs(1.267 - 1.262)

write.csv(df_rmse_300, file='pump_rmse_300.csv', row.names=F)
write.csv(df_rmse_3000, file='pump_rmse_3000.csv', row.names=F)