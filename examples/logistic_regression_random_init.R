# Logistic regression example (3 parameters)

# Load the packages ---------------------------------------
library("MCMCmaxlik")
library(nimble)
nimbleOptions(experimentalEnableDerivs = TRUE)
library('parallel')

# Model specification -------------------------------------
code <- nimbleCode({
  beta0 ~ dnorm(0, sd = 10000)
  beta1 ~ dnorm(0, sd = 10000)
  sigma_RE ~ dunif(0, 1000)
  for(i in 1:N) {
    beta2[i] ~ dnorm(0, sd = sigma_RE)
    logit(p[i]) <- beta0 + beta1 * x[i] + beta2[i]
    r[i] ~ dbin(p[i], n[i])
  }
})

## constants, data, and initial values
constants <- list(N = 10)

data <- list(
  r = c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46),
  n = c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79),
  x = c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
)

inits <- list(beta0 = 0, beta1 = 0, sigma_RE = 1)

logreg <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check = FALSE)

paramNodesLogreg <- logreg$getNodeNames(topOnly = T)
Clogreg <- compileNimble(logreg)


# Testing Algorithms ------------------------------------------------------
# Compile the necessary functions.
compiledFunsLogreg <- buildMCMCmaxlik(logreg, paramNodesLogreg)

init <- c(0, 0, 1)
boundary <- list(c(-40, 40), c(-40, 40), c(0.05, 10))
numMCMCSamples <- 20

set.seed(1500)
df_list_20 <- lapply(1:30, function(b) {
  init <- c(runif(1, min=-10, max=10), 
            runif(1, min=-10, max=10),
            runif(1, min=0.05, max=15))
  resultsLogregFixed <- computeMLE(logreg, paramNodesLogreg,
                                   method="fixed", paramInit=init,
                                   stepsize=0.05,
                                   compiledFuns=compiledFunsLogreg,
                                   numMCMCSamples=numMCMCSamples,
                                   maxIter=300,
                                   boundary=boundary)
  resultsLogregSmallFixed <- computeMLE(logreg, paramNodesLogreg,
                                        method="fixed", paramInit=init,
                                        stepsize=0.005,
                                        compiledFuns=compiledFunsLogreg,
                                        numMCMCSamples=numMCMCSamples,
                                        maxIter=300,
                                        boundary=boundary)
  #resultsLogregNR <- computeMLE(logreg, paramNodesLogreg,
  #                              method="NR", paramInit=init,
  #                         compiledFuns=compiledFunsLogreg,
  #                             numMCMCSamples=numMCMCSamples,
  #                              maxIter=300,
  #                              boundary=boundary)
  
  resultsLogregAdam <- computeMLE(logreg, paramNodesLogreg,
                                  method="adam", paramInit=init,
                                  compiledFuns=compiledFunsLogreg,
                                  numMCMCSamples=numMCMCSamples,
                                  maxIter=300,
                                  boundary=boundary)
  
  resultsLogreg1D <- computeMLE(logreg, paramNodesLogreg,
                                method="ga1D", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300)
  
  df <- data.frame(do.call(rbind, list(resultsLogregFixed$MLE,
                                       resultsLogregSmallFixed$MLE,
                                      # resultsLogregNR$MLE,
                                       resultsLogregAdam$MLE,
                                       resultsLogreg1D$MLE)))
  df$method <- c('Fixed (0.05)', 'Fixed (0.005)', 
                 #'Newton-Raphson', 
                 'Adam', '1D sampling')
  
  names(df) <- c("beta0", "beta1", "sigma_RE", "method")
  df
})

df_combined_20 <- rbindlist(df_list_20)
df_rmse_20 <- df_combined_20 %>% group_by(method) %>%
  summarize(RMSE_beta0=sqrt(mean((beta0 - (-0.519))^2)),
            RMSE_beta1=sqrt(mean((beta1 - 1.019)^2)),
            RMSE_sigma_RE=sqrt(mean((sigma_RE - 0.307)^2)))

set.seed(1600)
numMCMCSamples <- 300
df_list_300 <- lapply(1:30, function(b) {
  init <- c(runif(1, min=-10, max=10), 
            runif(1, min=-10, max=10),
            runif(1, min=0.05, max=15))
  resultsLogregFixed <- computeMLE(logreg, paramNodesLogreg,
                                   method="fixed", paramInit=init,
                                   stepsize=0.05,
                                   compiledFuns=compiledFunsLogreg,
                                   numMCMCSamples=numMCMCSamples,
                                   maxIter=300,
                                   boundary=boundary)
  resultsLogregSmallFixed <- computeMLE(logreg, paramNodesLogreg,
                                        method="fixed", paramInit=init,
                                        stepsize=0.005,
                                        compiledFuns=compiledFunsLogreg,
                                        numMCMCSamples=numMCMCSamples,
                                        maxIter=300,
                                        boundary=boundary)
  #resultsLogregNR <- computeMLE(logreg, paramNodesLogreg,
  #                              method="NR", paramInit=init,
  #                         compiledFuns=compiledFunsLogreg,
  #                             numMCMCSamples=numMCMCSamples,
  #                              maxIter=300,
  #                              boundary=boundary)
  
  resultsLogregAdam <- computeMLE(logreg, paramNodesLogreg,
                                  method="adam", paramInit=init,
                                  compiledFuns=compiledFunsLogreg,
                                  numMCMCSamples=numMCMCSamples,
                                  maxIter=300,
                                  boundary=boundary)
  
  resultsLogreg1D <- computeMLE(logreg, paramNodesLogreg,
                                method="ga1D", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300)
  
  df <- data.frame(do.call(rbind, list(resultsLogregFixed$MLE,
                                       resultsLogregSmallFixed$MLE,
                                       # resultsLogregNR$MLE,
                                       resultsLogregAdam$MLE,
                                       resultsLogreg1D$MLE)))
  df$method <- c('Fixed (0.05)', 'Fixed (0.005)', 
                 #'Newton-Raphson', 
                 'Adam', '1D sampling')
  
  names(df) <- c("beta0", "beta1", "sigma_RE", "method")
  df
})

df_combined_300 <- rbindlist(df_list_300)
df_rmse_300 <- df_combined_300 %>% group_by(method) %>%
  summarize(RMSE_beta0=sqrt(mean((beta0 - (-0.519))^2)),
            RMSE_beta1=sqrt(mean((beta1 - 1.019)^2)),
            RMSE_sigma_RE=sqrt(mean((sigma_RE - 0.307)^2)))

write.csv(df_rmse_20, file='logistic_regression_rmse_20.csv', row.names=F)
write.csv(df_rmse_300, file='logistic_regression_rmse_300.csv', row.names=F)


LogregMCEM <- buildMCEM_AD(model=logreg, 
                          latentNodes = 'beta2',
                          C = 0.01,
                          boxConstraints = list( 
                            list(c('sigma_RE'), 
                                 limits = c(0, 1000) ) ))
LogregMCEM_time <- proc.time()
out <- LogregMCEM$run()
LogregMCEM_time <- proc.time() - LogregMCEM_time

LogregMCEM_MLEs <- matrix(nrow=3, ncol=30)

for (i in 1:30) {
  init <- c(runif(1, min=-10, max=10), 
            runif(1, min=-10, max=10),
            runif(1, min=0.05, max=15))
  LogregMCEM_MLEs[, i] <- LogregMCEM$run(thetaInit = init)
}

LogregMCEM_rmse <- apply(
  apply(LogregMCEM_MLEs, 2, function(x) {(x - c(-0.519, 1.019, 0.307))^2}),
  1,
  function(y) {sqrt(mean(y))}
)

df_rmse_20 <- rbind(
  df_rmse_20, 
  data.frame(method='MCEM', RMSE_beta0=LogregMCEM_rmse[1], RMSE_beta1=LogregMCEM_rmse[2],
             RMSE_sigma_RE=LogregMCEM_rmse[3])
)

df_rmse_300 <- rbind(
  df_rmse_300, 
  data.frame(method='MCEM', RMSE_beta0=LogregMCEM_rmse[1], RMSE_beta1=LogregMCEM_rmse[2],
             RMSE_sigma_RE=LogregMCEM_rmse[3])
)

write.csv(df_rmse_20, file='logistic_regression_rmse_20.csv', row.names=F)
write.csv(df_rmse_300, file='logistic_regression_rmse_300.csv', row.names=F)