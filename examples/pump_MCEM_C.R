# Pump example (2 parameters)

# Load the packages ---------------------------------------
library(nimble)
nimbleOptions(experimentalEnableDerivs = TRUE)
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


# Testing Algorithms ------------------------------------------------------

# Compile the necessary functions.
#compiledFunsPump <- buildMCMCmaxlik(pump, paramNodesPump)

source("~/Desktop/MCMCmaxlik-dev/packages/MCMCmaxlik/R/MCEM_AD_build.R")


 C_opts = c(0.001, 0.01, 0.1)
 
 ## can only run one without R session being aborted
 ## have to restart R before every run, is there a way to update C without re-compiling?
 
i=3
   pumpMCEM <- buildMCEM_AD(model = pump, latentNodes = 'theta',
                            boxConstraints = list( list(c('alpha', 'beta'),
                                                        limits = c(0, Inf) ) ),
                            C = C_opts[i])
   
   MCEM_time <- proc.time()
   out <- pumpMCEM$run()
   MCEM_time <- proc.time() - MCEM_time
   

#alpha      beta 
#0.8235263 1.2630080 
#0.8265624 1.2462309 
#0.8260041 1.2546285 

#user  system elapsed 
#3.707   0.057   3.759 
#0.171   0.003   0.173 
#0.129   0.002   0.131 
 
