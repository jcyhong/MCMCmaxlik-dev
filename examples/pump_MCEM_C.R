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

#### C ####
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
   
### gamma ###

#probability of deciding that the algorithm has converged, that is, that the difference between two Q functions is less than C, when in fact it has not. Default is 0.05.
   
gamma_opts = c(0.01, 0.05, 0.1) ## don't have to run 2, same as 1 above

i=3
pumpMCEM <- buildMCEM_AD(model = pump, latentNodes = 'theta',
                         boxConstraints = list( list(c('alpha', 'beta'),
                                                     limits = c(0, Inf) ) ),
                         C = 0.001, gamma = gamma_opts[i])

MCEM_time <- proc.time()
out <- pumpMCEM$run()
MCEM_time <- proc.time() - MCEM_time

#alpha      beta 
#0.8241349 1.2679344 
#0.8235263 1.2630080 
#0.8294857 1.2861582 

#user  system elapsed 
#1.150   0.019   1.163 
#3.707   0.057   3.759 
# 0.518   0.009   0.522 

### alpha ### 

#probability of a type one error - here, the probability of accepting a parameter estimate that does not increase the likelihood. Default is 0.25.

alpha_opts = c(0.125, 0.25, 0.5) ## don't have to run 2, same as 1 above

i=3
pumpMCEM <- buildMCEM_AD(model = pump, latentNodes = 'theta',
                         boxConstraints = list( list(c('alpha', 'beta'),
                                                     limits = c(0, Inf) ) ),
                         C = 0.001, alpha = alpha_opts[i])

MCEM_time <- proc.time()
out <- pumpMCEM$run()
MCEM_time <- proc.time() - MCEM_time

#alpha      beta 
#0.8230365 1.2624168    
#0.8235263 1.2630080 
#0.8194606 1.2442112 

#user  system elapsed 
#2.665   0.057   2.758 
#3.707   0.057   3.759 
#0.539   0.013   0.548 

### beta ###
   
#probability of a type two error - here, the probability of rejecting a parameter estimate that does increase the likelihood. Default is 0.25.
    
beta_opts = c(0.125, 0.25, 0.5) ## don't have to run 2, same as 1 above
 
i=3
pumpMCEM <- buildMCEM_AD(model = pump, latentNodes = 'theta',
                         boxConstraints = list( list(c('alpha', 'beta'),
                                                     limits = c(0, Inf) ) ),
                         C = 0.001, beta = beta_opts[i])

MCEM_time <- proc.time()
out <- pumpMCEM$run()
MCEM_time <- proc.time() - MCEM_time

#alpha      beta 
#0.8251551 1.2615388 
#0.8235263 1.2630080 
#0.8298914 1.2805786 

#user  system elapsed 
#1.732   0.014   1.741 
#3.707   0.057   3.759 
#1.111   0.014   1.119 
