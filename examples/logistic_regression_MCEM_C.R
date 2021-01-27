# Logistic regression example (3 parameters)

# Load the packages ---------------------------------------
library("MCMCmaxlik")
library(nimble)
nimbleOptions(experimentalEnableDerivs = TRUE)

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

source("~/Desktop/MCMCmaxlik-dev/packages/MCMCmaxlik/R/MCEM_AD_build.R")
C_opts = c(0.001, 0.01, 0.1)


i <- 3
LogregMCEM <- buildMCEM_AD(
  model = logreg, latentNodes = "beta2",
  boxConstraints = list(
    list(c("sigma_RE"),
      limits = c(0, 1000)
    )
  ),
  C = C_opts[i]
)

MCEM_time <- proc.time()
out <- LogregMCEM$run()
MCEM_time <- proc.time() - MCEM_time


#beta0      beta1   sigma_RE 
#-0.5479683  1.3180001  0.2498723 
#-0.5413166  1.3136490  0.2616243 
#-0.5525307  1.3124394  0.2842736 

#user  system elapsed 
#12.136   0.081  12.227 
# 2.109   0.038   2.136 
# 0.763   0.020   0.779 