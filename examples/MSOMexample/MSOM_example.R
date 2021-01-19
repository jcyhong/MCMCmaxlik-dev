# From Zipkin et al. (2010) and then Ponisio et al. (2020)

library(nimble)
nimbleOptions(experimentalEnableDerivs = TRUE)
library(reshape)
source("setup.R")

survey.data <- read.csv("data/occupancy_data.csv")
species.groups <- read.csv("data/species_groups.csv")
survey.dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

model.input <- prepMutiSpData(survey.data,
                              survey.dates,
                              species.groups,
                              habitat,
                              n.zeroes =0, ## don't augment data
                              remove.zs = FALSE,
                              hyper.param = TRUE)

model.input$inits <- c(model.input$inits,
                       list(cato.occ.mean = 0.5,
                            sigma.ucato = 1,
                            fcw.occ.mean = 0.5,
                            sigma.ufcw = 1,
                            cato.det.mean = 0.5,
                            sigma.vcato = 1,
                            fcw.det.mean = 0.5,
                            sigma.vfcw = 1,
                            mu.a1 = 0,
                            sigma.a1 = 1,
                            mu.a2 = 0,
                            sigma.a2 = 1,
                            mu.a3 = 0,
                            sigma.a3 = 1,
                            mu.a4 = 0,
                            sigma.a4 = 1,
                            mu.b1 = 0,
                            sigma.b1 = 1,
                            mu.b2 = 0,
                            sigma.b2 = 1))
                            
ms.ss.occ <- nimbleCode({
  ## Define prior distributions for community-level model
  ## parameters
  cato.occ.mean ~ dunif(0,1)
  mu.ucato <- log(cato.occ.mean) - log(1-cato.occ.mean)
  sigma.ucato ~ dunif(0, 100)
  tau.ucato <-  1/(sigma.ucato*sigma.ucato)

  fcw.occ.mean ~ dunif(0,1)
  mu.ufcw <- log(fcw.occ.mean) - log(1-fcw.occ.mean)
  sigma.ufcw ~ dunif(0, 100)
  tau.ufcw <-  1/(sigma.ufcw*sigma.ufcw)

  cato.det.mean ~ dunif(0,1)
  mu.vcato <- log(cato.det.mean) - log(1-cato.det.mean)
  sigma.vcato ~ dunif(0, 100)
  tau.vcato <-  1/(sigma.vcato*sigma.vcato)

  fcw.det.mean ~ dunif(0,1)
  mu.vfcw <- log(fcw.det.mean) - log(1-fcw.det.mean)
  sigma.vfcw ~ dunif(0, 100)
  tau.vfcw <-  1/(sigma.vfcw*sigma.vfcw)

  ## random effects
  mu.a1 ~ dnorm(0, 0.001)
  sigma.a1 ~ dunif(0, 100)
  tau.a1 <-  1/(sigma.a1*sigma.a1)
  mu.a2 ~ dnorm(0, 0.001)
  sigma.a2 ~ dunif(0, 100)
  tau.a2 <-  1/(sigma.a2*sigma.a2)
  mu.a3 ~ dnorm(0, 0.001)
  sigma.a3 ~ dunif(0, 100)
  tau.a3 <-  1/(sigma.a3*sigma.a3)
  mu.a4 ~ dnorm(0, 0.001)
  sigma.a4 ~ dunif(0, 100)
  tau.a4 <-  1/(sigma.a4*sigma.a4)
  mu.b1 ~ dnorm(0, 0.001)
  sigma.b1 ~ dunif(0, 100)
  tau.b1 <-  1/(sigma.b1*sigma.b1)
  mu.b2 ~ dnorm(0, 0.001)
  sigma.b2 ~ dunif(0, 100)
  tau.b2 <-  1/(sigma.b2*sigma.b2)

  for (i in 1:(num.species)) {
    ## Create priors for species i from the community level prior
    ## distributions
    u.cato[i] ~ dnorm(mu.ucato, tau.ucato)
    u.fcw[i] ~ dnorm(mu.ufcw, tau.ufcw)
    a1[i] ~ dnorm(mu.a1, tau.a1)
    a2[i] ~ dnorm(mu.a2, tau.a2)
    a3[i] ~ dnorm(mu.a3, tau.a3)
    a4[i] ~ dnorm(mu.a4, tau.a4)
    v.cato[i] ~ dnorm(mu.vcato, tau.vcato)
    v.fcw[i] ~ dnorm(mu.vfcw, tau.vfcw)
    b1[i] ~ dnorm(mu.b1, tau.b1)
    b2[i] ~ dnorm(mu.b2, tau.b2)
    ## Create a loop to estimate the Z matrix (true occurrence for
    ## species i at point j).

    for (j in 1:num.points) {
      logit(psi[j,i]) <- u.cato[i]*(1-habitat.ind[j]) +
        u.fcw[i]*habitat.ind[j] +
        a1[i]*ufc.linear[j] +
        a2[i]*ufc.quadratic[j] +
        a3[i]*ba.linear[j] +
        a4[i]*ba.quadratic[j]
      mu.psi[j,i] <- psi[j,i]
      Z[j,i] ~ dbern(mu.psi[j,i])
      ## Create a loop to estimate detection for species i at point k
      ## during sampling period k.
      for (k in 1:num.reps[j]) {
        logit(p[j,k,i]) <-  v.cato[i]*(1-habitat.ind[j]) +
          v.fcw[i]*habitat.ind[j] +
          b1[i]*date.linear[j,k] +
          b2[i]*date.quadratic[j,k]
        mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
        X[j,k,i] ~ dbern(mu.p[j,k,i])
      }
    }
  }
})

MSOMmodel <- nimbleModel(ms.ss.occ, data = model.input$data, constants = model.input$constants, inits = model.input$inits)
cMSOMmodel <- compileNimble(MSOMmodel)

## Unfortunately the setup of Ponisio et al did not fully initialize this model
## so we use the automated initializer.
MSOMinit <- initializeModel(MSOMmodel)
cMSOMinit <- compileNimble(MSOMinit, project = MSOMmodel)
cMSOMinit$run()
cMSOMmodel$calculate()

latentNodes <- MSOMmodel$getNodeNames(latentOnly = TRUE, stochOnly = TRUE, includeData = FALSE)

source("MCEM_AD_build.R")


MSOMmcem <- buildMCEM_AD(MSOMmodel, latentNodes = latentNodes,
                         C = 0.05, alpha = .01, beta = .01, gamma = .01)
ptm <- proc.time()
MSOMmle <- MSOMmcem$run(initM = 1000)
time_taken <- proc.time() - ptm
time_taken 

