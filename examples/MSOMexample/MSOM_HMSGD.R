# MSOM example (no Adadelta)

# From Zipkin et al. (2010) and then Ponisio et al. (2020)

setwd("~/Desktop/research/maxlik/MCMCmaxlik-dev/examples/MSOMexample/")

library(nimble)
library(dplyr)
nimbleOptions(experimentalEnableDerivs = TRUE)
library(reshape)
source("setup.R")
library("MCMCmaxlik")

survey.data <- read.csv("data/occupancy_data.csv")
species.groups <- read.csv("data/species_groups.csv")
survey.dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

model.input <- prepMutiSpData(survey.data,
                              survey.dates,
                              species.groups,
                              habitat,
                              n.zeroes=0, ## don't augment data
                              remove.zs = FALSE,
                              hyper.param = TRUE)

model.input$inits <- c(model.input$inits,
                       list(mu.ucato = 0,
                            log.sigma.ucato = 0,
                            mu.ufcw = 0,
                            log.sigma.ufcw = 0,
                            mu.vcato = 0,
                            log.sigma.vcato = 0,
                            mu.vfcw = 0,
                            log.sigma.vfcw = 0,
                            mu.a1 = 0,
                            log.sigma.a1 = 0,
                            mu.a2 = 0,
                            log.sigma.a2 = 0,
                            mu.a3 = 0,
                            log.sigma.a3 = 0,
                            mu.a4 = 0,
                            log.sigma.a4 = 0,
                            mu.b1 = 0,
                            log.sigma.b1 = 0,
                            mu.b2 = 0,
                            log.sigma.b2 = 0))

ms.ss.occ <- nimbleCode({
  ## Define prior distributions for community-level model
  ## parameters

  # Parmetrized by mu instead.
  # cato.occ.mean ~ dunif(0,1)
  #mu.ucato <- log(cato.occ.mean) - log(1-cato.occ.mean)
  mu.ucato ~ dnorm(0, 0.001)
  log.sigma.ucato ~ dnorm(0, 0.001)
  sigma.ucato <- exp(log.sigma.ucato)
  tau.ucato <-  1/(sigma.ucato*sigma.ucato)

  #fcw.occ.mean ~ dunif(0,1)
  #mu.ufcw <- log(fcw.occ.mean) - log(1-fcw.occ.mean)
  mu.ufcw ~ dnorm(0, 0.001)
  log.sigma.ufcw ~ dnorm(0, 0.001)
  sigma.ufcw <- exp(log.sigma.ufcw)
  tau.ufcw <-  1/(sigma.ufcw*sigma.ufcw)

  #cato.det.mean ~ dunif(0,1)
  #mu.vcato <- log(cato.det.mean) - log(1-cato.det.mean)
  mu.vcato ~ dnorm(0, 0.001)
  log.sigma.vcato ~ dnorm(0, 0.001)
  sigma.vcato <- exp(log.sigma.vcato)
  tau.vcato <-  1/(sigma.vcato*sigma.vcato)

  #fcw.det.mean ~ dunif(0,1)
  #mu.vfcw <- log(fcw.det.mean) - log(1-fcw.det.mean)
  mu.vfcw ~ dnorm(0, 0.001)
  log.sigma.vfcw ~ dnorm(0, 0.001)
  sigma.vfcw <- exp(log.sigma.vfcw)
  tau.vfcw <-  1/(sigma.vfcw*sigma.vfcw)

  ## random effects
  mu.a1 ~ dnorm(0, 0.001)
  log.sigma.a1 ~ dnorm(0, 0.001)
  sigma.a1 <- exp(log.sigma.a1)
  tau.a1 <-  1/(sigma.a1*sigma.a1)
  mu.a2 ~ dnorm(0, 0.001)
  log.sigma.a2 ~ dnorm(0, 0.001)
  sigma.a2 <- exp(log.sigma.a2)
  tau.a2 <-  1/(sigma.a2*sigma.a2)
  mu.a3 ~ dnorm(0, 0.001)
  log.sigma.a3 ~ dnorm(0, 0.001)
  sigma.a3 <- exp(log.sigma.a3)
  tau.a3 <-  1/(sigma.a3*sigma.a3)
  mu.a4 ~ dnorm(0, 0.001)
  log.sigma.a4 ~ dnorm(0, 0.001)
  sigma.a4 <- exp(log.sigma.a4)
  tau.a4 <-  1/(sigma.a4*sigma.a4)
  mu.b1 ~ dnorm(0, 0.001)
  log.sigma.b1 ~ dnorm(0, 0.001)
  sigma.b1 <- exp(log.sigma.b1)
  tau.b1 <-  1/(sigma.b1*sigma.b1)
  mu.b2 ~ dnorm(0, 0.001)
  log.sigma.b2 ~ dnorm(0, 0.001)
  sigma.b2 <- exp(log.sigma.b2)
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

MSOMmodel <- nimbleModel(ms.ss.occ, data = model.input$data, 
                         constants = model.input$constants, inits = model.input$inits)
cMSOMmodel <- compileNimble(MSOMmodel)

## Unfortunately the setup of Ponisio et al did not fully initialize this model
## so we use the automated initializer.
MSOMinit <- initializeModel(MSOMmodel)
cMSOMinit <- compileNimble(MSOMinit, project = MSOMmodel)
cMSOMinit$run()
cMSOMmodel$calculate()

latentNodes <- MSOMmodel$getNodeNames(latentOnly = TRUE, stochOnly = TRUE, includeData = FALSE)
MSOMParamNodes <- MSOMmodel$getNodeNames(topOnly = TRUE, includeData = FALSE)
compiledFunsMSOM <- buildMCMCmaxlik(MSOMmodel, MSOMParamNodes)

init <- unlist(list(mu.ucato = 0,
                    log.sigma.ucato = 0,
                    mu.ufcw = 0,
                    log.sigma.ufcw = 0,
                    mu.vcato = 0,
                    log.sigma.vcato = 0,
                    mu.vfcw = 0,
                    log.sigma.vfcw = 0,
                    mu.a1 = 0,
                    log.sigma.a1 = 0,
                    mu.a2 = 0,
                    log.sigma.a2 = 0,
                    mu.a3 = 0,
                    log.sigma.a3 = 0,
                    mu.a4 = 0,
                    log.sigma.a4 = 0,
                    mu.b1 = 0,
                    log.sigma.b1 = 0,
                    mu.b2 = 0,
                    log.sigma.b2 = 0))
boundary <- list(c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40))

numMCMCSamples <- 300
maxIter <- 300

set.seed(1500)

getMSOMResults <- function(init, 
                           boundary=boundary, 
                           numMCMCSamples) {
  # 1. Fixed step size ----------------------------------------
  resultsMSOMFixed <- computeMLE(MSOMmodel, MSOMParamNodes,
                                 method="fixed", paramInit=init,
                                 stepsize=0.05,
                                 compiledFuns=compiledFunsMSOM,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=maxIter,
                                 boundary=boundary)
  
  # 2. Small fixed step size ----------------------------------------
  resultsMSOMFixedSmall <- computeMLE(MSOMmodel, MSOMParamNodes,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsMSOM,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=maxIter,
                                      boundary=boundary)
  
  # 4. Adam ----------------------------------------
  resultsMSOMAdam <- computeMLE(MSOMmodel, MSOMParamNodes,
                            method="adam", paramInit=init,
                            compiledFuns=compiledFunsMSOM,
                            numMCMCSamples=numMCMCSamples,
                            maxIter=maxIter,
                            boundary=boundary)
  
  # 5. Newton-Raphson ----------------------------------------
  resultsMSOMNR <- computeMLE(MSOMmodel, MSOMParamNodes,
                          method="NR", paramInit=init,
                          compiledFuns=compiledFunsMSOM,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=maxIter,
                          boundary=boundary)
  
  # 6. 1-D sampling ----------------------------------------
  # resultsMSOM1D <- computeMLE(MSOMmodel, MSOMParamNodes,
  #                         method="ga1D", paramInit=init,
  #                         compiledFuns=compiledFunsMSOM,
  #                         numMCMCSamples=numMCMCSamples,
  #                         maxIter=maxIter, 
  #                         boundary=boundary)
  
  list(resultsMSOMFixed, resultsMSOMFixedSmall,
       resultsMSOMAdam, resultsMSOMNR#, resultsMSOM1D
       )
}

resultsMSOM <- list(resultsMSOMFixed, resultsMSOMFixedSmall,
     resultsMSOMAdam, resultsMSOMNR#, resultsMSOM1D
)

saveRDS(resultsMSOM, "resultsMSOM_300")

scale_conversion <- function(x) {
  sapply(1:length(x), function(i) {
    ifelse(i %% 2, ifelse(i >= 9, x[i], exp(x[i]) / (1 + exp(x[i]))), exp(x[i]))
  })
}

resultsMSOM <- readRDS("resultsMSOM_300")

df_MSOM <- data.frame(t(sapply(1:length(resultsMSOM), function(i) {
  scale_conversion(resultsMSOM[[i]]$MLE)
})))
names(df_MSOM) <- c("cato.occ.mean", "sigma.ucato", "fcw.occ.mean", "sigma.ufcw",
                    "cato.det.mean", "sigma.vcato", "fcw.det.mean", "sigma.vfcw",
                    "mu.a1", "sigma.a1", "mu.a2", "sigma.a2",
                    "mu.a3", "sigma.a3", "mu.a4", "sigma.a4",
                    "mu.b1", "sigma.b1", "mu.b2", "sigma.b2")
df_MSOM$method <- c('Fixed (0.05)', 'Fixed (0.005)', 'Adam', "NR")
df_ae <- df_MSOM %>% group_by(method) %>%
  summarize(RMSE.cato.occ.mean=sqrt(mean((cato.occ.mean - 0.353384378)^2)),
            RMSE.sigma.ucato=sqrt(mean((sigma.ucato - 3.258728858)^2)),
            RMSE.fcw.occ.mean=sqrt(mean((fcw.occ.mean - 0.464779842)^2)),
            RMSE.sigma.ufcw=sqrt(mean((sigma.ufcw - 2.727642539)^2)),
            RMSE.cato.det.mean=sqrt(mean((cato.det.mean - 0.186099275)^2)),
            RMSE.sigma.vcato=sqrt(mean((sigma.vcato - 1.328024673)^2)),
            RMSE.fcw.det.mean=sqrt(mean((fcw.det.mean - 0.248226209)^2)),
            RMSE.sigma.vfcw=sqrt(mean((sigma.vfcw - 1.192028288)^2)),
            RMSE.mu.a1=sqrt(mean((mu.a1 - 0.456863874)^2)),
            RMSE.sigma.a1=sqrt(mean((sigma.a1 - 0.637116252)^2)),
            RMSE.mu.a2=sqrt(mean((mu.a2 - 0.005105438)^2)),
            RMSE.sigma.a2=sqrt(mean((sigma.a2 - 0.242693905)^2)),
            RMSE.mu.a3=sqrt(mean((mu.a3 - (-0.154347721))^2)),
            RMSE.sigma.a3=sqrt(mean((sigma.a3 - 0.171261355)^2)),
            RMSE.mu.a4=sqrt(mean((mu.a4 - 0.127042312)^2)),
            RMSE.sigma.a4=sqrt(mean((sigma.a4 - 0.251311617)^2)),
            RMSE.mu.b1=sqrt(mean((mu.b1 - (-0.132541853))^2)),
            RMSE.sigma.b1=sqrt(mean((sigma.b1 - 0.200357512)^2)),
            RMSE.mu.b2=sqrt(mean((mu.b2 - 0.085964990)^2)),
            RMSE.sigma.b2=sqrt(mean((sigma.b2 - 0.047909323)^2))
  )