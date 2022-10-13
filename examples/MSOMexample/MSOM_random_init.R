# Supplementary materials for "Sampling-Based Approaches to Maximum 
# Likelihood Estimation for Latent Variable Models":
# Code for the GLMM example (6 parameters)
#
# Description:
# The following code contains the numerical experiments for the GLMM example: 
# fixed step size, small fixed step size, adadelta, adam, Newton-Raphson, 
# and 1-D sampling. MCEM (from the R package NIMBLE) is used as a benchmark.
#
# Background/motivation:
# Excerpts from http://www.jstor.org/stable/pdf/2532317.pdf:
# "The scientific question addressed in the study is whether geographically 
# isolated populations of salamanders develop barriers to successful mating. 
# Females from two distinct populations called whiteside (W) and rough-butt 
# (R) were paired for mating to males from their own and from the other 
# population...  A second question is whether there exists heterogeneity 
# among individuals in mating success and, if so, whether it is greater for 
# males or females."
# 
#
##########################################################################

# Load the packages ---------------------------------------
library("MCMCmaxlik")
require(nimble)
library(parallel)
library(data.table)
library(reshape)
library(dplyr)
this_cluster <- makeCluster(6)

# Fix the seed for reproducibility.
set.seed(12301)

survey.data <- read.csv("data/occupancy_data.csv")
species.groups <- read.csv("data/species_groups.csv")
survey.dates <- read.csv("data/survey_dates.csv")
habitat <- read.csv("data/habitat.csv")

run_allcode <- function(b, survey.data, species.groups, survey.dates, habitat) {
  # MSOM example
  
  # From Zipkin et al. (2010) and then Ponisio et al. (2020)
  
  library(nimble)
  library("MCMCmaxlik")
  nimbleOptions(experimentalEnableDerivs = TRUE)
  library(reshape)
  createEncounterHistory <- function(survey.data){
    ## The detection/non-detection data is reshaped into a three
    ## dimensional array X where the first dimension, j, is the point;
    ## the second dimension, k, is the rep; and the last dimension, i, is
    ## the species.
    survey.data$occupancy <- rep(1, dim(survey.data)[1])
    
    X = melt(survey.data,
             id.var = c("species", "point", "repetition"),
             measure.var = "occupancy")
    
    X = cast(X, point ~ repetition ~ species)
    
    ## Change abundance to presence/absence
    ## Set counts 1 and greater to indicators
    X[which(X > 0)] <- 1
    
    ## Create all zero encounter histories to add to the detection array X
    ## as part of the data augmentation to account for additional
    ## species (beyond the n observed species).
    X.zero = matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
    
    return(list(X=X,X.zero=X.zero))
  }
  
  addMissingData <- function(histories, survey.dates){
    ## 'Add' missing data:set X and X.zero for the unsurveyed
    ## repetitions to NA
    X <- histories$X
    X.zero <- histories$X.zero
    for(point in 1:length(unique(survey.data$point))){
      point.index <- which(survey.dates$point == row.names(X)[point])
      missing <- is.na(survey.dates[point.index,][,-(1:1)])
      X[point, missing,] <- NA
      X.zero[point, missing] <- NA
    }
    return(list(X=X,X.zero=X.zero))
  }
  
  zeroAugment <- function(histories, n.zeroes){
    ## X.aug is the augmented version of X.  The first n species were
    ## actually observed and the n+1 through n.zeroes species are all
    ## zero encounter histories create an empty 3D array with n + n.zero
    ## species
    X <- histories$X
    X.zero <- histories$X.zero
    X.dim <- dim(X)
    X.dim[3] <- X.dim[3] + n.zeroes
    X.aug <- array( NA, dim = X.dim)
    
    ## fill in the array with the occurrence data
    X.aug[,,1:dim(X)[3]] <-  X
    
    ## fill the zero histories
    X.aug[,,-(1:dim(X)[3])] <- rep(X.zero, n.zeroes)
    return(X.aug)
  }
  
  siteLevelStandardized <- function(parameter){
    m <- mean(parameter, na.rm = TRUE)
    sd <- sd(parameter, na.rm = TRUE)
    linear <- (parameter - m)/sd
    quadratic <- linear * linear
    return(list(linear=linear,quadratic=quadratic))
  }
  
  surveyLevelStandardized <- function(parameter){
    parameter <- as.matrix(parameter)
    m <- mean(parameter, na.rm = TRUE)
    sd <- sd(parameter, na.rm = TRUE)
    linear <- (parameter - m) /  sd
    quadratic <- linear * linear
    linear <- as.matrix(linear)
    quadratic <- as.matrix(quadratic)
    return(list(linear=linear,quadratic=quadratic))
  }
  
  
  
  reformatData <- function(survey.data,
                           survey.dates,
                           species.groups,
                           habitat,
                           n.zeroes){
    
    manipulateData <- function(survey.data,
                               survey.dates,
                               species.groups,
                               habitat, n.zeroes){
      histories <- createEncounterHistory(survey.data)
      histories <- addMissingData(histories, survey.dates)
      X.aug <- zeroAugment(histories, n.zeroes)
      
      num.species <- length(unique(survey.data$species))
      num.points <- length(unique(survey.data$point))
      num.reps <- apply(survey.dates[,-(1:1)],
                        1,
                        function(x) length(which(!is.na(x))))
      
      ## Create an indicator vector for each assemblage (ground, mid-story)
      ground <- mid <- rep(0, dim(species.groups)[1])
      ground[which(species.groups$group == 1)] <- 1
      mid[which(species.groups$group == 2)] <- 1
      
      ## Create a vector to indicate which habitat type each point is in
      ## (CATO = 1; FCW  = 0)
      habitat.ind <- rep(0, num.points)
      habitat.ind[grep("CAT", row.names(histories$X))] <- 1
      
      ## Standardize variables
      ufc <-  siteLevelStandardized(habitat$ufc)
      ufc.linear <- ufc$linear
      ufc.quadratic <- ufc$quadratic
      
      ba <- siteLevelStandardized(habitat$ba)
      ba.linear <- ba$linear
      ba.quadratic <- ba$quadratic
      
      date <- surveyLevelStandardized(survey.dates[,-(1:1)])
      date.linear <- date$linear
      date.quadratic <- date$quadratic
      
      return(list(X.aug = X.aug, num.species = num.species,
                  num.points = num.points,
                  num.reps = num.reps, ground = ground, mid = mid,
                  habitat.ind = habitat.ind, ufc.linear = ufc.linear,
                  ufc.quadratic = ufc.quadratic, ba.linear = ba.linear,
                  ba.quadratic = ba.quadratic, date.linear = date.linear,
                  date.quadratic = date.quadratic))
    }
    
    
    return(manipulateData(survey.data,
                          survey.dates,
                          species.groups,
                          habitat,
                          n.zeroes))
  }
  
  
  ## prep data for model
  prepMutiSpData <- function(survey.data,
                             survey.dates,
                             species.groups,
                             habitat,
                             n.zeroes,
                             remove.zs=TRUE,
                             vectorized=TRUE,
                             hyper.param=TRUE){
    ## reformat data
    data <- reformatData(survey.data,
                         survey.dates,
                         species.groups,
                         habitat,
                         n.zeroes)
    
    num.species <- data$num.species
    num.points <- data$num.points
    num.reps <- data$num.reps
    
    ## Z data for whether or not a species was ever observed
    ## zs with 1s as 1s and 0s as NAs
    zs <- apply(data$X, c(1, 3), max, na.rm=TRUE)
    zs[zs == 0] <- NA
    zs[!is.finite(zs)] <- NA
    
    model.data <- list(Z = zs,
                       X = data$X.aug,
                       ground = data$ground,
                       mid = data$mid,
                       habitat.ind = data$habitat.ind,
                       ufc.linear = data$ufc.linear,
                       ufc.quadratic = data$ufc.quadratic,
                       ba.linear = data$ba.linear,
                       ba.quadratic = data$ba.quadratic,
                       date.linear = data$date.linear,
                       date.quadratic = data$date.quadratic)
    
    psi.mean.draw <- runif(1, 0.25, 1)
    
    ## initial values
    omega.draw <- runif(1, num.species/(num.species + n.zeroes), 1)
    
    ## inital conditions. 1 should be NA, NA should be a 0 or 1
    zinits <- zs
    zinits[zinits == 1] <- 2
    ## zinits[is.na(zinits)] <- 1
    zinits[is.na(zinits)] <- sample(0:1, sum(is.na(zinits)),
                                    replace=TRUE)
    zinits[zinits == 2] <- NA
    
    if(hyper.param){
      inits <-list(Z=zinits,
                   omega = omega.draw,
                   w = c(rep(1, num.species),
                         rbinom(n.zeroes, size = 1, prob = omega.draw)),
                   u.cato = rnorm(num.species + n.zeroes),
                   v.cato = rnorm(num.species + n.zeroes),
                   u.fcw = rnorm(num.species + n.zeroes) ,
                   v.fcw = rnorm(num.species + n.zeroes),
                   a1 = rnorm(num.species + n.zeroes),
                   a2 = rnorm(num.species + n.zeroes),
                   a3 = rnorm(num.species + n.zeroes),
                   a4 = rnorm(num.species + n.zeroes),
                   b1 = rnorm(num.species + n.zeroes),
                   b2 = rnorm(num.species + n.zeroes))
    } else{
      inits <-list(Z=zinits,
                   omega = omega.draw,
                   w = c(rep(1, num.species),
                         rbinom(n.zeroes, size = 1, prob = omega.draw)),
                   u.cato = rnorm(1),
                   v.cato = rnorm(1),
                   u.fcw = rnorm(1) ,
                   v.fcw = rnorm(1),
                   a1 = rnorm(1),
                   a2 = rnorm(1),
                   a3 = rnorm(1),
                   a4 = rnorm(1),
                   b1 = rnorm(1),
                   b2 = rnorm(1))
    }
    
    ## constants
    constants <- list(num.species = num.species,
                      num.points = num.points,
                      num.reps = num.reps,
                      n.zeroes = n.zeroes)
    
    ## for the non-data augmented case
    if(n.zeroes == 0){
      inits[c("w", "omega")] <- NULL
      constants[c("n.zeroes")] <- NULL
      
    }
    
    ## these were removed from the model
    model.data[["ground"]] <- NULL
    model.data[["mid"]] <- NULL
    
    if(vectorized){
      ## Since calculations with date_linear and date_quadratic are now
      ## vectorized, we'll set the NAs to 0
      model.data$date.linear[is.na(model.data$date.linear)] <- 0
      model.data$date.quadratic[is.na(model.data$date.quadratic)] <- 0
    }
    if(remove.zs) {
      ## zs are removed from these models
      model.data[["Z"]] <- NULL
      inits[["Z"]] <- NULL
      ## additional constants and dats for models where z is removed
      constants$max.num.reps <- max(constants$num.reps)
      model.data$onesRow <- matrix(rep(1, constants$max.num.reps),
                                   nrow=1)
    }
    return(list(constants=constants,
                inits=inits,
                data=model.data))
  }
  
  ## new function for generating initial values for use in nimbleModel(),
  ## using posterior summary statistics
  ## -DT
  genInits <- function(summary) {
    e <- new.env()
    for(n in names(summary))    eval(parse(text = paste0('e$', n, ' <- summary[\'', n, '\']'))[[1]])
    new_inits <- list()
    for(n in ls(e)) new_inits[[n]] <- e[[n]]
    return(new_inits)
  }
  
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
  
  MSOMmodel <- nimbleModel(ms.ss.occ, data = model.input$data, constants = model.input$constants, inits = model.input$inits)
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
  
  init <- runif(20, -1, 1)

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
  maxIter <- 600
  
  # Fixed step size ----------------------------------------
  resultsMSOMFixed <- computeMLE(MSOMmodel, MSOMParamNodes,
                                 method="fixed", paramInit=init,
                                 stepsize=0.05,
                                 compiledFuns=compiledFunsMSOM,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=maxIter,
                                 boundary=boundary)
  
  # Small fixed step size ----------------------------------------
  resultsMSOMFixedSmall <- computeMLE(MSOMmodel, MSOMParamNodes,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsMSOM,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=maxIter,
                                      boundary=boundary)
  
  # Adam ----------------------------------------
  resultsMSOMAdam <- computeMLE(MSOMmodel, MSOMParamNodes,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsMSOM,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=maxIter,
                                boundary=boundary)
  
  scale_conversion <- function(x) {
    sapply(1:length(x), function(i) {
      ifelse(i %% 2, ifelse(i >= 9, x[i], exp(x[i]) / (1 + exp(x[i]))), exp(x[i]))
    })
  }
  
  df <- data.frame(do.call(rbind, list(scale_conversion(resultsMSOMFixed$MLE),
                                       scale_conversion(resultsMSOMFixedSmall$MLE),
                                       scale_conversion(resultsMSOMAdam$MLE))))
  df$method <- c('Fixed (0.05)', 'Fixed (0.005)', 'Adam')
  
  names(df) <- c("cato.occ.mean", "sigma.ucato", "fcw.occ.mean", "sigma.ufcw",
                 "cato.det.mean", "sigma.vcato", "fcw.det.mean", "sigma.vfcw",
                 "mu.a1", "sigma.a1", "mu.a2", "sigma.a2",
                 "mu.a3", "sigma.a3", "mu.a4", "sigma.a4",
                 "mu.b1", "sigma.b1", "mu.b2", "sigma.b2",
                 "method")
  list(df, init, resultsMSOMFixed, resultsMSOMFixedSmall, resultsMSOMAdam)
}


all_results <- parLapply(cl = this_cluster, X = 1:30, 
                         fun = run_allcode, 
                         survey.data = survey.data,
                         species.groups = species.groups, 
                         survey.dates = survey.dates, 
                         habitat = habitat)

stopCluster(this_cluster)

df_list_all <- lapply(all_results, function(result) {result[[1]]})
df_combined_all <- rbindlist(df_list_all)
df_rmse <- df_combined_all %>% group_by(method) %>%
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

write.csv(df_rmse, file='MSOM_rmse_300_600.csv', row.names=F)


# MCEM results ----

run_MCEM <- function(b, survey.data, species.groups, survey.dates, habitat) {
  # MSOM example
  
  # From Zipkin et al. (2010) and then Ponisio et al. (2020)
  
  library(nimble)
  library("MCMCmaxlik")
  nimbleOptions(experimentalEnableDerivs = TRUE)
  library(reshape)
  createEncounterHistory <- function(survey.data){
    ## The detection/non-detection data is reshaped into a three
    ## dimensional array X where the first dimension, j, is the point;
    ## the second dimension, k, is the rep; and the last dimension, i, is
    ## the species.
    survey.data$occupancy <- rep(1, dim(survey.data)[1])
    
    X = melt(survey.data,
             id.var = c("species", "point", "repetition"),
             measure.var = "occupancy")
    
    X = cast(X, point ~ repetition ~ species)
    
    ## Change abundance to presence/absence
    ## Set counts 1 and greater to indicators
    X[which(X > 0)] <- 1
    
    ## Create all zero encounter histories to add to the detection array X
    ## as part of the data augmentation to account for additional
    ## species (beyond the n observed species).
    X.zero = matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
    
    return(list(X=X,X.zero=X.zero))
  }
  
  addMissingData <- function(histories, survey.dates){
    ## 'Add' missing data:set X and X.zero for the unsurveyed
    ## repetitions to NA
    X <- histories$X
    X.zero <- histories$X.zero
    for(point in 1:length(unique(survey.data$point))){
      point.index <- which(survey.dates$point == row.names(X)[point])
      missing <- is.na(survey.dates[point.index,][,-(1:1)])
      X[point, missing,] <- NA
      X.zero[point, missing] <- NA
    }
    return(list(X=X,X.zero=X.zero))
  }
  
  zeroAugment <- function(histories, n.zeroes){
    ## X.aug is the augmented version of X.  The first n species were
    ## actually observed and the n+1 through n.zeroes species are all
    ## zero encounter histories create an empty 3D array with n + n.zero
    ## species
    X <- histories$X
    X.zero <- histories$X.zero
    X.dim <- dim(X)
    X.dim[3] <- X.dim[3] + n.zeroes
    X.aug <- array( NA, dim = X.dim)
    
    ## fill in the array with the occurrence data
    X.aug[,,1:dim(X)[3]] <-  X
    
    ## fill the zero histories
    X.aug[,,-(1:dim(X)[3])] <- rep(X.zero, n.zeroes)
    return(X.aug)
  }
  
  siteLevelStandardized <- function(parameter){
    m <- mean(parameter, na.rm = TRUE)
    sd <- sd(parameter, na.rm = TRUE)
    linear <- (parameter - m)/sd
    quadratic <- linear * linear
    return(list(linear=linear,quadratic=quadratic))
  }
  
  surveyLevelStandardized <- function(parameter){
    parameter <- as.matrix(parameter)
    m <- mean(parameter, na.rm = TRUE)
    sd <- sd(parameter, na.rm = TRUE)
    linear <- (parameter - m) /  sd
    quadratic <- linear * linear
    linear <- as.matrix(linear)
    quadratic <- as.matrix(quadratic)
    return(list(linear=linear,quadratic=quadratic))
  }
  
  
  
  reformatData <- function(survey.data,
                           survey.dates,
                           species.groups,
                           habitat,
                           n.zeroes){
    
    manipulateData <- function(survey.data,
                               survey.dates,
                               species.groups,
                               habitat, n.zeroes){
      histories <- createEncounterHistory(survey.data)
      histories <- addMissingData(histories, survey.dates)
      X.aug <- zeroAugment(histories, n.zeroes)
      
      num.species <- length(unique(survey.data$species))
      num.points <- length(unique(survey.data$point))
      num.reps <- apply(survey.dates[,-(1:1)],
                        1,
                        function(x) length(which(!is.na(x))))
      
      ## Create an indicator vector for each assemblage (ground, mid-story)
      ground <- mid <- rep(0, dim(species.groups)[1])
      ground[which(species.groups$group == 1)] <- 1
      mid[which(species.groups$group == 2)] <- 1
      
      ## Create a vector to indicate which habitat type each point is in
      ## (CATO = 1; FCW  = 0)
      habitat.ind <- rep(0, num.points)
      habitat.ind[grep("CAT", row.names(histories$X))] <- 1
      
      ## Standardize variables
      ufc <-  siteLevelStandardized(habitat$ufc)
      ufc.linear <- ufc$linear
      ufc.quadratic <- ufc$quadratic
      
      ba <- siteLevelStandardized(habitat$ba)
      ba.linear <- ba$linear
      ba.quadratic <- ba$quadratic
      
      date <- surveyLevelStandardized(survey.dates[,-(1:1)])
      date.linear <- date$linear
      date.quadratic <- date$quadratic
      
      return(list(X.aug = X.aug, num.species = num.species,
                  num.points = num.points,
                  num.reps = num.reps, ground = ground, mid = mid,
                  habitat.ind = habitat.ind, ufc.linear = ufc.linear,
                  ufc.quadratic = ufc.quadratic, ba.linear = ba.linear,
                  ba.quadratic = ba.quadratic, date.linear = date.linear,
                  date.quadratic = date.quadratic))
    }
    
    
    return(manipulateData(survey.data,
                          survey.dates,
                          species.groups,
                          habitat,
                          n.zeroes))
  }
  
  
  ## prep data for model
  prepMutiSpData <- function(survey.data,
                             survey.dates,
                             species.groups,
                             habitat,
                             n.zeroes,
                             remove.zs=TRUE,
                             vectorized=TRUE,
                             hyper.param=TRUE){
    ## reformat data
    data <- reformatData(survey.data,
                         survey.dates,
                         species.groups,
                         habitat,
                         n.zeroes)
    
    num.species <- data$num.species
    num.points <- data$num.points
    num.reps <- data$num.reps
    
    ## Z data for whether or not a species was ever observed
    ## zs with 1s as 1s and 0s as NAs
    zs <- apply(data$X, c(1, 3), max, na.rm=TRUE)
    zs[zs == 0] <- NA
    zs[!is.finite(zs)] <- NA
    
    model.data <- list(Z = zs,
                       X = data$X.aug,
                       ground = data$ground,
                       mid = data$mid,
                       habitat.ind = data$habitat.ind,
                       ufc.linear = data$ufc.linear,
                       ufc.quadratic = data$ufc.quadratic,
                       ba.linear = data$ba.linear,
                       ba.quadratic = data$ba.quadratic,
                       date.linear = data$date.linear,
                       date.quadratic = data$date.quadratic)
    
    psi.mean.draw <- runif(1, 0.25, 1)
    
    ## initial values
    omega.draw <- runif(1, num.species/(num.species + n.zeroes), 1)
    
    ## inital conditions. 1 should be NA, NA should be a 0 or 1
    zinits <- zs
    zinits[zinits == 1] <- 2
    ## zinits[is.na(zinits)] <- 1
    zinits[is.na(zinits)] <- sample(0:1, sum(is.na(zinits)),
                                    replace=TRUE)
    zinits[zinits == 2] <- NA
    
    if(hyper.param){
      inits <-list(Z=zinits,
                   omega = omega.draw,
                   w = c(rep(1, num.species),
                         rbinom(n.zeroes, size = 1, prob = omega.draw)),
                   u.cato = rnorm(num.species + n.zeroes),
                   v.cato = rnorm(num.species + n.zeroes),
                   u.fcw = rnorm(num.species + n.zeroes) ,
                   v.fcw = rnorm(num.species + n.zeroes),
                   a1 = rnorm(num.species + n.zeroes),
                   a2 = rnorm(num.species + n.zeroes),
                   a3 = rnorm(num.species + n.zeroes),
                   a4 = rnorm(num.species + n.zeroes),
                   b1 = rnorm(num.species + n.zeroes),
                   b2 = rnorm(num.species + n.zeroes))
    } else{
      inits <-list(Z=zinits,
                   omega = omega.draw,
                   w = c(rep(1, num.species),
                         rbinom(n.zeroes, size = 1, prob = omega.draw)),
                   u.cato = rnorm(1),
                   v.cato = rnorm(1),
                   u.fcw = rnorm(1) ,
                   v.fcw = rnorm(1),
                   a1 = rnorm(1),
                   a2 = rnorm(1),
                   a3 = rnorm(1),
                   a4 = rnorm(1),
                   b1 = rnorm(1),
                   b2 = rnorm(1))
    }
    
    ## constants
    constants <- list(num.species = num.species,
                      num.points = num.points,
                      num.reps = num.reps,
                      n.zeroes = n.zeroes)
    
    ## for the non-data augmented case
    if(n.zeroes == 0){
      inits[c("w", "omega")] <- NULL
      constants[c("n.zeroes")] <- NULL
      
    }
    
    ## these were removed from the model
    model.data[["ground"]] <- NULL
    model.data[["mid"]] <- NULL
    
    if(vectorized){
      ## Since calculations with date_linear and date_quadratic are now
      ## vectorized, we'll set the NAs to 0
      model.data$date.linear[is.na(model.data$date.linear)] <- 0
      model.data$date.quadratic[is.na(model.data$date.quadratic)] <- 0
    }
    if(remove.zs) {
      ## zs are removed from these models
      model.data[["Z"]] <- NULL
      inits[["Z"]] <- NULL
      ## additional constants and dats for models where z is removed
      constants$max.num.reps <- max(constants$num.reps)
      model.data$onesRow <- matrix(rep(1, constants$max.num.reps),
                                   nrow=1)
    }
    return(list(constants=constants,
                inits=inits,
                data=model.data))
  }
  
  ## new function for generating initial values for use in nimbleModel(),
  ## using posterior summary statistics
  ## -DT
  genInits <- function(summary) {
    e <- new.env()
    for(n in names(summary))    eval(parse(text = paste0('e$', n, ' <- summary[\'', n, '\']'))[[1]])
    new_inits <- list()
    for(n in ls(e)) new_inits[[n]] <- e[[n]]
    return(new_inits)
  }
  
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
  
  MSOMmodel <- nimbleModel(ms.ss.occ, data = model.input$data, constants = model.input$constants, inits = model.input$inits)
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
  
  ## Working on update to MCEM so Johnny can use it with AD
  
  ## Notes for future reference
  ## getMCEMRanges should not be a nimbleFunction
  ## We could check that maxNodes are actually connected to latentNodes instead of including all non-data stochastic nodes that are not latentNodes
  ## We could use the parameterTransform system to avoid box constraints.
  ## mcmc_Latent_Conf should not have to monitor all var names
  ## burnIn can be a run-time parameter
  ## In calc_E_llk_gen: paramValues are repeatedly put in model even when diff is FALSE.
  ##                    fixedCalcNodes are never used
  ##                    need for paramDepDetermNodes_latent is not clear
  ## The MBB is really wasteful in re-calculating values.
  ## 1/nSamples should be 1/(nSamples - burnin), but the factor is different only be a constant so might not matter.
  ##    It will however matter for the std. err. 
  ## Pathological behavior: if you start at the MLE, it spins the Monte Carlo sample size up enormously to try to get a confident improvement in likelihood, when there is none possible.
  
  # Calculates Q function if diff = 0, calculates difference in Q functions if diff = 1.
  calc_E_llk_gen_AD = nimbleFunction(
    name = 'calc_E_llk_gen_AD',
    setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0){
      ## fixedCalcNodes <- model$getDependencies(fixedNodes)	
      latentCalcNodes <- model$getDependencies(sampledNodes)
      lengthFixedNodes <- length(model$expandNodeNames(fixedNodes, returnScalarComponents = TRUE))
      paramDepDetermNodes_fixed <- model$getDependencies(fixedNodes, determOnly = TRUE)
      allCalcNodes <- model$topologicallySortNodes(c(paramDepDetermNodes_fixed, latentCalcNodes))
      paramDepDetermNodes_latent <- model$getDependencies(latentCalcNodes, determOnly = TRUE) ## to remove
      areFixedDetermNodes <- length(paramDepDetermNodes_fixed) > 0
      areLatentDetermNodes <- length(paramDepDetermNodes_latent) >0 ## to remove
      logLik <- as.numeric(0)
      gradientLik <- rep(as.numeric(0), lengthFixedNodes)
      lastParamValues <- rep(-Inf, lengthFixedNodes)
    },
    run = function(paramValues = double(1)) {
      burnin <- burnIn
      run_internal(paramValues, burnin)
      lastParamValues <<- paramValues
      return(logLik)
      returnType(double(0))
    },
    methods = list(
      resetForNewSample = function() {
        lastParamValues <<- rep(-Inf, lengthFixedNodes)
      },
      run_internal = function(paramValues = double(1), burnin = integer(0, default = 0)) {
        nSamples <- getsize(mvSample)
        sum_LL <- 0
        sum_gradient <- numeric(value = 0, length = lengthFixedNodes)
        values(model, fixedNodes) <<- paramValues
        for(i in (burnin+1):nSamples){
          nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
          sample_derivs <- nimDerivs(model$calculate(allCalcNodes), wrt = fixedNodes, order = c(0, 1))
          sample_LL <- sample_derivs$value[1]
          if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
            stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')' to see if node 'y[3]' is the cause of the problem). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
          sum_LL <- sum_LL + sample_LL
          sum_gradient <- sum_gradient + sample_derivs$jacobian[1,]
        }
        logLik <<- sum_LL / (nSamples - burnIn)
        if(is.nan(logLik))
          logLik <<- -Inf
        gradientLik <<- sum_gradient / (nSamples - burnIn)
      },
      sample_logLiks = function(paramValues = double(1)) {
        burnin <- burnIn
        nSamples <- getsize(mvSample)
        values(model, fixedNodes) <<- paramValues
        logLiks <- numeric(length = nSamples - burnin, init = FALSE)
        for(i in (burnin+1):nSamples){
          nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
          logLiks[i - burnin] <- model$calculate(allCalcNodes)
        }
        return(logLiks)
        returnType(double(1))
      },
      gradient = function(paramValues = double(1)) {
        burnin <- burnIn
        if(all(paramValues == lastParamValues))
          return(gradientLik)
        run_internal(paramValues, burnin)
        lastParamValues <<- paramValues
        return(gradientLik)
        returnType(double(1))
      },
      ## old run
      old_run = function(paramValues = double(1), oldParamValues = double(1), diff = integer(0)){
        nSamples = getsize(mvSample)
        mean_LL <- 0
        
        for(i in (burnIn+1):nSamples){
          nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
          values(model, fixedNodes) <<- paramValues  #first use new params, then old ones
          if(areFixedDetermNodes){
            simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
          }
          if(areLatentDetermNodes){
            simulate(model, paramDepDetermNodes_latent)	#	Fills in the deterministic nodes
          }
          sample_LL = calculate(model, latentCalcNodes)
          if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
            stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')' to see if node 'y[3]' is the cause of the problem). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
          mean_LL = mean_LL + sample_LL
          if(diff == 1){
            values(model, fixedNodes) <<- oldParamValues #now old params
            if(areFixedDetermNodes){
              simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
            }
            if(areLatentDetermNodes){
              simulate(model, paramDepDetermNodes_latent)  #	Fills in the deterministic nodes
            }
            sample_LL = calculate(model, latentCalcNodes)
            if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
              stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')'). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
            mean_LL = mean_LL - sample_LL
          }
        }
        values(model, fixedNodes) <<- paramValues  #first use new params, then old ones
        if(areFixedDetermNodes){
          simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
        }
        mean_LL <- mean_LL / nSamples
        if(is.nan(mean_LL)){
          mean_LL = -Inf	
        }
        returnType(double())
        return(mean_LL)
      }
    )
  )
  
  getMCEMRanges_AD <- function(model, maxNodes, buffer){
    low_limits = rep(-Inf, length(maxNodes) ) 
    hi_limits  = rep(Inf,  length(maxNodes) )
    nodes <- model$expandNodeNames(maxNodes)
    for(i in seq_along(nodes)) {
      low_limits[i] = getBound(model, nodes[i], 'lower') + abs(buffer)
      hi_limits[i]  = getBound(model, nodes[i], 'upper')  - abs(buffer)
    }
    return(list(low_limits, hi_limits))
  }
  
  buildMCEM_AD <- function(model,
                           latentNodes,
                           burnIn = 500 ,
                           mcmcControl = list(adaptInterval = 100),
                           boxConstraints = list(),
                           buffer = 10^-6,
                           alpha = 0.25,
                           beta = 0.25, 
                           gamma = 0.05,
                           C = 0.001,
                           numReps = 300,
                           forceNoConstraints = FALSE,
                           verbose = TRUE) {
    require(mcmcse)
    latentNodes = model$expandNodeNames(latentNodes)
    latentNodes <- intersect(latentNodes, model$getNodeNames(stochOnly = TRUE))
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    allStochNonDataNodes = model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    if(buffer == 0)
      warning("'buffer' is zero. This can cause problems if the likelihood function is degenerate on boundary")
    if(buffer < 0)
      stop("'buffer' must be non-negative.")
    
    if(length(setdiff(latentNodes, allStochNonDataNodes) ) != 0 )
      stop('latentNodes provided not found in model')
    maxNodes <- model$expandNodeNames(setdiff(allStochNonDataNodes, latentNodes),
                                      returnScalarComponents = TRUE)
    if(any(model$isDiscrete(maxNodes)))
      stop(paste0("MCEM cannot optimize over discrete top-level parameters. The following top-level parameters in your model are discrete: ",
                  paste0(maxNodes[model$isDiscrete(maxNodes)], collapse = ', ')))
    
    limits <- getMCEMRanges_AD(model, maxNodes, buffer)
    low_limits = limits[[1]]
    hi_limits  = limits[[2]]
    
    ## will be assign()'ed from within run() for use in getAsymptoticCov
    paramMLEs <- NA
    mcmcIters <- NA
    
    constraintNames = list()
    for(i in seq_along(boxConstraints) )
      constraintNames[[i]] = model$expandNodeNames(boxConstraints[[i]][[1]])
    for(i in seq_along(constraintNames) ) {
      limits = boxConstraints[[i]][[2]]
      inds = which(maxNodes %in% constraintNames[[i]])
      if(length(inds) == 0)
        stop(paste("warning: provided a constraint for nodes", constraintNames[[i]], ", but those nodes do not exist in the model!"))
      tooLowNodes <- which(limits[1] + abs(buffer) < low_limits[inds])
      tooHighNodes <- which(limits[2] - abs(buffer) > hi_limits[inds])
      if(length(tooLowNodes) > 0) warning(paste0("User-specified lower bound for ", constraintNames[[i]][tooLowNodes],
                                                 " is below lower bound detected by NIMBLE.  "))
      if(length(tooHighNodes) > 0) warning(paste0("User-specified upper bound for ", constraintNames[[i]][tooHighNodes],
                                                  " is above upper bound detected by NIMBLE.  "))
      low_limits[inds] = limits[1] + abs(buffer)
      hi_limits[inds] = limits[2] - abs(buffer)
    }
    if(any(low_limits>=hi_limits))
      stop('lower limits greater than or equal to upper limits!')
    if(forceNoConstraints ||
       identical(low_limits, rep(-Inf, length(low_limits))) &&
       identical(hi_limits, rep(Inf, length(hi_limits))))
      optimMethod = "BFGS"
    else 
      optimMethod = "L-BFGS-B"
    
    if(length(latentNodes) == 0)
      stop('no latentNodes')
    
    if(length(maxNodes) == 0)
      stop('no nodes to be maximized over')
    resetFunctions <- FALSE
    if(is(model, "RmodelBaseClass") ){
      Rmodel = model
      if(is(model$CobjectInterface, "uninitializedField")){
        cModel <- compileNimble(model)
      }
      else{
        cModel = model$CobjectInterface
        resetFunctions <- TRUE
      }
    }
    else{
      cModel <- model
      Rmodel <- model$Rmodel
      resetFunctions <- TRUE
    }
    
    zAlpha <- qnorm(alpha, 0, 1, lower.tail=FALSE)
    zBeta <- qnorm(beta, 0, 1, lower.tail=FALSE)
    zGamma <- qnorm(gamma, 0, 1, lower.tail=FALSE)
    
    mcmc_Latent_Conf <- configureMCMC(Rmodel, nodes = latentNodes, monitors = model$getVarNames(), control = mcmcControl, print = FALSE)
    Rmcmc_Latent <- buildMCMC(mcmc_Latent_Conf)
    sampledMV <- Rmcmc_Latent$mvSamples
    mvBlock <- modelValues(Rmodel)
    Rcalc_E_llk <- calc_E_llk_gen_AD(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)
    #  RvarCalc <- calc_asympVar(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvBlock, mvSample = sampledMV, numReps = numReps)
    #  RgetCov <- bootstrapGetCov(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)
    
    cmcmc_Latent = compileNimble(Rmcmc_Latent, project = Rmodel, resetFunctions = resetFunctions)
    #  cGetCov = compileNimble(RgetCov, project = Rmodel)  
    #  cvarCalc <- compileNimble(RvarCalc, project = Rmodel)
    cCalc_E_llk = compileNimble(Rcalc_E_llk, project = Rmodel)  
    nParams = length(maxNodes)
    run <- function(initM = 1000, thetaInit = NULL){
      if(burnIn >= initM)
        stop('mcem quitting: burnIn > initial m value')
      cmcmc_Latent$run(1, reset = TRUE)	# To get valid initial values 
      if (is.null(thetaInit)) {
        theta <- values(cModel, maxNodes)
      } else {
        theta <- thetaInit
      }
      if(optimMethod == "L-BFGS-B"){
        for(i in seq_along(maxNodes) ) {  # check that initial values satisfy constraints
          if(identical(low_limits[i], -Inf) && (hi_limits[i] < Inf)){
            if(theta[i] > hi_limits[i]){
              theta[i] <- hi_limits[i] - 1
            }
          }
          else if(identical(hi_limits[i], Inf) && (low_limits[i] > -Inf)){
            if(theta[i] < low_limits[i]){
              theta[i] <- low_limits[i] + 1
            }
          }
          else if((low_limits[i] > -Inf) && (hi_limits[i] < Inf)){
            if(!(theta[i] >= low_limits[i] & theta[i] <= hi_limits[i])){
              theta[i] = (low_limits[i] + hi_limits[i])/2
            }
          }
        }
        values(cModel, maxNodes) <<- theta
        cModel$simulate(cModel$getDependencies(maxNodes, self = FALSE))
      }
      m <- initM 
      endCrit <- C+1 #ensure that first iteration runs
      sigSq <-0 #use initM as m value for first step
      diff <- 1 # any nonzero value can be used here, gets overwritten quickly in algo
      itNum <- 0
      while(endCrit > C){ 
        acceptCrit <- 0
        #starting sample size calculation for this iteration
        m <- burnIn + ceiling(max(m - burnIn, sigSq*((zAlpha + zBeta)^2)/((diff)^2)))
        cmcmc_Latent$run(m, reset = TRUE)   #initial mcmc run of size m
        thetaPrev <- theta  #store previous theta value
        itNum <- itNum + 1
        while(acceptCrit == 0){
          cCalc_E_llk$resetForNewSample()
          if(optimMethod == "L-BFGS-B")
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, gr = cCalc_E_llk$gradient,
                                control = list(fnscale = -1), method = 'L-BFGS-B', 
                                lower = low_limits, upper = hi_limits)
          if(optimMethod == "BFGS")
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, gr = cCalc_E_llk$gradient,
                                control = list(fnscale = -1), method = 'BFGS')
          
          theta = optimOutput$par
          sample_logLiks_prev <- cCalc_E_llk$sample_logLiks(thetaPrev)
          sample_logLiks_new <- cCalc_E_llk$sample_logLiks(theta)
          mcse_result <- mcmcse::mcse(sample_logLiks_new - sample_logLiks_prev,
                                      size = ceiling(min(1000, (m - burnIn)/20)), ## backward compatible
                                      method = "obm",
                                      r = 1) ## What is the lugsail?
          ##        sigSq <- cvarCalc$run(m, theta, thetaPrev) 
          ## ase <- sqrt(sigSq) #asymptotic std. error
          ase <- mcse_result$se
          diff <- mean(sample_logLiks_new) - mean(sample_logLiks_prev)
          ## Should have: mean(sample_logLiks_new) == mcse_result$est == optimOutput$value
          ## diff <- cCalc_E_llk$run(theta, thetaPrev, 1)
          if((diff - zAlpha*ase)<0){ #swamped by mc error
            cat("Monte Carlo error too big: increasing MCMC sample size.\n")
            mAdd <- ceiling((m-burnIn)/2)  #from section 2.3, additional mcmc samples will be taken if difference is not great enough
            cmcmc_Latent$run(mAdd, reset = FALSE)
            m <- m + mAdd
          }
          else{
            acceptCrit <- 1
            endCrit <- diff + zGamma*ase #evaluate ending criterion
            if(itNum == 1)
              endCrit <- C+1 #ensure that at least two iterations are run
            
            if(verbose == T){
              cat("Iteration Number: ", itNum, ".\n", sep = "")
              cat("Current number of MCMC iterations: ", m, ".\n", sep = "")
              output = optimOutput$par
              names(output) = maxNodes
              cat("Parameter Estimates: \n", sep = "")
              print(output)
              cat("Convergence Criterion: ", endCrit, ".\n", sep = "")
            }
          }
        }
      }
      output <- optimOutput$par
      assign('paramMLEs', output, envir = parent.env(environment()))
      assign('mcmcIters', m, envir = parent.env(environment()))
      names(output) <- maxNodes
      return(output)
    }
    
    ## This has not been updated
    estimateCov = function(MLEs = NA, useExistingSamples = FALSE){
      delta <- .0001
      if(!(length(MLEs) == 1 && is.na(MLEs))){
        if(!is.numeric(MLEs)) stop("MLEs argument must be numeric.")
        if(!identical(sort(names(MLEs)), sort(maxNodes))){
          stop(paste('MLEs argument must be a named vector with MLEs for all of the following parameters: ', paste(maxNodes, collapse = ", ")))
        }
        covMLEs <- unname(MLEs[maxNodes])
      }
      else{
        if(length(paramMLEs) == 1 && is.na(paramMLEs)){
          stop(paste('No MLEs argument was provided, and the run() method has not been called yet.  Please call the run() method first or provide a named vector of MLEs.'))
        }
        else{
          covMLEs <- unname(paramMLEs)
        }
      }
      if(dim(as.matrix(cmcmc_Latent$mvSamples))[1]<2){
        if(useExistingSamples == TRUE){
          warning('MCMC over latent states has not been run yet, cannot have useExistingSamples = TRUE')
          useExistingSamples <- FALSE
        }
      }
      values(cModel, maxNodes) <- covMLEs
      calculate(cModel, cModel$getDependencies(maxNodes))
      if(!useExistingSamples){
        if(is.na(mcmcIters)) mcmcIters <- 20000
        cmcmc_Latent$run(mcmcIters)
      }
      FIM <- cGetCov$run(covMLEs, delta)
      cov <- solve(FIM)
      colnames(cov) <- maxNodes
      rownames(cov) <- maxNodes
      return(cov)
    }
    return(list(run = run,
                estimateCov = estimateCov))
  }
  
  ## This has not been updated
  ## bootstrapGetCov <- nimbleFunction(
  ##     name = 'bootstrapGetCov',
  ##     setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0){
  ##       fixedCalcNodes <- model$getDependencies(fixedNodes)	
  ##       latentCalcNodes <- model$getDependencies(sampledNodes)
  ##       paramDepDetermNodes_fixed <- model$getDependencies(fixedNodes, determOnly = TRUE) 
  ##       paramDepDetermNodes_latent <- model$getDependencies(latentCalcNodes, determOnly = TRUE) 
  ##       areFixedDetermNodes <- length(paramDepDetermNodes_fixed) > 0
  ##       areLatentDetermNodes <- length(paramDepDetermNodes_latent) > 0
  ##     },
  ##   run = function(theta = double(1), delta = double(0)){
  ##     nSamples <- getsize(mvSample)
  ##     paramLengths <- length(theta)
  ##     fxph <- numeric(paramLengths)
  ##     fxmh <- numeric(paramLengths)
  ##     grad <- numeric(paramLengths)
  ##     derivxy <- matrix(0, nrow = paramLengths, ncol = paramLengths)
  ##     meanGrad <- numeric(paramLengths)
  ##     meanGradGrad <-  matrix(0, nrow = paramLengths, ncol = paramLengths)
  ##     meanDerivxy <- matrix(0, nrow = paramLengths, ncol = paramLengths)
  ##     for(i in (burnIn+1):nSamples){
  ##       values(model, fixedNodes) <<- theta
  ##       if(areFixedDetermNodes){
  ##         simulate(model, paramDepDetermNodes_fixed)  
  ##       }
  ##       if(areLatentDetermNodes){
  ##         simulate(model, paramDepDetermNodes_latent)
  ##       }
  ##       nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
  ##       if(areFixedDetermNodes){
  ##         simulate(model, paramDepDetermNodes_fixed) 
  ##       }
  ##       if(areLatentDetermNodes){
  ##         simulate(model, paramDepDetermNodes_latent)	
  ##       }
  ##       origValue <- calculate(model, latentCalcNodes)
  ##       for(iNode in 1:paramLengths){
  ##         theta[iNode] <- theta[iNode] + delta
  ##         values(model, fixedNodes) <<- theta 
  ##         if(areFixedDetermNodes){
  ##           simulate(model, paramDepDetermNodes_fixed)  
  ##         }
  ##         if(areLatentDetermNodes){
  ##           simulate(model, paramDepDetermNodes_latent)	
  ##         }
  ##         fxph[iNode] <- calculate(model, latentCalcNodes)
  ##         theta[iNode] <- theta[iNode] - 2*delta
  ##         values(model, fixedNodes) <<- theta 
  ##         if(areFixedDetermNodes){
  ##           simulate(model, paramDepDetermNodes_fixed)  
  ##         }
  ##         if(areLatentDetermNodes){
  ##           simulate(model, paramDepDetermNodes_latent)
  ##         }
  ##         fxmh[iNode] <- calculate(model, latentCalcNodes)
  ##         grad[iNode] <- (fxph[iNode] - fxmh[iNode])/(2*delta)
  ##         theta[iNode] <- theta[iNode] + delta
  ##         derivxy[iNode, iNode] <- (fxph[iNode] -2*origValue + fxmh[iNode])/(delta^2)
  ##       }
  ##       for(iNode in 1:paramLengths){
  ##         if(iNode != paramLengths){
  ##           for(jNode in (iNode+1):paramLengths){
  ##             theta[iNode] <- theta[iNode] + delta
  ##             theta[jNode] <- theta[jNode] + delta
  ##             values(model, fixedNodes) <<- theta 
  ##             if(areFixedDetermNodes){
  ##               simulate(model, paramDepDetermNodes_fixed) 
  ##             }
  ##             if(areLatentDetermNodes){
  ##               simulate(model, paramDepDetermNodes_latent)	
  ##             }
  ##             fxyph <- calculate(model, latentCalcNodes)
  ##             theta[iNode] <- theta[iNode] - 2*delta
  ##             theta[jNode] <- theta[jNode] - 2*delta
  ##             values(model, fixedNodes) <<- theta
  ##             if(areFixedDetermNodes){
  ##               simulate(model, paramDepDetermNodes_fixed) 
  ##             }
  ##             if(areLatentDetermNodes){
  ##               simulate(model, paramDepDetermNodes_latent)
  ##             }
  ##             fxymh <- calculate(model, latentCalcNodes)
  ##             derivxy[iNode, jNode] <- (fxyph - fxph[iNode] - fxph[jNode] + 2*origValue - fxmh[iNode] - fxmh[jNode] + fxymh)/(2*delta^2)
  ##             theta[iNode] <- theta[iNode] + delta
  ##             theta[jNode] <- theta[jNode] + delta
  ##           }
  ##         }
  ##       }
  ##       meanGrad <- meanGrad + grad
  ##       meanDerivxy <- meanDerivxy + derivxy
  ##       meanGradGrad <- meanGradGrad + (grad%*%t(grad))
  ##     }
  ##     meanGrad <- meanGrad / (nSamples - burnIn)
  ##     meanDerivxy <- meanDerivxy /  (nSamples - burnIn)
  ##     meanGradGrad <- meanGradGrad /  (nSamples - burnIn)
  ##     for(iNode in 1:paramLengths){
  ##       if(iNode != paramLengths){
  ##         for(jNode in (iNode+1):paramLengths){
  ##           meanDerivxy[jNode, iNode] <- meanDerivxy[iNode, jNode] ## smarter way to do this with matrix mult?
  ##         }
  ##       }
  ##     }
  ##     returnType(double(2))
  ##     returnMat <- -meanDerivxy - meanGradGrad +(meanGrad%*%t(meanGrad))
  ##     return(returnMat)
  ##   },where = getLoadingNamespace())
  
  init <- runif(20, -1, 1)
  
  MSOMmcem <- buildMCEM_AD(MSOMmodel, latentNodes = latentNodes,
                           C = 0.05, alpha = .01, beta = .01, gamma = .01)
  
  scale_conversion <- function(x) {
    sapply(1:length(x), function(i) {
      ifelse(i %% 2, ifelse(i >= 9, x[i], exp(x[i]) / (1 + exp(x[i]))), exp(x[i]))
    })
  }
  
  df <- data.frame(matrix(nrow = 5, ncol=20))
  
  for (j in 1:5) {
    init <- runif(20, -1, 1)
    df[j, ] <- scale_conversion(MSOMmcem$run(initM = 1000, thetaInit = init))
  }
  
  names(df) <- c("cato.occ.mean", "sigma.ucato", "fcw.occ.mean", "sigma.ufcw",
                 "cato.det.mean", "sigma.vcato", "fcw.det.mean", "sigma.vfcw",
                 "mu.a1", "sigma.a1", "mu.a2", "sigma.a2",
                 "mu.a3", "sigma.a3", "mu.a4", "sigma.a4",
                 "mu.b1", "sigma.b1", "mu.b2", "sigma.b2")
  df
}


MCEM_results <- parLapply(cl = this_cluster, X = 1:6, 
                         fun = run_MCEM, 
                         survey.data = survey.data,
                         species.groups = species.groups, 
                         survey.dates = survey.dates, 
                         habitat = habitat)


