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
library(dplyr)
this_cluster <- makeCluster(4)

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

