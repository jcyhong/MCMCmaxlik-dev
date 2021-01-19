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

# Fix the seed for reproducibility.
set.seed(21745)

# Load the packages ---------------------------------------
library("MCMCmaxlik")
require(glmm)
require(nimble)
library(parallel)
library(data.table)
this_cluster <- makeCluster(4)

set.seed(12301)

run_allcode <- function(b) {
  library(nimble)
  require(glmm)
  library("MCMCmaxlik")
  nimbleOptions(experimentalEnableDerivs = TRUE)
  # Model specification -------------------------------------
  data(salamander)
  names(salamander)
  nrow(salamander) 
  #### organize data ####
  orgF=cbind.data.frame(Findex=1:length(unique(salamander$Female)),idF=unique(salamander$Female))
  
  orgM=cbind.data.frame(Mindex=1:length(unique(salamander$Male)),idM=unique(salamander$Male))
  
  salamander2=merge(salamander,orgF,by.x="Female",by.y="idF")
  salamander3=merge(salamander2,orgM,by.x="Male",by.y="idM")
  
  isRR=ifelse(salamander3$Cross=="R/R",1,0)
  #head(isRR)
  
  ## indicator for type of pairing
  isRW=ifelse(salamander3$Cross=="R/W",1,0)
  #head(isRW)
  
  isWR=ifelse(salamander3$Cross=="W/R",1,0)
  #head(isWR,20)
  
  isWW=ifelse(salamander3$Cross=="W/W",1,0)
  
  
  #### nimble setup ####
  glmmCode <- nimbleCode({
    for(i in 1:N){
      logit(theta[i])<-beta1*isRR[i]+beta2*isRW[i]+beta3*isWR[i]+beta4*isWW[i]+REF[Findex[i]]+REM[Mindex[i]]
      y[i]~dbin(theta[i],1)
    }
    for(i in 1:numFemales){
      REF[i]~dnorm(0, sd=sqrt(varF))
      
    }
    for(i in 1:numMales){
      REM[i]~dnorm(0, sd=sqrt(varM))
      
      
    }
    
    beta1 ~ dnorm(0, sd = 10000)
    beta2 ~ dnorm(0, sd = 10000)
    beta3 ~ dnorm(0, sd = 10000)
    beta4 ~ dnorm(0, sd = 10000)
    varF ~ dunif(0, 1000)
    varM ~ dunif(0, 1000)
    
  })
  
  glmmConstants <- list(N = nrow(salamander3),numFemales=length(unique(salamander3$Female)),numMales=length(unique(salamander3$Male)),Findex=salamander3$Findex,Mindex=salamander3$Mindex)
  
  glmmData <- list(
    isRR=isRR,isWR=isWR,isRW=isRW,isWW=isWW,y=salamander3$Mate
  )
  
  glmmInits <- list(beta1 = 1, beta2 = 1, beta3=1,beta4=1,varF=5,varM=5)
  
  glmmMod <- nimbleModel(code=glmmCode, constants=glmmConstants, data=glmmData, inits=glmmInits, check = FALSE)
  
  glmmParamNodes <- glmmMod$getNodeNames(topOnly = T,stochOnly=T)
  
  Cglmm <- compileNimble(glmmMod)
  
  compiledFunsGlmm <- buildMCMCmaxlik(glmmMod, glmmParamNodes)
  
  init <- c(runif(1, min=-10, max=10), 
            runif(1, min=-10, max=10),
            runif(1, min=-10, max=10), 
            runif(1, min=-10, max=10),
            runif(1, min=0.05, max=15),
            runif(1, min=0.05, max=15))
  boundary <- list(c(-40, 40), c(-40, 40), c(-40, 40),
                   c(-40, 40), c(0.01, 100), c(0.01, 100))
  numMCMCSamples <- 300
  maxIter <- 300
  
  # 1. Fixed step size ----------------------------------------
  resultsGLMMFixed <- computeMLE(glmmMod, glmmParamNodes,
                                 method="fixed", paramInit=init,
                                 stepsize=0.05,
                                 compiledFuns=compiledFunsGlmm,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=maxIter,
                                 boundary=boundary)
  
  # 2. Small fixed step size ----------------------------------------
  resultsGLMMFixedSmall <- computeMLE(glmmMod, glmmParamNodes,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsGlmm,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=maxIter,
                                      boundary=boundary)
  
  # 4. Adam ----------------------------------------
  resultsGLMMAdam <- computeMLE(glmmMod, glmmParamNodes,
                            method="adam", paramInit=init,
                            compiledFuns=compiledFunsGlmm,
                            numMCMCSamples=numMCMCSamples,
                            maxIter=maxIter,
                            boundary=boundary)
  
  # 5. Newton-Raphson ----------------------------------------
  # resultsGLMMNR <- computeMLE(glmmMod, glmmParamNodes,
  #                         method="NR", paramInit=init,
  #                         compiledFuns=compiledFunsGlmm,
  #                         numMCMCSamples=numMCMCSamples,
  #                         tol=1e-20,
  #                         maxIter=maxIter,
  #                         boundary=boundary)
  
  # 6. 1-D sampling ----------------------------------------
  resultsGLMM1D <- computeMLE(glmmMod, glmmParamNodes,
                          method="ga1D", paramInit=init,
                          compiledFuns=compiledFunsGlmm,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=maxIter, 
                          boundary=boundary)  
  
  df <- data.frame(do.call(rbind, list(resultsGLMMFixed$MLE,
                                       resultsGLMMFixedSmall$MLE,
                                       # resultsGLMMNR$MLE,
                                       resultsGLMMAdam$MLE,
                                       resultsGLMM1D$MLE)))
  df$method <- c('Fixed (0.05)', 'Fixed (0.005)', 
                 #'Newton-Raphson', 
                 'Adam', '1D sampling')
  
  names(df) <- c("beta1", "beta2", "beta3", "beta4", "varF", "varM", "method")
  df
}

df_list_all <- parLapply(cl = this_cluster, X = 1:30, 
          fun = run_allcode)

stopCluster(this_cluster)

df_combined_all <- rbindlist(df_list_all)
df_rmse <- df_combined_all %>% group_by(method) %>%
  summarize(RMSE_beta1=sqrt(mean((beta1 - 1.008)^2)),
            RMSE_beta2=sqrt(mean((beta2 - 0.306)^2)),
            RMSE_beta3=sqrt(mean((beta3 - (-1.896))^2)),
            RMSE_beta4=sqrt(mean((beta4 - 0.990)^2)),
            RMSE_varF=sqrt(mean((varF - 1.174)^2)),
            RMSE_varM=sqrt(mean((varM - 1.041)^2)))

write.csv(df_rmse, file='glmm_rmse.csv', row.names=F)
