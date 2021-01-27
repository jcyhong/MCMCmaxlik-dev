# Load the packages ---------------------------------------
library("MCMCmaxlik")
require(glmm)
require(nimble)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
require(lme4)
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

## checks
range(salamander3$Findex)
length(unique(salamander$Female))
range(salamander3$Mindex)
length(unique(salamander$Male))


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


#source("MCEM_with_output.R")
glmmNew <- glmmMod$newModel()
latentNodes <- glmmNew$getNodeNames(latentOnly=TRUE, stochOnly=TRUE)


source("~/Desktop/MCMCmaxlik-dev/packages/MCMCmaxlik/R/MCEM_AD_build.R")
C_opts = c(0.001, 0.01, 0.1)

i <- 3
glmmMCEM <- buildMCEM_AD(
  model = glmmNew, latentNodes = latentNodes,
  C = C_opts[i]
)

MCEM_time <- proc.time()
out <- glmmMCEM$run()
MCEM_time <- proc.time() - MCEM_time


#beta1      beta2      beta3      beta4       varF       varM 
#1.0137778  0.3236878 -1.9439178  0.9937086  1.3875624  1.2328804 
#1.0007060  0.3036434 -1.9475255  0.9860114  1.3951910  1.2247479 
#1.1400793  0.4181533 -2.0556457  1.0989191  1.8081469  1.5473603 

#user  system elapsed 
#261.607   1.555 265.566 
# 51.052   0.278  51.060 
# 11.190   0.259  10.105 