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

compiledFunsGlmm <- buildMCMCmaxlik(glmmMod, glmmParamNodes)


init <- rep(2, 6)
boundary <- list(c(-40, 40), c(-40, 40), c(-40, 40),
                 c(-40, 40), c(0.01, 100), c(0.01, 100))
numMCMCSamples <- 300
maxIter <- 300

getGLMMResults <- function(init, 
                           boundary=list(c(-40, 40), c(-40, 40), c(-40, 40),
                                         c(-40, 40), c(0.01, 100), c(0.01, 100)), 
                           numMCMCSamples) {
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
  resultsAdam <- computeMLE(glmmMod, glmmParamNodes,
                            method="adam", paramInit=init,
                            compiledFuns=compiledFunsGlmm,
                            numMCMCSamples=numMCMCSamples,
                            maxIter=maxIter,
                            boundary=boundary)
  
  # 5. Newton-Raphson ----------------------------------------
  resultsNR <- computeMLE(glmmMod, glmmParamNodes,
                          method="NR", paramInit=init,
                          compiledFuns=compiledFunsGlmm,
                          numMCMCSamples=numMCMCSamples,
                          tol=1e-20,
                          maxIter=maxIter,
                          boundary=boundary)
  
  # 6. 1-D sampling ----------------------------------------
  results1D <- computeMLE(glmmMod, glmmParamNodes,
                          method="ga1D", paramInit=init,
                          compiledFuns=compiledFunsGlmm,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=maxIter, 
                          boundary=boundary)
  
  list(resultsGLMMFixed, resultsGLMMFixedSmall,
       resultsAdam, resultsNR, results1D)
}

sal2=glmer(Mate~ -1 + Cross + (1|Female) + (1|Male),
           data=salamander, family=binomial)
resultsGLMMGlmer <- c(coef(summary(sal2))[, 'Estimate'], 
                      attr(summary(sal2)$varcor$Female, 'stddev')^2,
                      attr(summary(sal2)$varcor$Male, 'stddev')^2)
gmreDevfunGLMM <- update(sal2, devFunOnly=TRUE)
# sigmaRE, beta0, beta1
# log-likelihood difference


getGLMMSummary <- function(resultsGLMM) {
  timesGLMM <- sapply(resultsGLMM,
                      function(results) {
                        exec.time <- results$execution.time[3]
                        conv.time <- results$convergence.time[3]
                        conv.iter <- results$convergence.iter
                        if(is.null(conv.time)) conv.time <- NA
                        if(is.null(conv.iter)) conv.iter <- NA
                        benchmarkMLE <- resultsGLMMGlmer
                        names(benchmarkMLE) <- c('beta1', 'beta2', 'beta3', 'beta4', 'varF', 'varM')
                        benchmarkMLE <- benchmarkMLE[c('varF', 'varM', 'beta1', 'beta2', 'beta3', 'beta4')]
                        benchmarkMLE[1:2] <- sqrt(benchmarkMLE[1:2])
                        MLE <- results$MLE
                        names(MLE) <- c('beta1', 'beta2', 'beta3', 'beta4', 'varF', 'varM')
                        MLE <- MLE[c('varF', 'varM', 'beta1', 'beta2', 'beta3', 'beta4')]
                        MLE[1:2] <- sqrt(MLE[1:2])
                        c(round(results$MLE, 3), 
                          round(exec.time, 3), 
                          round(conv.time, 3), 
                          conv.iter,
                          round((gmreDevfunGLMM(benchmarkMLE) - gmreDevfunGLMM(MLE)) / -2, 5))
                      })
  rownames(timesGLMM) <- c(glmmParamNodes, 'exec.time', 'conv.time', 'conv.iter', 'loglik_diff')
  colnames(timesGLMM) <- c('fixed', 'small_fixed', 'Adam', 'NR', '1D')
  t(timesGLMM)
}



set.seed(1000)
resultsGLMM_2_300 <- getGLMMResults(init=rep(2, 6), numMCMCSamples=300)
timesGLMM_2_300 <- getGLMMSummary(resultsGLMM_2_300)

resultsGLMM_4_300 <- getGLMMResults(init=rep(4, 6), numMCMCSamples=300)
timesGLMM_4_300 <- getGLMMSummary(resultsGLMM_4_300)

rbind(cbind(timesGLMM_2_300, init=2),
      cbind(timesGLMM_4_300, init=4))

write.csv(rbind(cbind(timesGLMM_2_300, init=2),
                cbind(timesGLMM_4_300, init=4)), 'times_glmm2.csv')

# MCEM
resultsGLMMMCEM_MLE <- c(1.0168445, 0.3181842, -1.9470261, 0.9946614, 1.3877689, 1.2448641)

benchmarkMLE <- resultsGLMMGlmer
names(benchmarkMLE) <- c('beta1', 'beta2', 'beta3', 'beta4', 'varF', 'varM')
benchmarkMLE <- benchmarkMLE[c('varF', 'varM', 'beta1', 'beta2', 'beta3', 'beta4')]
benchmarkMLE[1:2] <- sqrt(benchmarkMLE[1:2])

names(resultsGLMMMCEM_MLE) <- c('beta1', 'beta2', 'beta3', 'beta4', 'varF', 'varM')
resultsGLMMMCEM_MLE <- resultsGLMMMCEM_MLE[c('varF', 'varM', 'beta1', 'beta2', 'beta3', 'beta4')]
resultsGLMMMCEM_MLE[1:2] <- sqrt(resultsGLMMMCEM_MLE[1:2])

round((gmreDevfunGLMM(benchmarkMLE) - gmreDevfunGLMM(resultsGLMMMCEM_MLE)) / -2, 5)



ptm <- proc.time()
sal <- glmm(Mate ~ 0 + Cross, random=list(~ 0 + Female,
                                          ~ 0 + Male), 
            varcomps.names=c("F", "M"), data=salamander,
            family.glmm=bernoulli.glmm, m=10^5, debug=TRUE)
glmmTime <- proc.time() - ptm 
glmmTime ## 813.043 
summary(sal)

# Call:
#   glmm(fixed = Mate ~ 0 + Cross, random = list(~0 + Female, ~0 + 
#                                                  Male), varcomps.names = c("F", "M"), data = salamander, family.glmm = bernoulli.glmm, 
#        m = 10^5, debug = TRUE)
# 
# 
# Link is: "logit (log odds)"
# 
# Fixed Effects:
#   Estimate Std. Error z value Pr(>|z|)    
# CrossR/R   0.9964     0.3871   2.574   0.0101 *  
#   CrossR/W   0.3410     0.4363   0.781   0.4345    
# CrossW/R  -1.8957     0.4335  -4.373 1.23e-05 ***
#   CrossW/W   1.0173     0.4860   2.093   0.0363 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# Variance Components for Random Effects (P-values are one-tailed):
#   Estimate Std. Error z value Pr(>|z|)/2   
# F   1.2832     0.5569   2.304    0.01061 * 
#   M   1.1525     0.4199   2.745    0.00303 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ptm <- proc.time()
sal2=glmer(Mate ~ -1 + Cross + (1|Female) + (1|Male),
           data=salamander, family=binomial)
glmerTime <- proc.time() - ptm 
glmerTime ## 0.404  
summary(sal2)

# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: Mate ~ -1 + Cross + (1 | Female) + (1 | Male)
# Data: salamander
# 
# AIC      BIC   logLik deviance df.resid 
# 430.6    453.9   -209.3    418.6      354 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.0508 -0.6156  0.2714  0.5973  2.5507 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# Female (Intercept) 1.174    1.084   
# Male   (Intercept) 1.041    1.020   
# Number of obs: 360, groups:  Female, 60; Male, 60
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# CrossR/R   1.0082     0.3937   2.561   0.0105 *  
#   CrossR/W   0.3062     0.3747   0.817   0.4138    
# CrossW/R  -1.8960     0.4460  -4.251 2.13e-05 ***
#   CrossW/W   0.9904     0.3912   2.532   0.0113 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   CrsR/R CrsR/W CrsW/R
# CrossR/W  0.280              
# CrossW/R  0.112 -0.019       
# CrossW/W  0.042  0.249  0.145

set.seed(32856)
#source("MCEM_with_output.R")
glmmNew <- glmmMod$newModel()
latentNodes <- glmmNew$getNodeNames(latentOnly=TRUE, stochOnly=TRUE)

#glmmMCEM <- buildMCEM(model=glmmNew, latentNodes=latentNodes,
#                      theta0=rep(2, 6)) 

glmmMCEM <- buildMCEM(model=glmmNew, latentNodes=latentNodes) 

ptm <- proc.time()
resultsMCEM <- glmmMCEM$run(initM = 1000)
timeGLMMMCEM <- proc.time() - ptm

save(resultsMCEM, "MCEM_GLMM.RData")

timeMCEM 


## killed it after an hour or so

# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 1.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta1      beta2      beta3      beta4       varF       varM 
# 0.7719557  0.3153691 -1.2421619  0.7727807  4.8870236  4.9359067 
# Convergence Criterion: 1.001.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 2.
# Current number of MCMC iterations: 1000.
# Parameter Estimates: 
#   beta1      beta2      beta3      beta4       varF       varM 
# 0.7560798  0.2426937 -1.4150895  0.7603411  4.8581485  4.9220263 
# Convergence Criterion: 0.1776775.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 3.
# Current number of MCMC iterations: 1250.
# Parameter Estimates: 
#   beta1      beta2      beta3      beta4       varF       varM 
# 0.7566014  0.2365373 -1.4296400  0.7556359  4.9064599  4.8710893 
# Convergence Criterion: 0.04049478.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 4.
# Current number of MCMC iterations: 2188.
# Parameter Estimates: 
#   beta1      beta2      beta3      beta4       varF       varM 
# 0.7537331  0.2449870 -1.4282083  0.7519991  4.8582776  4.8518666 
# Convergence Criterion: 0.02039753.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 5.
# Current number of MCMC iterations: 2188.
# Parameter Estimates: 
#   beta1      beta2      beta3      beta4       varF       varM 
# 0.7605328  0.2403929 -1.4287206  0.7614988  4.8453013  4.8852084 
# Convergence Criterion: 0.01239757.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 6.
# Current number of MCMC iterations: 4298.
# Parameter Estimates: 
#   beta1      beta2      beta3      beta4       varF       varM 
# 0.7569682  0.2437938 -1.4269022  0.7564405  4.8582642  4.8921791 
# Convergence Criterion: 0.003558387.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Iteration Number: 7.
# Current number of MCMC iterations: 6197.
# Parameter Estimates: 
#   beta1      beta2      beta3      beta4       varF       varM 
# 0.7539376  0.2429794 -1.4295914  0.7562341  4.8676279  4.9083186 
# Convergence Criterion: 0.003545863.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|
#   Monte Carlo error too big: increasing MCMC sample size.
# |-------------|-------------|-------------|-------------|
#   |-------------------------------------------------------|



#save(resultsMCEM, file="MCEM_GLMM.RData")

#timeInfo=cbind(c(timeGLMMFixed, timeGLMMSmallFixed,
# timeGLMMAdadelta, timeGLMMAdam, timeGLMM1D, timeGLMMMCEM),
# c("fixed","smallFixed","adadelta","adam","1D","mcem"))


#write.csv(timeInfo,"timeGLMM.csv",row.names=F)

load(file="GLMM1D.RData")
load(file="GLMMAdam.RData")
#load(file="GLMMAdadelta.RData")
load(file="GLMMFixedSmall.RData")
load(file="GLMMFixed.RData")

iteratesFixed <- data.frame(
  iter=1:nrow(resultsGLMMFixed$param),
  resultsGLMMFixed$param,
  method=rep("Fixed (0.05)", nrow(resultsGLMMFixed$param))
)
iteratesSmallFixed <- data.frame(
  iter=1:nrow(resultsGLMMFixedSmall$param),
  resultsGLMMFixedSmall$param,
  method=rep("Fixed (0.005)", nrow(resultsGLMMFixedSmall$param))
)
iteratesNR <- data.frame(
  iter=1:nrow(resultsNR$param),
  resultsNR$param,
  method=rep("Newton-Raphson", nrow(resultsNR$param))
)
iterates1D <- data.frame(
  iter=1:nrow(results1D$param),
  results1D$param,
  method=rep("1D sampling", nrow(results1D$param))
)
iteratesAdam <- data.frame(
  iter=1:nrow(resultsAdam$param),
  resultsAdam$param, 
  method=rep("Adam", nrow(resultsAdam$param))
)

iterates <- rbind(iteratesFixed, iteratesSmallFixed,
                  iteratesNR,
                  iteratesAdam,
                  iterates1D) 
names(iterates) <- c("iter", glmmParamNodes, "method")

write.csv(iterates, "glmmResultsAll.csv", row.names=F)

iteratesMelted <- melt(iterates, id.vars=c("iter", "method"))
levels(iteratesMelted$variable) <- c(expression(beta[1]),
                                     expression(beta[2]),
                                     expression(beta[3]),
                                     expression(beta[4]),
                                     expression(sigma[F]^2), 
                                     expression(sigma[M]^2))
library(RColorBrewer)
cbPalette <- brewer.pal(8, "Dark2")


trajectoryPlot <- ggplot(iteratesMelted, aes(x=iter, y=value, colour=method, linetype=method)) 
trajectoryPlot <- trajectoryPlot + 
  geom_line(alpha=0.8, lwd=0.8)
trajectoryPlot<- trajectoryPlot + 
  facet_grid(variable ~ ., labeller=label_parsed, scales="free_y")
trajectoryPlot <- trajectoryPlot + theme(legend.position="none")
trajectoryPlot <- trajectoryPlot + 
  scale_color_manual(values=cbPalette)
trajectoryPlot <- trajectoryPlot + xlab("number of iterations")
trajectoryPlot<- trajectoryPlot + 
  theme(axis.title=element_text(face="bold", size=20),
        axis.text=element_text(size=15),
        strip.text.y=element_text(size=15))  +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank(), legend.text=element_text(size=15))

pdf(paste0("glmm_", 300, ".pdf"), width = 7.3)
trajectoryPlot
dev.off()

cbPalette <- brewer.pal(8, "Dark2")
trajectoryPlot <- ggplot(iteratesLogreg,
                         aes(x=iter, y=value, colour=method)) 
trajectoryPlot <- trajectoryPlot + geom_line(alpha=0.8) + 
  geom_point(data=iteratesMelted, mapping=aes(x=iter, y=value, shape=method))
trajectoryPlot <- trajectoryPlot + 
  facet_grid(variable ~ .,
             labeller=label_parsed,
             scales="free_y")
trajectoryPlot<- trajectoryPlot
trajectoryPlot <- trajectoryPlot + 
  scale_color_manual(values=cbPalette,labels=unique(iterates$method)) + 
  scale_shape_manual(values=c(0, 5, 6, 15,8),labels=unique(iterates$method))
trajectoryPlot<- trajectoryPlot + xlab("number of iterations")
trajectoryPlot <- trajectoryPlot + 
  theme(axis.title=element_text(face="bold", size=20),
        axis.text=element_text(size=15),
        strip.text.y=element_text(size=15)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank(), legend.text=element_text(size=15))

## below doesn't work but we don't really need the legend seperate

# Extract Legend 
g_legend <- function(a.gplot) { 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}

pdf(paste0("legend_glmm_", 300, ".pdf"), width=12)
legend <- g_legend(trajectoryPlot) 
grid.draw(legend) 
dev.off()
