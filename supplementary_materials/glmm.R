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

# Model specification -------------------------------------
data(salamander)
names(salamander)
nrow(salamander) 

isRR <- ifelse(salamander$Cross == "R/R", 1, 0)
head(isRR)

isRW <- ifelse(salamander$Cross == "R/W", 1, 0)
head(isRW)

isWR <- ifelse(salamander$Cross == "W/R", 1, 0)
head(isWR, 20)

isWW <- ifelse(salamander$Cross == "W/W", 1, 0)
head(isWW, 20)

code <- nimbleCode({
  for(i in 1:N) {
    logit(theta[i]) <- beta1 * isRR[i] + beta2 * isRW[i] + 
      beta3 * isWR[i] + beta4 * isWW[i] + REF[i] + REM[i]
    y[i] ~ dbin(theta[i], 1)
    REF[i] ~ dnorm(0, varF)
    REM[i] ~ dnorm(0, varM)
  }
  beta1 ~ dnorm(0, sd=10000)
  beta2 ~ dnorm(0, sd=10000)
  beta3 ~ dnorm(0, sd=10000)
  beta4 ~ dnorm(0, sd=10000)
  varF ~ dunif(0, 1000)
  varM ~ dunif(0, 1000)
})

constants <- list(N=nrow(salamander))

data <- list(
  isRR=isRR, isWR=isWR, isRW=isRW, isWW=isWW, y=salamander$Mate
)

inits <- list(beta1=1, beta2=1, beta3=1, beta4=1, varF=5, varM=5)

glmmMod <- nimbleModel(code=code, constants=constants, data=data, 
                       inits=inits, check=FALSE)

paramNodes <- glmmMod$getNodeNames(topOnly=T, stochOnly=T) 
## "beta1" "beta2" "beta3" "beta4" "varF"  "varM"

Cglmm <- compileNimble(glmmMod)

compiledFunsglmm <- buildMCMCmaxlik(Cglmm, paramNodes)

init <- rep(2, 6)
boundary <- list(c(-40, 40), c(-40, 40), c(-40, 40),
                 c(-40, 40), c(0.01, 100), c(0.01, 100))
numMCMCSamples <- 300

# 1. Fixed step size ----------------------------------------
ptm <- proc.time()
resultsGLMMFixed <- computeMLE(glmmMod, paramNodes,
                               method="fixed", paramInit=init,
                               stepsize=0.05,
                               compiledFuns=compiledFunsglmm,
                               numMCMCSamples=numMCMCSamples,
                               maxIter=300,
                               boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeGLMMFixed <-proc.time() - ptm 
timeGLMMFixed ## 102.840 
save(resultsGLMMFixed, file="GLMMFixed.RData")

apply(tail(resultsGLMMFixed$param, 20), 2, mean, trim=0.2)
#0.8371804  0.2690307 -1.5762730  0.8357959  1.7182524  2.2435778

ptm <- proc.time()
resultsGLMMFixed2 <- computeMLE(glmmMod, paramNodes,
                               method="fixed", paramInit=init,
                               stepsize=0.05,
                               compiledFuns=compiledFunsglmm,
                               numMCMCSamples=numMCMCSamples,
                               maxIter=300,
                               boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeGLMMFixed2 <- proc.time() - ptm  
timeGLMMFixed2 ##  102.760 
resultsGLMMFixed2$iter ## 300
save(resultsGLMMFixed2, file="GLMMFixedCC.RData")

apply(tail(resultsGLMMFixed$param, 20), 2, mean, trim=0.2)
#0.8371804  0.2690307 -1.5762730  0.8357959  1.7182524  2.2435778

# 2. Small fixed step size ----------------------------------------
ptm <- proc.time()
resultsGLMMFixedSmall <- computeMLE(glmmMod, paramNodes,
                                    method="fixed", paramInit=init,
                                    stepsize=0.005,
                                    compiledFuns=compiledFunsglmm,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=3000,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeGLMMFixedSmall <- proc.time() - ptm  
timeGLMMFixedSmall ## 1019.076 
save(resultsGLMMFixedSmall, file="GLMMFixedSmall3000.RData")

apply(tail(resultsGLMMFixedSmall$param, 20), 2, mean, trim=0.2)
##0.8816105  0.2855453 -1.6549337  0.8874643  1.5867583  1.3504009

ptm <- proc.time()
resultsGLMMFixedSmall2 <- computeMLE(glmmMod, paramNodes,
                                    method="fixed", paramInit=init,
                                    stepsize=0.005,
                                    compiledFuns=compiledFunsglmm,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=3000,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeGLMMFixedSmall2 <- proc.time() - ptm 
timeGLMMFixedSmall2 ## 327.310 
resultsGLMMFixedSmall2$iter ## 971
save(resultsGLMMFixedSmall2, file="GLMMFixedSmallCC3000.RData")

apply(tail(resultsGLMMFixedSmall2$param, 20), 2, mean, trim=0.2)
#0.8603595  0.2776129 -1.6058155  0.8543346  1.8159546  1.5992581

ptm <- proc.time()
resultsGLMMFixedSmall3 <- computeMLE(glmmMod, paramNodes,
                                    method="fixed", paramInit=init,
                                    stepsize=0.005,
                                    compiledFuns=compiledFunsglmm,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeGLMMFixedSmall3 <- proc.time() - ptm  
timeGLMMFixedSmall3 ## 101.730 
save(resultsGLMMFixedSmall3, file="GLMMFixedSmall.RData")

apply(tail(resultsGLMMFixedSmall3$param, 20), 2, mean, trim=0.2)
## 0.8526400  0.2687585 -1.6110869  0.8586432  1.6645518  1.7643731

ptm <- proc.time()
resultsGLMMFixedSmall4 <- computeMLE(glmmMod, paramNodes,
                                     method="fixed", paramInit=init,
                                     stepsize=0.005,
                                     compiledFuns=compiledFunsglmm,
                                     numMCMCSamples=numMCMCSamples,
                                     maxIter=300,
                                     boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeGLMMFixedSmall4 <- proc.time() - ptm  
timeGLMMFixedSmall4 ## 100.038 
resultsGLMMFixedSmall4$iter ## 300
save(resultsGLMMFixedSmall4, file="GLMMFixedSmallCC.RData")

apply(tail(resultsGLMMFixedSmall4$param, 20), 2, mean, trim=0.2)
#0.8482261  0.2738440 -1.5956676  0.8475491  1.7945180  1.8507742

# 3. Adadelta ----------------------------------------
ptm <- proc.time()
resultsAdadelta <- computeMLE(glmmMod, paramNodes,
                              method="adadelta", paramInit=init,
                              compiledFuns=compiledFunsglmm,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeGLMMAdadelta <- proc.time() - ptm 
timeGLMMAdadelta ## 101.764 
save(resultsAdadelta, file="GLMMAdadelta.RData")

apply(tail(resultsAdadelta$param, 20), 2, mean, trim=0.2)
# 1.0754894  0.3623355 -1.9918989  1.0855375  1.0004103  0.4887503

ptm <- proc.time()
resultsAdadelta2 <- computeMLE(glmmMod, paramNodes,
                              method="adadelta", paramInit=init,
                              compiledFuns=compiledFunsglmm,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeGLMMAdadelta2 <- proc.time() - ptm  
timeGLMMAdadelta2 ## 97.762 

resultsAdadelta2$iter ## 285
save(resultsAdadelta2, file="GLMMAdadeltaCC.RData")

apply(tail(resultsAdadelta2$param, 20), 2, mean, trim=0.2)
#1.0380292  0.3586400 -1.9700764  1.0516159  0.8552225  0.6292048

# 4. Adam ----------------------------------------
ptm <- proc.time()
resultsAdam <- computeMLE(glmmMod, paramNodes,
                          method="adam", paramInit=init,
                          compiledFuns=compiledFunsglmm,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=300,
                          boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeGLMMAdam <- proc.time() - ptm   
timeGLMMAdam ## 103.482 
save(resultsAdam, file="GLMMAdam.RData")

apply(tail(resultsAdam$param, 20), 2, mean, trim=0.2)
#0.8357023  0.2757600 -1.5892911  0.8443124  2.8834021  1.4572084

ptm <- proc.time()
resultsAdam2 <- computeMLE(glmmMod, paramNodes,
                          method="adam", paramInit=init,
                          compiledFuns=compiledFunsglmm,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=300,
                          boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeGLMMAdam2 <- proc.time() - ptm 
timeGLMMAdam2 ##  102.562 
resultsAdam2$iter ## 300
save(resultsAdam2, file="GLMMAdamCC.RData")

apply(tail(resultsAdam2$param, 20), 2, mean, trim=0.2)
#0.8814696  0.2821036 -1.6574235  0.8878936  1.4115316  1.4581775

# 5. Newton-Raphson ----------------------------------------
ptm <- proc.time()
resultsNR <- computeMLE(glmmMod, paramNodes=paramNodes,
                        method="NR", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        tol=1e-20,
                        maxIter=300,
                        boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)
timeGLMMNR <- proc.time() - ptm 
## error


ptm <- proc.time()
resultsNR2 <- computeMLE(glmmMod, paramNodes=paramNodes,
                        method="NR", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        tol=1e-20,
                        maxIter=300,
                        boundary=boundary,trackEffSizeGrad=F,skipConvCheck=F)
timeGLMMNR2 <- proc.time() - ptm 
## error

# 6. 1-D sampling ----------------------------------------
ptm <- proc.time()
results1D <- computeMLE(glmmMod, paramNodes,
                        method="ga1D", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        maxIter=300,skipConvCheck=T)
timeGLMM1D <- proc.time() - ptm 
timeGLMM1D ## 550.234 
save(results1D, file="GLMM1D.Rdata")

apply(tail(results1D$param, 20), 2, mean, trim=0.2)
#0.8918794  0.2916108 -1.6772665  0.8583025  0.8958334  5.1852529

ptm <- proc.time()
results1D2 <- computeMLE(glmmMod, paramNodes,
                        method="ga1D", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        maxIter=300,skipConvCheck=T)
timeGLMM1D2 <- proc.time() - ptm 
timeGLMM1D2 ## 545.236 
results1D2$iter ## 300
save(results1D2, file="GLMM1DCC.Rdata")

apply(tail(results1D2$param, 20), 2, mean, trim=0.2)
## 0.7640408  0.2326564 -1.4693103  0.7794910  3.8002812  4.3891180


ptm <- proc.time()
sal <- glmm(Mate ~ 0 + Cross, random=list(~ 0 + Female,
                                          ~ 0 + Male), 
            varcomps.names=c("F", "M"), data=salamander,
            family.glmm=bernoulli.glmm, m=10^5, debug=TRUE)
glmmTime <- proc.time() - ptm 
glmmTime ## 813.043 
summary(sal)
save(sal, file="glmmMod.RData")

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
sal2=glmer(Mate~ -1 + Cross + (1|Female) + (1|Male),
           data=salamander, family=binomial)
glmerTime <- proc.time() - ptm 
glmerTime ## 0.404  
summary(sal2)
save(sal2, file="glmerMod.RData")

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
load(file="GLMMAdadelta.RData")
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
# iteratesLogregNR <- data.frame(
#   iter=1:nrow(resultsLogregNR$param),
#   resultsLogregNR$param,
#   method=rep("Newton-Raphson", nrow(resultsLogregNR$param))
# )
iterates1D <- data.frame(
  iter=1:nrow(results1D$param),
  results1D$param,
  method=rep("1D sampling", nrow(results1D$param))
)
iteratesAdadelta <- data.frame(
  iter=1:nrow(resultsAdadelta$param),
  resultsAdadelta$param,
  method=rep("Adadelta", nrow(resultsAdadelta$param))
)
iteratesAdam <- data.frame(
  iter=1:nrow(resultsAdam$param),
  resultsAdam$param, 
  method=rep("Adam", nrow(resultsAdam$param))
)

iterates <- rbind(iteratesFixed, iteratesSmallFixed,
                  #iteratesLogregNR,
                  iteratesAdadelta, iteratesAdam,
                  iterates1D) 
names(iterates) <- c("iter", paramNodes, "method")

write.csv(iterates, "glmmResultsAll.csv", row.names=F)

iteratesMelted <- melt(iterates, id.vars=c("iter", "method"))
levels(iteratesMelted$variable) <- c(expression(beta[1]),
                                     expression(beta[2]),
                                     expression(beta[3]),
                                     expression(beta[4]),
                                     "varM", "varF")
library(RColorBrewer)
cbPalette <- brewer.pal(8, "Dark2")


trajectoryPlot <- ggplot(iteratesMelted, aes(x=iter, y=value, colour=method)) 
trajectoryPlot <- trajectoryPlot + 
  geom_line(alpha=0.8) + 
  geom_point(data=iteratesMelted, mapping=aes(x=iter, y=value, shape=method))
trajectoryPlot<- trajectoryPlot + 
  facet_grid(variable ~ ., labeller=label_parsed, scales="free_y")
trajectoryPlot <- trajectoryPlot + theme(legend.position="none")
trajectoryPlot <- trajectoryPlot + 
  scale_color_manual(values=cbPalette)
trajectoryPlot <- trajectoryPlot + xlab("number of iterations")
trajectoryPlot<- trajectoryPlot + 
  theme(axis.title=element_text(face="bold", size=20),
        axis.text=element_text(size=15),
        strip.text.y=element_text(size=15))

pdf(paste0("glmm_", 300, ".pdf"))
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
