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
                               boundary=boundary)
timeGLMMFixed <- proc.time() - ptm  ##  103.498 
save(resultsGLMMFixed, file="GLMMFixed.RData")

apply(tail(resultsGLMMFixed$param, 100), 2, mean, trim=0.2)

# 2. Small fixed step size ----------------------------------------
ptm <- proc.time()
resultsGLMMFixedSmall <- computeMLE(glmmMod, paramNodes,
                                    method="fixed", paramInit=init,
                                    stepsize=0.005,
                                    compiledFuns=compiledFunsglmm,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary)
timeGLMMFixedSmall <- proc.time() - ptm  ## 104.788 
save(resultsGLMMFixedSmall, file="GLMMFixedSmall.RData")

apply(tail(resultsGLMMFixedSmall$param, 100), 2, mean, trim=0.2)

# 3. Adadelta ----------------------------------------
ptm <- proc.time()
resultsAdadelta <- computeMLE(glmmMod, paramNodes,
                              method="adadelta", paramInit=init,
                              compiledFuns=compiledFunsglmm,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary)
timeGLMMAdadelta <- proc.time() - ptm  ## 107.366 
save(resultsAdadelta, file="GLMMAdadelta.RData")

apply(tail(resultsAdadelta$param, 100), 2, mean, trim=0.2)

# 4. Adam ----------------------------------------
ptm <- proc.time()
resultsAdam <- computeMLE(glmmMod, paramNodes,
                          method="adam", paramInit=init,
                          compiledFuns=compiledFunsglmm,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=300,
                          boundary=boundary)
timeGLMMAdam <- proc.time() - ptm   ## 107.962 
save(resultsAdam, file="GLMMAdam.RData")

apply(tail(resultsAdam$param, 100), 2, mean, trim=0.2)

# 5. Newton-Raphson ----------------------------------------
ptm <- proc.time()
resultsNR <- computeMLE(glmmMod, paramNodes=paramNodes,
                        method="NR", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        tol=1e-20,
                        maxIter=300,
                        boundary=boundary)
timeGLMMNR <- proc.time() - ptm 
## error

# 6. 1-D sampling ----------------------------------------
ptm <- proc.time()
results1D <- computeMLE(glmmMod, paramNodes,
                        method="ga1D", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        maxIter=300)
timeGLMM1D <- proc.time() - ptm ## 511.821 
save(results1D, file="GLMM1D.Rdata")

mean(results1D$param[201:300,1], trim=0.2)  ## 0.8590177
mean(results1D$param[201:300,2], trim=0.2)  ## 0.2767505
mean(results1D$param[201:300,3], trim=0.2)  ## -1.606323
mean(results1D$param[201:300,4], trim=0.2)  ## 0.8699418
mean(results1D$param[201:300,5], trim=0.2)  ##  4.032243
mean(results1D$param[201:300,6], trim=0.2)  ## 1.394891

ptm <- proc.time()
sal <- glmm(Mate ~ 0 + Cross, random=list(~ 0 + Female,
                                          ~ 0 + Male), 
            varcomps.names=c("F", "M"), data=salamander,
            family.glmm=bernoulli.glmm, m=10^5, debug=TRUE)
glmmTime <- proc.time() - ptm 
glmmTime ## 1181.431 
summary(sal)
save(sal, file="glmmMod.RData")

# Call:
#   glmm(fixed=Mate ~ 0 + Cross, random=list(~0 + Female, ~0 + 
#        Male), varcomps.names=c("F", "M"), data=salamander, 
#        family.glmm=bernoulli.glmm, m=10^5, debug=TRUE)
# 
# 
# Link is: "logit (log odds)"
# 
# Fixed Effects:
#   Estimate Std. Error z value Pr(>|z|)    
# CrossR/R   0.9880     0.3624   2.726  0.00641 ** 
#   CrossR/W   0.3268     0.4103   0.797  0.42564    
# CrossW/R  -1.9791     0.4855  -4.077 4.57e-05 ***
#   CrossW/W   1.0599     0.4054   2.614  0.00894 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# Variance Components for Random Effects (P-values are one-tailed):
#   Estimate Std. Error z value Pr(>|z|)/2    
# F   1.3082     0.3665   3.570   0.000179 ***
#   M   1.4909     0.6052   2.464   0.006875 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ptm <- proc.time()
sal2=glmer(Mate~ -1 + Cross + (1|Female) + (1|Male),
           data=salamander, family=binomial)
glmerTime <- proc.time() - ptm 
glmerTime ## 0.466 
summary(sal2)
save(sal2, file="glmerMod.RData")

# Generalized linear mixed model fit by maximum likelihood 
# (Laplace Approximation) ['glmerMod']
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
# 
# 

set.seed(32856)
source("MCEM_with_output.R")
glmmNew <- glmmMod$newModel()
latentNodes <- glmmNew$getNodeNames(latentOnly=TRUE, stochOnly=TRUE)

glmmMCEM <- buildMCEM(model=glmmNew, latentNodes=latentNodes,
                      theta0=rep(2, 6)) 

ptm <- proc.time()
resultsMCEM <- glmmMCEM()
timeGLMMMCEM <- proc.time() - ptm

save(resultsMCEM, "MCEM_GLMM.RData")

timeMCEM 


## killed it after
## iter 42
# beta1      beta2      beta3      beta4       varF       varM 
# 0.8439873  0.2744783 -1.5834728  0.8501573  1.8215415  1.8495618 
# Convergence Criterion: 0.008104164.
# user   system  elapsed 
# 9347.784   22.710 9393.257

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


