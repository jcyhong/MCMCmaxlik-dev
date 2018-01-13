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
setwd("~/Desktop/2018-01-12_johnny")
alg=c("fixed","fixed small","adadelta","adam","NR","1D")
execution.iter=execution.time=convergence.time=convergence.iter=rep(NA,length(alg))
glmmTimeIterInfo=cbind.data.frame(alg,execution.iter,execution.time,convergence.iter,convergence.time)
glmmTimeIterInfo$alg=as.character(glmmTimeIterInfo$alg)
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
resultsGLMMFixed <- computeMLE(glmmMod, paramNodes,
                               method="fixed", paramInit=init,
                               stepsize=0.05,
                               compiledFuns=compiledFunsglmm,
                               numMCMCSamples=numMCMCSamples,
                               maxIter=300,
                               boundary=boundary)

save(resultsGLMMFixed, file="GLMMFixed.RData")

apply(tail(resultsGLMMFixed$param, 20), 2, mean, trim=0.2)
#0.8371804  0.2690307 -1.5762730  0.8357959  1.7182524  2.2435778

glmmTimeIterInfo$execution.iter[1]=resultsGLMMFixed$execution.iter
glmmTimeIterInfo$execution.time[1]=resultsGLMMFixed$execution.time[3]
glmmTimeIterInfo$convergence.iter[1]=ifelse(is.null(resultsGLMMFixed$convergence.iter),NA,resultsGLMMFixed$convergence.iter)
glmmTimeIterInfo$convergence.time[1]=ifelse(is.null(resultsGLMMFixed$convergence.time),NA,resultsGLMMFixed$convergence.time[3])


# 2. Small fixed step size ----------------------------------------

resultsGLMMFixedSmall <- computeMLE(glmmMod, paramNodes,
                                    method="fixed", paramInit=init,
                                    stepsize=0.005,
                                    compiledFuns=compiledFunsglmm,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary)

save(resultsGLMMFixedSmall, file="GLMMFixedSmall.RData")

apply(tail(resultsGLMMFixedSmall$param, 20), 2, mean, trim=0.2)
##0.8816105  0.2855453 -1.6549337  0.8874643  1.5867583  1.3504009

glmmTimeIterInfo$execution.iter[2]=resultsGLMMFixedSmall$execution.iter
glmmTimeIterInfo$execution.time[2]=resultsGLMMFixedSmall$execution.time[3]
glmmTimeIterInfo$convergence.iter[2]=ifelse(is.null(resultsGLMMFixedSmall$convergence.iter),NA,resultsGLMMFixedSmall$convergence.iter)
glmmTimeIterInfo$convergence.time[2]=ifelse(is.null(resultsGLMMFixedSmall$convergence.time),NA,resultsGLMMFixedSmall$convergence.time[3])


# 3. Adadelta ----------------------------------------
resultsAdadelta <- computeMLE(glmmMod, paramNodes,
                              method="adadelta", paramInit=init,
                              compiledFuns=compiledFunsglmm,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary)

save(resultsAdadelta, file="GLMMAdadelta.RData")

apply(tail(resultsAdadelta$param, 20), 2, mean, trim=0.2)
# 1.0754894  0.3623355 -1.9918989  1.0855375  1.0004103  0.4887503

glmmTimeIterInfo$execution.iter[3]=resultsAdadelta$execution.iter
glmmTimeIterInfo$execution.time[3]=resultsAdadelta$execution.time[3]
glmmTimeIterInfo$convergence.iter[3]=ifelse(is.null(resultsAdadelta$convergence.iter),NA,resultsAdadelta$convergence.iter)
glmmTimeIterInfo$convergence.time[3]=ifelse(is.null(resultsAdadelta$convergence.time),NA,resultsAdadelta$convergence.time[3])

# 4. Adam ----------------------------------------

resultsAdam <- computeMLE(glmmMod, paramNodes,
                          method="adam", paramInit=init,
                          compiledFuns=compiledFunsglmm,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=300,
                          boundary=boundary)

save(resultsAdam, file="GLMMAdam.RData")

apply(tail(resultsAdam$param, 20), 2, mean, trim=0.2)
#0.8357023  0.2757600 -1.5892911  0.8443124  2.8834021  1.4572084

glmmTimeIterInfo$execution.iter[4]=resultsAdam$execution.iter
glmmTimeIterInfo$execution.time[4]=resultsAdam$execution.time[3]
glmmTimeIterInfo$convergence.iter[4]=ifelse(is.null(resultsAdam$convergence.iter),NA,resultsAdam$convergence.iter)
glmmTimeIterInfo$convergence.time[4]=ifelse(is.null(resultsAdam$convergence.time),NA,resultsAdam$convergence.time[3])

# 5. Newton-Raphson ----------------------------------------

resultsNR <- computeMLE(glmmMod, paramNodes=paramNodes,
                        method="NR", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        tol=1e-20,
                        maxIter=300,
                        boundary=boundary)

glmmTimeIterInfo$execution.iter[5]=resultsNR$execution.iter
glmmTimeIterInfo$execution.time[5]=resultsNR$execution.time[3]
glmmTimeIterInfo$convergence.iter[5]=ifelse(is.null(resultsNR$convergence.iter),NA,resultsNR$convergence.iter)
glmmTimeIterInfo$convergence.time[5]=ifelse(is.null(resultsNR$convergence.time),NA,resultsNR$convergence.time[3])


# 6. 1-D sampling ----------------------------------------

results1D <- computeMLE(glmmMod, paramNodes,
                        method="ga1D", paramInit=init,
                        compiledFuns=compiledFunsglmm,
                        numMCMCSamples=numMCMCSamples,
                        maxIter=300)

save(results1D, file="GLMM1D.Rdata")

apply(tail(results1D$param, 20), 2, mean, trim=0.2)
#0.8918794  0.2916108 -1.6772665  0.8583025  0.8958334  5.1852529

glmmTimeIterInfo$execution.iter[6]=results1D$execution.iter
glmmTimeIterInfo$execution.time[6]=results1D$execution.time[3]
glmmTimeIterInfo$convergence.iter[6]=ifelse(is.null(results1D$convergence.iter),NA,results1D$convergence.iter)
glmmTimeIterInfo$convergence.time[6]=ifelse(is.null(results1D$convergence.time),NA,results1D$convergence.time[3])

write.csv(glmmTimeIterInfo,"glmmTimeIterInfo.csv",row.names=F)

glmmTimeIterInfo$convergence.time=rep(0,nrow(glmmTimeIterInfo))

require(ggplot2)
require(gridExtra)


g1<-ggplot(data=glmmTimeIterInfo,aes(x=alg,y=execution.time))+geom_bar(stat="identity")+
  xlab("algorithm")+ylab("execution time (s)")+ggtitle("GLMM Example \n Maximum Number of Iterations: 300")

g2<-ggplot(data=glmmTimeIterInfo,aes(x=alg,y=convergence.time))+geom_bar(stat="identity")+
  xlab("algorithm")+ylab("convergence time (s)")+ylim(0,1)
#+ggtitle("Pump Example \n Maximum Number of Iterations: 300")


#grid.arrange(g1,g2,g3,g4,ncol=2)
grid.arrange(g1,g2,ncol=2)

par(mfrow=c(3,1))
plot(resultsGLMMFixed$param[,1],xlim=c(0,350),ylim=c(0,2),main="GLMM Example",xlab="iter",ylab=expression(beta[1]),type="l")
lines(resultsGLMMFixedSmall$param[,1],col="red",lwd=2)
lines(resultsAdadelta$param[,1],col="blue")
lines(resultsAdam$param[,1],col="forestgreen",lwd=3)
lines(results1D$param[,1],col="goldenrod",lwd=2)
lines(resultsNR$param[,1],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)

plot(resultsGLMMFixed$param[,2],xlim=c(0,350),ylim=c(-0.5,2),main="GLMM Example",xlab="iter",ylab=expression(beta[2]),type="l")
lines(resultsGLMMFixedSmall$param[,2],col="red",lwd=2)
lines(resultsAdadelta$param[,2],col="blue")
lines(resultsAdam$param[,2],col="forestgreen",lwd=3)
lines(results1D$param[,2],col="goldenrod",lwd=2)
lines(resultsNR$param[,2],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)

plot(resultsGLMMFixed$param[,3],xlim=c(0,350),ylim=c(-2.5,2),main="GLMM Example",xlab="iter",ylab=expression(beta[3]),type="l")
lines(resultsGLMMFixedSmall$param[,3],col="red",lwd=2)
lines(resultsAdadelta$param[,3],col="blue")
lines(resultsAdam$param[,3],col="forestgreen",lwd=3)
lines(results1D$param[,3],col="goldenrod",lwd=2)
lines(resultsNR$param[,3],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)


plot(resultsGLMMFixed$param[,4],xlim=c(0,350),ylim=c(0,2),main="GLMM Example",xlab="iter",ylab=expression(beta[4]),type="l")
lines(resultsGLMMFixedSmall$param[,4],col="red",lwd=2)
lines(resultsAdadelta$param[,4],col="blue")
lines(resultsAdam$param[,4],col="forestgreen",lwd=3)
lines(results1D$param[,4],col="goldenrod",lwd=2)
lines(resultsNR$param[,4],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)

plot(resultsGLMMFixed$param[,5],xlim=c(0,350),ylim=c(0,3),main="GLMM Example",xlab="iter",ylab=expression(sigma^2[F]),type="l")
lines(resultsGLMMFixedSmall$param[,5],col="red",lwd=2)
lines(resultsAdadelta$param[,5],col="blue")
lines(resultsAdam$param[,5],col="forestgreen",lwd=3)
lines(results1D$param[,5],col="goldenrod",lwd=2)
lines(resultsNR$param[,5],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)

plot(resultsGLMMFixed$param[,6],xlim=c(0,350),ylim=c(0,3),main="GLMM Example",xlab="iter",ylab=expression(sigma^2[M]),type="l")
lines(resultsGLMMFixedSmall$param[,6],col="red",lwd=2)
lines(resultsAdadelta$param[,6],col="blue")
lines(resultsAdam$param[,6],col="forestgreen",lwd=3)
lines(results1D$param[,6],col="goldenrod",lwd=2)
lines(resultsNR$param[,6],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)


