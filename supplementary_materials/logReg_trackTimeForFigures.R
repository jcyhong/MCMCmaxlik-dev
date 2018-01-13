# Supplementary materials for "Sampling-Based Approaches to Maximum 
# Likelihood Estimation for Latent Variable Models":
# Code for the logistic regression example (3 parameters)
#
# Description:
# The following code contains the numerical experiments for the logistic 
# regression model: fixed step size, small fixed step size, adadelta, adam, 
# Newton-Raphson, and 1-D sampling. MCEM (from the R package NIMBLE) is used
# as a benchmark.
#
##########################################################################
setwd("~/Desktop/2018-01-12_johnny")
alg=c("fixed","fixed small","adadelta","adam","NR","1D")
execution.iter=execution.time=convergence.time=convergence.iter=rep(NA,length(alg))
logRegTimeIterInfo=cbind.data.frame(alg,execution.iter,execution.time,convergence.iter,convergence.time)
logRegTimeIterInfo$alg=as.character(logRegTimeIterInfo$alg)
# Fix the seed for reproducibility.
set.seed(21745)

# Load the packages ---------------------------------------
library("MCMCmaxlik")

# Model specification -------------------------------------
code <- nimbleCode({
  beta0 ~ dnorm(0, sd=10000)
  beta1 ~ dnorm(0, sd=10000)
  sigma_RE ~ dunif(0, 1000)
  for(i in 1:N) {
    beta2[i] ~ dnorm(0, sd=sigma_RE)
    logit(p[i]) <- beta0 + beta1 * x[i] + beta2[i]
    r[i] ~ dbin(p[i], n[i])
  }
})

## constants, data, and initial values
constants <- list(N=10)

data <- list(
  r=c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46),
  n=c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79),
  x=c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
)

inits <- list(beta0=0, beta1=0, sigma_RE=1)

logreg <- nimbleModel(code=code, constants=constants, 
                      data=data, inits=inits, check=FALSE)

paramNodesLogreg <- logreg$getNodeNames(topOnly=T)
Clogreg <- compileNimble(logreg)


# Testing Algorithms ------------------------------------------------------
# Compile the necessary functions.
compiledFunsLogreg <- buildMCMCmaxlik(logreg, paramNodesLogreg)

init <- c(0, 0, 1)
boundary <- list(c(-40, 40), c(-40, 40), c(0.05, 10))
numMCMCSamples <- 300

# 1. Fixed step size ----------------------------------------

resultsLogregFixed <- computeMLE(logreg, paramNodesLogreg,
                                 method="fixed", paramInit=init,
                                 stepsize=0.05,
                                 compiledFuns=compiledFunsLogreg,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=300,
                                 boundary=boundary)

mean(tail(resultsLogregFixed$param[,1],20), trim=.2)  ## -0.856524
mean(tail(resultsLogregFixed$param[,2],20), trim=.2)  ## 1.183854
mean(tail(resultsLogregFixed$param[,3],20), trim=.2)  ## 1.175493
save(resultsLogregFixed, file="logRegFixed.RData")

logRegTimeIterInfo$execution.iter[1]=resultsLogregFixed$execution.iter
logRegTimeIterInfo$execution.time[1]=resultsLogregFixed$execution.time[3]
logRegTimeIterInfo$convergence.iter[1]=ifelse(is.null(resultsLogregFixed$convergence.iter),NA,resultsLogregFixed$convergence.iter)
logRegTimeIterInfo$convergence.time[1]=ifelse(is.null(resultsLogregFixed$convergence.time),NA,resultsLogregFixed$convergence.time[3])



# 2. Small fixed step size ----------------------------------------

resultsLogregSmallFixed <- computeMLE(logreg, paramNodesLogreg,
                                      method="fixed", paramInit=init,
                                      stepsize=0.005,
                                      compiledFuns=compiledFunsLogreg,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary)

mean(tail(resultsLogregSmallFixed$param[,1],20), trim=.2)  ## -0.5468611
mean(tail(resultsLogregSmallFixed$param[,2],20), trim=.2)  ##  1.304298
mean(tail(resultsLogregSmallFixed$param[,3],20), trim=.2)  ## 0.2528262
save(resultsLogregSmallFixed, file="logRegSmallFixed.RData")

logRegTimeIterInfo$execution.iter[2]=resultsLogregSmallFixed$execution.iter
logRegTimeIterInfo$execution.time[2]=resultsLogregSmallFixed$execution.time[3]
logRegTimeIterInfo$convergence.iter[2]=ifelse(is.null(resultsLogregSmallFixed$convergence.iter),NA,resultsLogregSmallFixed$convergence.iter)
logRegTimeIterInfo$convergence.time[2]=ifelse(is.null(resultsLogregSmallFixed$convergence.time),NA,resultsLogregSmallFixed$convergence.time[3])


# 3. Adadelta ----------------------------------------

resultsLogregAdadelta <- computeMLE(logreg, paramNodesLogreg,
                                    method="adadelta", paramInit=init,
                                    compiledFuns=compiledFunsLogreg,
                                    numMCMCSamples=numMCMCSamples,
                                    maxIter=300,
                                    boundary=boundary)

mean(tail(resultsLogregAdadelta$param[, 1], 20), trim=.2)  ## -0.575974
mean(tail(resultsLogregAdadelta$param[, 2], 20), trim=.2)  ## 1.3102
mean(tail(resultsLogregAdadelta$param[, 3], 20), trim=.2)  ##  0.3415709
save(resultsLogregAdadelta, file="logRegAdadelta.RData")


logRegTimeIterInfo$execution.iter[3]=resultsLogregAdadelta$execution.iter
logRegTimeIterInfo$execution.time[3]=resultsLogregAdadelta$execution.time[3]
logRegTimeIterInfo$convergence.iter[3]=ifelse(is.null(resultsLogregAdadelta$convergence.iter),NA,resultsLogregAdadelta$convergence.iter)
logRegTimeIterInfo$convergence.time[3]=ifelse(is.null(resultsLogregAdadelta$convergence.time),NA,resultsLogregAdadelta$convergence.time[3])


# 4. Adam ----------------------------------------
ptm <- proc.time()
resultsLogregAdam <- computeMLE(logreg, paramNodesLogreg,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                stepsize=0.2,
                                eps=1e-4,
                                maxIter=300,
                                boundary=boundary)
timeLogRegAdam <- proc.time() - ptm 
timeLogRegAdam ## 3.005 
mean(tail(resultsLogregAdam$param[, 1], 20), trim=.2) ## -0.5575068
mean(tail(resultsLogregAdam$param[, 2], 20), trim=.2) ##  1.313344
mean(tail(resultsLogregAdam$param[, 3], 20), trim=.2) ##  0.2433034
save(resultsLogregAdam, file="logRegAdam.RData")

logRegTimeIterInfo$execution.iter[4]=resultsLogregAdam$execution.iter
logRegTimeIterInfo$execution.time[4]=resultsLogregAdam$execution.time[3]
logRegTimeIterInfo$convergence.iter[4]=ifelse(is.null(resultsLogregAdam$convergence.iter),NA,resultsLogregAdam$convergence.iter)
logRegTimeIterInfo$convergence.time[4]=ifelse(is.null(resultsLogregAdam$convergence.time),NA,resultsLogregAdam$convergence.time[3])




# 5. Newton-Raphson ----------------------------------------

resultsLogregNR <- computeMLE(logreg, paramNodesLogreg,
                              method="NR", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=numMCMCSamples,
                              maxIter=300,
                              boundary=boundary)
#Error in solve.default(approxHessian, gradCurr) : 
#  system is computationally singular: reciprocal condition number = 1.16616e-22

mean(tail(resultsLogregNR$param[,1],20), trim=.2) ##  -6.389591
mean(tail(resultsLogregNR$param[,2],20), trim=.2) ## -10.43822
mean(tail(resultsLogregNR$param[,3],20), trim=.2) ## 10 (boundary)

save(resultsLogregNR, file="logRegNR.RData")

logRegTimeIterInfo$execution.iter[5]=resultsLogregNR$execution.iter
logRegTimeIterInfo$execution.time[5]=resultsLogregNR$execution.time[3]
logRegTimeIterInfo$convergence.iter[5]=ifelse(is.null(resultsLogregNR$convergence.iter),NA,resultsLogregNR$convergence.iter)
logRegTimeIterInfo$convergence.time[5]=ifelse(is.null(resultsLogregNR$convergence.time),NA,resultsLogregNR$convergence.time[3])



# 6. 1-D sampling ----------------------------------------

resultsLogreg1D <- computeMLE(logreg, paramNodesLogreg,
                              method="ga1D", paramInit=init,
                              compiledFuns=compiledFunsLogreg,
                              numMCMCSamples=300, numMCMCSamples1D=300, 
                              maxIter=300)

mean(tail(resultsLogreg1D$param[, 1]), trim=.2)  ## -0.5133427
mean(tail(resultsLogreg1D$param[, 2]), trim=.2)  ## 1.278134
mean(tail(resultsLogreg1D$param[, 3]), trim=.2)  ## 0.2491426
save(resultsLogreg1D, file="logReg1D.RData")

logRegTimeIterInfo$execution.iter[6]=resultsLogreg1D$execution.iter
logRegTimeIterInfo$execution.time[6]=resultsLogreg1D$execution.time[3]
logRegTimeIterInfo$convergence.iter[6]=ifelse(is.null(resultsLogreg1D$convergence.iter),NA,resultsLogreg1D$convergence.iter)
logRegTimeIterInfo$convergence.time[6]=ifelse(is.null(resultsLogreg1D$convergence.time),NA,resultsLogreg1D$convergence.time[3])


write.csv(logRegTimeIterInfo,"logRegTimeIterInfo.csv",row.names=F)


require(ggplot2)
require(gridExtra)


g1<-ggplot(data=logRegTimeIterInfo,aes(x=alg,y=execution.time))+geom_bar(stat="identity")+
  xlab("algorithm")+ylab("execution time (s)")+ggtitle("Logistic Regression Example \n Maximum Number of Iterations: 300")

g2<-ggplot(data=logRegTimeIterInfo,aes(x=alg,y=convergence.time))+geom_bar(stat="identity")+
  xlab("algorithm")+ylab("convergence time (s)")
#+ggtitle("Pump Example \n Maximum Number of Iterations: 300")


#grid.arrange(g1,g2,g3,g4,ncol=2)
grid.arrange(g1,g2,ncol=2)


par(mfrow=c(3,1))
plot(resultsLogregFixed$param[,1],xlim=c(0,350),main="Logistic Regression Example",xlab="iter",ylab=expression(beta[0]),type="l")
lines(resultsLogregSmallFixed$param[,1],col="red",lwd=2)
lines(resultsLogregAdadelta$param[,1],col="blue")
lines(resultsLogregAdam$param[,1],col="forestgreen",lwd=3)
lines(resultsLogreg1D$param[,1],col="goldenrod",lwd=2)
lines(resultsLogregNR$param[,1],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)

plot(resultsLogregFixed$param[,2],xlim=c(0,350),main="Logistic Regression Example",xlab="iter",ylab=expression(beta[1]),type="l")
lines(resultsLogregSmallFixed$param[,2],col="red",lwd=2)
lines(resultsLogregAdadelta$param[,2],col="blue")
lines(resultsLogregAdam$param[,2],col="forestgreen",lwd=3)
lines(resultsLogreg1D$param[,2],col="goldenrod",lwd=2)
lines(resultsLogregNR$param[,2],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)

plot(resultsLogregFixed$param[,3],xlim=c(0,350),main="Logistic Regression Example",xlab="iter",ylab=expression(sigma[RE]),type="l")
lines(resultsLogregSmallFixed$param[,3],col="red",lwd=2)
lines(resultsLogregAdadelta$param[,3],col="blue")
lines(resultsLogregAdam$param[,3],col="forestgreen",lwd=3)
lines(resultsLogreg1D$param[,3],col="goldenrod",lwd=2)
lines(resultsLogregNR$param[,3],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod","magenta"),
       c("fixed","fixed small","adadelta","adam","1D","NR"),lty=1,lwd=2)




