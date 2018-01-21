setwd("~/Desktop/2018-01-12_johnny")
alg=c("fixed","fixed small","adadelta","adam","NR","1D")
execution.iter=execution.time=convergence.time=convergence.iter=rep(NA,length(alg))
stateSpaceTimeIterInfo=cbind.data.frame(alg,execution.iter,execution.time,convergence.iter,convergence.time)
stateSpaceTimeIterInfo$alg=as.character(stateSpaceTimeIterInfo$alg)
library("MCMCmaxlik")


set.seed(931154)
rho=0.5
sigPN=5
sigOE=5
x0=5
TT=100
x=rep(rnorm(1,rho*x0,sd=sqrt(sigPN^2/(1-rho*rho))),TT)
y=rep(rnorm(1,x[1],sd=sigOE),TT)
for(i in 2:TT){
  x[i] =rnorm(1,rho * x[i - 1] , sd = sigPN)
  y[i] =rnorm(1,x[i], sd = sigOE)
}


ssCode <- nimbleCode({
  rho ~ dunif(-0.9999, 0.9999)
  #mu ~ dnorm(0, sd = 1000) # 0
  sigPN ~ dunif(1e-04, 100)
  sigOE ~ dunif(1e-04, 100)
  x0 ~dnorm(0, sd = 1000)
  x[1] ~ dnorm(rho*x0, sd = sqrt(sigPN^2/(1-rho*rho)))
  #x[1] ~ dnorm(rho*x0, sd = sigPN)
  #x[1] ~ dnorm(mu, sd = sqrt(sigPN^2/(1-rho*rho)))
  y[1] ~ dnorm(x[1], sd = sigOE)
  for (i in 2:TT) {
    #x[i] ~ dnorm(rho * (x[i - 1] - mu) + mu, sd = sigPN)
    x[i] ~ dnorm(rho * x[i - 1] , sd = sigPN)
    y[i] ~ dnorm(x[i], sd = sigOE)
  }
}) ## take out mu

Consts <- list(TT = 100)

data <- list(
  x=x,
  y=y
)


## build the model object
ssModel <- nimbleModel(ssCode, constants = Consts,data=data)



write.csv(cbind(x,y),"statespaceModData.csv")

paramNodesSS <- ssModel$getNodeNames(topOnly=T)
CssModel <- compileNimble(ssModel)




# 
# Consts <- list(TT = 100)
# 
# ssModel <- nimbleModel(ssCode, constants = Consts)
# 
# ## put values in the top-level parameters
# ssModel$rho <- 0.5
# ssModel$x0 <- 5.0
# ssModel$sigPN <- 5
# ssModel$sigOE <- 5
# 
# ## calculate any dependencies that may be needed to simulate states and data
# ssModel$calculate(ssModel$getDependencies( c('rho','x0','sigPN','sigOE'), determOnly = TRUE))
# 
# 
# ## use the model to simulate states and data
# set.seed(12345)
# ssModel$simulate( ssModel$getDependencies(c('x','y')))
# ## Look at what we got:
# ssModel$x
# ssModel$y
# ssModel$setData('y')
# 
# 
# write.csv(cbind(ssModel$x,ssModel$y),"statespaceModData.csv")
# 
# paramNodesSS <- ssModel$getNodeNames(topOnly=T)
# CssModel <- compileNimble(ssModel)





compiledFunsSSmodel <- buildMCMCmaxlik(ssModel, paramNodesSS)



# Testing Algorithms


# 1. Fixed step size ----------------------------------------
# Compile the necessary functions.

init <- c(0.75,10,10,10)
boundary <- list(c(-0.9999, 0.9999),c(1e-04, 100),c(1e-04, 100),c(-100000,100000))
numMCMCSamples <- 300

resultsSSModFixed <- computeMLE(ssModel, paramNodesSS,
                                method="fixed", paramInit=init,
                                compiledFuns=compiledFunsSSmodel,
                                stepsize=0.05,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300,
                                boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)

mean(tail(resultsSSModFixed$param[, 1], 20), trim=.2)  #  0.5076828
mean(tail(resultsSSModFixed$param[, 2], 20), trim=.2)   # 4.086307
mean(tail(resultsSSModFixed$param[, 3], 20), trim=.2) #  6.391043
mean(tail(resultsSSModFixed$param[, 4], 20), trim=.2)  # 10.14549

save(resultsSSModFixed, file="SSModFixed.RData")

stateSpaceTimeIterInfo$execution.iter[1]=resultsSSModFixed$execution.iter
stateSpaceTimeIterInfo$execution.time[1]=resultsSSModFixed$execution.time[3]
stateSpaceTimeIterInfo$convergence.iter[1]=ifelse(is.null(resultsSSModFixed$convergence.iter),NA,resultsSSModFixed$convergence.iter)
stateSpaceTimeIterInfo$convergence.time[1]=ifelse(is.null(resultsSSModFixed$convergence.time),NA,resultsSSModFixed$convergence.time[3])


# 2. Small fixed step size ----------------------------------------
resultsSSModSmallFixed <- computeMLE(ssModel,paramNodesSS,
                                     method="fixed", paramInit=init,
                                     compiledFuns=compiledFunsSSmodel,
                                     stepsize=0.005,
                                     numMCMCSamples=numMCMCSamples,
                                     maxIter=300,
                                     boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)

mean(tail(resultsSSModSmallFixed$param[, 1], 20), trim=.2)  ## 0.5219547
mean(tail(resultsSSModSmallFixed$param[, 2], 20), trim=.2) ## 4.650584
mean(tail(resultsSSModSmallFixed$param[, 3], 20), trim=.2)  ##  5.891227
mean(tail(resultsSSModSmallFixed$param[, 4], 20), trim=.2) ## 10.35237
save(resultsSSModSmallFixed, file="SSModSmallFixed.RData")

stateSpaceTimeIterInfo$execution.iter[2]=resultsSSModSmallFixed$execution.iter
stateSpaceTimeIterInfo$execution.time[2]=resultsSSModSmallFixed$execution.time[3]
stateSpaceTimeIterInfo$convergence.iter[2]=ifelse(is.null(resultsSSModSmallFixed$convergence.iter),NA,resultsSSModSmallFixed$convergence.iter)
stateSpaceTimeIterInfo$convergence.time[2]=ifelse(is.null(resultsSSModSmallFixed$convergence.time),NA,resultsSSModSmallFixed$convergence.time[3])



# 3. Adadelta ----------------------------------------

resultsSSModAdadelta <- computeMLE(ssModel, paramNodesSS,
                                   method="adadelta", paramInit=init,
                                   compiledFuns=compiledFunsSSmodel,
                                   numMCMCSamples=numMCMCSamples,
                                   maxIter=300,
                                   boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)

mean(tail(resultsSSModAdadelta$param[, 1], 20), trim=.2)  ##  0.5224348
mean(tail(resultsSSModAdadelta$param[, 2], 20), trim=.2)  ##4.203184
mean(tail(resultsSSModAdadelta$param[, 3], 20), trim=.2)  ## 6.277823
mean(tail(resultsSSModAdadelta$param[, 4], 20), trim=.2)  ## 12.17862
save(resultsSSModAdadelta, file="SSModAdadelta.RData")

stateSpaceTimeIterInfo$execution.iter[3]=resultsSSModAdadelta$execution.iter
stateSpaceTimeIterInfo$execution.time[3]=resultsSSModAdadelta$execution.time[3]
stateSpaceTimeIterInfo$convergence.iter[3]=ifelse(is.null(resultsSSModAdadelta$convergence.iter),NA,resultsSSModAdadelta$convergence.iter)
stateSpaceTimeIterInfo$convergence.time[3]=ifelse(is.null(resultsSSModAdadelta$convergence.time),NA,resultsSSModAdadelta$convergence.time[3])





# 4. Adam ----------------------------------------

resultsSSModAdam <- computeMLE(ssModel, paramNodesSS,
                               method="adam", paramInit=init,
                               compiledFuns=compiledFunsSSmodel,
                               numMCMCSamples=numMCMCSamples,
                               maxIter=300,
                               #maxIter=500,
                               boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)

mean(tail(resultsSSModAdam$param[, 1], 20), trim=.2)  ##0.4713427
mean(tail(resultsSSModAdam$param[, 2], 20), trim=.2)  ## 4.898543
mean(tail(resultsSSModAdam$param[, 3], 20), trim=.2)  ## 5.740359
mean(tail(resultsSSModAdam$param[, 4], 20), trim=.2) ## 17.39691
save(resultsSSModAdam, file="SSModAdam.RData")

stateSpaceTimeIterInfo$execution.iter[4]=resultsSSModAdam$execution.iter
stateSpaceTimeIterInfo$execution.time[4]=resultsSSModAdam$execution.time[3]
stateSpaceTimeIterInfo$convergence.iter[4]=ifelse(is.null(resultsSSModAdam$convergence.iter),NA,resultsSSModAdam$convergence.iter)
stateSpaceTimeIterInfo$convergence.time[4]=ifelse(is.null(resultsSSModAdam$convergence.time),NA,resultsSSModAdam$convergence.time[3])



# 5. Newton-Raphson ----------------------------------------

resultsSSModNR <- computeMLE(ssModel, paramNodes=paramNodesSS,
                             method="NR", paramInit=init,
                             compiledFuns=compiledFunsSSmodel,
                             numMCMCSamples=numMCMCSamples,
                             tol=1e-20,
                             maxIter=300,
                             boundary=boundary,trackEffSizeGrad=F,skipConvCheck=T)

mean(tail(resultsSSModNR$param[, 1], 20), trim=.2) ## -0.7385711
mean(tail(resultsSSModNR$param[, 2], 20), trim=.2)  ## 100
mean(tail(resultsSSModNR$param[, 3], 20), trim=.2) ##  100
mean(tail(resultsSSModNR$param[, 4], 20), trim=.2)  ##  2631.368
save(resultsSSModNR, file="SSModNR.RData")

stateSpaceTimeIterInfo$execution.iter[5]=resultsSSModNR$execution.iter
stateSpaceTimeIterInfo$execution.time[5]=resultsSSModNR$execution.time[3]
stateSpaceTimeIterInfo$convergence.iter[5]=ifelse(is.null(resultsSSModNR$convergence.iter),NA,resultsSSModNR$convergence.iter)
stateSpaceTimeIterInfo$convergence.time[5]=ifelse(is.null(resultsSSModNR$convergence.time),NA,resultsSSModNR$convergence.time[3])




# 6. 1-D sampling ----------------------------------------
resultsSSMod1D <- computeMLE(ssModel, paramNodesSS,
                             method="ga1D", paramInit=init,
                             compiledFuns=compiledFunsSSmodel,
                             numMCMCSamples=300, numMCMCSamples1D=300, 
                             maxIter=300,skipConvCheck=T)

mean(tail(resultsSSMod1D$param[, 1], 20), trim=.2)  ## 0.5128038
mean(tail(resultsSSMod1D$param[, 2], 20), trim=.2)   ## 4.515337
mean(tail(resultsSSMod1D$param[, 3], 20), trim=.2) ## 6.117298
mean(tail(resultsSSMod1D$param[, 4], 20), trim=.2)  ##  10.16848
save(resultsSSMod1D, file="SSMod1D.RData")

stateSpaceTimeIterInfo$execution.iter[6]=resultsSSMod1D$execution.iter
stateSpaceTimeIterInfo$execution.time[6]=resultsSSMod1D$execution.time[3]
stateSpaceTimeIterInfo$convergence.iter[6]=ifelse(is.null(resultsSSMod1D$convergence.iter),NA,resultsSSMod1D$convergence.iter)
stateSpaceTimeIterInfo$convergence.time[6]=ifelse(is.null(resultsSSMod1D$convergence.time),NA,resultsSSMod1D$convergence.time[3])



write.csv(stateSpaceTimeIterInfo,"stateSpaceTimeIterInfo.csv",row.names=F)

setwd("/Users/Sara/Desktop/2018-01-12_johnny")
stateSpaceTimeIterInfo<-read.csv("stateSpaceTimeIterInfo.csv")

require(ggplot2)
require(gridExtra)

stateSpaceTimeIterInfo$convergence.time=rep(0,6)

#g1<-
ggplot(data=subset(stateSpaceTimeIterInfo,!is.na(execution.time)),aes(x=alg,y=execution.time))+geom_bar(stat="identity")+
  xlab("algorithm")+ylab("execution time (s)")+ggtitle("State Space Model Example \n Maximum Number of Iterations: 300")

#g2<-ggplot(data=stateSpaceTimeIterInfo,aes(x=alg,y=convergence.time))+geom_bar(stat="identity")+
#  xlab("algorithm")+ylab("convergence time (s)")
#+ggtitle("Pump Example \n Maximum Number of Iterations: 300")


#grid.arrange(g1,g2,g3,g4,ncol=2)
#grid.arrange(g1,g2,ncol=2)

load("SSModFixed.RData")
load("SSModSmallFixed.RData")
load("SSModAdadelta.RData")
load("SSModAdam.RData")
load("SSMod1D.RData")

par(mfrow=c(2,1))
plot(resultsSSModFixed$param[,1],xlim=c(0,400),ylim=c(0,2),main="State Space Model Example",xlab="iter",ylab=expression(rho),type="l")
lines(resultsSSModSmallFixed$param[,1],col="red",lwd=2)
lines(resultsSSModAdadelta$param[,1],col="blue")
lines(resultsSSModAdam$param[,1],col="forestgreen",lwd=3)
lines(resultsSSMod1D$param[,1],col="goldenrod",lwd=2)
#lines(resultsNR$param[,1],col="magenta",lwd=2)
legend("topright",col=c("black","red","blue","forestgreen","goldenrod"#,"magenta"
                           ),
       c("fixed","fixed small","adadelta","adam","1D"#,"NR"
         ),lty=1,lwd=2)

plot(resultsSSModFixed$param[,2],xlim=c(0,400),ylim=c(0,15),main="State Space Model Example",xlab="iter",ylab=expression(sigma[PN]),type="l")
lines(resultsSSModSmallFixed$param[,2],col="red",lwd=2)
lines(resultsSSModAdadelta$param[,2],col="blue")
lines(resultsSSModAdam$param[,2],col="forestgreen",lwd=3)
lines(resultsSSMod1D$param[,2],col="goldenrod",lwd=2)
#lines(resultsNR$param[,2],col="magenta",lwd=2)
legend("topright",col=c("black","red","blue","forestgreen","goldenrod"#,"magenta"
                           ),
       c("fixed","fixed small","adadelta","adam","1D"#,"NR"
         ),lty=1,lwd=2)

plot(resultsSSModFixed$param[,3],xlim=c(0,400),ylim=c(4,10),main="State Space Model Example",xlab="iter",ylab=expression(sigma[OE]),type="l")
lines(resultsSSModSmallFixed$param[,3],col="red",lwd=2)
lines(resultsSSModAdadelta$param[,3],col="blue")
lines(resultsSSModAdam$param[,3],col="forestgreen",lwd=3)
lines(resultsSSMod1D$param[,3],col="goldenrod",lwd=2)
#lines(resultsNR$param[,3],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod"#,"magenta"
                           ),
       c("fixed","fixed small","adadelta","adam","1D"#,"NR"
         ),lty=1,lwd=2)


plot(resultsSSModFixed$param[,4],xlim=c(0,400),ylim=c(-5,15),main="State Space Model Example",xlab="iter",ylab=expression(x[0]),type="l")
lines(resultsSSModSmallFixed$param[,4],col="red",lwd=2)
lines(resultsSSModAdadelta$param[,4],col="blue")
lines(resultsSSModAdam$param[,4],col="forestgreen",lwd=3)

lines(resultsSSMod1D$param[,4],col="goldenrod",lwd=2)
#lines(resultsNR$param[,4],col="magenta",lwd=2)
legend("bottomright",col=c("black","red","blue","forestgreen","goldenrod"#,"magenta"
                           ),
       c("fixed","fixed small","adadelta","adam","1D"#,"NR"
         ),lty=1,lwd=2)
