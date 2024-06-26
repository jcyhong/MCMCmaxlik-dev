# Simulation studies for different MCMC sample sizes and scientific sample sizes.
# Model: salamander (glmm with cross random effects)

# Load the packages ---------------------------------------
library("MCMCmaxlik")
require(glmm)
require(nimble)
library(reshape2) 
library(ggplot2)
library(gridExtra)
library(grid)
require(lme4)
library(dplyr)
nimbleOptions(experimentalEnableDerivs = TRUE)

# Model specification -------------------------------------
salamander3 = read.csv('salamander_new/first30F_n180.csv.csv')
N = nrow(salamander3) 

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

glmmConstants <- list(
  N = nrow(salamander3),
  numFemales=length(unique(salamander3$Findex)),
  numMales=length(unique(salamander3$Mindex)),
  Findex=salamander3$Findex,
  Mindex=salamander3$Mindex
)

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
maxIter <- 100

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
  
  # 3. Adam ----------------------------------------
  resultsAdam <- computeMLE(glmmMod, glmmParamNodes,
                            method="adam", paramInit=init,
                            compiledFuns=compiledFunsGlmm,
                            numMCMCSamples=numMCMCSamples,
                            maxIter=maxIter,
                            boundary=boundary)
  
  # 4. Newton-Raphson ----------------------------------------
  resultsNR <- computeMLE(glmmMod, glmmParamNodes,
                          method="NR", paramInit=init,
                          compiledFuns=compiledFunsGlmm,
                          numMCMCSamples=numMCMCSamples,
                          tol=1e-20,
                          maxIter=maxIter,
                          boundary=boundary)
  
  # 5. 1-D sampling ----------------------------------------
  results1D <- computeMLE(glmmMod, glmmParamNodes,
                          method="ga1D", paramInit=init,
                          compiledFuns=compiledFunsGlmm,
                          numMCMCSamples=numMCMCSamples,
                          numMCMCSamples1D=numMCMCSamples,
                          maxIter=maxIter, 
                          boundary=boundary)
  
  list(resultsGLMMFixed, resultsGLMMFixedSmall,
       resultsAdam, resultsNR, results1D)
}

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

# Fix the seed for reproducibility.
set.seed(21745)

resultsGLMMList <- list()
numMCMCSamplesGrid <- c(50, 100, 200, 300)
for (i in 1:length(numMCMCSamplesGrid)) {
  resultsGLMMList[[i]] <- getGLMMResults(
    init=rep(2, 6), numMCMCSamples=numMCMCSamplesGrid[i])
}

timesGLMMList <- lapply(resultsGLMMList, getGLMMSummary)

flattenTimesList <- function(timesList, methodNames, numMCMCSamples) {
  do.call(rbind, lapply(1:length(timesList), function(i) {
    data.frame(
      method=factor(methodNames, levels=methodNames), 
      MCMCSampleSize=numMCMCSamples[i], 
      exec.time=timesList[[i]][, 'exec.time']
    )
  }))
}

flattenParam <- function(paramMatrix, paramNames) {
  numRows <- nrow(paramMatrix)
  numCols <- ncol(paramMatrix)
  iteration <- rep(0:(numRows - 1), numCols)
  parameter <- factor(rep(paramNames, each=numRows), levels=paramNames)
  data.frame(iteration=iteration, parameter=parameter, estimate=as.vector(paramMatrix))
}

flattenResult <- function(result, paramNames, methodNames) {
  do.call(rbind, lapply(1:length(result), function(i) {
    df <- flattenParam(result[[i]]$param, paramNames=paramNames)
    df$method <- methodNames[i]
    df
  })
  ) %>% mutate(method=factor(method, levels=methodNames))
}

flattenResultList <- function(resultList, paramNames, methodNames, numMCMCSamples) {
  do.call(rbind, lapply(1:length(resultList), function(i) {
    df <- flattenResult(resultList[[i]], paramNames=paramNames, methodNames=methodNames)
    df$MCMCSampleSize <- numMCMCSamples[i]
    df
  })
  ) %>% mutate(MCMCSampleSize=factor(MCMCSampleSize, levels=numMCMCSamples))
}

methodNames <- c('Fixed (0.05)', 'Fixed (0.005)', 'Adam', 'Newton-Raphson', '1D Sampling')
resultsGLMMListFlattened <- flattenResultList(
  resultsGLMMList, paramNames=glmmParamNodes, 
  methodNames=methodNames, numMCMCSamples=numMCMCSamplesGrid
)
resultsGLMMTimesFlattened <- flattenTimesList(
  timesGLMMList, 
  methodNames=methodNames,
  numMCMCSamples=numMCMCSamplesGrid)

resultPlot <- ggplot(
  resultsGLMMListFlattened %>% 
    mutate(
      MCMCSampleSize=factor(MCMCSampleSize, levels=numMCMCSamplesGrid)
    ), 
  aes(x=iteration, y=estimate, col=MCMCSampleSize, lty=MCMCSampleSize)) +
  facet_grid(parameter~method, scales="free_y") +
  geom_line() +
  ggtitle(paste0('N = ', N)) +
  theme(legend.position = "top") +
  theme(
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    legend.text=element_text(size=13)
  ) +
  labs(color='MCMC sample size', lty='MCMC sample size') 

timesPlot <- ggplot(
  resultsGLMMTimesFlattened,
  aes(x=MCMCSampleSize, y=exec.time, col=method, lty=method)) +
  geom_line(size=1.2) +
  ggtitle(paste0('N = ', N)) +
  xlab('MCMC sample size') +
  ylab('execution time (seconds)')
timesPlot

write.csv(resultsGLMMListFlattened, paste0("salamander_results/results_n", N, ".csv"), row.names = F)
write.csv(resultsGLMMTimesFlattened, paste0("salamander_results/times_n", N, ".csv"), row.names = F)

pdf(paste0("glmm_N_", N, ".pdf"), width = 12)
resultPlot
dev.off()

pdf(paste0("glmm_N_", N, "_times.pdf"), width = 7)
timesPlot
dev.off()


# Combining the exec time results ----
timesCombined <- rbind(
  cbind(read.csv("salamander_results/times_n180.csv"), ScientificSampleSize=180),
  cbind(read.csv("salamander_results/times_n360.csv"), ScientificSampleSize=360),
  cbind(read.csv("salamander_results/times_n720.csv"), ScientificSampleSize=720),
  cbind(read.csv("salamander_results/times_n1440.csv"), ScientificSampleSize=1440)
)

facetLabelAppender <- function(string, prefix = "N_mcmc = ") paste0(prefix, string)

plotTime1 <- ggplot(timesCombined, 
       aes(x=ScientificSampleSize, y=exec.time, col=method, lty=method)
       ) + 
  geom_line(size=1.1) +
  xlab("MCMC sample size (N_mcmc)") +
  ylab('execution time (seconds)') +
  facet_grid(~MCMCSampleSize, labeller=as_labeller(facetLabelAppender))

plotTime2 <- ggplot(timesCombined, 
       aes(x=MCMCSampleSize, y=exec.time, col=method, lty=method)
) + 
  geom_line(size=1.1) +
  xlab("scientific sample size (N_obs)") +
  ylab('execution time (seconds)') +
  facet_grid(~ScientificSampleSize, 
             labeller=as_labeller(function(x) facetLabelAppender(x, 'N_obs = ')))

pdf(paste0("glmm_combined_times.pdf"), width = 7)
grid.arrange(
  plotTime1, plotTime2, nrow=2,
  top="Execution time for varying MCMC sample sizes (N_mcmc) and scientific sample sizes (N_obs)"
)
dev.off()


# Combining the trajectory results ----
trajectoryCombined <- rbind(
  cbind(read.csv("salamander_results/results_n180.csv"), ScientificSampleSize=180),
  cbind(read.csv("salamander_results/results_n360.csv"), ScientificSampleSize=360),
  cbind(read.csv("salamander_results/results_n720.csv"), ScientificSampleSize=720),
  cbind(read.csv("salamander_results/results_n1440.csv"), ScientificSampleSize=1440)
)

trajectoryPlot <- function(N) {
  ggplot(
    trajectoryCombined %>% 
      filter(ScientificSampleSize == N) %>% 
      mutate(
        MCMCSampleSize=factor(MCMCSampleSize, levels=numMCMCSamplesGrid),
        method=factor(method, levels=methodNames)
      ), 
    aes(x=iteration, y=estimate, col=MCMCSampleSize, lty=MCMCSampleSize)) +
    facet_grid(parameter~method, scales="free_y") +
    geom_line() +
    ggtitle(paste0('N_obs = ', N)) +
    theme(legend.position = "top") +
    theme(
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 20),
      legend.title=element_text(size=20),
      legend.text=element_text(size=20),
      plot.title = element_text(size = 30)
    ) +
    labs(color='MCMC sample size', lty='MCMC sample size') 
}

t180 <- trajectoryPlot(N = 180)
t360 <- trajectoryPlot(N = 360)
t720 <- trajectoryPlot(N = 720)
t1440 <- trajectoryPlot(N = 1440)

pdf("glmm_combined_trajectories.pdf", width=20, height=20)
grid.arrange(t180, t360, t720, t1440, nrow=2, 
  top=textGrob("Estimate trajectories", gp=gpar(fontsize=35,font=8)))
dev.off()


