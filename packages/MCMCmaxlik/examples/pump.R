# Pump example (2 parameters)

# Load the packages ---------------------------------------
library(nimble)
nimbleOptions(experimentalEnableDerivs = TRUE)
library("MCMCmaxlik")

# Model specification -------------------------------------

pumpCode <- nimbleCode({
  for (i in 1:N) {
    theta[i] ~ dgamma(alpha, beta)
    lambda[i] <- theta[i] * t[i]
    x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1, 1.0)
})
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                         31.4, 1.05, 1.05, 2.1, 10.5))
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))

# Build the model.
pump <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

paramNodesPump <- pump$getNodeNames(topOnly = T)
compileNimble(pump)

# Contour plot prep work ---------------

fr2 <- function(alpha, beta) {
  x <- pumpData$x
  t <- pumpConsts$t
  first_term <- (-alpha - beta) + (0.1 - 1) * log(beta)
  second_term <- 0
  for (i in 1:length(x)) {
    second_term <- second_term + lgamma(alpha + x[i]) - 
      (alpha + x[i]) * log(beta + t[i])
  }
  third_term <- length(x) * (alpha * log(beta) - lgamma(alpha))
  #-(first_term + second_term + third_term)
  -(second_term + third_term)
}

fr2withprior <- function(alpha, beta) {
  x <- pumpData$x
  t <- pumpConsts$t
  first_term <- (-alpha - beta) + (0.1 - 1) * log(beta)
  second_term <- 0
  for (i in 1:length(x)) {
    second_term <- second_term + lgamma(alpha + x[i]) - 
      (alpha + x[i]) * log(beta + t[i])
  }
  third_term <- length(x) * (alpha * log(beta) - lgamma(alpha))
  -(first_term + second_term + third_term)
}

# Generate matrix for contour plot.
alphaSeq <- seq(.01, 10, length=1000)
betaSeq <- seq(.01, 16, length=1000)
mat1 <- outer(alphaSeq, betaSeq, fr2)
contourdf <- data.frame(cbind(expand.grid(alphaSeq, betaSeq), as.vector(mat1)))
names(contourdf) <- c("alpha", "beta", "density")

# Testing Algorithms ------------------------------------------------------

# Compile the necessary functions.
compiledFunsPump <- buildMCMCmaxlik(pump, paramNodesPump)

init <- c(10, 10)
boundary <- list(c(0.05, 10000), c(0.05, 10000))
numMCMCSamples <- 300

getPumpResults <- function(init, 
                           boundary=list(c(0.05, 10000), c(0.05, 10000)), 
                           numMCMCSamples) {
  resultsPumpFixed <- computeMLE(pump, paramNodesPump,
                                 method="fixed", paramInit=init,
                                 compiledFuns=compiledFunsPump,
                                 stepsize=0.05,
                                 numMCMCSamples=numMCMCSamples,
                                 maxIter=300,
                                 boundary=boundary)
  
  resultsPumpSmallFixed <- computeMLE(pump, paramNodesPump,
                                      method="fixed", paramInit=init,
                                      compiledFuns=compiledFunsPump,
                                      stepsize=0.005,
                                      numMCMCSamples=numMCMCSamples,
                                      maxIter=300,
                                      boundary=boundary)
  
  resultsPumpNR <- computeMLE(pump, paramNodes=paramNodesPump,
                              method="NR", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              tol=1e-20,
                              maxIter=300,
                              boundary=boundary)
  
  resultsPumpAdam <- computeMLE(pump, paramNodesPump,
                                method="adam", paramInit=init,
                                compiledFuns=compiledFunsPump,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300)
  
  resultsPump1D <- computeMLE(pump, paramNodesPump,
                              method="ga1D", paramInit=init,
                              compiledFuns=compiledFunsPump,
                              numMCMCSamples=numMCMCSamples,
                              numMCMCSamples1D=numMCMCSamples,
                              maxIter=300)
  return(list(resultsPumpFixed, resultsPumpSmallFixed,
              resultsPumpAdam, resultsPumpNR, resultsPump1D))
}

getPumpSummary <- function(resultsPump) {
  times <- sapply(resultsPump,
                  function(results) {
                    exec.time <- results$execution.time[3]
                    conv.time <- results$convergence.time[3]
                    conv.iter <- results$convergence.iter
                    if(is.null(conv.time)) conv.time <- NA
                    if(is.null(conv.iter)) conv.iter <- NA
                    c(round(results$MLE, 3), 
                      round(exec.time, 3), 
                      round(conv.time, 3), 
                      conv.iter,
                      round(fr2(results$MLE[1], results$MLE[2]) - fr2(0.823, 1.262), 5))
                  })
  rownames(times) <- c('alpha', 'beta', 'exec.time', 'conv.time', 'conv.iter', 'loglik_diff')
  colnames(times) <- c('fixed', 'small_fixed', 'Adam', 'NR',  '1D')
  t(times)
}

set.seed(1000)
resultsPump_10_10_300 <- getPumpResults(init=c(10, 10), numMCMCSamples=300)
timesPump_10_10_300 <- getPumpSummary(resultsPump_10_10_300)
timesPump_10_10_300

set.seed(2000)
resultsPump_10_2_300 <- getPumpResults(init=c(10, 2), numMCMCSamples=300)
timesPump_10_2_300 <- getPumpSummary(resultsPump_10_2_300)
timesPump_10_2_300

set.seed(3000)
resultsPump_10_10_3000 <- getPumpResults(init=c(10, 10), numMCMCSamples=3000)
timesPump_10_10_3000 <- getPumpSummary(resultsPump_10_10_3000)
timesPump_10_10_3000

set.seed(4000)
resultsPump_10_2_3000 <- getPumpResults(init=c(10, 2), numMCMCSamples=3000)
timesPump_10_2_3000 <- getPumpSummary(resultsPump_10_2_3000)
timesPump_10_2_3000

timesPumpAll <- 
  rbind(cbind(timesPump_10_10_300, alpha_init=10, beta_init=10, numMCMCSamples=300),
        cbind(timesPump_10_2_300, alpha_init=10, beta_init=2, numMCMCSamples=300),
        cbind(timesPump_10_10_300, alpha_init=10, beta_init=10, numMCMCSamples=300),
        cbind(timesPump_10_2_3000, alpha_init=10, beta_init=2, numMCMCSamples=3000))

write.csv(timesPumpAll, 'times_pump.csv')

# MCEM time: 0.425
resultsPumpMCEM_MLE <- c(0.8190847, 1.2537176)
round(fr2(resultsPumpMCEM_MLE[1], resultsPumpMCEM_MLE[2]) - fr2(0.823, 1.262), 5)


# pumpMCEM <- buildMCEM_AD(model = pump, latentNodes = 'theta',
#                          boxConstraints = list( list(c('alpha', 'beta'),  
#                                                      limits = c(0, Inf) ) ))
# MCEM_time <- proc.time()
# out <- pumpMCEM$run()
# MCEM_time <- proc.time() - MCEM_time

# pumpMCEM <- buildMCEM(model = pump, latentNodes = 'theta',
#                       boxConstraints = list( list( c('alpha', 'beta'),
#                                                    limits = c(0, Inf) ) ))
# MCEM_time <- proc.time()
# out <- pumpMCEM$run()
# MCEM_time <- proc.time() - MCEM_time

# Results --------------------------------
iteratesFixed <- data.frame(resultsPumpFixed$param,
                            method=rep("Fixed (0.05)",
                                       nrow(resultsPumpFixed$param)))
iteratesSmallFixed <- data.frame(resultsPumpSmallFixed$param,
                                 method=rep("Fixed (0.005)",
                                            nrow(resultsPumpSmallFixed$param)))
iteratesNR <- data.frame(resultsPumpNR$param,
                         method=rep("Newton-Raphson", nrow(resultsPumpNR$param)))
iteratesAdam <- data.frame(resultsPumpAdam$param, 
                           method=rep("Adam", nrow(resultsPumpAdam$param)))
iterates1D <- data.frame(resultsPump1D$param,
                         method=rep("1D sampling", nrow(resultsPump1D$param)))
iterates <- rbind(iteratesFixed, 
                  iteratesSmallFixed, 
                  iteratesNR,
                  iteratesAdam,
                  iterates1D)
names(iterates) <- c("alpha", "beta", "method")


library(RColorBrewer)
cbPalette <- brewer.pal(8, "Dark2")

pdf(paste0("pump_all_", numMCMCSamples, "_",
           paste0(init, collapse="_"),
           ".pdf"), width=7.5)
par(mar=c(5, 5, 1, 1) + 0.1, xpd=T)
contour(alphaSeq, betaSeq, mat1, nlevels = 30, 
        xlab=expression(alpha),
        ylab=expression(beta),
        cex.lab=2,
        cex.axis=2)
for (method in 1:length(levels(iterates$method))) {
  with(iterates[iterates$method==levels(iterates$method)[method], ],
       lines(alpha, beta, col=cbPalette[method], lty=(1:5)[method],
             lwd=3))
  # with(iterates[iterates$method==levels(iterates$method)[method], ],
  #      points(alpha, beta, col=cbPalette[method]))
}
dev.off()

pdf(paste0("pump_all_", numMCMCSamples, "_",
           paste0(init, collapse="_"), '_zoom_in',
           ".pdf"), width=7.5)
par(mar=c(5, 5, 1, 1) + 0.1, xpd=F)
contour(alphaSeq, betaSeq, mat1, nlevels = 30, 
        xlab=expression(alpha),
        ylab=expression(beta),
        cex.lab=2,
        cex.axis=2,
        xlim=c(0, 2),
        ylim=c(0, 2))
for (method in 1:length(levels(iterates$method))) {
  with(iterates[iterates$method==levels(iterates$method)[method], ],
       lines(alpha, beta, col=cbPalette[method], lty=(1:6)[method],
             lwd=3))
  # with(iterates[iterates$method==levels(iterates$method)[method], ],
  #      points(alpha, beta, col=cbPalette[method]))
}
dev.off()

pdf("pump_legend.pdf", width = 12, height=3.5)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, legend=levels(iterates$method),
       lty=1:6, lwd=3, col=cbPalette[1:length(levels(iterates$method))], cex = 1.4,
       xjust=0.5, yjust=0.5, ncol=3)
dev.off()

#title(paste0("pump log marginal likelihood\nMCMC samples: ", numMCMCSamples),
#     cex.main=2)
