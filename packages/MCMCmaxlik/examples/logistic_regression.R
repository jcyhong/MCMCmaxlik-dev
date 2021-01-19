# Logistic regression example (3 parameters)

# Load the packages ---------------------------------------
library("MCMCmaxlik")
library(nimble)
nimbleOptions(experimentalEnableDerivs = TRUE)

# Model specification -------------------------------------
code <- nimbleCode({
  beta0 ~ dnorm(0, sd = 10000)
  beta1 ~ dnorm(0, sd = 10000)
  sigma_RE ~ dunif(0, 1000)
  for(i in 1:N) {
    beta2[i] ~ dnorm(0, sd = sigma_RE)
    logit(p[i]) <- beta0 + beta1 * x[i] + beta2[i]
    r[i] ~ dbin(p[i], n[i])
  }
})

## constants, data, and initial values
constants <- list(N = 10)

data <- list(
  r = c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46),
  n = c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79),
  x = c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
)

inits <- list(beta0 = 0, beta1 = 0, sigma_RE = 1)

logreg <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check = FALSE)

paramNodesLogreg <- logreg$getNodeNames(topOnly = T)
Clogreg <- compileNimble(logreg)


# Testing Algorithms ------------------------------------------------------
# Compile the necessary functions.
compiledFunsLogreg <- buildMCMCmaxlik(logreg, paramNodesLogreg)

boundary <- list(c(-40, 40), c(-40, 40), c(0.05, 10))

getLogregResults <- function(init, 
                             boundary=list(c(-40, 40), c(-40, 40), c(0.05, 10)), 
                             numMCMCSamples) {
  
  resultsLogregFixed <- computeMLE(logreg, paramNodesLogreg,
                                   method="fixed", paramInit=init,
                                   stepsize=0.05,
                                   compiledFuns=compiledFunsLogreg,
                                   numMCMCSamples=numMCMCSamples,
                                   maxIter=300,
                                   boundary=boundary)
  resultsLogregSmallFixed <- computeMLE(logreg, paramNodesLogreg,
                                        method="fixed", paramInit=init,
                                        stepsize=0.005,
                                        compiledFuns=compiledFunsLogreg,
                                        numMCMCSamples=numMCMCSamples,
                                        maxIter=300,
                                        boundary=boundary)
  resultsLogregNR <- computeMLE(logreg, paramNodesLogreg,
                                method="NR", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300,
                                boundary=boundary)
  
  resultsLogregAdam <- computeMLE(logreg, paramNodesLogreg,
                                  method="adam", paramInit=init,
                                  compiledFuns=compiledFunsLogreg,
                                  numMCMCSamples=numMCMCSamples,
                                  maxIter=300,
                                  boundary=boundary)
  
  resultsLogreg1D <- computeMLE(logreg, paramNodesLogreg,
                                method="ga1D", paramInit=init,
                                compiledFuns=compiledFunsLogreg,
                                numMCMCSamples=numMCMCSamples,
                                maxIter=300)
  list(resultsLogregFixed, resultsLogregSmallFixed,
       resultsLogregAdam, resultsLogregNR, resultsLogreg1D)
}

require(lme4)
r=data$r
n=data$n
x=data$x


xlong<-c()
rlong<-c()
indiv<-c()

for(i in 1:length(r)){
  xlong=c(xlong,rep(x[i],n[i]))
  rlong<-c(rlong,rep(1,r[i]),rep(0,n[i]-r[i]))
  indiv<-c(indiv,rep(i,n[i]))
}

data=cbind.data.frame(xlong,rlong,indiv)
gmre=glmer(rlong~xlong+(1|indiv),family=binomial(),data=data)
gmreDevfun <- update(gmre, devFunOnly=TRUE)
# sigmaRE, beta0, beta1
# log-likelihood difference
resultsLogregGlmer <- c(coef(summary(gmre))[, 'Estimate'], 
                        attr(summary(gmre)$varcor$indiv, 'stddev'))
names(resultsLogregGlmer) <- paramNodesLogreg
loglik_diff_logreg <- sapply(1:nrow(logregMLEs), function(i) {
  (gmreDevfun(logregMLEs[nrow(logregMLEs), c('sigma_RE', 'beta0', 'beta1')]) - 
     gmreDevfun(logregMLEs[i, c('sigma_RE', 'beta0', 'beta1')])) / -2
})


getLogregSummary <- function(resultsLogreg) {
  timesLogreg <- sapply(resultsLogreg,
                        function(results) {
                          exec.time <- results$execution.time[1]
                          conv.time <- results$convergence.time[1]
                          conv.iter <- results$convergence.iter
                          if(is.null(conv.time)) conv.time <- NA
                          if(is.null(conv.iter)) conv.iter <- NA
                          names(results$MLE) <- paramNodesLogreg
                          c(round(results$MLE, 3), 
                            round(exec.time, 3), 
                            round(conv.time, 3), 
                            conv.iter,
                            round((gmreDevfun(resultsLogregGlmer[c('sigma_RE', 'beta0', 'beta1')]) - 
                               gmreDevfun(results$MLE[c('sigma_RE', 'beta0', 'beta1')])) / -2, 5))
                        })
  rownames(timesLogreg) <- c('beta_0', 'beta_1', 'sigma', 'exec.time', 'conv.time', 'conv.iter',
                             'loglik_diff')
  colnames(timesLogreg) <- c('fixed', 'small_fixed', 'Adam', 'NR', '1D')
  t(timesLogreg)
}

set.seed(1000)
resultsLogreg_0_0_1_20 <- getLogregResults(init=c(0, 0, 1), numMCMCSamples=20)
timesLogreg_0_0_1_20 <- getLogregSummary(resultsLogreg_0_0_1_20)
set.seed(1800)
resultsLogreg_1_1_4_20 <- getLogregResults(init=c(-1, -1, 4), numMCMCSamples=20)
timesLogreg_1_1_4_20 <- getLogregSummary(resultsLogreg_1_1_4_20)
set.seed(2000)
resultsLogreg_0_0_1_300 <- getLogregResults(init=c(0, 0, 1), numMCMCSamples=300)
timesLogreg_0_0_1_300 <- getLogregSummary(resultsLogreg_0_0_1_300)
set.seed(3000)
resultsLogreg_1_1_4_300 <- getLogregResults(init=c(-1, -1, 4), numMCMCSamples=300)
timesLogreg_1_1_4_300 <- getLogregSummary(resultsLogreg_1_1_4_300)

write.csv(t(timesLogreg), 'times_glmm1.csv')

resultsLogregMCEM_MLE <- c(-0.5444048, 1.3102369, 0.2514697)
round((gmreDevfun(resultsLogregGlmer[c('sigma_RE', 'beta0', 'beta1')]) - 
         gmreDevfun(resultsLogregMCEM_MLE[c(3, 1, 2)])) / -2, 5)


# LogregMCEM <- buildMCEM(model=logreg, 
#                         latentNodes = 'beta2',
#                         boxConstraints = list( 
#                           list(c('sigma_RE'), 
#                                limits = c(0, 1000) ) ))
# LogregMCEM_time <- proc.time()
# out <- LogregMCEM$run()
# LogregMCEM_time <- proc.time() - LogregMCEM_time
# 
# 
# # lme4
# r=c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46)
# n=c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79)
# x=c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
# 
# xlong<-c()
# rlong<-c()
# indiv<-c()
# 
# for(i in 1:length(r)){
#   xlong=c(xlong,rep(x[i],n[i]))
#   rlong<-c(rlong,rep(1,r[i]),rep(0,n[i]-r[i]))
#   indiv<-c(indiv,rep(i,n[i]))
# }
# 
# data=cbind.data.frame(xlong,rlong,indiv)
# 
# require(lme4)
# ptm <- proc.time()
# gmre=glmer(rlong~xlong+(1|indiv),family=binomial(),data=data)
# glmmTime <- proc.time() - ptm ##  0.093 
# 
# summary(gmre)


# Results -----------
iteratesLogregFixed <- data.frame(
  iter=1:nrow(resultsLogregFixed$param),
  resultsLogregFixed$param,
  method=rep("Fixed (0.05)", nrow(resultsLogregFixed$param))
)
iteratesLogregSmallFixed <- data.frame(
  iter=1:nrow(resultsLogregSmallFixed$param),
  resultsLogregSmallFixed$param,
  method=rep("Fixed (0.005)", nrow(resultsLogregSmallFixed$param))
)
iteratesLogregNR <- data.frame(
  iter=1:nrow(resultsLogregNR$param),
  resultsLogregNR$param,
  method=rep("Newton-Raphson", nrow(resultsLogregNR$param))
)
iteratesLogreg1D <- data.frame(
  iter=1:nrow(resultsLogreg1D$param),
  resultsLogreg1D$param,
  method=rep("1D sampling", nrow(resultsLogreg1D$param))
)
iteratesLogregAdam <- data.frame(
  iter=1:nrow(resultsLogregAdam$param),
  resultsLogregAdam$param, 
  method=rep("Adam", nrow(resultsLogregAdam$param))
)

if (numMCMCSamples > 100) {
  iteratesLogreg <- rbind(iteratesLogregFixed, iteratesLogregSmallFixed,
                          iteratesLogregNR,
                          iteratesLogregAdam,
                          iteratesLogreg1D) 
} else {
  iteratesLogreg <- rbind(iteratesLogregFixed, iteratesLogregSmallFixed,
                          iteratesLogregNR, iteratesLogregAdam,
                          iteratesLogreg1D) 
}
names(iteratesLogreg) <- c("iter", paramNodesLogreg, "method")

library(ggplot2)
library(reshape2)
iteratesLogregMelted <- melt(iteratesLogreg, id.vars=c("iter", "method"))
levels(iteratesLogregMelted$variable) <- c(expression(beta[0]),
                                           expression(beta[1]),
                                           expression(sigma["RE"]))

library(RColorBrewer)
cbPalette <- brewer.pal(8, "Dark2")
# trajectoryPlotLogreg <- ggplot(iteratesLogregMelted, aes(x=iter, y=value, colour=method)) 
# trajectoryPlotLogreg <- trajectoryPlotLogreg + geom_line(alpha=0.8)
# trajectoryPlotLogreg <- trajectoryPlotLogreg + facet_grid(variable ~ .,
#                                                           labeller=label_parsed,
#                                                           scales = "free_y")
# trajectoryPlotLogreg <- trajectoryPlotLogreg 
# trajectoryPlotLogreg <- trajectoryPlotLogreg + 
#   scale_color_manual(values=cbPalette)
# trajectoryPlotLogreg <- trajectoryPlotLogreg + ggtitle(
#   paste0("logistic regression\n start: (",
#          paste(init, collapse=", "),
#          "); MCMC samples: ",
#          numMCMCSamples)
# )
# trajectoryPlotLogreg <- trajectoryPlotLogreg + xlab("number of iterations")
# trajectoryPlotLogreg <- trajectoryPlotLogreg + 
#   theme(plot.title = element_text(face = "bold", size = 20),
#         axis.title = element_text(face = "bold", size = 20),
#         axis.text = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         legend.title = element_text(face = "bold", size = 20))
# pdf(paste0("logistic_regression_", numMCMCSamples, ".pdf"))
# trajectoryPlotLogreg

cbPalette <- brewer.pal(8, "Dark2")
trajectoryPlotLogreg <- ggplot(iteratesLogregMelted, aes(x=iter, y=value, colour=method, linetype=method)) 
trajectoryPlotLogreg <- trajectoryPlotLogreg + geom_line(alpha=0.8, lwd=0.8)
trajectoryPlotLogreg <- trajectoryPlotLogreg + facet_grid(variable ~ .,
                                                          labeller=label_parsed,
                                                          scales = "free_y")
trajectoryPlotLogreg <- trajectoryPlotLogreg + theme(legend.position="none")
trajectoryPlotLogreg <- trajectoryPlotLogreg + 
  scale_color_manual(values=cbPalette)
trajectoryPlotLogreg <- trajectoryPlotLogreg + xlab("number of iterations")
trajectoryPlotLogreg <- trajectoryPlotLogreg + 
  theme(axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(size = 15),
        strip.text.y = element_text(size = 15))
pdf(paste0("logistic_regression_", numMCMCSamples, ".pdf"))
trajectoryPlotLogreg
dev.off()


trajectoryPlotLogreg_zoom <- ggplot(iteratesLogregMelted[iteratesLogregMelted$method != 'Fixed (0.05)' &
                                                           iteratesLogregMelted$iter >= 200, ], 
                                    aes(x=iter, y=value, colour=method, linetype=method))
trajectoryPlotLogreg_zoom <- trajectoryPlotLogreg_zoom + geom_line(alpha=0.8, lwd=0.8)
trajectoryPlotLogreg_zoom <- trajectoryPlotLogreg_zoom + facet_grid(variable ~ .,
                                                          labeller=label_parsed,
                                                          scales = "free_y")
trajectoryPlotLogreg_zoom <- trajectoryPlotLogreg_zoom + theme(legend.position="none")
trajectoryPlotLogreg_zoom <- trajectoryPlotLogreg_zoom + 
  scale_color_manual(values=cbPalette[-1])
trajectoryPlotLogreg_zoom <- trajectoryPlotLogreg_zoom + xlab("number of iterations")
trajectoryPlotLogreg_zoom <- trajectoryPlotLogreg_zoom + 
  theme(axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(size = 15),
        strip.text.y = element_text(size = 15))
trajectoryPlotLogreg_zoom
pdf(paste0("logistic_regression_", numMCMCSamples, "_zoom_in.pdf"))
trajectoryPlotLogreg_zoom
dev.off()


#-------------------------------

# Plot the legend only.

library(ggplot2); library(gridExtra); library(grid); library(dplyr)

cbPalette <- brewer.pal(8, "Dark2")
trajectoryPlotLogreg <- ggplot(iteratesLogregMelted, aes(x=iter, y=value, colour=method, linetype=method)) 
trajectoryPlotLogreg <- trajectoryPlotLogreg + geom_line(alpha=0.8, lwd=0.8)
trajectoryPlotLogreg <- trajectoryPlotLogreg + facet_grid(variable ~ .,
                                                          labeller=label_parsed,
                                                          scales = "free_y")
trajectoryPlotLogreg <- trajectoryPlotLogreg
trajectoryPlotLogreg <- trajectoryPlotLogreg + 
  scale_color_manual(values=cbPalette)
trajectoryPlotLogreg <- trajectoryPlotLogreg + xlab("number of iterations")
trajectoryPlotLogreg <- trajectoryPlotLogreg + 
  theme(axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(size = 15),
        strip.text.y = element_text(size = 15)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank(), legend.text = element_text(size = 15))

# Extract Legend 
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

pdf(paste0("legend_logistic_regression_", numMCMCSamples, ".pdf"),
    width=12)
legend <- g_legend(trajectoryPlotLogreg) 
grid.draw(legend) 
dev.off()
