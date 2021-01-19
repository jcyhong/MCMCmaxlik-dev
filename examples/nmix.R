library("MCMCmaxlik")
library(unmarked)
require(nimble)
library(parallel)
library(data.table)
install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
library(INLA)
nimbleOptions(experimentalEnableDerivs = TRUE)

set.seed(12345)
sim.nmix <- function(n.sites = 72, # number of study sites
                     n.surveys = 3, # short term replicates
                     n.years = 9, # number of years
                     b0 = 2.0, # intercept log(lambda)
                     b1 = 2.0, # x1 slope log(lambda)
                     b2 = -3.0, # x2 slope log(lambda)
                     b3 = 1.0, # x3 slope log(lambda)
                     a0 = 1.0, # intercept logit(p)
                     a1 = -2.0, # x1 slope logit(p)
                     a4 = 1.0, # x4 slope logit(p)
                     th = 3.0 # overdisperison parameter
                     ){
  # make empty N and Y arrays
  if(n.years %% 2 == 0) {n.years <- n.years + 1}
  N.tr <- array(dim = c(n.sites, n.years))
  Y <- array(dim = c(n.sites, n.surveys, n.years))
  Y.m <- array(dim = c(n.sites, n.surveys, n.years))

  # create abundance covariate values
  x1 <- array(as.numeric(scale(runif(n = n.sites, -0.5, 0.5), scale = F)),
              dim = c(n.sites, n.years))
  x2 <- array(as.numeric(scale(runif(n = n.sites, -0.5, 0.5), scale = F)),
              dim = c(n.sites, n.years))
  yrs <- 1:n.years; yrs <- (yrs - mean(yrs)) / (max(yrs - mean(yrs))) / 2
  x3 <- array(rep(yrs, each = n.sites), dim = c(n.sites, n.years))
  
  # fill true N array
  lam.tr <- exp(b0 + b1 * x1 + b2 * x2 + b3 * x3)
  for(i in 1:n.sites){
    for(k in 1:n.years){
      N.tr[i, k] <- rnbinom(n = 1, mu = lam.tr[i, k], size = th)
    }}
  
  # create detection covariate values
  x1.p <- array(x1[,1], dim = c(n.sites, n.surveys, n.years))
  x4 <- array(as.numeric(scale(runif(n = n.sites * n.surveys * n.years,
                                     -0.5, 0.5), scale = F)), dim = c(n.sites, n.surveys, n.years))
  
  # average x4 per site-year for example 1
  x4.m <- apply(x4, c(1, 3), mean, na.rm = F)
  out1 <- c()
  for(k in 1:n.years){
    chunk1 <- x4.m[ , k]
    chunk2 <- rep(chunk1, n.surveys)
    out1 <- c(out1, chunk2)
    }
  x4.m.arr <- array(out1, dim = c(n.sites, n.surveys, n.years))
  
  # fill Y.m count array using x4.m for example 1
  p.tr1 <- plogis(a0 + a1 * x1.p + a4 * x4.m.arr)
  for (i in 1:n.sites){
    for (k in 1:n.years){
      for (j in 1:n.surveys){
        Y.m[i, j, k] <- rbinom(1, size = N.tr[i, k], prob = p.tr1[i, j, k])
      }}}
  
  # fill Y count array using x4 for example 2
  p.tr2 <- plogis(a0 + a1 * x1.p + a4 * x4)

   for (i in 1:n.sites){
     for (k in 1:n.years){
       for (j in 1:n.surveys){
         Y[i, j, k] <- rbinom(1, size = N.tr[i, k], prob = p.tr2[i, j, k])
       }}}
  
  # format Y.m for data frame output for inla and unmarked
  Y.m.df <- Y.m[ , , 1]
  for(i in 2:n.years){
    y.chunk <- Y.m[ , , i]
    Y.m.df <- rbind(Y.m.df, y.chunk)
  }
  
  # format covariates for data frame output for inla and unmarked
  x1.df <- rep(x1[ , 1], n.years)
  x2.df <- rep(x2[ , 1], n.years)
  x3.df <- rep(x3[1, ], each = n.sites)
  x1.p.df <- rep(x1.p[ , 1, 1], n.years)
  x4.df <- c(x4.m)
  
  # put together data frames for inla and unmarked
  inla.df <- unmk.df <- data.frame(y1 = Y.m.df[ , 1], y2 = Y.m.df[ , 2],
                                   y3 = Y.m.df[ , 3], x1 = x1.df, x2 = x2.df, x3 = x3.df,
                                   x1.p = x1.p.df, x4.m = x4.df)
  
  # return all necessary data for examples 1 and 2
  return(list(inla.df = inla.df, unmk.df = unmk.df, n.sites = n.sites,
              n.surveys = n.surveys, n.years = n.years, x1 = x1[ , 1],
              x2 = x2[ , 1], x3 = x3[1, ], x4 = x4, x4.m = x4.m, x4.m.arr = x4.m.arr,
              Y = Y, Y.m = Y.m, lam.tr = lam.tr, N.tr = N.tr, x1.p = x1.p[ , 1, 1]
              ))
  } # end sim.nmix function
sim.data <- sim.nmix()
  
nmixCode <- nimbleCode({
  a0 ~ dnorm(0, sd=100) 
  a1 ~ dnorm(0, sd=100)
  a4 ~ dnorm(0, sd=100)
  b0 ~ dnorm(0, sd=100)
  b1 ~ dnorm(0, sd=100)
  b2 ~ dnorm(0, sd=100)
  b3 ~ dnorm(0, sd=100)
  th ~ dunif(0, 5)
  for (k in 1:n.years){
    for (i in 1:n.sites){
      N[i, k] ~ dnegbin(prob=prob[i, k], size=th)
      prob[i, k] <- th / (th + lambda[i, k])
      log(lambda[i, k]) <- b0 + (b1 * x1[i]) + (b2 * x2[i]) + (b3 * x3[k])
      for (j in 1:n.surveys){
        Y.m[i + (j - 1) * n.sites, k] ~ dbin(p[i + (j - 1) * n.sites, k], N[i,k])
        p[i + (j - 1) * n.sites, k] <- exp(lp[i + (j - 1) * n.sites, k]) / (1 + exp(lp[i + (j - 1) * n.sites, k]))
        lp[i + (j - 1) * n.sites, k] <- a0 + (a1 * x1.p[i]) + (a4 * x4.m[i, k])
        }}}
})


nmixConstants <- list(n.years=sim.data$n.years,
                      n.sites=sim.data$n.sites,
                      n.surveys=sim.data$n.surveys)
       
nmixData <- list(
  Y.m=rbind(sim.data$Y.m[, 1, ], sim.data$Y.m[, 2, ], sim.data$Y.m[, 3, ]),
  #Y.m=sim.data$Y.m,
  x1.p=sim.data$x1.p,
  x1=sim.data$x1,
  x2=sim.data$x2,
  x3=sim.data$x3,
  x4.m=sim.data$x4.m
)

nmixInits <- list(a0=1, a1=-2, a4=1, b0=2, b1=2, b2=-3, b3=1, th=3)

nmixMod <- nimbleModel(code=nmixCode, constants=nmixConstants,
                       data=nmixData, inits=nmixInits, check = FALSE)

nmixParamNodes <- nmixMod$getNodeNames(topOnly = T,stochOnly=T)

Cnmix <- compileNimble(nmixMod)

compiledFunsnmix <- buildMCMCmaxlik(nmixMod, nmixParamNodes)


init <- c(runif(1, min=-5, max=5), 
          runif(1, min=-5, max=5),
          runif(1, min=-5, max=5), 
          runif(1, min=-5, max=5),
          runif(1, min=-5, max=5),
          runif(1, min=-5, max=5),
          runif(1, min=-5, max=5),
          runif(1, min=0.05, max=5))
boundary <- list(c(-40, 40), c(-40, 40), c(-40, 40),
                 c(-40, 40), c(-40, 40), c(-40, 40),
                 c(-40, 40), c(0.01, 100))
numMCMCSamples <- 300
maxIter <- 300


resultsAdam <- computeMLE(nmixMod, nmixParamNodes,
                          method="adam", paramInit=init,
                          compiledFuns=compiledFunsnmix,
                          numMCMCSamples=numMCMCSamples,
                          maxIter=maxIter,
                          boundary=boundary)

nmixMCMC <- buildMCMC(nmixMod)

nmixMCMC$run(1)



nmixMod$simulate("N")
nmixMod$N
