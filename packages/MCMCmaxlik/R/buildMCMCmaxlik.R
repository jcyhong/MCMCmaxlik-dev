#' Construct compiled functions for finding MLE.
#'
#' This function builds a list of compiled functions for finding the MLE.
#' @param model a nimble model object
#' @param paramNodes a character vector indicating the top-level parameters
#' @keywords MLE
#' @export
#' @examples
#' buildMCMCmaxlik()

buildMCMCmaxlik <- function(model, paramNodes){
  sampler_RW_rotated <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      ## control list extraction
      logScale      <- if(!is.null(control$log))           control$log           else FALSE
      reflective    <- if(!is.null(control$reflective))    control$reflective    else FALSE
      adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
      adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200
      scale         <- if(!is.null(control$scale))         control$scale         else 1
      computeGrad        <- control$computeGrad
      decideIncludePrior <- control$decideIncludePrior
      D <- length(model$expandNodeNames(target))
      upperBound         <- 10^8
      lowerBound         <- -10^8
      ## node list generation
      calcNodes  <- model$getDependencies(target)
      calcNodesWithoutSelf <- model$getDependencies(target, self=FALSE)
      ## numeric value generation
      scaleOriginal <- scale
      timesRan      <- 0
      timesAccepted <- 0
      timesAdapted  <- 0
      
      optimalAR     <- 0.44
      gamma1        <- 0
      grad          <- numeric(2)
      ## checks
      if(logScale & reflective)      stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    },
    run = function() {
      currentValue <- values(model, target)
      propLogScale <- 0
      grad <<- computeGrad$grad / sqrt(sum(computeGrad$grad ^ 2))
      includePrior <- decideIncludePrior$includePrior
      propValue <- currentValue + rnorm(1, mean = 0,  sd = scale) * grad
      
      values(model, target) <<- propValue
      if (includePrior) logMHR <- calculateDiff(model, calcNodes) + propLogScale 
      else {
        logMHR <- calculateDiff(model, calcNodesWithoutSelf) + propLogScale
        for (i in 1:length(propValue)) {
          if (propValue[i] > upperBound & propValue[i] < lowerBound) logMHR <- -Inf 
        }
      }
      jump <- decide(logMHR)
      if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
      if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
      adaptiveProcedure = function(jump = logical()) {
        timesRan <<- timesRan + 1
        if(jump)     timesAccepted <<- timesAccepted + 1
        
        if(timesRan %% adaptInterval == 0) {
          acceptanceRate <- timesAccepted / timesRan
          timesAdapted <<- timesAdapted + 1
          gamma1 <<- 1/((timesAdapted + 3)^0.8)
          gamma2 <- 10 * gamma1
          adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
          scale <<- scale * adaptFactor
          timesRan <<- 0
          timesAccepted <<- 0
        }
      },
      reset = function() {
        scale <<- scaleOriginal
        timesRan      <<- 0
        timesAccepted <<- 0
        timesAdapted  <<- 0
        
        gamma1 <<- 0
      }
    ), where = getLoadingNamespace()
  )
  
  nfSetParams <- nimbleFunction(
    setup = function(model, paramNodes) {
      deps <- model$getDependencies(paramNodes, determOnly = TRUE)
    },
    run = function(P = double(1)) {
      values(model, paramNodes) <<- P
      model$calculate(deps)
    }
  )
  
  nfSetLatent <- nimbleFunction(
    setup = function(model, latentNodes) {
      deps <- model$getDependencies(latentNodes, determOnly = TRUE)
    },
    run = function(P = double(1)) {
      values(model, latentNodes) <<- P
      model$calculate(deps)
    }
  )
  
  setParams <- nfSetParams(model, paramNodes)
  latentNodes <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE)
  setLatent <- nfSetLatent(model, latentNodes)
  mcmcConf <- configureMCMC(model, nodes = latentNodes, monitors = latentNodes)
  MCMC <- buildMCMC(mcmcConf)
  LLR <- nfLLR(model, paramNodes, latentNodes, MCMC$mvSamples)
  
  nfApproxGrad <- nimbleFunction(
    setup = function(model, paramNodes, stateNodes, mvSamples, nfLLR) {
      paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      calcNodes <- model$getDependencies(stateNodes)
      D <- length(model$expandNodeNames(paramNodes))
      grad <- numeric(2)
      myLLR <- nfLLR(model, paramNodes, stateNodes, mvSamples)
    },
    run = function(delta = double(0), includePrior = logical(0),
                   burninFrac = double(0)) {
      P <- values(model, paramNodes)
      m <- getsize(mvSamples)
      lr <- numeric(D)
      for (j in 1:D) {
        inc <- numeric(D, value = 0)
        inc[j] <- delta
        lr[j] <- myLLR$run(P, inc, includePrior, burninFrac)
      }
      grad <<- (lr - 1) / delta
    }
  )
  
  nfApproxHess <- nimbleFunction(
    setup = function(model, paramNodes, stateNodes, mvSamples, nfLLR) {
      paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      calcNodes <- model$getDependencies(stateNodes)
      D <- length(model$expandNodeNames(paramNodes))
      myLLR2 <- nfLLR(model, paramNodes, stateNodes, mvSamples)
    },
    run = function(delta = double(0), includePrior = logical(0), burninFrac = double(0)) {
      P <- values(model, paramNodes)
      m <- getsize(mvSamples)
      lr1 <- numeric(D)
      for (j in 1:D) {
        inc <- numeric(D, value = 0)
        inc[j] <- delta
        lr1[j] <- myLLR2$run(P, inc, includePrior, burninFrac)
      }
      grad <- (lr1 - 1) / delta
      
      hess <- matrix(0, nrow = D, ncol = D)
      for (i in 1:D) {
        for (j in i:D) {
          inc <- numeric(D, value = 0)
          inc[i] <- delta
          inc[j] <- inc[j] + delta 
          lr <- myLLR2$run(P, inc, includePrior, burninFrac)
          gradij <- (lr - 1) / delta
          hess[i, j] <- - grad[i] * grad[j] +
            1 / delta * (gradij - grad[i] - grad[j])
          if (i != j) hess[j, i] <- hess[i, j]
        }
      }
      returnType(double(2))
      return(hess)
    }
  )
  
  nfIncludePrior <- nimbleFunction(
    setup = function() {
      includePrior <- TRUE
    },
    run = function(withPrior = logical(0)) {
      includePrior <<- withPrior
    }
  )
  
  computeGrad <- nfApproxGrad(model, paramNodes, latentNodes, MCMC$mvSamples, nfLLR)
  computeHess <- nfApproxHess(model, paramNodes, latentNodes, MCMC$mvSamples, nfLLR)
  decideIncludePrior <- nfIncludePrior()
  mcmc1DConf <- configureMCMC(model)
  mcmc1DConf$removeSamplers(paramNodes) 
  mcmc1DConf$addSampler(target=paramNodes, type=sampler_RW_rotated,
                        control=list(computeGrad=computeGrad,
                                     decideIncludePrior=decideIncludePrior))
  MCMC1D <- buildMCMC(mcmc1DConf)
  compiledFuns <- compileNimble(setParams, LLR, setLatent,
                                computeGrad, 
                                computeHess,
                                decideIncludePrior,
                                MCMC, 
                                MCMC1D,
                                project = model, resetFunctions = TRUE) 
  return(compiledFuns)
}

nfLLR <- function(model, paramNodes, stateNodes, mvSamples) {
  nf <- nimbleFunction(
    setup = function(model, paramNodes, stateNodes, mvSamples){
      ## What to calculate when we replace parameters
      paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      ## What to calculate for each row of the MC sample
      calcNodes <- model$getDependencies(stateNodes)
    },
    run = function(P = double(1), delta = double(1), includePrior = logical(0),
                   burninFrac = double(0)) {
      ans <- 0
      m <- getsize(mvSamples)
      lower <- ceiling(burninFrac) + 1
      for(i in lower:m) {
        ## put MCMC output into model
        copy(from = mvSamples, to = model, row = i, logProb = FALSE)
        ## Put P + delta parameters in the model and get numerator
        values(model, paramNodes) <<- P + delta
        model$calculate(paramDeps)
        logNumerator <- model$calculate(calcNodes)
        ## Put P in the model and get denominator
        values(model, paramNodes) <<- P
        model$calculate(paramDeps)
        logDenominator <- model$calculate(calcNodes)
        ## tally result
        ans <- ans + exp(logNumerator - logDenominator)
      }
      ans <- ans / (m - lower + 1)
      if (includePrior) {
        values(model, paramNodes) <<- P + delta
        model$calculate(paramNodes)
        priorRatioNum <- exp(model$getLogProb(paramNodes))
        values(model, paramNodes) <<- P
        model$calculate(paramNodes)
        priorRatioDenom <- exp(model$getLogProb(paramNodes))
        ans <- ans * priorRatioNum / priorRatioDenom
      }
      return(ans)
      returnType(double())
    })
  return(nf(model, paramNodes, stateNodes, mvSamples))
}
