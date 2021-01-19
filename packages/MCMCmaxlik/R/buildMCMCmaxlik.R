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
  
  ## nfADgrad <- nimbleFunction(
  ##   setup = function(model, paramNodes, stateNodes, mvSamples) {
  ##     # paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
  ##     calcNodes <- c(model$getDependencies(stateNodes), 
  ##                    model$getDependencies(paramNodes, determOnly = TRUE))
  ##     D <- length(model$expandNodeNames(paramNodes))
  ##     grad <- numeric(D)
  ##   },
  ##   run = function(includePrior = logical(0), burninFrac = double(0)) {
  ##     m <- getsize(mvSamples)
  ##     ans <- numeric(D)
  ##     lower <- ceiling(burninFrac * m) + 1
  ##     for (i in lower:m) {
  ##       copy(from=mvSamples, to=model, row=i, logProb=FALSE)
  ##       temp <- nimDerivs(model$calculate(calcNodes), wrt=paramNodes, order = c(1,2))
  ##       ans <- ans + temp$jacobian[1, ]
  ##     }
  ##     grad <<- ans / (m - lower + 1)
  ##   }
  ## )
  
  nfADGradHess <- nimbleFunction(
    setup = function(model, paramNodes, stateNodes, mvSamples) {
      # paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      calcNodes <- c(model$getDependencies(stateNodes), 
                     model$getDependencies(paramNodes, determOnly = TRUE))
      calcNodes <- model$topologicallySortNodes(calcNodes)
      D <- length(model$expandNodeNames(paramNodes))
      grad <- numeric(D)
      hess <- nimMatrix(0, nrow = D, ncol = D)
    },
    run = function(includePrior = logical(0), burninFrac = double(0), gradient = logical(0), hessian = logical(0)) {
      m <- getsize(mvSamples)
      order <- numeric(length = gradient + hessian)
      if(gradient) {
        ansGrad <- numeric(D)
        order[1] <- 1
      }
      if(hessian) {
        ansHess <- nimMatrix(0, nrow = D, ncol = D)
        if(gradient)
          order[2] <- 2
        else
          order[1] <- 2
      }
      lower <- ceiling(burninFrac * m) + 1
      for (i in lower:m) {
        copy(from=mvSamples, to=model, row=i, logProb=FALSE)
        temp <- nimDerivs(model$calculate(calcNodes), wrt=paramNodes, 
                          order = order)
        if(gradient)
          ansGrad <- ansGrad + temp$jacobian[1, ]
        if(hessian) {
          for(j in 1:D) {
            for(k in 1:D) {
              ansHess[j,k] <- ansHess[j,k] + temp$hessian[j, k, 1]
            }
          }
        }
      }
      if(hessian)
        hess <<- ansHess / (m - lower + 1)
      if(gradient)
        grad <<- ansGrad / (m - lower + 1)
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
  
  setParamsInternal <- nfSetParams(model, paramNodes)
  # setParam and setLatent will become wrappers to the compiled nimbleFunctions.
  # This will allow them to be used without modification of any calling code.
  # To do so, we make a list with an elemet named "run".
  setParams <- list(run = function(p, checkNames = FALSE) {
    # If a matrix was provided, use the last row, in a way that preserved column names
    if(is.matrix(p)) p <- p[ dim(p)[1], ]
    # If there are no names, use directly, assuming p is ordered to match paramNodes
    if(is.null(names(p))) return(compiledFuns$setParamsInternal$run(p))
    # If there are names, use them to order p correctly to match paramNodes.
    # There is by default error checking on correctness or completeness of names of p.  They must match!
    if(checkNames) {
      pTrial <- p[paramNodes]
      if(any(is.na(names(pTrial)))) {
        iBad <- which(is.na(names(pTrial)))
        warning(paste("Problem with names in setParams. Names of paramNodes that were expected but not found include", paste(paramNodes[iBad], sep=", ", collapse = ", ")))        
      }
      return(compiledFuns$setParamsInternal$run(pTrial))
    }
    compiledFuns$setParamsInternal$run(p[paramNodes])
  })
  latentNodes <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE)
  setLatentInternal <- nfSetLatent(model, latentNodes)
  setLatent <- list(run = function(p, checkNames = FALSE) {
    # If a matrix was provided, use the last row, in a way that preserved column names
    if(is.matrix(p)) p <- p[ dim(p)[1], ]
    # If there are no names, use directly, assuming p is ordered to match paramNodes
    if(is.null(names(p))) return(compiledFuns$setLatentInternal$run(p))
    # If there are names, use them to order p correctly to match paramNodes.
    # There is by default error checking on correctness or completeness of names of p.  They must match!
    if(checkNames) {
      pTrial <- p[latentNodes]
      if(any(is.na(names(pTrial)))) {
        iBad <- which(is.na(names(pTrial)))
        warning(paste("Problem with names in setLatent. Names of latentNodes that were expected but not found include", paste(paramNodes[iBad], sep=", ", collapse = ", ")))
      }
      return(compiledFuns$setLatentInternal$run(pTrial))
    }
    compiledFuns$setLatentInternal$run(p[latentNodes])
  })
  mcmcConf <- configureMCMC(model, nodes = latentNodes, monitors = latentNodes)
  MCMC <- buildMCMC(mcmcConf)
  
  #  computeGrad <- nfADgrad(model, paramNodes, latentNodes, MCMC$mvSamples)
  computeGradHess <- nfADGradHess(model, paramNodes, latentNodes, MCMC$mvSamples)
  decideIncludePrior <- nfIncludePrior()
  mcmc1DConf <- configureMCMC(model)
  mcmc1DConf$removeSamplers(paramNodes) 
  mcmc1DConf$addSampler(target=paramNodes, type=sampler_RW_rotated,
                        control=list(computeGrad=computeGradHess,
                                     decideIncludePrior=decideIncludePrior))
  MCMC1D <- buildMCMC(mcmc1DConf)
  compiledFuns <- compileNimble(setParamsInternal, setLatentInternal, #computeGrad, 
                                computeGradHess,
                                decideIncludePrior, MCMC, MCMC1D,
                                project = model, resetFunctions = TRUE)
  result <- list(setParams = setParams,
                 setLatent = setLatent,
                 computeGradHess = compiledFuns$computeGradHess,
                 decideIncludePrior = compiledFuns$decideIncludePrior,
                 MCMC = compiledFuns$MCMC,
                 MCMC1D = compiledFuns$MCMC1D)
  #  return(compiledFuns)
  return(result)
}
