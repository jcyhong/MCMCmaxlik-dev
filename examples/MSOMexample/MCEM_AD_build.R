## Working on update to MCEM so Johnny can use it with AD

## Notes for future reference
## getMCEMRanges should not be a nimbleFunction
## We could check that maxNodes are actually connected to latentNodes instead of including all non-data stochastic nodes that are not latentNodes
## We could use the parameterTransform system to avoid box constraints.
## mcmc_Latent_Conf should not have to monitor all var names
## burnIn can be a run-time parameter
## In calc_E_llk_gen: paramValues are repeatedly put in model even when diff is FALSE.
##                    fixedCalcNodes are never used
##                    need for paramDepDetermNodes_latent is not clear
## The MBB is really wasteful in re-calculating values.
## 1/nSamples should be 1/(nSamples - burnin), but the factor is different only be a constant so might not matter.
##    It will however matter for the std. err. 
## Pathological behavior: if you start at the MLE, it spins the Monte Carlo sample size up enormously to try to get a confident improvement in likelihood, when there is none possible.

# Calculates Q function if diff = 0, calculates difference in Q functions if diff = 1.
calc_E_llk_gen_AD = nimbleFunction(
    name = 'calc_E_llk_gen_AD',
  setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0){
   ## fixedCalcNodes <- model$getDependencies(fixedNodes)	
    latentCalcNodes <- model$getDependencies(sampledNodes)
    lengthFixedNodes <- length(model$expandNodeNames(fixedNodes, returnScalarComponents = TRUE))
    paramDepDetermNodes_fixed <- model$getDependencies(fixedNodes, determOnly = TRUE)
    allCalcNodes <- model$topologicallySortNodes(c(paramDepDetermNodes_fixed, latentCalcNodes))
    paramDepDetermNodes_latent <- model$getDependencies(latentCalcNodes, determOnly = TRUE) ## to remove
    areFixedDetermNodes <- length(paramDepDetermNodes_fixed) > 0
    areLatentDetermNodes <- length(paramDepDetermNodes_latent) >0 ## to remove
    logLik <- as.numeric(0)
    gradientLik <- rep(as.numeric(0), lengthFixedNodes)
    lastParamValues <- rep(-Inf, lengthFixedNodes)
  },
  run = function(paramValues = double(1)) {
    burnin <- burnIn
    run_internal(paramValues, burnin)
    lastParamValues <<- paramValues
    return(logLik)
    returnType(double(0))
  },
  methods = list(
    resetForNewSample = function() {
      lastParamValues <<- rep(-Inf, lengthFixedNodes)
    },
    run_internal = function(paramValues = double(1), burnin = integer(0, default = 0)) {
      nSamples <- getsize(mvSample)
      sum_LL <- 0
      sum_gradient <- numeric(value = 0, length = lengthFixedNodes)
      values(model, fixedNodes) <<- paramValues
      for(i in (burnin+1):nSamples){
        nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
        sample_derivs <- nimDerivs(model$calculate(allCalcNodes), wrt = fixedNodes, order = c(0, 1))
        sample_LL <- sample_derivs$value[1]
        if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
          stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')' to see if node 'y[3]' is the cause of the problem). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
        sum_LL <- sum_LL + sample_LL
        sum_gradient <- sum_gradient + sample_derivs$jacobian[1,]
      }
      logLik <<- sum_LL / (nSamples - burnIn)
      if(is.nan(logLik))
        logLik <<- -Inf
      gradientLik <<- sum_gradient / (nSamples - burnIn)
    },
    sample_logLiks = function(paramValues = double(1)) {
      burnin <- burnIn
      nSamples <- getsize(mvSample)
      values(model, fixedNodes) <<- paramValues
      logLiks <- numeric(length = nSamples - burnin, init = FALSE)
      for(i in (burnin+1):nSamples){
        nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
        logLiks[i - burnin] <- model$calculate(allCalcNodes)
      }
      return(logLiks)
      returnType(double(1))
    },
    gradient = function(paramValues = double(1)) {
      burnin <- burnIn
      if(all(paramValues == lastParamValues))
        return(gradientLik)
      run_internal(paramValues, burnin)
      lastParamValues <<- paramValues
      return(gradientLik)
      returnType(double(1))
    },
    ## old run
    old_run = function(paramValues = double(1), oldParamValues = double(1), diff = integer(0)){
      nSamples = getsize(mvSample)
      mean_LL <- 0
      
      for(i in (burnIn+1):nSamples){
        nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
        values(model, fixedNodes) <<- paramValues  #first use new params, then old ones
        if(areFixedDetermNodes){
          simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
        }
        if(areLatentDetermNodes){
          simulate(model, paramDepDetermNodes_latent)	#	Fills in the deterministic nodes
        }
        sample_LL = calculate(model, latentCalcNodes)
        if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
          stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')' to see if node 'y[3]' is the cause of the problem). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
        mean_LL = mean_LL + sample_LL
        if(diff == 1){
          values(model, fixedNodes) <<- oldParamValues #now old params
          if(areFixedDetermNodes){
            simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
          }
          if(areLatentDetermNodes){
            simulate(model, paramDepDetermNodes_latent)  #	Fills in the deterministic nodes
          }
          sample_LL = calculate(model, latentCalcNodes)
          if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
            stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')'). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
          mean_LL = mean_LL - sample_LL
        }
      }
      values(model, fixedNodes) <<- paramValues  #first use new params, then old ones
      if(areFixedDetermNodes){
        simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
      }
      mean_LL <- mean_LL / nSamples
      if(is.nan(mean_LL)){
        mean_LL = -Inf	
      }
      returnType(double())
      return(mean_LL)
    }
  )
)

getMCEMRanges_AD <- function(model, maxNodes, buffer){
  low_limits = rep(-Inf, length(maxNodes) ) 
  hi_limits  = rep(Inf,  length(maxNodes) )
  nodes <- model$expandNodeNames(maxNodes)
  for(i in seq_along(nodes)) {
    low_limits[i] = getBound(model, nodes[i], 'lower') + abs(buffer)
    hi_limits[i]  = getBound(model, nodes[i], 'upper')  - abs(buffer)
  }
  return(list(low_limits, hi_limits))
}

buildMCEM_AD <- function(model,
                      latentNodes,
                      burnIn = 500 ,
                      mcmcControl = list(adaptInterval = 100),
                      boxConstraints = list(),
                      buffer = 10^-6,
                      alpha = 0.25,
                      beta = 0.25, 
                      gamma = 0.05,
                      C = 0.001,
                      numReps = 300,
                      forceNoConstraints = FALSE,
                      verbose = TRUE) {
  require(mcmcse)
  latentNodes = model$expandNodeNames(latentNodes)
  latentNodes <- intersect(latentNodes, model$getNodeNames(stochOnly = TRUE))
  dataNodes <- model$getNodeNames(dataOnly = TRUE)
  allStochNonDataNodes = model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
  if(buffer == 0)
    warning("'buffer' is zero. This can cause problems if the likelihood function is degenerate on boundary")
  if(buffer < 0)
    stop("'buffer' must be non-negative.")
  
  if(length(setdiff(latentNodes, allStochNonDataNodes) ) != 0 )
    stop('latentNodes provided not found in model')
  maxNodes <- model$expandNodeNames(setdiff(allStochNonDataNodes, latentNodes),
                                    returnScalarComponents = TRUE)
  if(any(model$isDiscrete(maxNodes)))
      stop(paste0("MCEM cannot optimize over discrete top-level parameters. The following top-level parameters in your model are discrete: ",
                  paste0(maxNodes[model$isDiscrete(maxNodes)], collapse = ', ')))
  
  limits <- getMCEMRanges_AD(model, maxNodes, buffer)
  low_limits = limits[[1]]
  hi_limits  = limits[[2]]
  
  ## will be assign()'ed from within run() for use in getAsymptoticCov
  paramMLEs <- NA
  mcmcIters <- NA
  
  constraintNames = list()
  for(i in seq_along(boxConstraints) )
    constraintNames[[i]] = model$expandNodeNames(boxConstraints[[i]][[1]])
  for(i in seq_along(constraintNames) ) {
    limits = boxConstraints[[i]][[2]]
    inds = which(maxNodes %in% constraintNames[[i]])
    if(length(inds) == 0)
      stop(paste("warning: provided a constraint for nodes", constraintNames[[i]], ", but those nodes do not exist in the model!"))
    tooLowNodes <- which(limits[1] + abs(buffer) < low_limits[inds])
    tooHighNodes <- which(limits[2] - abs(buffer) > hi_limits[inds])
    if(length(tooLowNodes) > 0) warning(paste0("User-specified lower bound for ", constraintNames[[i]][tooLowNodes],
                                               " is below lower bound detected by NIMBLE.  "))
    if(length(tooHighNodes) > 0) warning(paste0("User-specified upper bound for ", constraintNames[[i]][tooHighNodes],
                                                " is above upper bound detected by NIMBLE.  "))
    low_limits[inds] = limits[1] + abs(buffer)
    hi_limits[inds] = limits[2] - abs(buffer)
  }
  if(any(low_limits>=hi_limits))
    stop('lower limits greater than or equal to upper limits!')
  if(forceNoConstraints ||
     identical(low_limits, rep(-Inf, length(low_limits))) &&
     identical(hi_limits, rep(Inf, length(hi_limits))))
    optimMethod = "BFGS"
  else 
    optimMethod = "L-BFGS-B"

  if(length(latentNodes) == 0)
    stop('no latentNodes')
  
  if(length(maxNodes) == 0)
    stop('no nodes to be maximized over')
  resetFunctions <- FALSE
  if(is(model, "RmodelBaseClass") ){
    Rmodel = model
    if(is(model$CobjectInterface, "uninitializedField")){
      cModel <- compileNimble(model)
    }
    else{
      cModel = model$CobjectInterface
      resetFunctions <- TRUE
    }
  }
  else{
    cModel <- model
    Rmodel <- model$Rmodel
    resetFunctions <- TRUE
  }
  
  zAlpha <- qnorm(alpha, 0, 1, lower.tail=FALSE)
  zBeta <- qnorm(beta, 0, 1, lower.tail=FALSE)
  zGamma <- qnorm(gamma, 0, 1, lower.tail=FALSE)
  
  mcmc_Latent_Conf <- configureMCMC(Rmodel, nodes = latentNodes, monitors = model$getVarNames(), control = mcmcControl, print = FALSE)
  Rmcmc_Latent <- buildMCMC(mcmc_Latent_Conf)
  sampledMV <- Rmcmc_Latent$mvSamples
  mvBlock <- modelValues(Rmodel)
  Rcalc_E_llk <- calc_E_llk_gen_AD(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)
#  RvarCalc <- calc_asympVar(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvBlock, mvSample = sampledMV, numReps = numReps)
#  RgetCov <- bootstrapGetCov(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)
  
  cmcmc_Latent = compileNimble(Rmcmc_Latent, project = Rmodel, resetFunctions = resetFunctions)
#  cGetCov = compileNimble(RgetCov, project = Rmodel)  
#  cvarCalc <- compileNimble(RvarCalc, project = Rmodel)
  cCalc_E_llk = compileNimble(Rcalc_E_llk, project = Rmodel)  
  nParams = length(maxNodes)
  run <- function(initM = 1000){
    if(burnIn >= initM)
      stop('mcem quitting: burnIn > initial m value')
    cmcmc_Latent$run(1, reset = TRUE)	# To get valid initial values 
    theta <- values(cModel, maxNodes)
    if(optimMethod == "L-BFGS-B"){
      for(i in seq_along(maxNodes) ) {  # check that initial values satisfy constraints
        if(identical(low_limits[i], -Inf) && (hi_limits[i] < Inf)){
          if(theta[i] > hi_limits[i]){
            theta[i] <- hi_limits[i] - 1
          }
        }
        else if(identical(hi_limits[i], Inf) && (low_limits[i] > -Inf)){
          if(theta[i] < low_limits[i]){
            theta[i] <- low_limits[i] + 1
          }
        }
        else if((low_limits[i] > -Inf) && (hi_limits[i] < Inf)){
          if(!(theta[i] >= low_limits[i] & theta[i] <= hi_limits[i])){
            theta[i] = (low_limits[i] + hi_limits[i])/2
          }
        }
      }
      values(cModel, maxNodes) <<- theta
      cModel$simulate(cModel$getDependencies(maxNodes, self = FALSE))
    }
    m <- initM 
    endCrit <- C+1 #ensure that first iteration runs
    sigSq <-0 #use initM as m value for first step
    diff <- 1 # any nonzero value can be used here, gets overwritten quickly in algo
    itNum <- 0
    while(endCrit > C){ 
      acceptCrit <- 0
      #starting sample size calculation for this iteration
      m <- burnIn + ceiling(max(m - burnIn, sigSq*((zAlpha + zBeta)^2)/((diff)^2)))
      cmcmc_Latent$run(m, reset = TRUE)   #initial mcmc run of size m
      thetaPrev <- theta  #store previous theta value
      itNum <- itNum + 1
      while(acceptCrit == 0){
        cCalc_E_llk$resetForNewSample()
        if(optimMethod == "L-BFGS-B")
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, gr = cCalc_E_llk$gradient,
                                control = list(fnscale = -1), method = 'L-BFGS-B', 
                                lower = low_limits, upper = hi_limits)
        if(optimMethod == "BFGS")
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, gr = cCalc_E_llk$gradient,
                                control = list(fnscale = -1), method = 'BFGS')
        
        theta = optimOutput$par
        sample_logLiks_prev <- cCalc_E_llk$sample_logLiks(thetaPrev)
        sample_logLiks_new <- cCalc_E_llk$sample_logLiks(theta)
        mcse_result <- mcse(sample_logLiks_new - sample_logLiks_prev,
                            size = ceiling(min(1000, (m - burnIn)/20)), ## backward compatible
                            method = "obm",
                            r = 1) ## What is the lugsail?
        ##        sigSq <- cvarCalc$run(m, theta, thetaPrev) 
        ## ase <- sqrt(sigSq) #asymptotic std. error
        ase <- mcse_result$se
        diff <- mean(sample_logLiks_new) - mean(sample_logLiks_prev)
        ## Should have: mean(sample_logLiks_new) == mcse_result$est == optimOutput$value
        ## diff <- cCalc_E_llk$run(theta, thetaPrev, 1)
        if((diff - zAlpha*ase)<0){ #swamped by mc error
          cat("Monte Carlo error too big: increasing MCMC sample size.\n")
          mAdd <- ceiling((m-burnIn)/2)  #from section 2.3, additional mcmc samples will be taken if difference is not great enough
          cmcmc_Latent$run(mAdd, reset = FALSE)
          m <- m + mAdd
        }
        else{
          acceptCrit <- 1
          endCrit <- diff + zGamma*ase #evaluate ending criterion
          if(itNum == 1)
            endCrit <- C+1 #ensure that at least two iterations are run
          
          if(verbose == T){
            cat("Iteration Number: ", itNum, ".\n", sep = "")
            cat("Current number of MCMC iterations: ", m, ".\n", sep = "")
            output = optimOutput$par
            names(output) = maxNodes
            cat("Parameter Estimates: \n", sep = "")
            print(output)
            cat("Convergence Criterion: ", endCrit, ".\n", sep = "")
          }
        }
      }
    }
    output <- optimOutput$par
    assign('paramMLEs', output, envir = parent.env(environment()))
    assign('mcmcIters', m, envir = parent.env(environment()))
    names(output) <- maxNodes
    return(output)
  }

  ## This has not been updated
  estimateCov = function(MLEs = NA, useExistingSamples = FALSE){
    delta <- .0001
    if(!(length(MLEs) == 1 && is.na(MLEs))){
      if(!is.numeric(MLEs)) stop("MLEs argument must be numeric.")
      if(!identical(sort(names(MLEs)), sort(maxNodes))){
        stop(paste('MLEs argument must be a named vector with MLEs for all of the following parameters: ', paste(maxNodes, collapse = ", ")))
      }
      covMLEs <- unname(MLEs[maxNodes])
    }
    else{
      if(length(paramMLEs) == 1 && is.na(paramMLEs)){
        stop(paste('No MLEs argument was provided, and the run() method has not been called yet.  Please call the run() method first or provide a named vector of MLEs.'))
      }
      else{
        covMLEs <- unname(paramMLEs)
      }
    }
    if(dim(as.matrix(cmcmc_Latent$mvSamples))[1]<2){
      if(useExistingSamples == TRUE){
        warning('MCMC over latent states has not been run yet, cannot have useExistingSamples = TRUE')
        useExistingSamples <- FALSE
      }
    }
    values(cModel, maxNodes) <- covMLEs
    calculate(cModel, cModel$getDependencies(maxNodes))
    if(!useExistingSamples){
      if(is.na(mcmcIters)) mcmcIters <- 20000
      cmcmc_Latent$run(mcmcIters)
    }
    FIM <- cGetCov$run(covMLEs, delta)
    cov <- solve(FIM)
    colnames(cov) <- maxNodes
    rownames(cov) <- maxNodes
    return(cov)
  }
  return(list(run = run,
              estimateCov = estimateCov))
}

## This has not been updated
## bootstrapGetCov <- nimbleFunction(
##     name = 'bootstrapGetCov',
##     setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0){
##       fixedCalcNodes <- model$getDependencies(fixedNodes)	
##       latentCalcNodes <- model$getDependencies(sampledNodes)
##       paramDepDetermNodes_fixed <- model$getDependencies(fixedNodes, determOnly = TRUE) 
##       paramDepDetermNodes_latent <- model$getDependencies(latentCalcNodes, determOnly = TRUE) 
##       areFixedDetermNodes <- length(paramDepDetermNodes_fixed) > 0
##       areLatentDetermNodes <- length(paramDepDetermNodes_latent) > 0
##     },
##   run = function(theta = double(1), delta = double(0)){
##     nSamples <- getsize(mvSample)
##     paramLengths <- length(theta)
##     fxph <- numeric(paramLengths)
##     fxmh <- numeric(paramLengths)
##     grad <- numeric(paramLengths)
##     derivxy <- matrix(0, nrow = paramLengths, ncol = paramLengths)
##     meanGrad <- numeric(paramLengths)
##     meanGradGrad <-  matrix(0, nrow = paramLengths, ncol = paramLengths)
##     meanDerivxy <- matrix(0, nrow = paramLengths, ncol = paramLengths)
##     for(i in (burnIn+1):nSamples){
##       values(model, fixedNodes) <<- theta
##       if(areFixedDetermNodes){
##         simulate(model, paramDepDetermNodes_fixed)  
##       }
##       if(areLatentDetermNodes){
##         simulate(model, paramDepDetermNodes_latent)
##       }
##       nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
##       if(areFixedDetermNodes){
##         simulate(model, paramDepDetermNodes_fixed) 
##       }
##       if(areLatentDetermNodes){
##         simulate(model, paramDepDetermNodes_latent)	
##       }
##       origValue <- calculate(model, latentCalcNodes)
##       for(iNode in 1:paramLengths){
##         theta[iNode] <- theta[iNode] + delta
##         values(model, fixedNodes) <<- theta 
##         if(areFixedDetermNodes){
##           simulate(model, paramDepDetermNodes_fixed)  
##         }
##         if(areLatentDetermNodes){
##           simulate(model, paramDepDetermNodes_latent)	
##         }
##         fxph[iNode] <- calculate(model, latentCalcNodes)
##         theta[iNode] <- theta[iNode] - 2*delta
##         values(model, fixedNodes) <<- theta 
##         if(areFixedDetermNodes){
##           simulate(model, paramDepDetermNodes_fixed)  
##         }
##         if(areLatentDetermNodes){
##           simulate(model, paramDepDetermNodes_latent)
##         }
##         fxmh[iNode] <- calculate(model, latentCalcNodes)
##         grad[iNode] <- (fxph[iNode] - fxmh[iNode])/(2*delta)
##         theta[iNode] <- theta[iNode] + delta
##         derivxy[iNode, iNode] <- (fxph[iNode] -2*origValue + fxmh[iNode])/(delta^2)
##       }
##       for(iNode in 1:paramLengths){
##         if(iNode != paramLengths){
##           for(jNode in (iNode+1):paramLengths){
##             theta[iNode] <- theta[iNode] + delta
##             theta[jNode] <- theta[jNode] + delta
##             values(model, fixedNodes) <<- theta 
##             if(areFixedDetermNodes){
##               simulate(model, paramDepDetermNodes_fixed) 
##             }
##             if(areLatentDetermNodes){
##               simulate(model, paramDepDetermNodes_latent)	
##             }
##             fxyph <- calculate(model, latentCalcNodes)
##             theta[iNode] <- theta[iNode] - 2*delta
##             theta[jNode] <- theta[jNode] - 2*delta
##             values(model, fixedNodes) <<- theta
##             if(areFixedDetermNodes){
##               simulate(model, paramDepDetermNodes_fixed) 
##             }
##             if(areLatentDetermNodes){
##               simulate(model, paramDepDetermNodes_latent)
##             }
##             fxymh <- calculate(model, latentCalcNodes)
##             derivxy[iNode, jNode] <- (fxyph - fxph[iNode] - fxph[jNode] + 2*origValue - fxmh[iNode] - fxmh[jNode] + fxymh)/(2*delta^2)
##             theta[iNode] <- theta[iNode] + delta
##             theta[jNode] <- theta[jNode] + delta
##           }
##         }
##       }
##       meanGrad <- meanGrad + grad
##       meanDerivxy <- meanDerivxy + derivxy
##       meanGradGrad <- meanGradGrad + (grad%*%t(grad))
##     }
##     meanGrad <- meanGrad / (nSamples - burnIn)
##     meanDerivxy <- meanDerivxy /  (nSamples - burnIn)
##     meanGradGrad <- meanGradGrad /  (nSamples - burnIn)
##     for(iNode in 1:paramLengths){
##       if(iNode != paramLengths){
##         for(jNode in (iNode+1):paramLengths){
##           meanDerivxy[jNode, iNode] <- meanDerivxy[iNode, jNode] ## smarter way to do this with matrix mult?
##         }
##       }
##     }
##     returnType(double(2))
##     returnMat <- -meanDerivxy - meanGradGrad +(meanGrad%*%t(meanGrad))
##     return(returnMat)
##   },where = getLoadingNamespace())
