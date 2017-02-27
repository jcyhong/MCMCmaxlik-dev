#' NIMBLE function approximating likelihood ratio using finite differences
#'
#' @param model a compiled NIMBLE model
#' @param paramNodes a character vector containing the names of 
#' parameters of interest (usually the top-level parameters).
#' @param stateNodes a character vector containing names of latent nodes
#' @param mvSamples modelValues object that contains the MCMC samples
#' @param P value to potentially add a delta to and evaluate the likelihood ratio at
#'   delta: vector of length equal to the lenth of paramNodes, delta for each parameter, 
#'   only one should be non-zero.
#' @export
#' @author Perry de Valpine

nfLLR <- function(model, paramNodes, stateNodes, mvSamples) {
  nf <- nimbleFunction(
    setup = function(model, paramNodes, stateNodes, mvSamples){
      ## What to calculate when we replace parameters
      paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      ## What to calculate for each row of the MC sample
      calcNodes <- model$getDependencies(stateNodes)
    },
    run = function(P = double(1), delta = double(1)) {
      ans <- 0
      m <- getsize(mvSamples)
      for(i in 1:m) {
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
      ans <- ans / m
      return(ans)
      returnType(double())
    }) 
  return(nf(model, paramNodes, stateNodes, mvSamples))
}