#' Approximating Posterior Gradient via MCMC and Finite Difference Approximation
#' 
#' @description
#' \code{approxGrad} returns an finite difference approximation of the 
#' gradient of the posterior density using conditional MCMC samples of
#' latent variables.
#' 
#' @param theta a numeric vector containing the parameter values at which
#'  the gradient is evaluated.
#' @param model a compiled NIMBLE model.
#' @param paramNodes a character vector containing the names of 
#' parameters of interest (usually the top-level parameters).
#' @param compiledFuns a list of functions, obtained by \code{setupMCGfuns()}.
#' @param numMCMCSamples a scalar indicating the MCMC sample size. Default to 10000.
#' @param delta a scalar indicating the finite difference. Default to 1e-4.
#' @param postMode a logical value indicating whether prior is factored into 
#' the calculation. Default to \code{TRUE}.
#'
#' @return vNew estimated marginal posterior gradient vector at the 
#' specified \code{theta}
#' @export

approxGrad <- function(theta, model, paramNodes, compiledFuns,
                       numMCMCSamples=10000, delta=1e-4,postMode=T) {
  D <- length(paramNodes)
  compiledFuns$setParams$run(theta)
  compiledFuns$MCMC$run(numMCMCSamples)
  v <- rep(NA, D)
  for (i in 1:D) {
    inc <- rep(0, D)
    inc[i] <- delta
    v[i]<-compiledFuns$MCLR$run(theta, inc)
  }
  
  if(postMode){
    vP=approxPrior(theta,model,paramNodes,delta)
    vNew=(v*vP-1)/delta
  }else{
    vNew=(v-1)/delta
  }
  return(vNew)
}