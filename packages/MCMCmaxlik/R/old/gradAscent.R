#' Gradient (MCMC and Finite-Difference Approximated) Ascent for
#' Maximizing Marginal Likelihood 
#' @param model a compiled NIMBLE model.
#' @param paramNodes a character vector containing the names of 
#' parameters of interest (usually the top-level parameters).
#' @param compiledFuns a list of functions, obtained by \code{setupMCGfuns()}.
#' @param theta0 a vector of length equal to the length of \code{paramNodes}, 
#' starting values for each of the parameters in paramNodes
#' @param postMode a logical value indicating whether prior is factored into 
#' the calculation. Default to \code{TRUE}.
#' @param stepsize a scalar of fixed stepsize. Default to 8e-4.
#' @param maxIter a scalar indicating the maximum number of iterations. 
#' Default to 100.
#' @param numMCMCSamples a scalar indicating the MCMC sample size. Default to 10000.
#' @param delta a scalar indicating the finite difference. Default to 1e-4.
#' @param tol a scalar indicating tolerance for current gradient 
#' (difference from zero). Default to 1e-4.
#' @param boundary a list of length equal to the length of paramNodes,
#' each element of the list is a vector of length 2
#' with the minimum and the maximum for the corresponding parameter in
#' paramNodes, -Inf and Inf are accepted.
#' @export
#' @return  A list containing a list and a scalar
#' \describe{
#'  \item{\code{theta}}{a list of the estimates for each of \code{paramNodes} 
#'  per iteration}
#'  \item{\code{iter}}{iteration that the algorithm was terminated at, 
#'  if equal to \code{maxIter}, consider increasing \code{maxIter}}
#' }

gradAscent <- function(model, paramNodes, compiledFuns,
                       theta0, postMode=T,stepsize=0.0008, maxIter=100,
                       numMCMCSamples=10000,
                       delta=1e-4, tol=1e-4,boundary) {
  
  # Initialization
  thetaCur <- theta0
  iter <- 1
  thr <- Inf
  thetaList <- list(thetaCur)
  vNew <- approxGrad(thetaCur, model, paramNodes, compiledFuns,
                     numMCMCSamples=numMCMCSamples,
                     delta=delta, postMode=postMode) ## NaN starting out
  
  #print(paste("approxGrad:",vNew,collapse=""))
  

  while (iter < maxIter & thr > tol) { 
    
    thetaNew <- thetaCur + stepsize * vNew
    #print(paste("thetaNew:", thetaNew,collapse=""))
    for(i in 1:length(paramNodes)){
      if(thetaNew[i]<boundary[[i]][1]){
        thetaNew[i]=boundary[[i]][1]
      }else if(thetaNew[i]>boundary[[i]][2]){
        thetaNew[i]=boundary[[i]][2]
      }else{
        thetaNew[i]=thetaNew[i]
      }
    }
    #print(paste("thetaNewAdj:", thetaNew,collapse=""))
    values(model,paramNodes)<-thetaNew
    model$calculate(paramNodes)
   
    
    thetaList <- c(thetaList, list(thetaNew))
    vNew <- approxGrad(thetaCur, model, paramNodes, compiledFuns,
                    numMCMCSamples=numMCMCSamples, delta=delta,postMode=postMode)
    
    #print(paste("approxGrad:",vNew,collapse=""))
    
    thr <- sum((vNew)^2)
    thetaCur <- thetaNew
    iter <- iter+1
    #print(paste("threshold:",thr,collapse=""))
    print(iter)
  }
  return(list(theta=thetaList,iter=iter))
}