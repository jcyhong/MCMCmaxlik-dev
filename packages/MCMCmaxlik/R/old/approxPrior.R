#' Computing the Prior Ratio of Finite Difference
#' 
#' @description
#' \code{approxPrior} returns the prior ratio of finite difference
#' (that is, the prior density at theta plus an increment divided 
#' by the prior density at theta). 
#' 
#' @param theta values of each parameter to evaluate the gradient at
#' @param model a compiled NIMBLE model
#' @param paramNodes a character vector containing the names of 
#' parameters of interest (usually the top-level parameters).
#' @param delta a scalar indicating the finite difference. Default to 1e-4.
#' @export
#' @return finite difference ratio of the prior evaluated at theta

approxPrior <- function(theta, model, paramNodes,
                        delta=1e-4) {
  D <- length(paramNodes)
  v <- rep(NA, D)
  values(model,paramNodes) <- theta
  model$calculate(paramNodes)
  vDenom <- model$getLogProb(paramNodes)
  vDenom <- exp(vDenom)
  for (i in 1:D) {
    inc <- rep(0, D)
    inc[i] <- delta
    values(model, paramNodes) <- theta + inc
    model$calculate(paramNodes)
    
    v[i] <- model$getLogProb(paramNodes)
    v[i] <- exp(v[i])
  }
  return(v / vDenom)
}
