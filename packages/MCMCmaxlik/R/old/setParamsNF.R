#' NIMBLE function to update the value of a node

#' @param model a NIMBLE model
#' @param paramNodes a character vector containing the names of 
#' parameters of interest (usually the top-level parameters).
#' @param P: vector of values to put into \code{paramNodes}.
#' @export
#' @return update the value of \code{paramNodes}.

setParamsNF <- function(model, paramNodes) {
  nf <- nimbleFunction(
    setup = function(model, paramNodes) {
      deps <- model$getDependencies(paramNodes, determOnly = TRUE)
    },
    run = function(P = double(1)) {
      values(model, paramNodes) <<- P
      model$calculate(deps)
    }
  )
  return(nf(model, paramNodes))
}
