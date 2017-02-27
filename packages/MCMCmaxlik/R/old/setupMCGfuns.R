#' Compiling All Relevant Functions
#'
#' @param model a compiled NIMBLE model
#' @param paramNodes a character vector containing the names of 
#' parameters of interest (usually the top-level parameters).
#' @export
#' @return list of compiled functions, containing the compiled version 
#' of setParamsNF, MCMC, and likelihood ratio estimation

setupMCGfuns <- function(model, paramNodes) {
  setParams <- setParamsNF(model, paramNodes)
  latentNodes <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE)
  mcmcConf <- configureMCMC(model, nodes = latentNodes, monitors = latentNodes)
  MCMC <- buildMCMC(mcmcConf)
  MCLR <- nfLLR(model, paramNodes, latentNodes, MCMC$mvSamples)
  compiledFuns <- compileNimble(setParams, MCMC, MCLR, project = model)
  return(compiledFuns)
}
