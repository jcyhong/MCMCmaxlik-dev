#' Finding the mode of a kernel density estimate of samples
#' 
#' @description
#' \code{getKernelMode} computes the mode of the the kernel density estimate 
#' of samples.
#' @param samples a vector of samples (1 dimensional samples) 
#' @param kern the name of the kernel. Default to \code{"gaussian"}.
#' @param bdwth bandwidth for kernel. Default to \code{"nrd0"} 
#' (rule-of-thumb for choosing the bandwidth of a Gaussian kernel 
#' density estimator).
#' @return mode of estimated density
#' @export
#' @examples 
#' set.seed(1)
#' test <- rnorm(100,0,1)
#' plot(density(test))
#' mm <- getKernelMode(test) ##0.3294585
#' abline(v=mm)
#' @references 
#' http://stackoverflow.com/questions/16255622/peak-of-the-kernel-density-estimation

getKernelMode <- function(samples, kern="gaussian", bdwth="nrd0") {
  d <- density(samples, bw=bdwth, kernel=kern)
  hdrResult <- hdrcde::hdr(den=d, prob=0)
  return(hdrResult$mode)
}
