library("dlm")

# Load the data.
dt <- read.csv("statespaceModData.csv")
df <- dt[, -1]
colnames(df) <- c("x", "y")

# Model setup
# Remark: dlm does not seem to support uncentered state transitions.
# The MLE of mu is the sample mean of y (I think).
mu <- mean(df$y)

# The data is then centered before fitting the model.
yt <- scale(df$y, center=TRUE, scale=FALSE)
# The parameters to be estimated are rho, sigPN, sigOE.
# Set parameter restrictions: rho in (-1, 1), sigPN > 0, sigOE > 0.
param_rest <- function(param){
  return( c(exp(param[1])/(1+exp(param[1])), exp(param[2]), exp(param[3])) ) 
}

# Set up the state-space model
ssm <- function(param){
  param <- param_rest(param)
  return( dlm(FF=1,V=param[3]^2,GG=param[1],W=param[2]^2,
              m0=0,C0=param[2]^2/(1-param[1]^2)) )
}

# MLE
# Remark: dlm assumes that y0 does not exist.
# -> need to drop the first element in the simulated data.
fit1 <- dlmMLE(y=yt[2:length(yt)], parm=c(0,1,1), build=ssm, hessian=T)
# rho, sigPN, sigOE
param_rest(fit1$par) # 0.3330513 5.2121850 5.3901920
# mu
mu # 7.332008