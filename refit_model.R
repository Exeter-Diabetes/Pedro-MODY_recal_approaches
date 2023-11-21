#:--------------------------------------------------------
#   This file contains the refit model of the MODY calculator
#
#:--------------------------------------------------------

library(nimble)
library(tidyverse)


### Load representative dataset
# representative_dataset <- readRDS("representative_dataset.rds")


### Fit Bayesian Refit nimble model
model_code <- nimbleCode({
  
  ## likelihood component
  for (j in 1:nind) {
    
    ## prediction from UNITED
    logit(p[j]) <- beta0 + inprod(beta[1:np], x[j, 1:np])
    
    ## likelihood
    M[j] ~ dbern(p[j])
    
  }
  
  ## regression priors
  beta0 ~ dnorm(0, 0.01)
  for(j in 1:np) {
    beta[j] ~ dnorm(0, 0.01)
  }
  
})


#:--- UNITED covariates
x <- as.matrix(select(representative_dataset, pardm, agerec, hba1c, agedx, sex))
M <- representative_dataset$mody


### set up data list for NIMBLE
data <- list(M = M)

### set up other components of model
consts <- list(nind = length(M),             # number of individuals in UNITED
               np = ncol(x),                 # number of covariates
               x = x)                        # covariate matrix for UNITED\

### set up initial values for parameters in the model
initFn <- function(M, x, np) {
  
  ## simulate from priors
  beta0 <- rnorm(1)
  beta <- rnorm(np)
  
  ## simulate missing information
  Z <- t(beta %*% t(x))
  M1 <- Z
  M1 <- exp(M1) / (1 + exp(M1))
  M1 <- rbinom(length(M1), 1, M1)
  
  M1[!is.na(M)] <- NA
  
  ## return initial values
  list(beta0 = beta0,
       beta = beta,
       M = M1
  )
}

# adding this so that the logProb isn't -Inf
logProb <- "-Inf"

while (logProb == "-Inf") {
  
  ### set list of initial values for NIMBLE
  inits <- initFn(data$M, consts$x, consts$np)
  
  ### define the model, data, inits and constants
  model <- nimbleModel(code = model_code, constants = consts, data = data, inits = inits)
  
  logProb <- model$calculate()
  
}


### compile the model
cModel <- compileNimble(model)

### configure MCMC and monitor parameters needed for prediction
config <- configureMCMC(cModel)

### build the model
built <- buildMCMC(config)
cBuilt <- compileNimble(built)

### run the model
run <- runMCMC(cBuilt,
               niter = 1000,
               nburnin = 1000,
               nchains = 2,
               progressBar = TRUE,
               summary = TRUE,
               samplesAsCodaMCMC = TRUE,
               thin = 1)
