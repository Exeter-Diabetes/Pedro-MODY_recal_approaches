#:--------------------------------------------------------
#   This file contains the MODY calculator equation using case-control data
#
#:--------------------------------------------------------

library(nimble)
library(tidyverse)


### Load control-case dataset
# case_control_dataset <- readRDS("case_control.rds")


### Fit Bayesian nimble model
model_code <- nimbleCode({
  
  ## case-control likelihood component
  for (i in 1:nCC) {
    
    ## regression mean
    logit(pCC[i]) <- beta0 + inprod(beta[1:np], xCC[i, 1:np])
    
    ## likelihood
    MCC[i] ~ dbern(pCC[i])
  }
  
  ## regression priors
  beta0 ~ dnorm(0, 0.01)
  for(j in 1:np) {
    beta[j] ~ dnorm(0, 0.01)
  }
  
})

#:--- Case control covariates
xCC <- as.matrix(select(case_control_dataset, pardm, agerec, hba1c, agedx, sex))
MCC <- as.numeric(as.character(case_control_dataset$mody))

### set up data list for NIMBLE
data <- list(MCC = MCC)

### set up other components of model
consts <- list(nCC = length(MCC),            # number of individuals in Case control
               np = ncol(xCC),               # number of covariates
               xCC = xCC)                    # covariate matrix for Case control

### set up initial values for parameters in the model
initFn <- function(np) {
  
  ## simulate from priors
  beta0 <- rnorm(1)
  beta <- rnorm(np)
  
  ## return initial values
  list(beta0 = beta0,
       beta = beta
  )
}


# adding this so that the logProb isn't -Inf
logProb <- "-Inf"

while (logProb == "-Inf") {
  
  ### set list of initial values for NIMBLE
  inits <- initFn(consts$np)
  
  ### define the model, data, inits and constants
  model <- nimbleModel(code = model_code, constants = consts, data = data, inits = inits)
  
  logProb <- model$calculate()
  
}

### compile the model
cModel <- compileNimble(model)

### configure MCMC and monitor parameters needed for prediction
config <- configureMCMC(cModel)
config$printMonitors()


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





