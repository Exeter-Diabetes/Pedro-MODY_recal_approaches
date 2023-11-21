#:--------------------------------------------------------
#   This file contains the recalibration model
#
#:--------------------------------------------------------

library(nimble)
library(tidyverse)


### Load control-case dataset
# case_control_dataset <- readRDS("case_control.rds")

### Load representative dataset
# representative_dataset <- readRDS("representative_dataset.rds")


### Fit Bayesian Recalibration nimble model
model_code <- nimbleCode({
  
  ## case-control likelihood component
  for (i in 1:nCC) {
    
    ## regression mean
    logit(pCC[i]) <- beta0 + inprod(beta[1:np], xCC[i, 1:np])
    
    ## likelihood
    MCC[i] ~ dbern(pCC[i])
  }
  
  ## likelihood component
  for (j in 1:nind) {
    
    ## prediction from UNITED
    Z[j] <- beta0 + inprod(beta[1:np], x[j, 1:np])
    
    ## regression mean - with shrinkage applied
    logit(p[j]) <- gamma0 + gamma1 * Z[j]
    
    ## likelihood
    M[j] ~ dbern(p[j])
    
  }
  
  ## regression priors
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  for(j in 1:np) {
    beta[j] ~ dnorm(0, 0.01)
  }
  
})

#:--- Case control covariates
xCC <- as.matrix(select(case_control_dataset, pardm, agerec, hba1c, agedx, sex))
MCC <- as.numeric(as.character(case_control_dataset$mody))

#:--- Representative covariates
x <- as.matrix(select(representative_dataset, pardm, agerec, hba1c, agedx, sex))
M <- representative_dataset$mody

### set up data list for NIMBLE
data <- list(M = M, MCC = MCC)

### set up other components of model
consts <- list(nind = length(M),             # number of individuals in UNITED
               nCC = length(MCC),            # number of individuals in Case control
               np = ncol(x),                 # number of covariates
               xCC = xCC,                    # covariate matrix for Case control
               x = x)                        # covariate matrix for UNITED\

### set up initial values for parameters in the model
initFn <- function(M, x, np) {
  
  ## simulate from priors
  gamma0 <- rnorm(1)
  gamma1 <- rnorm(1)
  beta0 <- rnorm(1)
  beta <- rnorm(np)
  
  ## simulate missing information
  Z <- t(beta %*% t(x))
  M1 <- gamma0 + gamma1 * Z
  M1 <- exp(M1) / (1 + exp(M1))
  M1 <- rbinom(length(M1), 1, M1)
  
  M1[!is.na(M)] <- NA
  
  ## return initial values
  list(gamma0 = gamma0,
       gamma1 = gamma1,
       beta0 = beta0,
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

