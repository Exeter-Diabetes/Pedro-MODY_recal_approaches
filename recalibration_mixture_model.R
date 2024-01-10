#:--------------------------------------------------------
#   This file contains the recalibration mixture model
#
#:--------------------------------------------------------

library(nimble)
library(tidyverse)


### Load control-case dataset
# case_control_dataset <- readRDS("case_control.rds")

### Load representative dataset
# representative_dataset <- readRDS("representative_dataset.rds") %>%
#   mutate(T = ifelse(C == 0 | A == 1, 1, 0)) # T is 1 if Cn or Ap


###### Recalibration Mixture model

## calculate the prior distribution
pM <- rbeta(100000, 2 + 1, 247 - 2 + 1)
pT_M <- rbeta(100000, 6 + 1, 282 - 6 + 1)
p <- pM * pT_M

## fit beta distribution to this
fn <- function(pars, psamp) {
  alpha <- pars[1]
  beta <- pars[2]
  if(alpha < 0 | beta < 0) return(NA)
  sum(dbeta(psamp, alpha, beta, log = TRUE))
}
pbeta <- optim(c(1, 1), fn, control = list(fnscale = -1), psamp = p)


### Define model structure
model_code <- nimbleCode({
  
  ## case-control likelihood component
  for (i in 1:nCC) {
    
    ## regression mean
    logit(pCC[i]) <- beta0 + inprod(beta[1:np], xCC[i, 1:np])
    
    ## likelihood
    MCC[i] ~ dbern(pCC[i])
  }
  
  ## real world likelihood
  for (j in 1:nind) { # iterate through patients from real world data
    
    ## prediction from UNITED
    Z[j] <- beta0 + inprod(beta[1:np], x[j, 1:np])
    
    ## regression of real world data - shrinkage applied
    logit(p_shrinkage[j]) <- gamma0 + gamma1 * Z[j]
    
    ## MODY modelling
    M[j] ~ dbern(p[j]) # probability of MODY
    p[j] <- T[j] * pMp_Cn_or_Ap + ( 1 - T[j] ) * p_shrinkage[j] # decomposition of MODY based on T
    T[j] ~ dbern(pT[j]) # explaining how T should be defined
    
    logit(pT[j]) <- beta_t0 + inprod(beta_t[1:np_T], xT[j, 1:np_T]) + inprod(beta_spline[1:np_spline], x_spline[j, 1:np_spline])
  }
  
  ## testing priors
  pMp_Cn_or_Ap ~ dbeta(pMp_Cn_or_Ap_alpha, pMp_Cn_or_Ap_beta)       # probability of MODY given C- or A+ is < 0.01
  
  beta0 ~ dnorm(0, sd = 10)
  for(j in 1:np) {
    beta[j] ~ dnorm(0, sd = 10)
  }
  
  beta_t0 ~ dnorm(0, sd = 10)
  for(j in 1:np_T) {
    beta_t[j] ~ dnorm(0, sd = 10)
  }
  for(j in 1:np_spline) {
    beta_spline[j] ~ dnorm(0, sd = 10)
  }
  
  ## regression priors
  gamma0 ~ dnorm(0, sd = 10)
  gamma1 ~ dnorm(1, sd = 10)
  
})

### Set up data for model
#:--- Case control covariates
xCC <- as.matrix(select(case_control_dataset, pardm, agerec, hba1c, agedx, sex))
MCC <- as.numeric(as.character(case_control_dataset$mody))

#:--- Representative covariates
x <- as.matrix(select(representative_dataset, pardm, agerec, hba1c, agedx, sex))
M <- representative_dataset$mody
T <- representative_dataset$T
xT <- as.matrix(select(representative_dataset, bmi, agedx, pardm, agerec))
x_spline <- representative_dataset %>%
  select(bmi, agedx, agerec) %>%
  mutate(bmi_spline = as.numeric(rms::rcs(bmi, 3)[, 2]),
         agedx_spline = as.numeric(rms::rcs(agedx, 3)[, 2]),
         agerec_spline = as.numeric(rms::rcs(agerec, 3)[, 2])) %>%
  select(-bmi, -agedx, -agerec) %>%
  set_names(c("bmi_spline", "agedx_spline", "agerec_spline")) %>%
  as.matrix()

# necessary for predictions after model being fitted
rcs_parms <- data.frame(
  bmi = attr(rms::rcs(representative_dataset$bmi, 3)[, 2], "parms"),
  agedx = attr(rms::rcs(representative_dataset$agedx, 3)[, 2], "parms"),
  agerec = attr(rms::rcs(representative_dataset$agerec, 3)[, 2], "parms")
)


### set up data list for NIMBLE
data <- list(M = M,
             MCC = MCC,
             T = T)

### set up other components of model
consts <- list(nind = length(M),             # number of individuals in UNITED
               nCC = length(MCC),            # number of individuals in Case control
               np = ncol(x),                 # number of covariates
               np_T = ncol(xT),                 # number of covariates T
               np_spline = ncol(x_spline),   # number of covariates splines
               xCC = xCC,                    # covariate matrix for Case control
               x = x,                        # covariate matrix for UNITED
               xT = xT,                        # covariate matrix for UNITED T
               x_spline = x_spline,          # covariate matrix for UNITED splines
               pMp_Cn_or_Ap_alpha = pbeta$par[1],
               pMp_Cn_or_Ap_beta = pbeta$par[2])

### set up initial values for parameters in the model
initFn <- function(M, x, T, np, np_T, np_spline, pbeta) {
  
  ## simulate from priors
  gamma0 <- rnorm(1)
  gamma1 <- rnorm(1)
  pMp_Cn_or_Ap <- rbeta(1, pbeta[1], pbeta[2])
  beta0 <- rnorm(1)
  beta <- rnorm(np)
  beta_t0 <- rnorm(1)
  beta_t <- rnorm(np_T)
  beta_spline <- rnorm(np_spline)
  
  ## simulate missing information for MODY
  Z <- t(beta %*% t(x))
  M1 <- gamma0 + gamma1 * Z
  M1 <- exp(M1) / (1 + exp(M1))
  M1 <- rbinom(length(M1), 1, M1)
  M1[!is.na(M)] <- NA
  
  # simulate latent missing variable
  T1 <- T
  T1[is.na(T)] <- 1
  T1[!is.na(T)] <- NA
  
  ## return initial values
  list(gamma0 = gamma0,
       gamma1 = gamma1,
       beta0 = beta0,
       beta = beta,
       beta_t0 = beta_t0,
       beta_t = beta_t,
       beta_spline = beta_spline,
       pMp_Cn_or_Ap = pMp_Cn_or_Ap,
       T = T1,
       M = M1
  )
}


# adding this so that the logProb isn't -Inf
logProb = "-Inf"

while (logProb == "-Inf") {
  
  ### set list of initial values for NIMBLE
  inits <- initFn(data$M, consts$x, data$T, consts$np, consts$np_T, consts$np_spline, pbeta$par)
  
  ### define the model, data, inits and constants
  model <- nimbleModel(code = model_code, constants = consts, data = data, inits = inits)
  
  logProb = model$calculate()
  
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
               inits = inits,
               progressBar = TRUE,
               summary = TRUE,
               samplesAsCodaMCMC = TRUE,
               thin = 1)
