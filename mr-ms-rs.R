# Author: Jeffrey W. Doser
# Code: analysis code to fit the data (mr-ms-rs-data.R) with the 
#       JAGS model (mr-ms-rs-jags.txt) in R

# Empty memory
# rm(list = ls())

# Working Directory, Data, Packages ---------------------------------------

# Set working directory ---------------
# setwd()

# Read in data ------------------------
load("mr-ms-rs-data.R")
  # Files: y: data, separated by time period of detection
         # eta: data, not separated by time period of detection
         # year.covs, percentForest, regen, basalArea: abundance covariates
         # day, time: detection covariates
         # guilds: guild membership for each species
         # J: number of sites in each region
         # max.J: max(J)
         # n: total number of species
         # n.regions: number of species observed in each park
         # n.years: number of years
         # R: number of parks
         # r.sp: species observed at each park
         # sp: list of all species observed across network
         # T: number of time intervals
         # year.by.regions: number of years monitored at each park
         # year.regions: specific years observed at each park

  # Parks: (1) ACAD (2) MABI (3) MIMA (4) MORR (5) ROVA 
         # (6) SAGA (7) SAGA (8) WEFA

# Install and load packages -----------
if(!require(jagsUI)) {install.packages("jagsUI");require(jagsUI)}

# Prep data for use in Jags -----------------------------------------------

# Bundle the data for Bayesian analysis using jags
bugs.data = list(n = n.regions, R = R, J = J, T = T, y = y, 
         	 years = year.by.regions, 
                 year.covs = year.covs, day = day, time = time, 
                 regen = regen, percentForest = percentForest, 
                 basalArea = basalArea)

# Parameters to monitor ---------------
parameters <- c('bp.1', 'mu.lambda', 'beta.1.lambda', 'beta.2.lambda', 
		'beta.3.lambda', 'beta.4.lambda', 'beta.1.lambda.r', 
		'beta.2.lambda.r', 'beta.3.lambda.r', 'beta.4.lambda.r', 
                'beta.5.lambda', 'beta.5.lambda.r', 'mu.lambda.r', 
		'mu.lambda.sp', 'beta.1.lambda.sp', 'beta.2.lambda.sp',
	       	'beta.3.lambda.sp', 'beta.4.lambda.sp', 'beta.5.lambda.sp', 
                'mu.p', 'beta.1.p', 'beta.2.p', 'beta.3.p', 'beta.4.p', 
                'beta.1.p.r', 'beta.2.p.r', 'beta.3.p.r', 'mu.p.r', 
                'beta.1.p.sp', 'beta.2.p.sp', 'beta.3.p.sp', 'mu.p.sp', 'N'
)

# Initial Values ----------------------
inits <- function() {
  list(mean.p = runif(1), beta.1.p = rnorm(1), beta.2.p = rnorm(1),
       beta.3.p = rnorm(1), tau.p.int = runif(1, 0.5, 1), 
       tau.p.1 = runif(1, 0.5, 1), tau.p.2 = runif(1, 0.5, 1), 
       tau.p.3 = runif(1, 0.5, 1), tau.p.4 = runif(1, 0.5, 1), 
       mu.lambda = rnorm(1), beta.1.lambda = rnorm(1), 
       tau.int.lambda = runif(1, 0.5, 1), tau.1.lambda = runif(1, 0.5, 1),
       beta.2.lambda = rnorm(1), tau.2.lambda = runif(1, 0.5, 1), 
       beta.3.lambda = rnorm(1), tau.3.lambda = runif(1, 0.5, 1),
       beta.4.lambda = rnorm(1), tau.4.lambda = runif(1, 0.5, 1),
       beta.5.lambda = rnorm(1), tau.5.lambda = runif(1, 0.5, 1)
  )
}

# Set the MCMC settings ---------------
# Takes a while
# n.burn <- 45000
# n.iter <- 50000
# n.thin <- 2
# n.chains <- 3
# Small run
n.burn <- 100
n.iter <- 500
n.thin <- 1
n.chains <- 1

# Fit the model -----------------------------------------------------------
out <- jags(bugs.data, inits, parameters, 'mr-ms-rs-jags.txt', 
            n.chains = n.chains, n.thin = n.thin, n.burnin = n.burn,
            n.iter = n.iter, parallel = TRUE)

# Save the jags fit -------------------
date <- Sys.Date()
# Potentially a better name for the file results
#file.name <- paste("jagsFit-", "mr-ms-rs-", n.iter, "-iter-", date, ".R", sep = "")
file.name <- 'mr-ms-rs-results.R'
save(out, file = file.name)
