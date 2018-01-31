# get.beta.dist.R
# given an ancestry estimate and standard deviation, get the beta distribution approximation
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory for Cancer Research
# SAIC-Frederick, Inc
# Created August 21, 2012
# Last Modified August 21, 2012

get.alpha <- function(p, sd)
{
                  # if proportion is estimated "near 0", alpha needs to be 1
  alpha <- ifelse(p <= 0.5 & p * (1 - p) < sd^2,
                  1,
                  
                  # if proportion is estimated "near 1", estimate alpha based on sd only (with beta == 1)
           ifelse(p > 0.5 & p * (1 - p) < sd^2,
                  NA, # pick this up later -- much faster
                  
                  # under "normal" circumstances, this will be the estimation we use
                  p * (p * (1 - p) / sd^2 - 1)))

  # pick up this bit -- much faster this way
  if(any(is.na(alpha)))
  {
    for(i in (1:length(alpha))[is.na(alpha)])
    {
      alpha[i] <- uniroot(function(a) a / ((a + 1)^2 * (a + 2)) - sd[i]^2, c(1, 1e7))$root
    }
  }

  return(alpha)
}

get.beta <- function(p, sd)
{
                 # if proportion is estimated "near 1", beta needs to be 1
  beta <- ifelse(p >= 0.5 & p * (1 - p) < sd^2,
                 1,

                 # if proportion is "near 0", estimate beta based on sd only (with alpha == 1)
          ifelse(p < 0.5 & p * (1 - p) < sd^2,
                 NA,
  
                 # under "normal" circumstances, this will be the estimation we use
                 (1 - p) * (p * (1 - p) / sd^2 - 1)))

  # pick up this bit -- much faster this way
  if(any(is.na(beta)))
  {
    for(i in (1:length(beta))[is.na(beta)])
    {
      beta[i] <- uniroot(function(b) b / ((b + 1)^2 * (b + 2)) - sd^2, c(1, 1e7))$root
    }
  }
     
  return(beta)
}
