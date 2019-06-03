library(glmmTMB)
library(bbmle)
library(ggplot2)

beta_ab2muphi <- function(a, b) {
  return(c(
    a/(a+b), 
    a*b/((a+b)^2 * (a+b+1))
  ))
}

beta_muphi2ab <- function(mu, phi) {
  return(c(
    mu^2 * (1 - mu) / phi - mu,
    mu / phi * (1 - mu)^2 + mu - 1
  ))
}

