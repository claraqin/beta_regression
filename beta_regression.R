library(glmmTMB)
library(bbmle)
library(ggplot2)

# For reparameterizing the beta distribution
# "Classical" parameterization
beta_ab2muphi <- function(a, b) {
  return(c(
    a/(a+b), 
    a*b/((a+b)^2 * (a+b+1))
  ))
}

# Center and shape parameterization
beta_muphi2ab <- function(mu, phi) {
  return(c(
    mu^2 * (1 - mu) / phi - mu,
    mu / phi * (1 - mu)^2 + mu - 1
  ))
}

# Function for plotting beta pdf
beta_pdf <- function(a, b, add=FALSE) {
  grid <- seq(0,1, by=0.01)
  if(add) {
    points(x = grid,
           y = dbeta(grid, a, b))
  } else {
    plot(x = grid,
           y = dbeta(grid, a, b))
  }
}

beta_pdf(1, 1, FALSE)
beta_pdf(0.5, 0.5, TRUE)
beta_pdf(0.9, 0.9, TRUE)
beta_pdf(1, 0.9, TRUE)

# Function for plotting beta pdf with center-shape parameterization
beta_pdf2 <- function(mu, phi, add=FALSE) {
  ab_parameters <- beta_muphi2ab(mu, phi)
  a <- ab_parameters[1]
  b <- ab_parameters[2]
  grid <- seq(0,1, by=0.01)
  if(add) {
    points(x = grid,
           y = dbeta(grid, a, b))
  } else {
    plot(x = grid,
         y = dbeta(grid, a, b))
  }
}

beta_pdf2(0.5, 0.25)


# Randomly generate data with bimodal beta-distributed residuals