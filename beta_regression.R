library(glmmTMB)
library(bbmle)
library(ggplot2)
library(dplyr)

# For reparameterizing the beta distribution
# From "classical" parameterization
beta_ab2muphi <- function(a, b) {
  return(c(
    a/(a + b), 
    a + b
  ))
}

# From mean-and-dispersion parameterization
beta_muphi2ab <- function(mu, phi) {
  return(c(
    mu * phi,
    phi * (1 - mu)
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
beta_pdf(5, 2)

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


# Randomly generate data with bistable beta-distributed residuals
set.seed(1)
x1 <- rnorm(250)
x2 <- rnorm(250)
beta1 <- -0.4
beta2 <- 5
y <- boot::inv.logit( beta1*x1 + beta2*x2 )
hist(y)
# This actually already looks like it has a 
# bimodal beta error distribution, when in fact, 
# it's beta-distributed *without* errors. 
# Fitting a model to it should yield a perfect fit:
data <- data.frame(y, x1, x2)
glm1 <- glmmTMB(y ~ x1 + x2,
                data = data,
                family = poisson(link = "logit"))
summary(glm1)
# Yes, it returns the original coefficients

# Actually, inv.logit(rnorm(250, 0, x)) produces a
# seemingly bimodal distribution bounded at 0 and 1
# for any value of x < ~2.
hist(boot::inv.logit(rnorm(250, 0, 1)))
hist(boot::inv.logit(rnorm(250, 0, 2)))
hist(boot::inv.logit(rnorm(250, 0, 3)))

# Suppose we still need to add bimodal 
# beta-distributed errors. Since response values 
# need to be bounded at 0 and 1, keep generating 
# new bimodal beta-distributed errors, adding 
# them to the raw response if within bounds, or 
# rejecting them if exceeding the bounds. This
# maintains the beta distribution of the errors.
epsilon_star <- rbeta(250, 0.5, 0.5)
y_star <- y + epsilon_star - 0.5
while(sum(y < 0 | y > 1)) { # untested
  y <- ifelse(y_star < 0 | y_star > 1, y, y_star)
  epsilon_star <- rbeta(250, 0.5, 0.5)
  y_star <- y + epsilon_star
}
