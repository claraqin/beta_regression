library(glmmTMB)
# library(bbmle)
library(ggplot2)
library(dplyr)

# For reparameterizing the beta distribution
# From "classical" parameterization
beta_ab2muphi <- function(a, b) {
  return(list(
    mu = a/(a + b), 
    phi = a + b
  ))
}

# From mean-and-precision parameterization
beta_muphi2ab <- function(mu, phi) {
  return(list(
    a = mu * phi,
    b = phi * (1 - mu)
  ))
}

# For generating random beta-distributed 
# observations with mean-precision
# parameterization
rbeta2 <- function(n, mu, phi) {
  if(length(mu) > 1 | length(phi) > 1) {
    stop("Error: mu and phi parameters must",
         " be constants, not vectors.")
  }
  params <- beta_muphi2ab(mu, phi)
  return(rbeta(n, params$a[1], params$b[1]))
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

# Function for plotting beta pdf with 
# mean-precision parameterization
beta_pdf2 <- function(mu, phi, add=FALSE) {
  ab_parameters <- beta_muphi2ab(mu, phi)
  a <- ab_parameters$a
  b <- ab_parameters$b
  grid <- seq(0,1, by=0.01)
  if(add) {
    points(x = grid,
           y = dbeta(grid, a, b))
  } else {
    plot(x = grid,
         y = dbeta(grid, a, b))
  }
}

beta_pdf2(0.5, 2)
beta_pdf2(0.5, 1, TRUE)
beta_pdf2(0.5, 0.5, TRUE)

# Need to transform 0-1 scale because its values
# are inclusive of 0 and 1. We need it to
# exclude those endpoints, in order for logit
# link to avoid Inf and -Inf.
exclude01 <- function(x) {
  n <- length(x)
  return((x * (n - 1) + 0.5) / n)
}

# Randomly generate data with bistable beta-distributed residuals

# We assume the following model:
#    Y ~ Beta(mu, phi)
#    logit(mu) = XB
# In the first step, we generate random X,
# specify coefficients B, and calculate
# the vector mu.
set.seed(1)
x1 <- rnorm(250)
x2 <- rnorm(250)
beta1 <- -0.4
beta2 <- 5
mu <- boot::inv.logit( beta1*x1 + beta2*x2 )
hist(mu)
# This actually already looks like it has a 
# bimodal beta error distribution, when in fact, 
# it's beta-distributed *without* errors. 
# Fitting a model to it should yield a perfect fit:
data <- data.frame(mu, x1, x2)
glm1 <- glmmTMB(mu ~ x1 + x2,
                data = data,
                family = beta_family(link = "logit"))
summary(glm1)
# Yes, it returns the original coefficients

# Actually, inv.logit(rnorm(250, 0, x)) produces a
# seemingly bimodal distribution bounded at 0 and 1
# for any value of x < ~2.
hist(boot::inv.logit(rnorm(250, 0, 1)))
hist(boot::inv.logit(rnorm(250, 0, 2)))
hist(boot::inv.logit(rnorm(250, 0, 3)))

# Now to simulate response values with bimodal
# beta-distributed errors, assume mu_i to be the
# mean of the beta distribution from which y_i
# was drawn.
set.seed(1001)
phi <- 1.5 # arbitrarily assume precision = 1.5
out <- lapply(mu, function(x) {rbeta2(1, x, phi)}) 
y <- exclude01(unlist(out))

# Now attempt to fit model to the new data
data <- data.frame(y, x1, x2)
glm_beta <- glmmTMB(y ~ x1 + x2,
                    data = data,
                    family = beta_family(link = "logit"))
summary(glm_beta)
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.00576    0.07133  -0.081   0.9356    
# x1          -0.12210    0.07420  -1.645   0.0999 .  
# x2           1.46665    0.11360  12.911   <2e-16 ***

# Unable to return original coefficients,
# though at least the signs and relative
# magnitudes are consistent.

# What if data was generated from a beta 
# distribution with a higher (but still < 2)
# precision parameter?
set.seed(1489321)
phi_precise <- 1.9
out <- lapply(mu, function(x) {rbeta2(1, x, phi_precise)})
y_precise <- exclude01(unlist(out))
data <- data.frame(y_precise, x1, x2)
glm_beta_precise <- glmmTMB(y_precise ~ x1 + x2,
                            data = data, 
                            family = beta_family(link = "logit"))
summary(glm_beta_precise)
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.04590    0.06979   0.658    0.511    
# x1          -0.08107    0.07209  -1.124    0.261    
# x2           1.58275    0.12130  13.049   <2e-16 ***

# Still unable to return original coefficients, 
# though AIC is lower (better) than before.

# What if there was more data?
set.seed(89051)
x1 <- rnorm(10000)
x2 <- rnorm(10000)
beta1 <- -0.4
beta2 <- 1
mu <- boot::inv.logit( beta1*x1 + beta2*x2 )
phi <- 100
set.seed(1)
out <- lapply(mu, function(x) {rbeta2(1, x, phi)})
y <- exclude01(unlist(out))
data <- data.frame(y, x1, x2)
glm_beta <- glmmTMB(y ~ x1 + x2,
                    data = data,
                    family = beta_family(link = "logit"))
summary(glm_beta)
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.001569   0.002227     0.7    0.481    
# x1          -0.400665   0.002305  -173.8   <2e-16 ***
# x2           0.999744   0.002648   377.5   <2e-16 ***

# (Note: The model retrieves the original 
# parameters only if the precision is high.
# However, any precision > 2 would not
# represent a bistable (0,1) distribution.)

# What if it was a much simpler model?
# Try an intercept-only model.
beta0 <- 0.5
mu <- boot::inv.logit( beta0 * rep(1, 10000) )
phi <- 1.5
set.seed(1)
out <- lapply(mu, function(x) {rbeta2(1, x, phi)})
y <- exclude01(unlist(out))
data <- data.frame(y)
glm_beta <- glmmTMB(y ~ 1,
                    data = data,
                    family = beta_family(link = "logit"))
summary(glm_beta)
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.49598    0.01218   40.71   <2e-16 ***

# This model returns the
# original parameter, mu = 0.5!


# What if we try to use Bayesian modeling?
# (and change the coefficients and increase
# the precision parameter?)

# R2jags version
library(rjags)
library(R2jags)
n <- 10000
set.seed(987453)
x1 <- rnorm(n)
x2 <- rnorm(n)
beta1_true <- -0.4
beta2_true <- 1
mu_true <- boot::inv.logit( beta1_true*x1 + beta2_true*x2 )
phi_true <- 100 # much higher than before
set.seed(1)
out <- lapply(mu_true, function(x) {rbeta2(1, x, phi_true)})
y <- exclude01(unlist(out))
jags_data <- list(
  "y" = y,
  "x1" = x1,
  "x2" = x2,
  "N" = length(y)
)
jags_mod <- function() {
  # likelihood
  for(i in 1:N) {
    y[i] ~ dbeta(shape1[i], shape2[i])
    shape1[i] <- phi * mu[i]
    shape2[i] <- phi * (1 - mu[i])
    odds[i] <- exp(beta1 * x1[i] + beta2 * x2[i])
    mu[i] <- odds[i] / (1 + odds[i])
  }
  
  # priors
  beta1 ~ dnorm(0,5)
  beta2 ~ dnorm(0,5)
  phi ~ dnorm(100, 5)
}

# parameters to monitor
jags_par <- c("beta1", "beta2", "phi")

jags_draws  <- jags.parallel(data = jags_data,
                             model.file = jags_mod,
                             parameters.to.save = jags_par,
                             n.chains = 4, n.iter = 10000)


