library(glmmTMB)
# library(bbmle)
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

# From mean-and-precision parameterization
beta_muphi2ab <- function(mu, phi) {
  return(c(
    mu * phi,
    phi * (1 - mu)
  ))
}

# For generating random beta-distributed 
# observations with mean-precision
# parameterization
rbeta2 <- function(n, mu, phi) {
  params <- beta_muphi2ab(mu, phi)
  return(rbeta(n, params[1], params[2]))
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

beta_pdf2(0.5, 2)
beta_pdf2(0.5, 1, TRUE)
beta_pdf2(0.5, 0.5, TRUE)

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
                family = poisson(link = "logit"))
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
set.seed(1)
phi <- 1.5 # arbitrarily assume precision = 1.5
out <- lapply(mu, function(x) {rbeta2(1, x, phi)}) 
y <- unlist(out)

# Now attempt to fit model to the new data
data <- data.frame(y, x1, x2)
glm_beta <- glmmTMB(y ~ x1 + x2,
                    data = data,
                    family = beta_family(link = "logit"))
summary(glm_beta)
# Unable to return original coefficients,
# though at least the signs and relative
# magnitudes are consistent

# What if data was generated from a beta 
# distribution with a higher (but still < 2)
# precision parameter? 
set.seed(1)
phi_precise <- 1.9
out <- lapply(mu, function(x) {rbeta2(1, x, phi_precise)})
y_precise <- unlist(out)
data <- data.frame(y_precise, x1, x2)
glm_beta_precise <- glmmTMB(y_precise ~ x1 + x2,
                            data = data, 
                            family = beta_family(link = "logit"))
summary(glm_beta_precise)
# Still unable to return original coefficients, 
# though AIC is lower (better) than before.

# What if there was more data?
set.seed(1)
x1 <- rnorm(10000)
x2 <- rnorm(10000)
beta1 <- -0.4
beta2 <- 1
mu <- boot::inv.logit( beta1*x1 + beta2*x2 )
phi <- 100
set.seed(1)
out <- lapply(mu, function(x) {rbeta2(1, x, phi)})
y <- unlist(out)
data <- data.frame(y, x1, x2)
# Need to transform y because its values
# are inclusive of 0 and 1. We need it to
# exclude those endpoints.
data$y <- (data$y * (nrow(data) - 1) + 0.5) / nrow(data)
glm_beta <- glmmTMB(y ~ x1 + x2,
                    data = data,
                    family = beta_family(link = "logit"))
summary(glm_beta)
# Estimates:
# Overdispersion parameter for beta family (): 99.7 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.0004018  0.0022100     0.2    0.856    
# x1          -0.3998361  0.0022493  -177.8   <2e-16 ***
#   x2           1.0022905  0.0026594   376.9   <2e-16 ***

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
y <- unlist(out)
data <- data.frame(y)
glm_beta <- glmmTMB(y ~ 1,
                    data = data,
                    family = beta_family(link = "logit"))
summary(glm_beta)
# Estimate:
# beta0 = 0.4984 (0.0122)

# After switching to the reparameterized 
# rbeta (rbeta2), this model returns the
# original parameter, mu = 0.5!


# What if we try to use Bayesian modeling?
# greta
library(greta)
library(tensorflow)

# data (use same data as before)
set.seed(987453)
x1 <- rnorm(10000)
x2 <- rnorm(10000)
beta1_true <- -0.4
beta2_true <- 1
mu_true <- boot::inv.logit( beta1_true*x1 + beta2_true*x2 )
phi_true <- 100 # much higher than before
set.seed(1)
out <- lapply(mu_true, function(x) {rbeta2(1, x, phi_true)})
y <- unlist(out)
data <- data.frame(y, x1, x2)
# Need to transform y because its values
# are inclusive of 0 and 1. We need it to
# exclude those endpoints.
data$y <- (data$y * (nrow(data) - 1) + 0.5) / nrow(data)

# variables, priors, operations
beta1 <- normal(0, 5)
beta2 <- normal(0, 5)
mu <- boot::inv.logit(beta1 * data$x1 + beta2 * data$x2)
phi <- exponential(0.01)

shapes <- beta_muphi2ab(mu = mu, phi = phi)
shape1 <- shapes[1]
shape2 <- shapes[2]

# likelihood
distribution(y) <- beta(shape1 = shape1, shape2 = shape2)

# defining the model
m <- model(beta1, beta2, phi)

# plotting
plot(m)

# sampling
draws <- mcmc(m, n_samples = 1000,
              # initial_values = initials(
              #   beta1 = -0.4,
              #   beta2 = 1,
              #   precision = 100
              # )
              )

summary(draws)

# Iterations = 1:1000
# Thinning interval = 1 
# Number of chains = 4 
# Sample size per chain = 1000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD Naive SE Time-series SE
# beta1 -0.04366 0.92405 0.014610      3.147e-07
# beta2 -0.09874 0.03168 0.000501      1.029e-08
# phi    1.02519 0.10254 0.001621      8.530e-08
# 
# 2. Quantiles for each variable:
#   
#   2.5%      25%     50%      75%    97.5%
# beta1 -1.6436 -0.06669  0.4817  0.50471  0.50564
# beta2 -0.1340 -0.12763 -0.1007 -0.07176 -0.05964
# phi    0.8861  0.97982  1.0198  1.06515  1.17513

library(bayesplot)
pdf('figures/bayesplot1.pdf')
bayesplot::mcmc_trace(draws)
dev.off()


# Try more informative priors

set.seed(1)
out <- lapply(mu, function(x) {rbeta2(1, x, phi)}) # arbitrarily assume precision = 1
y <- unlist(out)
y <- as_data(y)

mean <- normal(0.5, 0.1, truncation = c(0,1))
precision <- normal(1.5, 0.1, truncation = c(0, Inf))

shapes <- beta_muphi2ab(mu = mean, phi = precision)
shape1 <- shapes[1]
shape2 <- shapes[2]

distribution(y) <- beta(shape1 = shape1, shape2 = shape2)

m <- model(mean, precision)

draws_informative <- mcmc(m, n_samples = 1000)

summary(draws_informative)

# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean       SD  Naive SE Time-series SE
# mean      0.2956 0.002458 3.887e-05      5.668e-05
# precision 2.0700 0.024781 3.918e-04      7.414e-04
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5%
# mean      0.2908 0.2939 0.2955 0.2972 0.3004
# precision 2.0208 2.0535 2.0700 2.0866 2.1186
