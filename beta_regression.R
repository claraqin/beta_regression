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

# greta
library(greta)
library(tensorflow)

# data
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
data <- data.frame(y, x1, x2)

# variables, priors, operations
beta1 <- normal(0, 5)
beta2 <- normal(0, 5)
mu <- boot::inv.logit(beta1 * data$x1 + beta2 * data$x2)
# phi <- exponential(0.01)
phi <- normal(100, 5, truncation = c(0,Inf)) # strong prior for phi

shapes <- beta_muphi2ab(mu = mu, phi = phi)
shape1 <- shapes$a
shape2 <- shapes$b

# likelihood
distribution(y) <- beta(shape1 = shape1, shape2 = shape2)

# defining the model
m <- model(beta1, beta2, phi)

# plotting
plot(m)

# sampling
draws <- mcmc(m, n_samples = 10000,
              initial_values = initials(
                beta1 = -0.4,
                beta2 = 1,
                phi = 100
              ))

summary(draws)

# Iterations = 1:10000
# Thinning interval = 1 
# Number of chains = 4 
# Sample size per chain = 10000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#        Mean        SD  Naive SE Time-series SE
# beta1 -0.40 1.836e-05 9.180e-08      3.242e-06
# beta2  1.00 3.452e-05 1.726e-07      6.139e-06
# phi   99.99 3.726e-03 1.863e-05      7.257e-04
# 
# 2. Quantiles for each variable:
#   
#          2.5%   25%   50%   75% 97.5%
# beta1 -0.4001 -0.40 -0.40 -0.40  -0.4
# beta2  1.0000  1.00  1.00  1.00   1.0
# phi   99.9828 99.99 99.99 99.99 100.0

# It only generates accurate estimates when
# (1) the precision parameter is high,
# (2) there is a strong prior for precision, and
# (3) the initial values are set to the true values.

# Even so, the chains don't converge within
# 10,000 samples.

library(bayesplot)
pdf('figures/bayesplot1.pdf')
bayesplot::mcmc_trace(draws)
dev.off()


# Let's try a simple one-covariate model,
# and decrease the size of the data,
# with everything else kept the same.
n <- 1000
x1 <- rnorm(n)
beta1_true <- -0.4
mu_true <- boot::inv.logit( beta1_true*x1 )
set.seed(213784921)
out <- lapply(mu_true, function(x) {rbeta2(1, x, phi_true)})
y <- as_data(exclude01(unlist(out)))
x1 <- as_data(x1)
beta1 <- normal(0, 5)
mu <- boot::inv.logit(beta1 * x1)
phi <- normal(100, 5, truncation = c(0,Inf)) # strong prior for phi
shapes <- beta_muphi2ab(mu = mu, phi = phi)
shape1 <- shapes$a
shape2 <- shapes$b
distribution(y) <- beta(shape1 = shape1, shape2 = shape2)
m2 <- model(beta1, phi)
plot(m2)
draws2 <- mcmc(m2, n_samples = 1000,
              # initial_values = initials(
              #   beta1 = -0.4,
              #   phi = 100
              # )
              )
summary(draws2)

pdf('figures/bayesplot2.pdf')
bayesplot::mcmc_trace(draws2)
dev.off()

# Let's try a simple one-covariate model,
# with wider range of x1, and a non-zero
# intercept.
n <- 1000
set.seed(135237)
x1 <- rnorm(n, 3, 5)
beta0_true <- 2
beta1_true <- -0.4
mu_true <- boot::inv.logit( beta0_true + beta1_true*x1 )
phi_true <- 100
hist(mu_true)
set.seed(213784921)
out <- lapply(mu_true, function(x) {rbeta2(1, x, phi_true)})
y <- as_data(exclude01(unlist(out)))
x1 <- as_data(x1)
beta0 <- normal(0, 5)
beta1 <- normal(0, 2)
mu <- boot::inv.logit(beta0 + beta1 * x1)
phi <- normal(100, 5, truncation = c(0,Inf)) # strong prior for phi
shapes <- beta_muphi2ab(mu = mu, phi = phi)
shape1 <- shapes$a
shape2 <- shapes$b
distribution(y) <- beta(shape1 = shape1, shape2 = shape2)
m3 <- model(beta0, beta1, phi)
plot(m3)
draws3 <- mcmc(m3, n_samples = 1000,
               # initial_values = initials(
               #   beta0 = 2,
               #   beta1 = -0.4,
               #   phi = 100
               # )
               )
summary(draws3)

pdf('figures/bayesplot3.pdf')
bayesplot::mcmc_trace(draws3)
dev.off()

# Let's try an even simpler model: a beta
# regression with only an intercept
n <- 1000
set.seed(135237)
beta0_true <- 2
mu_true <- boot::inv.logit( rep(beta0_true, n) )
phi_true <- 100
hist(mu_true)
set.seed(213784921)
out <- lapply(mu_true, function(x) {rbeta2(1, x, phi_true)})
y <- as_data(exclude01(unlist(out)))
beta0 <- normal(0, 5)
mu <- boot::inv.logit(beta0)
phi <- normal(100, 5, truncation = c(0,Inf)) # strong prior for phi
shapes <- beta_muphi2ab(mu = mu, phi = phi)
shape1 <- shapes$a
shape2 <- shapes$b
distribution(y) <- beta(shape1 = shape1, shape2 = shape2)
m4 <- model(beta0, phi)
plot(m4)
draws4 <- mcmc(m4, n_samples = 1000,
               # initial_values = initials(
               #   beta0 = 2,
               #   beta1 = -0.4,
               #   phi = 100
               # )
)
summary(draws4)

# Iterations = 1:1000
# Thinning interval = 1 
# Number of chains = 4 
# Sample size per chain = 1000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean        SD  Naive SE Time-series SE
# beta0     -2.115     14.55    0.2301         0.1623
# phi   107039.206 149443.85 2362.9147      8903.0707
# 
# 2. Quantiles for each variable:
#   
#   2.5%       25%    50%       75%     97.5%
# beta0 -16.88532 -16.59434 -4.022     12.49     14.84
# phi     0.01743   0.01761 85.840 240060.94 354822.85

pdf('figures/bayesplot4.pdf')
bayesplot::mcmc_trace(draws4)
dev.off()


# Let's try a simple linear Bayesian regression:
# data
n <- 100
set.seed(543258)
x1 <- rnorm(n, 3, 5)
beta0_true <- 2
beta1_true <- -0.4
y <- beta0_true + beta1_true * x1 + rnorm(100, 0, 0.5)
x1 <- as_data(x1)
y <- as_data(y)
#variables and priors
beta0 <- normal(0,5)
beta1 <- normal(0,5)
sd <- student(3, 0, 1, truncation = c(0, Inf))
# operations
mean <- beta0 + beta1 * x1
# likelihood
distribution(y) <- normal(mean, sd)
# model
m5 <- model(beta0, beta1, sd)
plot(m5)
# MCMC samples
draws5 <- mcmc(m5, n_samples = 1000)

summary(draws5)
bayesplot::mcmc_trace(draws5)
