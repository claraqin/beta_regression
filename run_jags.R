# Bayesian beta regression using JAGS: running the model

library(rjags)
library(R2jags)

# read in simulated data
jags_data <- readRDS("jags_data.rds")

# hierarchical model structure
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

# run JAGS
system.time(
  jags_draws  <- jags.parallel(data = jags_data,
                               model.file = jags_mod,
                               parameters.to.save = jags_par,
                               n.chains = 4, n.iter = 10000,
                               n.burnin = 1000)
)

# results
print(jags_draws)

plot(jags_draws)

traceplot(jags_draws)