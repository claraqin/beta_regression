# Bayesian beta regression using JAGS: data simulation

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

saveRDS(jags_data, file="jags_data.rds")
