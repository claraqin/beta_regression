# beta_regression

Beta regression demo using `glmmTMB` and `rjags`

`glmmTMB` uses maximum likelihood to estimate beta regression parameters.

For a demo of `glmmTMB`, see `glmmTMB_demo.R`.

`rjags` (based on JAGS) uses Bayesian hierarchical modeling and MCMC sampling to estimate beta regression parameters.

For a demo of `rjags`, first run `simulate_jags_data.R` to generate data (`jags_data.rds`), and then see `run_jags.R`.