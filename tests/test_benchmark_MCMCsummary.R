# test MCMCsummary using data objects in tests/testdata folder

testVis <- function(input_data, input_data2) {

  load(file = paste0("tests/testdata/", input_data, ".rda"), envir = .GlobalEnv)

  print(MCMCsummary(input_data2))
  print(MCMCsummary(input_data2, round = 2))
  print(MCMCsummary(input_data2, digits = 2))
  print(MCMCsummary(input_data2, HPD = TRUE, probs = .9))
  print(MCMCsummary(input_data2, HPD = TRUE, probs = .9, round = 2))
  print(MCMCsummary(input_data2, HPD = TRUE, probs = .9, digits = 2))
  print(MCMCsummary(input_data2, probs = c(.1, .5, .9)))
  print(MCMCsummary(input_data2, probs = c(.1, .5, .9), round = 2))
  print(MCMCsummary(input_data2, probs = c(.1, .5, .9), digits = 2))
  
}

testVis("MCMC_data", MCMC_data)
testVis("MCMC_data2", MCMC_data2)
testVis("jags_data", jags_data)
testVis("jagsparallel_data", jagsparallel_data)
testVis("jagsUI_data", jagsUI_data)
testVis("R2jags_data", R2jags_data)
testVis("matrix_data", matrix_data)
testVis("stan_data", stan_data)

# benchmark MCMCsummary and compare to master branch runtimes

library(microbenchmark)
library(ggplot2)

set.seed(1)

mbm <- microbenchmark(
  "jags large default" = {MCMCsummary(MCMC_data)}, 
  "jags large HPD" = {MCMCsummary(MCMC_data, HPD = TRUE, probs = .9) },
  "jags large custom ET" = {MCMCsummary(MCMC_data,  probs = c(.1, .5, .9))},
  "jags small default" = {MCMCsummary(jags_data)}, 
  "jags small HPD" = {MCMCsummary(jags_data, HPD = TRUE, probs = .9) },
  "jags small custom ET" = {MCMCsummary(jags_data,  probs = c(.1, .5, .9))},
  "stan default" = {MCMCsummary(stan_data)},
  "stan HPD" = {MCMCsummary(stan_data, HPD = TRUE, probs = .9, round = 2)},
  "stan custome ET" = {MCMCsummary(stan_data, probs = c(.1, .5, .9), digits = 2)})

autoplot(mbm)

