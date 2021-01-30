library(MCMCvis)

context("test_master")



# load data ---------------------------------------------------------------

load('../testdata/jags_data.rda')
load('../testdata/R2jags_data.rda')
load('../testdata/jagsparallel_data.rda')
load('../testdata/jagsUI_data.rda')
load('../testdata/stan_data.rda')
load('../testdata/matrix_data.rda')
load('../testdata/jagssamps_data.rda')
load('../testdata/threed_data.rda')
load('../testdata/nimble_mF_ch4.rda')



# run tests ---------------------------------------------------------------

testthat::test_that('MCMCsummary returns output for all supported object types',
          {
            #mcmc.list
            testthat::expect_equal(NROW(MCMCsummary(jags_data)), 1)
            #R2jags
            testthat::expect_equal(NROW(MCMCsummary(R2jags_data)), 2)
            #jags.parallel
            testthat::expect_equal(NROW(MCMCsummary(jagsparallel_data)), 2)
            #jagsUI
            testthat::expect_equal(NROW(MCMCsummary(jagsUI_data)), 2)
            #stan.fit
            testthat::expect_equal(NROW(MCMCsummary(stan_data)), 2)
            #matrix
            #testthat::expect_equal(NROW(MCMCsummary(matrix_data, Rhat = FALSE)), 3)
            #NIMBLE (list of matrices output)
            testthat::expect_equal(NROW(MCMCsummary(nimble_mF_ch4)), 1)
            #jags.samples - expect warning
            testthat::expect_error(MCMCsummary(jagssamps_data))
          })


testthat::test_that('MCMCpstr displays dimensions correctly for all object types',
          {
            #mcmc.list
            testthat::expect_output(str(MCMCpstr(MCMC_data)), 'List of 2')
            testthat::expect_equal(length(MCMCpstr(MCMC_data)$alpha), 6)
            #mcmc.list - 3d
            testthat::expect_output(str(MCMCpstr(threed_data)), 'List of 1')
            testthat::expect_equal(dim(MCMCpstr(threed_data)$alpha), c(2,2,2))
            #R2jags
            testthat::expect_output(str(MCMCpstr(R2jags_data)), 'List of 2')
            testthat::expect_equal(length(MCMCpstr(R2jags_data)$mu), 1)
            #jags.parallel
            testthat::expect_output(str(MCMCpstr(jagsparallel_data)), 'List of 2')
            testthat::expect_equal(length(MCMCpstr(jagsparallel_data)$mu), 1)
            #jagsUI
            testthat::expect_output(str(MCMCpstr(jagsUI_data)), 'List of 2')
            testthat::expect_equal(length(MCMCpstr(jagsUI_data)$mu), 1)
            #stan.fit
            testthat::expect_output(str(MCMCpstr(stan_data)), 'List of 2')
            testthat::expect_equal(length(MCMCpstr(stan_data)$mu), 1)
            #matrix
            testthat::expect_output(str(MCMCpstr(matrix_data)), 'List of 3')
            testthat::expect_equal(length(MCMCpstr(matrix_data)$alpha), 1)
            #NIMBLE (list of matrices output)
            testthat::expect_output(str(MCMCpstr(nimble_mF_ch4)), 'List of 1')
            testthat::expect_equal(length(MCMCpstr(nimble_mF_ch4)$mu), 1)
            #jags.samples - expect warning
            testthat::expect_error(MCMCpstr(jagssamps_data))
          })


# Insert test to make sure MCMpstr output matches MCMCsummary output (to make sure values are being placed in the correct position)

testthat::test_that('MCMCchains converts all supported object types to mcmc.list',
          {
            #mcmc.list
            testthat::expect_is(MCMCchains(MCMC_data, mcmc.list = TRUE), 'mcmc.list')
            #R2jags
            testthat::expect_is(MCMCchains(R2jags_data, mcmc.list = TRUE), 'mcmc.list')
            #jags.parallel
            testthat::expect_is(MCMCchains(jagsparallel_data, mcmc.list = TRUE), 'mcmc.list')
            #jagsUI
            testthat::expect_is(MCMCchains(jagsUI_data, mcmc.list = TRUE), 'mcmc.list')
            #stan.fit
            testthat::expect_is(MCMCchains(stan_data, mcmc.list = TRUE), 'mcmc.list')
            #matrix
            testthat::expect_error(MCMCchains(matrix_data, mcmc.list = TRUE))
            #NIMBLE (list of matrices)
            testthat::expect_is(MCMCchains(nimble_mF_ch4, mcmc.list = TRUE), 'mcmc.list')
            #jags.samples - expect warning
            testthat::expect_error(MCMCchains(jagssamps_data, mcmc.list = TRUE))
          })


testthat::test_that('MCMCsummary values agree with manual values derived from posterior chains', { 

  # mcmc.list - mean
  testthat::expect_equal(as.numeric(MCMCsummary(MCMC_data, param = 'alpha[1]', 
                                      ISB = FALSE, exact = TRUE, round = 2)[1]),
    round(mean(MCMCchains(MCMC_data, param = 'alpha[1]', ISB = FALSE, exact = TRUE)), 2))
  
  # mcmc.list - sd
  testthat::expect_equal(as.numeric(MCMCsummary(MCMC_data, param = 'alpha[1]', ISB = FALSE, 
                                      exact = TRUE, round = 2)[2]),
    round(sd(MCMCchains(MCMC_data, param = 'alpha[1]', ISB = FALSE, exact = TRUE)), 2))

  # mcmc.list - 2.5%
  testthat::expect_equal(as.numeric(MCMCsummary(MCMC_data, param = 'alpha[1]', 
                                      ISB = FALSE, exact = TRUE, round = 2)[3]),
    round(quantile(MCMCchains(MCMC_data,param = 'alpha[1]', ISB = FALSE, 
                              exact = TRUE), probs = 0.025)[[1]], 2))

  # mcmc.list - 50%
  testthat::expect_equal(as.numeric(MCMCsummary(MCMC_data, param = 'alpha[1]', 
                                      ISB = FALSE, exact = TRUE, round = 2)[4]),
    round(quantile(MCMCchains(MCMC_data,param = 'alpha[1]', 
                              ISB = FALSE, exact = TRUE), probs = 0.5)[[1]], 2))

  # mcmc.list - 97.5%
  testthat::expect_equal(as.numeric(MCMCsummary(MCMC_data, param = 'alpha[1]', 
                                      ISB = FALSE, exact = TRUE, round = 2)[5]),
    round(quantile(MCMCchains(MCMC_data,param = 'alpha[1]', 
                              ISB = FALSE, exact = TRUE), probs = 0.975)[[1]], 2))
  
  # mcmc.list - rhat
  testthat::expect_equal(as.numeric(MCMCsummary(MCMC_data, param = 'alpha[1]', 
                                      ISB = FALSE, exact = TRUE, round = 2)[6]),
    round(coda::gelman.diag(MCMCchains(MCMC_data, param = 'alpha[1]', 
                                       ISB = FALSE, exact = TRUE, mcmc.list = TRUE))$psrf[,1], 2))

  # stanfit - func = mean
  testthat::expect_equal(as.numeric(MCMCsummary(stan_data, param = 'mu', 
                                      ISB = FALSE, exact = TRUE, round = 2, func = mean)$func),
               as.numeric(MCMCsummary(stan_data, param = 'mu', 
                                      ISB = FALSE, exact = TRUE, round = 2)['mean']))
    
  })




testthat::test_that('MCMCsummary returns no errors for default and non-default specifications', { 

  # MCMC_data
  testthat::expect_error(MCMCsummary(MCMC_data), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, round = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, digits = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data, probs = c(.1, .5, .9), digits = 2), NA)
  
  # MCMC_data2
  testthat::expect_error(MCMCsummary(MCMC_data2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, round = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, digits = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(MCMC_data2, probs = c(.1, .5, .9), digits = 2), NA)

  # jags_data
  testthat::expect_error(MCMCsummary(jags_data), NA)
  testthat::expect_error(MCMCsummary(jags_data, round = 2), NA)
  testthat::expect_error(MCMCsummary(jags_data, digits = 2), NA)
  testthat::expect_error(MCMCsummary(jags_data, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(jags_data, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(jags_data, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(jags_data, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(jags_data, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(jags_data, probs = c(.1, .5, .9), digits = 2), NA)
  
  # jagsparallel_data
  testthat::expect_error(MCMCsummary(jagsparallel_data), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, round = 2), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, digits = 2), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(jagsparallel_data, probs = c(.1, .5, .9), digits = 2), NA)
  
  # jagsUI_data
  testthat::expect_error(MCMCsummary(jagsUI_data), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, round = 2), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, digits = 2), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(jagsUI_data, probs = c(.1, .5, .9), digits = 2), NA)
  
  # R2jags_data
  testthat::expect_error(MCMCsummary(R2jags_data), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, round = 2), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, digits = 2), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(R2jags_data, probs = c(.1, .5, .9), digits = 2), NA)
  
  # matrix_data
  # testthat::expect_error(MCMCsummary(matrix_data), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, round = 2), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, digits = 2), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, HPD = TRUE, prob = .9), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, HPD = TRUE, prob = .9, round = 2), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, HPD = TRUE, prob = .9, digits = 2), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, probs = c(.1, .5, .9)), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, probs = c(.1, .5, .9), round = 2), NA)
  # testthat::expect_error(MCMCsummary(matrix_data, probs = c(.1, .5, .9), digits = 2), NA)
  
  # stan_data
  testthat::expect_error(MCMCsummary(stan_data), NA)
  testthat::expect_error(MCMCsummary(stan_data, round = 2), NA)
  testthat::expect_error(MCMCsummary(stan_data, digits = 2), NA)
  testthat::expect_error(MCMCsummary(stan_data, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(stan_data, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(stan_data, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(stan_data, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(stan_data, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(stan_data, probs = c(.1, .5, .9), digits = 2), NA)
  
  # nimble_mF_ch4 - NIMBLE (list of matrices)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, round = 2), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, digits = 2), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, HPD = TRUE, prob = .9), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, HPD = TRUE, prob = .9, round = 2), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, HPD = TRUE, prob = .9, digits = 2), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, probs = c(.1, .5, .9)), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, probs = c(.1, .5, .9), round = 2), NA)
  testthat::expect_error(MCMCsummary(nimble_mF_ch4, probs = c(.1, .5, .9), digits = 2), NA)

})

# Add test to make sure colnames and rownames are correct for each object type (MCMCsummary, MCMCchains, MCMCpstr?)
