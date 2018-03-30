library(MCMCvis)

context("test_master")

## MCMCsummary
test_that('MCMCsummary returns output for all supported object types',
          {
            #mcmc.list
            expect_equal(NROW(MCMCsummary(MCMC_data)), 30)
            #R2jags
            expect_equal(NROW(MCMCsummary(R2jags_data)), 30)
            #jags.parallel
            expect_equal(NROW(MCMCsummary(jagsparallel_data)), 30)
            #jagsUI
            expect_equal(NROW(MCMCsummary(jagsUI_data)), 30)
            #stan.fit
            expect_equal(NROW(MCMCsummary(stan_data)), 11)
            #matrix
            expect_equal(NROW(MCMCsummary(matrix_data, Rhat = FALSE)), 30)
            #jags.samples - expect warning
            expect_error(MCMCsummary(jagssamps_data))
          })


## MCMCpstr
#threed_data
test_that('MCMCpstr displays dimensions correctly for all object types',
          {
            #mcmc.list
            expect_output(str(MCMCpstr(MCMC_data)), 'List of 3')
            expect_output(str(MCMCpstr(threed_data)), 'List of 1')
            expect_equal(dim(MCMCpstr(threed_data)$alpha), c(2,2,2))
            #R2jags
            # expect_equal(NROW(MCMCsummary(R2jags_data)), 30)
            # #jags.parallel
            # expect_equal(NROW(MCMCsummary(jagsparallel_data)), 30)
            # #jagsUI
            # expect_equal(NROW(MCMCsummary(jagsUI_data)), 30)
            # #stan.fit
            # expect_equal(NROW(MCMCsummary(stan_data)), 11)
            # #matrix
            # expect_equal(NROW(MCMCsummary(matrix_data, Rhat = FALSE)), 30)
            # #jags.samples - expect warning
            # expect_error(MCMCsummary(jagssamps_data))
          })


## MCMCchains
test_that('MCMCchains converts all supported object types to mcmc.list',
          {
            #mcmc.list
            expect_is(MCMCchains(MCMC_data, mcmc.list = TRUE), 'mcmc.list')
            #R2jags
            expect_is(MCMCchains(R2jags_data, mcmc.list = TRUE), 'mcmc.list')
            #jags.parallel
            expect_is(MCMCchains(jagsparallel_data, mcmc.list = TRUE), 'mcmc.list')
            #jagsUI
            expect_is(MCMCchains(jagsUI_data, mcmc.list = TRUE), 'mcmc.list')
            #stan.fit
            expect_is(MCMCchains(stan_data, mcmc.list = TRUE), 'mcmc.list')
            #matrix
            expect_error(MCMCchains(matrix_data, mcmc.list = TRUE))
            #jags.samples - expect warning
            expect_error(MCMCchains(jagssamps_data, mcmc.list = TRUE))
          })

