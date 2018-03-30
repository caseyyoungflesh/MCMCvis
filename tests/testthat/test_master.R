library(MCMCvis)

context("test_master")

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


test_that('MCMCpstr displays dimensions correctly for all object types',
          {
            #mcmc.list
            expect_output(str(MCMCpstr(MCMC_data)), 'List of 3')
            expect_equal(length(MCMCpstr(MCMC_data)$alpha), 10)
            #mcmc.list - 3d
            expect_output(str(MCMCpstr(threed_data)), 'List of 1')
            expect_equal(dim(MCMCpstr(threed_data)$alpha), c(2,2,2))
            #R2jags
            expect_output(str(MCMCpstr(R2jags_data)), 'List of 4')
            expect_equal(length(MCMCpstr(R2jags_data)$alpha), 1)
            #jags.parallel
            expect_output(str(MCMCpstr(jagsparallel_data)), 'List of 4')
            expect_equal(length(MCMCpstr(jagsparallel_data)$mu), 27)
            #jagsUI
            expect_output(str(MCMCpstr(jagsUI_data)), 'List of 4')
            expect_equal(length(MCMCpstr(jagsUI_data)$mu), 27)
            #stan.fit
            expect_output(str(MCMCpstr(stan_data)), 'List of 2')
            expect_equal(length(MCMCpstr(stan_data)$mu), 10)
            #matrix
            expect_output(str(MCMCpstr(matrix_data)), 'List of 3')
            expect_equal(length(MCMCpstr(matrix_data)$alpha), 10)
            #jags.samples - expect warning
            expect_error(MCMCpstr(jagssamps_data))
          })


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


test_that('MCMCsummary values agree with manual values derived from posterior chains',
          {
            #mcmc.list - mean
            expect_equal(MCMCsummary(MCMC_data,
                                      param = 'alpha\\[1\\]',
                                      ISB = FALSE)[1],
                         round(mean(MCMCchains(MCMC_data,
                                    param = 'alpha\\[1\\]',
                                    ISB = FALSE)), 2))
            #mcmc.list - sd
            expect_equal(MCMCsummary(MCMC_data,
                                     param = 'alpha\\[1\\]',
                                     ISB = FALSE)[2],
                         round(sd(MCMCchains(MCMC_data,
                                               param = 'alpha\\[1\\]',
                                               ISB = FALSE)), 2))
            #mcmc.list - 2.5%
            expect_equal(MCMCsummary(MCMC_data,
                                     param = 'alpha\\[1\\]',
                                     ISB = FALSE)[3],
                         round(quantile(MCMCchains(MCMC_data,
                                               param = 'alpha\\[1\\]',
                                               ISB = FALSE), probs = 0.025)[[1]], 2))
            #mcmc.list - 50%
            expect_equal(MCMCsummary(MCMC_data,
                                     param = 'alpha\\[1\\]',
                                     ISB = FALSE)[4],
                         round(quantile(MCMCchains(MCMC_data,
                                                   param = 'alpha\\[1\\]',
                                                   ISB = FALSE), probs = 0.5)[[1]], 2))
            #mcmc.list - 97.5%
            expect_equal(MCMCsummary(MCMC_data,
                                     param = 'alpha\\[1\\]',
                                     ISB = FALSE)[5],
                         round(quantile(MCMCchains(MCMC_data,
                                                   param = 'alpha\\[1\\]',
                                                   ISB = FALSE), probs = 0.975)[[1]], 2))
            #mcmc.list - rhat
            expect_equal(MCMCsummary(MCMC_data,
                                     param = 'alpha\\[1\\]',
                                     ISB = FALSE)[6],
                         round(coda::gelman.diag(MCMCchains(MCMC_data,
                                               param = 'alpha\\[1\\]',
                                               ISB = FALSE, mcmc.list = TRUE))$psrf[,1], 2))
          })

