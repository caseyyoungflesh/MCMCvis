library(MCMCvis)

context("test_master")




# create data -------------------------------------------------------------


library(rjags)
library(R2jags)
library(jagsUI)
library(rstan)

#create JAGS model
mf <- "
model {
for (i in 1:10)
{
  y[i] ~ dnorm(mu, 0.01);
}
mu ~ dnorm(0, 0.01)
}
"

data <- list(y = rnorm(10))
jm <- rjags::jags.model(textConnection(mf),
                 data = data,
                 n.chains = 3)
#rjags
jags_data <- rjags::coda.samples(jm,
                         variable.names = 'mu',
                         n.iter = 100)

#jags.samples
jagssamps_data <- rjags::jags.samples(jm,
                               variable.names = 'mu',
                               n.iter = 100)
#R2jags
R2jags_data <- R2jags::jags(data = data,
                    n.chains = 3,
                    model.file = textConnection(mf),
                    parameters.to.save = 'mu',
                    n.iter = 100)

#jags.parallel
jagsparallel_data <- R2jags::jags(data = data,
                    n.chains = 3,
                    model.file = textConnection(mf),
                    parameters.to.save = 'mu',
                    n.iter = 100)

#jagsUI
jagsUI_data <- jagsUI::autojags(data = data,
                                parameters.to.save = 'mu',
                                model.file = textConnection(mf),
                                n.chains = 3,
                                iter.increment = 10,
                                Rhat.limit=1.2,
                                max.iter= 100)

#Stan data

sm <- "
data {
real y[10];
}
parameters {
real mu;
}
model {
for (i in 1:10)
{
  y[i] ~ normal(mu, 10);
}
mu ~ normal(0, 10);
}
"

stan_data <- stan(model_code = sm,
                 data = data,
                 iter = 10)

#matrix data
matrix_data <- cbind(rnorm(100), rnorm(100), rnorm(100))
colnames(matrix_data) <- c('mu', 'alpha', 'beta')

#3d data
pnames <- rep(NA, 8)
id <- 0
for (i in 1:2)
{
  for (j in 1:2)
  {
    for (k in 1:2)
    {
      id <- sum(id, 1)
      pnames[id] <- c(paste0('alpha[', i, ',', j, ',', k, ']'))
    }
  }
}
nc <- 8
nr <- 10
means <- c(rnorm(nc, -10, 3))
sds <- abs(rnorm(nc, 5, 3))
threed_data <- coda::as.mcmc.list(
  lapply(1:3,
         function(i) coda::mcmc(matrix(rnorm(nc*nr, rep(means,each=nr), rep(sds, each=nr)),
                                       nrow=nr, dimnames=list(NULL,pnames)))))







# run tests ---------------------------------------------------------------



test_that('MCMCsummary returns output for all supported object types',
          {
            #mcmc.list
            expect_equal(NROW(MCMCsummary(jags_data)), 1)
            #R2jags
            expect_equal(NROW(MCMCsummary(R2jags_data)), 2)
            #jags.parallel
            expect_equal(NROW(MCMCsummary(jagsparallel_data)), 2)
            #jagsUI
            expect_equal(NROW(MCMCsummary(jagsUI_data)), 2)
            #stan.fit
            expect_equal(NROW(MCMCsummary(stan_data)), 2)
            #matrix
            expect_equal(NROW(MCMCsummary(matrix_data, Rhat = FALSE)), 3)
            #jags.samples - expect warning
            expect_error(MCMCsummary(jagssamps_data))
          })


test_that('MCMCpstr displays dimensions correctly for all object types',
          {
            #mcmc.list
            expect_output(str(MCMCpstr(MCMC_data)), 'List of 2')
            expect_equal(length(MCMCpstr(MCMC_data)$alpha), 6)
            #mcmc.list - 3d
            expect_output(str(MCMCpstr(threed_data)), 'List of 1')
            expect_equal(dim(MCMCpstr(threed_data)$alpha), c(2,2,2))
            #R2jags
            expect_output(str(MCMCpstr(R2jags_data)), 'List of 2')
            expect_equal(length(MCMCpstr(R2jags_data)$mu), 1)
            #jags.parallel
            expect_output(str(MCMCpstr(jagsparallel_data)), 'List of 2')
            expect_equal(length(MCMCpstr(jagsparallel_data)$mu), 1)
            #jagsUI
            expect_output(str(MCMCpstr(jagsUI_data)), 'List of 2')
            expect_equal(length(MCMCpstr(jagsUI_data)$mu), 1)
            #stan.fit
            expect_output(str(MCMCpstr(stan_data)), 'List of 2')
            expect_equal(length(MCMCpstr(stan_data)$mu), 1)
            #matrix
            expect_output(str(MCMCpstr(matrix_data)), 'List of 3')
            expect_equal(length(MCMCpstr(matrix_data)$alpha), 1)
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
