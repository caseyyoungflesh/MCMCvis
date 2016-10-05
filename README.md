MCMCvis
====

`MCMCvis` is an R package used to visualize, manipulate, and summarize MCMC output. MCMC output may be derived from Bayesian model output fit with Stan, JAGS, or other MCMC samplers.

The package contains four functions:

- `MCMCsummary` - summarize MCMC output for particular parameters of interest
- `MCMCtrace` - create trace and density plots of MCMC chains for particular parameters of interest
- `MCMCchains` - easily extract posterior chains from MCMC output for particular parameters of interest
- `MCMCplot` - create caterpillar plots from MCMC output for particular parameters of interest

`MCMCvis` was designed to perform key functions for MCMC analysis using minimal code, in order to free up time/brainpower for interpretation of analysis results. Functions support simple and straightforward subsetting of model parameters within the calls, and produce presentable and 'publication-ready' output.

Installation
------------

You can install the latest version with:
```{r}
install.packages('devtools')
devtools::install_github('caseyyoungflesh/MCMCvis', build_vignettes = TRUE)
```

Vignette
--------

The vignette for this package can be run using:
```{r}
vignette('MCMCvis')
```
