potools
====

`potools` is an R package used to manipulate and summarize MCMC output. MCMC output may be derived from Bayesian model ouput fit with JAGS or other MCMC samplers.

The package currently contains four functions:

- `posummary` - summarize MCMC output for particular parameters of interest
- `potrace` - plot MCMC chains for particular parameters of interest to check for convergence
- `pochains` - easily extract posterior chains from MCMC output for particular parameters of interest
- `poplot` - plot posterior distributions from MCMC output for particular parameters of interest

While several packages currently exist to evaluate MCMC output, they do not support simple subsetting of model parameters and/or do not provide an option to properly visualize MCMC output, particularly with large numbers of parameters. `potools` was designed to perform key functions for MCMC analysis using minimal code, in order to free up time/brainpower for interpretation of analysis results. 

Installation
------------

You can install the latest version with:
```{r}
install.packages('devtools')
devtools::install_github('caseyyoungflesh/potools', build_vignettes = TRUE)
```

Vignette
--------

The vignette for this package can be run using:
```{r}
vignette('potools')
```
