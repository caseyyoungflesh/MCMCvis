potools
====

**potools** is an R package used to manipulate and summarize MCMC output. MCMC output may be derived from Bayesian model ouput fit with JAGS or other MCMC samplers.

The package currently contains four functions:

- `posummary` - summarize MCMC output for particular parameters of interest
- `potrace` - plot MCMC chains for particular parameters of interest to check for convergence
- `pochains` - easily extract posterior chains from MCMC output for particular parameters of interest
- `poplot` - plot posterior distributions from MCMC output for particular parameters of interest


Installation
------------

You can install the latest version with:
```{r}
install.packages('devtools')
devtools::install_github('caseyyoungflesh/potools')
```
