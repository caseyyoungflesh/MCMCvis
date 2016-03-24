post
====

**post** is an R package used to summarize posterior MCMC chains (including parameter estimate extraction, chains extraction, and visualizing posterior distributions), obtained from Bayesian model output.


The package currently contains three functions:

- posummary
- pochains
- poplot

Workflow
--------

These functions are designed to be used in a workflow.

1) `posummary` - the summary information is extracted from the jags object

2) `pochains` - the posterior chains of interest are extracted from the jags object

3) `poplot` - the posterior chains of interest are plotted using density bars


```{r}
posummary(jags_object, params)
chains <- pochains(jags_object, params)
poplot(chains)
```

Installation
------------

You can install the latest version with:
```{r}
install.packages('devtools')
devtools::install_github('caseyyoungflesh/post')
```
