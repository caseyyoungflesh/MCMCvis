# MCMCvis <a href='https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html'><img src='figs/MCMC_vis_logo_v3_bw.png' align="right" height="245" /></a>

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MCMCvis)](http://cran.r-project.org/package=MCMCvis)
![devel 0.16.5](https://img.shields.io/badge/devel-0.16.5-red.svg)
[![Downloads](http://cranlogs.r-pkg.org/badges/MCMCvis?color=blue)](http://cran.rstudio.com/package=MCMCvis)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00640/status.svg)](https://doi.org/10.21105/joss.00640)


`MCMCvis` is an R package used to visualize, manipulate, and summarize MCMC output. MCMC output may be derived from Bayesian models fit with Stan, NIMBLE, JAGS, and other software. This includes (but is not limited to) output from the `rstan`, `cmdstanr`, `rstanarm`, `brms`, `rjags`, `nimble`, `R2jags`, and `jagsUI` packages - object types are recognized automatically.

The package contains six functions:

- `MCMCsummary` - summarize MCMC output for particular parameters of interest
- `MCMCpstr` - summarize MCMC output for particular parameters of interest while preserving parameter structure
- `MCMCtrace` - create trace and density plots of MCMC chains for particular parameters of interest
- `MCMCchains` - easily extract posterior chains from MCMC output for particular parameters of interest
- `MCMCplot` - create caterpillar plots from MCMC output for particular parameters of interest
- `MCMCdiag` - create a .txt file and save specified objects that summarize model inputs, outputs, and diagnostics

`MCMCvis` was designed to perform key functions for MCMC analysis using minimal code, in order to free up time/brainpower for interpretation of analysis results. Functions support simple and straightforward subsetting of model parameters within the calls, and produce presentable and 'publication-ready' output.

**This package can be cited as:**

Youngflesh, C. (2018) MCMCvis: Tools to visualize, manipulate, and summarize MCMC output. *Journal of Open Source Software*, 3(24), 640, https://doi.org/10.21105/joss.00640

## Installation

You can install the released version on CRAN with:
```{r}
install.packages('MCMCvis')
```

Or the latest, development version from Github with:
```{r}
install.packages('remotes')
remotes::install_github('caseyyoungflesh/MCMCvis', build_vignettes = TRUE)
```

## Vignette

The vignette for this package can be run using:
```{r}
vignette('MCMCvis')
```

## Examples

#### Summarize

```{r}
data(MCMC_data)
MCMCsummary(MCMC_data, params = 'beta', round = 2)

#>           mean    sd   2.5%    50%  97.5% Rhat n.eff
#> beta[1]  -4.62  6.54 -17.19  -4.62   8.38    1 10411
#> beta[2] -14.17  6.63 -27.15 -14.08  -1.41    1 10500
#> beta[3] -35.94  8.42 -52.60 -36.00 -19.26    1 10884
#> beta[4]   6.17 10.72 -14.67   6.11  27.27    1 10500
#> beta[5]   8.42  3.46   1.63   8.45  15.13    1 10500
#> beta[6] -12.05  2.34 -16.66 -12.05  -7.54    1 10500
```

```{r}
MCMCdiag(MCMC_data, file_name = 'model_summary.txt', 
          object_name = 'model-fit.rds', save_object = TRUE,
          mkdir = 'results-YYYY-MM-DD', round = 3)

#saves a text file summarizing model inputs, outputs, and diagnostics (and optionally saves model output and other objects of interest as '.rds' files into a new user-specified directory)
```

#### Evaluate

```{r}
PR <- rnorm(15000, 0, 32)
MCMCvis::MCMCtrace(MCMC_data, params = 'beta[1]', 
                   ISB = FALSE, exact = TRUE, priors = PR, ind = TRUE,
                   Rhat = TRUE, n.eff = TRUE, pdf = FALSE)
```
![](figs/Evaluate_ex.png)


#### Manipulate

```{r}
just_betas_mcmc_obj <- MCMCchains(MCMC_data, params = 'beta', mcmc.list = TRUE)
```

#### Visualize

```{r}
MCMCplot(object = MCMC_data, object2 = MCMC_data2,
         params = 'beta', rank = TRUE, offset = 0.14, ref_ovl = TRUE)
```
![](figs/Viz_ex.png)
