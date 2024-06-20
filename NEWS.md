NEWS
====

# 0.16.4

- `MCMCdiag` now accepts `params`, `excl`, `ISB`, and `exact` arguments to summarize and save out a subset of parameters (useful for model runs with many params where object may exceed memory limit)


# 0.16.3

- `MCMCdiag` fix spacing issue for min neff in summary file for `cmdstanr` objects
- `MCMCdiag` change behavior so that if specified `dir` does exist, save to working dir
- `MCMCdiag` if `mkdir` name exists, append '_1' to filename and create


# 0.16.2

- Fix CRAN documentation issue


# 0.16.1

- Fix bug in `MCMCsummary` that prevented `n.eff` from being displayed for `jagsUI` objects


# 0.16.0

- All functions now accept model objects fit with `cmdstanr` package
- Add Table of Contents to vignette


# 0.15.5

- `MCMCplot` fix bug associated with `ref_ovl = TRUE` and the specification of multiple objects. Function also now automatically lightens colors specified (instead of defaulting to gray) when `ref_ovl = TRUE`


# 0.15.4

- `MCMCsummary` fix bug associated with `pg0 = TRUE` when input was derived from Stan


# 0.15.3

- `MCMCsummary` now takes argument `pg0`. When `pg0 = TRUE` the proportion of the posterior that is greater than 0 is returned as a column
- `MCMCdiag` now takes any argument from `MCMCsummary` to modify printed summary output
- Clean up documentation


# 0.15.2:

- `MCMCdiag` add more descriptive message about output from `nimble` objects
- `MCMCtrace` fix rdev matrix warning


# 0.15.1:

- `MCMCplot` fix bug that prevented plotting with a single parameter


# 0.15.0:

- `MCMCdiag` function added to create .txt file summarizing model inputs and outputs and save model object as .rds file
- All functions now take `exact` as an argument to parse desired parameter. `ISB = TRUE` no longer uses regex matching (`exact` should be used to specify this).
- Add support for objects output from the `nimble` package
- `MCMCpstr` fix bug where the incorrect number of dimensions would be returned when subsetting a param with >1 dim
- `MCMCplot` now takes the argument `ci` which can be used to specify the credible interval displayed in the caterpillar plot
- `MCMCplot` now takes the argument `HPD` which can be used to specify whether highest posterior density intervals or equal-tailed intervals should be displayed in the caterpillar plot


# 0.14.3:

- `MCMCplot` fix bug where error would result when trying to plot a single parameter


# 0.14.2:

- `MCMCplot` can now plot models that have NA for posterior estimates. Useful when plotting two models side-by-side is desired and one model has parameters that the other does not (use MCMCchains to create matrix of draws, add NA-filled cols with missing parameter names, and plot)


# 0.14.1:

- `MCMCtrace` fix undesired behavior when `post_zm` is specified as `FALSE`


# 0.14.0:

- `MCMCsummary` for model objects fit with `jagsUI`, the function now returns Rhat and n.eff values that coincide with values calculated in that package (e.g., split-chain Rhat), rather than calculating these with the `coda` package. These values are fed into `MCMCtrace` when `Rhat` and/or `n.eff` are specified as `TRUE`.
- `MCMCsummary` fix bug where an error was returned when `func_name` was not specified


# 0.13.6:

- `MCMCtrace` now accepts expression input (for Greek characters etc.) from R for `main_den` and `main_tr` args


# 0.13.5:

- All functions now avoid using function `class` (and use `is` instead) as `matrix` objects will also be `array` as of R 4.0.0


# 0.13.4:

- `MCMCtrace` no longer changes working dir when `wd` is specified - just puts pdf there


# 0.13.3:

- `MCMCpstr` fix bug where dimensions were not displayed properly for Stan models
- `MCMCpstr` fix warning with objects derived from `rstanarm` and `brms`
- `MCMCsummary` change variable names to reflect what is displayed by `brms` package (rather than raw Stan output)


# 0.13.2:

- `MCMCsummary` fix bug thrown when one only one chain run for Stan models
- `MCMCsummary` fix Rhat label issue (from r_hat to Rhat)


# 0.13.1:

- `MCMCsummary` fix bug that caused stanfit parameter names to be numbers


# 0.13.0:

- `MCMCsummary` now takes `quantiles` as an argument, allowing user to specify which quantiles to return
- `MCMCsummary` now takes `HPD` as an argument to calculate highest posterior density intervals
- `MCMCsummary` now outputs as a data.frame rather than a matrix
- `MCMCchains` fix bug that output different parameters names for `stanreg` and `brms` objects
- `MCMCplot` fix bug that caused guidelines to plot over appropriate bounds


# 0.12.6:

- `MCMCplot` now has `object2`, `col2`, and `offset` arguments. Allows user to plot output from two separate models, as long as both model outputs have identical parameter names.


# 0.12.5:

- `MCMCplot` now has `guide_lines` argument, to plot lines to help reference which posterior corresponds to each parameter name
- Remove deprecated arguments


# 0.12.4:

- Add support for objects output from `rstanarm` and `brms` packages
- `MCMCplot` add support for plotting multiple colors on the same plot


# 0.12.3:

- `MCMCsummary` for `n.eff` the default is now `TRUE` (displays number of effective samples by default)
- `MCMCtrace` now has the option `plot`. When `FALSE` no plot is output. Used in conjunction with `PPO_out = TRUE` to to calculate PPO without plotting trace plots.


# 0.12.2:

- `MCMCsummary` fix bug where Stan input was not being sorted by parameter index


# 0.12.1:

- `MCMCtrace` fix spacing issue for Rhat and n.eff text when `Rhat = TRUE` and `n.eff = TRUE`
- `MCMCplot` add deprecation warnings for arguments


# 0.12.0:

- `MCMCsummary` for `stanfit` objects (model output derived from Stan) Rhat and n_eff are calculated using the `rstan` package. Note that `rstan` calculates Rhat and n_eff values slightly differently (more conservatively) than the `coda` package (commonly used to summarize model output derived from JAGS).
- `MCMCplot` the following argument names have been changed: 
    orig: `labels_sz`, new: `sz_labels`
    orig: `med_sz`, new: `sz_med`
    orig: `thick_sz`, new: `sz_thick`
    orig: `thin_sz`, new: `sz_thin`
    orig: `axis_text_sz`, new: `sz_ax_txt`
    orig: `ax_sz`, new: `sz_ax`
    orig: `tick_text_sz`, new: `sz_tick_txt`
    orig: `main_text_sz`, new: `sz_main_txt`
    orig: `tick_pos`, new: `pos_tick`
- `MCMCtrace` add ability to specify xlim and ylim for density plots
- `MCMCtrace` add ability to specify xlab and ylab for density and trace plots
- `MCMCtrace` add ability to specify title for trace and density plot
- `MCMCtrace` add ability to specify line width and line type for density and prior lines on density plots
- `MCMCtrace` add ability to specify color for density and priors lines on density plots
- `MCMCtrace` add ability to specify size and color of text when priors specified, and the position of ticks for density and trace plots
- `MCMCtrace` add ability to specify size of tick labels, axes labels, and thickness of axes
- `MCMCtrace` clean up plotting of trace plots when only two plots are plotted in window
- `MCMCtrace` add ability to include Rhat and number of effective samples on trace plots


# 0.11.3:

- `MCMCchains` fix error associated with coda::mcmc and `rjags` objects


# 0.11.2:

- `MCMCplot` now takes `guide_axis` argument. If `TRUE`, a second axis (x-axis if `HORIZ = TRUE`, y-axis if `HORIZ = FALSE`) is plotted to help interpret values on plot.


# 0.11.1:

- `MCMCchains` now assigns arbitrary names to columns (parameters) when input type is matrix (along with a warning that it is doing so)
- `MCMCtrace` now modifies the layout of trace plots when < 6 plots are generated
- `MCMCpstr` fix bug that prevented scalars when type = 'chains'


# 0.11.0:

- `MCMCsummary` `digits` argument is now NULL by default (all computed digits are returned upon default - any rounding must be explicitly specified)
- `MCMCsummary` `digits` argument uses `signif` rather than `round` for rounding (in other words, `digits` specifies number of significant digits rather than number of decimal places)
- `MCMCsummary` now takes `round` argument to round output to specified number of decimal places
- `MCMCpstr` no longer has the option to restrict the number of digits output (returns all digits)
- `MCMCpstr` now takes `type` as an argument. When `type = 'summary'` (default), values calculated with the `func` argument (default `mean`) are returned. When `type = 'chains'`, posterior chains are returned while preserving the parameter structure. Posterior samples are stored in the last dimension of the array for each element of the output list. In this way vector parameters are output in matrix format, matrix parameters are output in three dimension array format (within each element of the output list - one parameter for each list element).
- `MCMCpstr` now accepts output greater than length 1 from argument `func`. If output is greater than length 1, function output are stored in the last dimension of the array (within each element of the output list - one parameter for each list element).


# 0.10.4:

- `MCMCchains` now takes `chain_num` as an argument. When specified, single posterior chains can be output for a particular parameter of interest. Useful for determining the last value in an MCMC chain for each parameter (to be used as initial values for a subsequent model run).


# 0.10.3:

- All functions now return a warning for missing params in`excl` and `params` arguments instead of an error. This means that output will be returned even when values specified do not exist. Change was made to `MCMCchains` code but impacts all functions.
- Fix typos
- Add contributors Che-Castaldo and Hardy


# 0.10.2:

- `MCMCtrace` now takes `post_zm` as an argument. When `post_zm = FALSE`, x- and y-limits of density plots are scaled so that both the prior and posterior can be visualized on a plot (rather than zoomed on the posterior).


# 0.10.1:

- `MCMCpstr` fix bug to institute rounding when parameters as scalars or vectors
- `MCMCtrace` now takes `PPO_out` as an argument. When `PPO_out = TRUE`, the percent overlap between prior and posterior for each paras will be returned as a data.frame.
- Add CITATION file


# 0.10.0:

- Fix warning when feeding jags.parallel object to `MCMCsummary`
- Fix bug that produced an error when using `stanfit` objects
- Add automated tests to check package functions
- `MCMCtrace` now takes `open_pdf` as an argument. When `open_pdf = FALSE`, the generated pdf will not be opened in a viewer automatically.
- `MCMCtrace` now takes `gvals` as an argument. When simulated data are used to fit a model, the generating values used to simulate the data (true parameter values) will be plotted as vertical lines on the density plots.
- `MCMCplot` `ref_ovl` argument now defaults to FALSE (one color is plotting for all parameter estimates)
- Change `MCMC_data` (example data) so that it's smaller (only 5k iterations and two parameters)


# 0.9.4:

- Fix white space issue in `MCMCplot` when many parameters are plotted and large fig dimensions are used
- Fix label alignment issue in `MCMCplot` when `horiz = FALSE` and large numbers of parameters are plotted


# 0.9.3:

- Fix bug that prevented parameters from being sorted when using matrix input for `MCMCtrace`
- Add support for objects produced with the jagsUI package


# 0.9.2:

- `MCMCtrace` now takes matrix input (as with the other functions). One chain is assumed when matrix input is used.


# 0.9.1:

- Fix bug that produced errors when using the `jags.parallel` function in the `R2jags` package.
- All functions - when `ISB = FALSE`, `params` argument now takes the form of regular expressions
- Examples for `MCMCtrace` no longer open up external programs (pdf viewer) per CRAN policy


# 0.9.0:

- `MCMCpstr` function now added. Function returns summary output for a specified function while preserving structure of parameters (i.e., scalar, vector, matrix, array).
- `MCMCtrace` now takes a `priors` argument to visualize prior/posterior overlap. If specified, the prior (user specified as this information is not contained within the MCMC output) for a specified parameter is plotted on the same plot as the posterior output. Percent overlap between posterior and prior is also calculated and displayed.
- Fix bug in `MCMCchains` that caused incorrect alphabetization of parameter names when output from R2jags was used.


# 0.8.2:

- `MCMCsummary` greatly speed up calculation of Rhat values for objects with large numbers of parameters
- `MCMCchains` now takes the argument `mcmc.list`. If specified, `mcmc.list` object returned rather than a matrix.


# 0.8.1:

- Fix bug in `MCMCsummary` that displayed the same result twice when selecting only a single output parameter of interest
- Fix bug in `MCMCplot` that displayed the axis label too close to tick labels when `horiz = FALSE` and tick labels were very long
- `MCMCsummary` Rhat values always round to 2 digits
- `MCMCsummary` output from `func` argument rounded to specified digits for rest of output
- `MCMCsummary` now takes a `func_name` argument. If specified, column displaying output from `func` will be labeled with this name. If not specified, column will be labeled 'func'.
- `MCMCplot` change argument `x_axis_text_sz` and `x_tick_text_sz` to `axis_text_sz` and `tick_text_sz` respectively


# 0.8.0:

- Specification of parameters of interest now works slightly differently. The argument `ISB` (Ignore Square Bracket) has now been added. By default `ISB = TRUE` - `params` and `excl` match exactly to parameter names by default (ignoring square brackets). When `ISB = FALSE`, square brackets will not be ignored, and will match on partial names (as when using `grep`). This applies to all functions.
- `MCMCsummary` now takes a `func` argument. If a function is specified, it will be evaluated for all specified parameters and specified in the `MCMCsummary` output.
- `MCMCsummary` speed greatly increased. Parameters of interest are now sorted before calculations are made. Rhat values are no longer masked, but rather not calculated when `Rhat = FALSE`. These changes result in dramatic speed ups for large objects.
- `MCMCsummary` bug fixed that caused function to fail when only one chain was run
- `MCMCsummary` standard deviation added to summary output for each parameter.
- `MCMCsummary` number of effective samples added to summary output for each parameter. Default is `n.eff = FALSE` (metric will not be calculated or displayed).
- `MCMCtrace` default is now to write trace plots to pdf. Default number of iterations changed to 5000 from 2000.
- `MCMCplot` y-axis labels now vertical when `horiz = FALSE` to improve readability.
- `MCMCplot` bug that resulted in poor plot dimension choices in some circumstances now fixed.
- `MCMC_data` now contains three chains with 6000 iterations each.
- Error message now added about functions not taking objects produced from `jags.samples` function in the `coda` package. `coda.samples` should be used instead.


# 0.7.1:

- Fix bug in `MCMCplot` which incorrectly shaded parameter estimates when plotted vertically
- `MCMCsummary` now displays estimates for deviance with MCMC output fits with R2jags


# 0.7.0:

- Fix bug in `MCMCsummary` to do with Cholesky decomposition and calculating Rhat.
- Speed up processing of `MCMCsummary` for certain object types.
- Fix minor documentation errors in help files for several functions.
- `MCMCplot` labels now start at top and go down (more intuitive).
- `MCMCtrace` now plots only the last 2000 iterations of the posterior chains by default. The argument is now number of iterations to be plotted, rather than the starting iteration to plot. As such, `MCMCtrace` argument `iter_st` (start iteration) changed to `iter` (number of iterations from end).
- Add `horiz` argument to `MCMCplot` - caterpillar plots can now be plotted to run vertically rather than horizontally. Parameters are plotted left to right when plotted vertically.
- Remove extended lines on axes for `MCMCplot` - axes lines only goes to the end of ticks now.


# 0.6.3:

- Initial release
