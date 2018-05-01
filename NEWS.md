NEWS
====


0.10.2:

- `MCMCtrace` now takes `post_zm` as an argument. When `post_zm = FALSE`, x- and y-limits of density plots are scaled so that both the prior and posterior can be visualized on a plot (rather than zoomed on the posterior).


0.10.1:

- `MCMCpstr` fix bug to institute rounding when parameters as scalars or vectors
- `MCMCtrace` now takes `PPO_out` as an argmuent. When `PPO_out = TRUE`, the percent overlap between prior and posterior for each paras will be returned as a dataframe.
- Add CITATION file


0.10.0:

- Fix warning when feeding jags.parallel object to `MCMCsummary`
- Fix bug that produced an error when using `stanfit` objects
- Add automated tests to check package functions
- `MCMCtrace` now takes `open_pdf` as an argument. When `open_pdf = FALSE`, the generated pdf will not be opened in a viewer automatically.
- `MCMCtrace` now takes `gvals` as an argument. When simulated data are used to fit a model, the generating values used to simulate the data (true parameter values) will be plotted as vertical lines on the density plots.
- `MCMCplot` `ref_ovl` argument now defaults to FALSE (one color is plotting for all parameter estimates)
- Change `MCMC_data` (exmaple data) so that it's smaller (only 5k iterations and two parameters)


0.9.4:

- Fix whitespace issue in `MCMCplot` when many parameters are plotted and large fig dimensions are used
- Fix label alginment issue in `MCMCplot` when `horiz = FALSE` and large numbers of parameters are plotted


0.9.3:

- Fix bug that prevented parameters from being sorted when using matrix input for `MCMCtrace`
- Add support for objects produced with the jagsUI package


0.9.2:

- `MCMCtrace` now takes matrix input (as with the other functions). One chain is assumed when matrix input is used.


0.9.1:

- Fix bug that produced errors when using the `jags.parallel` function in the `R2jags` package.
- All functions - when `ISB = FALSE`, `params` argument now takes the form of regular expressions
- Examples for `MCMCtrace` no longer open up external programs (pdf viewer) per CRAN policy


0.9.0:

- `MCMCpstr` function now added. Function returns summary output for a specified function while preserving structure of parameters (i.e., scalar, vector, matrix, array).
- `MCMCtrace` now takes a `priors` argument to visualize prior/posterior overlap. If specified, the prior (user specified as this information is not contained within the MCMC output) for a specified parameter is plotted on the same plot as the posterior output. Percent overlap between posterior and prior is also calculated and displayed.
- Fix bug in `MCMCchains` that caused incorrect alphabetization of parameter names when output from R2jags was used.


0.8.2:

- `MCMCsummary` greatly speed up calculation of Rhat values for objects with large numbers of parameters
- `MCMCchains` now takes the argument `mcmc.list`. If specified, `mcmc.list` object returned rather than a matrix.


0.8.1:

- Fix bug in `MCMCsummary` that displayed the same result twice when selecting only a single output parameter of interest
- Fix bug in `MCMCplot` that displayed the axis label too close to tick labels when `horiz = FALSE` and tick labels were very long
- `MCMCsummary` Rhat values always round to 2 digits
- `MCMCsummary` output from `func` argument rounded to specified digits for rest of output
- `MCMCsummary` now takes a `func_name` argument. If specified, column displaying output from `func` will be labeled with this name. If not specified, column will be labeled 'func'.
- `MCMCplot` change argument `x_axis_text_sz` and `x_tick_text_sz` to `axis_text_sz` and `tick_text_sz` respectively


0.8.0:

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


0.7.1:

- Fix bug in `MCMCplot` which incorrectly shaded parameter estimates when plotted vertically
- `MCMCsummary` now displays estimates for deviance with MCMC output fits with R2jags


0.7.0:

- Fix bug in `MCMCsummary` to do with Cholesky decomposition and calculating Rhat.
- Speed up processing of `MCMCsummary` for certain object types.
- Fix minor documentation errors in help files for several functions.
- `MCMCplot` labels now start at top and go down (more intuitive).
- `MCMCtrace` now plots only the last 2000 iterations of the posterior chains by default. The argument is now number of iterations to be plotted, rather than the starting iteration to plot. As such, `MCMCtrace` argument `iter_st` (start iteration) changed to `iter` (number of iterations from end).
- Add `horiz` argument to `MCMCplot` - caterpillar plots can now be plotted to run vertically rather than horizontally. Parameters are plotted left to right when plotted vertically.
- Remove extended lines on axes for `MCMCplot` - axes lines only goes to the end of ticks now.


0.6.3:

- Initial release
