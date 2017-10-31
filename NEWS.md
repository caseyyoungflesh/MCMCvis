NEWS
====

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
- `MCMCsummary` now displays estimates for deviance


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
