NEWS
====

0.7.0:

*Fix bug in MCMCsummary to do with cholesky decomposition and calculating rhat.
*Speed up processing of MCMCsummary for certain object types.
*Fix documentation errors in function help.
*MCMCplot labels now start at top and go down (more intuitive).
*MCMCtrace now plots only the last 2000 iterations of the posterior chains by default. The argument is now number of iterations to be plotted, rather than the starting iteration to plot. As such, MCMCtrace argument 'iter_st' (start iteration) changed to 'iter' (number of interations from end).
*Add 'horiz' argument to MCMCplot - caterpillar plots can now be plotted to run vertically rather than horizontaly. Parameters are plotted left to right when plotted vertically.
*Remove extended lines on axes for MCMCplot - axes lines only goes to the end of ticks now.



0.6.3:

Initial release