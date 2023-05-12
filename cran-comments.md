CRAN comments
====

## Test environments
* local macOS install - R 4.2.1
* Ubuntu 20.04.1 (on rhub) - R-release
* Fedora (on rhub) - R-devel
* Windows Server 2022 (on rhub) - R-devel


## R CMD check results

There is one ERROR on Fedora via rhub, due to a [known issue](https://github.com/r-hub/rhub/issues/540) with `V8`, which is required for an `rstan` install (which is in `MCMCvis` Imports).

There are no WARNINGs.

There are currently 4 NOTES, all of which are expected.
#1: all systems
#2,3: only on Windows Server 2022, R-devel
#4: only Ubuntu Linux 20.04.1

1) The `cmdstanr` package is in Suggests and is not on CRAN:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Casey Youngflesh <caseyyoungflesh@gmail.com>'

Suggests or Enhances not in mainstream repositories:
     cmdstanr
Availability using Additional_repositories specification:
     cmdstanr   yes   https://mc-stan.org/r-packages/
```

2) A [known issue](https://github.com/r-hub/rhub/issues/503) with Rhub due to a bug/crash in MiKTeX:

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```

3) Similar to issue #2, only on Windows Server 2022.
```
* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  ''NULL''
```

4) A [known issue](https://github.com/r-hub/rhub/issues/548) with Rhub:

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
```


## Downstream dependencies

I have run R CMD check on downstream dependencies of MCMCvis. All packages passed.
