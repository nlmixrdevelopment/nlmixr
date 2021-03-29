# nlmixr 2.0.3

## User interface changes
 - Detection of initial conditions was rewritten to enable additional features
   in the initial conditions (#322). The most important user-facing change is
   that now arbitrary R expressions can be used when setting initial conditions
   such as `tvCL <- log(c(2,3,4))` (#253) instead of simply `tvCL <- log(3)`

 - The function as.nlmixrBounds() now supports adding the columns that are
   missing into the input data.frame.

 - omega definitions can be correlation matrices (#338)

 - Can specify `keep=` and `drop=` in the nlmixr function to keep and
   drop columns in nlmixr output.  Can also specify
   `control=list(keep=,drop=)` or `nlmixr(...,keep=,drop=)` to
   keep/drop columns (#260)

## `focei` changes:
 - Uses RxODE to re-arrange the problem so it does not include
   `if/else` in the model (ie. un-branched code). This allows
   sensitivities to be calculated in one pass saving time for multiple
   endpoint models and models with `if/else` in them.

- `linCmt()` now uses solved systems instead of translating to ODEs.
  - Uses `RxODE`/`stan`'s math headers to calculate the sensitivities
    of the super-positioned `linCmt()` solutions.
  - This uses the `advan` solutions and hence supports
    support time-varying covariates.

- `focei` now supports censoring in the same way `monolix` does, with
  `cens` and `limit` columns

- `focei` now allows `eta`s on dose-related modeled events like
  `alag`, `f`, etc by finite difference sensitivities.

- `focei` now supports 2 combined additive + proportional error
  models;
  - `combined1`: `trans(y) = trans(f) + (a+b*f^c)*err`
  - `combined2`: `trans(y) = trans(f) + sqrt(a^2+b^2*f^(2c))*err`

- `focei` `etaNudge` parameters were changed to use quadrature points
  covering 95% percent of a standard normal.

- With zero gradients, Gill differences are recomputed to try to find
  a non-zero gradient.

- Now when running if a zero gradient is detected, reset the problem
  (theta reset) and re-estimated with `outerOpt="bobyqa"`

- Now when running a model where the last objective function is not
  the minimum objective function, issue a warning and skip the
  covariance step. (See Issue #403)

- `focei` proportional and power models are more tolerant of 0
  predictions in your data


## SAEM changes

 - `saem` fits now gracefully fall back to the `focei` likelihood when
   they support files are no longer on the loaded disk

 - `saem` phi pile is now saved in the `RxODE::rxTempDir()` which can
   be customized to allow the `phi` file to remain after R has exited

 - `saem` fits now can add in `fo`, `foce` and `focei` likelihood

 - `saem` fits now use `liblsoda` by default and are multi-threaded when
   running (controlled by `RxODE`)

 - `saem` now supports time-varying covariates (like clock-time)

 - `saem` now supports 2 combined additive + proportional error models:
    - `combined1`: `trans(y) = trans(f) + (a+b*f^c)*err`
	- `combined2`: `trans(y) = trans(f) + sqrt(a^2+b^2*f^(2c))*err`

 - `saem` proportional and power models are more tolerant of 0
    predictions in your data

 - `saem` now supports censoring a similar way as `monolix` does, with
  `cens` and `limit` columns

 - Added the low-discrepancy "vandercorput" normal sequence from
   SITMO/`RxODE`; In theory this will cover the domain quickly and evenly
   possibly allowing a better solution.  This can be called with
   `saemControl(normal="vandercorput")`.

 - The default of `saem` additive + proportional error has been
   switched to `combined2`, which was the `focei` default, but you can
   change this back with `saemControl(addProp="combined2")`.  The
   table results will likely be different because in the last release
   the `saem` calculated `combined1` and then used these coefficients
   in the `combined2` focei problem.

## nlme changes

- `nlme` will now support 2 combined additive + proportional error models (if the patched version of nlme is used)
    - `combined1`: `y = f + (a+b*f)*err`
	- `combined2`: `y = f + sqrt(a^2+b^2*f^2)*err`
	- See https://github.com/nlmixrdevelopment/nlmixr/issues/428
	- Thanks to Johannes Ranke (@jranke) for the nlme patch and the catch

- Can switch with `nlmeControl(addProp="combined1")` to use the combined1 type of error model

## New Utilities

 - `bootstrapFit` now calculates the bootstrap confidence bands and
   (optionally) will compare with the theoretical chi-squared
   distribution to help assess their adequacy.

## General Changes

 - Added auto-completion of `nlmixr` object properties accessed by
   `$`. This works for major editors including `Rstudio`, `ESS`, and
   Base R itself.

 - Changed the way that Rstudio notebooks display `nlmixr` objects; It
   should be more legible in Rstudio.

 - Graphics have been revamped to show censoring (including adding
   ggplot stat/geom `geom_cens`) as well as use `RxODE`'s ggplot theme
   (`rxTheme()`).  Additionally time after dose is calculated as `tad`
   for all `nlmixr` models

 - Tables generation has been refactored; `npde` uses the `arma` and
   `RxODE` random number generators which may change results.  Also
   the default of `ties=TRUE` has been changed to `ties=FALSE`.
   `npde` calculations have been threaded with `OpenMP` to speed up
   the calculation as well.  This refactoring was required to have the
   `dv` imputation between `cwres` and `npde` use the same method.
   The `npde` option now calculates the decorrelated `npd` as well, (which is
   the recommended weighted residual; see Nguyen 2017)

## Bug Fixes

 - Aligned `saem` and `focei` additive + proportional error models, so
   `saem` `additive+proportional` outputs will be different using the
   correct `focei` method

Note this includes all the RxODE changes *including* dropping python.
