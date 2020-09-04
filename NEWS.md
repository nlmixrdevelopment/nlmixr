# Development nlmixr version

## User interface changes
 - Detection of initial conditions was rewritten to enable additional features
   in the initial conditions (#322). The most important user-facing change is
   that now arbitrary R expressions can be used when setting initial conditions
   such as `tvCL <- log(c(2,3,4))` (#253) instead of simply `tvCL <- log(3)`
 - The function as.nlmixrBounds() now supports adding the columns that are
   missing into the input data.frame.
 - omega definitions can be correlation matrices (#338)

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

## SAEM changes

 - `saem` fits now gracefully fall back to the `focei` likelihood when
   they support files are no longer on the loaded disk
 - `saem` fits now can add in `fo`, `foce` and `focei` likelihood
 - `saem` fits now use `liblsoda` by default and are multi-threaded when
   running (controlled by `RxODE`)
 - `saem` now supports time-varying covariates (like clock-time)

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

Note this includes all the RxODE changes *including* dropping python.
