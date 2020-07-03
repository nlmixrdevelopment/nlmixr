# nlmixr 1.1.1-9+

 - Detection of initial conditions was rewritten to enable additional features
   in the initial conditions (#322). The most important user-facing change is
   that now arbitrary R expressions can be used when setting initial conditions
   such as tvCL <- log(3) (#253).
 - The function as.nlmixrBounds() now supports adding the columns that are
   missing into the input data.frame.

# Before nlmixr 1.1.1-9

## `focei` changes:
 - Uses RxODE to re-arrange the problem so it does not include
   `if/else` in the model (ie. un-branched code). This allows
   sensitivities to be calculated in one pass saving time for multiple
   endpoint models and models with `if/else` in them.

- `linCmt()` now uses solved systems instead of translating to ODEs.
  - Uses `RxODE`/`stan`'s math headers to calculate the sensitivities
    of the super-positioned `linCmt()` solutions.
  - Since this currently uses the super-positioning, this does not
    support time-varying covariates.

## General Changes
 - Added auto-completion of `nlmixr` object properties accessed by
   `$`. This works for major editors including `Rstudio`, `ESS`, and
   Base R itself.

 - Changed the way that Rstudio notebooks display `nlmixr` objects; It
   should be more legible in Rstudio.

Note this includes all the RxODE changes *including* dropping python.
