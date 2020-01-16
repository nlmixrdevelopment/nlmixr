# nlmixr XXX

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

Note this includes all the RxODE changes *including* dropping python.
