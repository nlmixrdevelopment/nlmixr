## 


* Add values

Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags in up to 42 .Rd files, e.g.:
      addCovMultiple.Rd: \value
      as.focei.dynmodel.Rd: \value
      bootstrapFit.Rd: \value
      configsaem.Rd: \value
      covarSearchAuto.Rd: \value
      dot-nmGetData.Rd: \value
      ...

Added values to non-data documentation
  

You are setting options(warn=-1) in your function. This is not allowed.
Please rather use suppressWarnings() if really needed.



## Test environments
* local R installation, R 4.0.4
* ubuntu 16.04 (on travis-ci), R 4.0.4
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
