# Add values

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

Added values to non-data documentation; The following documentation
files still do not have \value statements since they document data or
the re-exported functions in nlmixr:


 [1] "Bolus_1CPT.Rd"    "Bolus_1CPTMM.Rd"  "Bolus_2CPT.Rd"    "Bolus_2CPTMM.Rd" 
 [5] "Infusion_1CPT.Rd" "invgaussian.Rd"   "metabolite.Rd"    "Oral_1CPT.Rd"    
 [9] "pheno_sd.Rd"      "pump.Rd"          "rats.Rd"          "reexports.Rd"    
[13] "theo_md.Rd"       "theo_sd.Rd"       "Wang2007.Rd"      "warfarin.Rd"     


# options(warn=-1) not allowed 

You are setting options(warn=-1) in your function. This is not allowed.
Please rather use suppressWarnings() if really needed.

Fixed, options(warn=-1) is no longer in the code base
