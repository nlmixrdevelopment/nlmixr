nlmixrTest(
{
  ## This is for long running examples
  ## These are not run by CRAN but should be checked
  for (t in c("dynmodel", "gnlmm", "foceiFit", "addCwres", "addNpde", "covarSearchAuto", "bootstrapFit", "nlmixr", "dynmodel.mcmc",
              "nmDocx", "configsaem")) {
    test_that(sprintf("%s test", t),
              eval(parse(text=paste0("expect_error(example(", t, ", ask=FALSE), NA)"))))
  }
}, test="examples")
