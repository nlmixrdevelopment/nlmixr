nlmixrTest(
{
  ## This is for long running examples
  ## These are not run by CRAN but should be checked
  for (t in c("dynmodel", "gnlmm", "foceiFit", "addCwres", "addNpde", "covarSearchAuto", "bootstrapFit", "saem_fit")) {
    test_that(sprintf("%s test", s),
              expect_error(example(t), NA))
  }
}, test="examples")
