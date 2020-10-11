test_that("saem issue #307", {

  model_onecmt <- function(){
    ini({
      tvka <- c(-10, -2, 10) # log(Absorption rate (1/hr))
      tvcl <- c(-10, -5, 10) # log(Clearance (L/hr))
      tvv <- c(-10, -3, 10) # log(Volume of Distribution (L/kg))
      tvf_sc <- c(0, 0.5, 1) # subcutaneous bioavailability
      eta_cl ~ 0.3 # IIV in clearance
      add.sd <- c(0, 0.1, 1) # Proportional error
    })
    model({
      cl = exp(tvcl + eta_cl)
      ka = exp(tvka)
      v = exp(tvv)
      kel = cl/v
      d/dt(SC) = -ka*SC
      d/dt(C) = ka*tvf_sc*SC - kel*C
      cp = C/v
      cp ~ add(add.sd)
    })
  }

  fit <- nlmixr(model_onecmt)

  expect_equal(fit$saem.theta.trans, c(2L, 1L, 3L, 4L))

})
