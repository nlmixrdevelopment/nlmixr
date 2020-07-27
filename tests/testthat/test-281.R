nlmixrTest({
  context("Covs in SAEM covs and model, Issue #281")
  test_that("saem building works with ", {
    skip_on_os("solaris") # gcc compiler has to be available and working

    Lesion7 <- function() {
      ini({
        temaxD <- -0.05 # typical value of drug emax
        tec50 <- -0.04 # typical value of ec50
        temaxT <- -0.9 # typical value of emax for placebo effect as a function of time
        tet50 <- 6.64 # typical value of et50
        tkin <- -7.8
        tslope <- -1 # typical value of growth parameter
        TRXslope <- -0.1
        tdelay <- 7.15

        eta.emaxD ~ 0.3
        eta.emaxT ~ 0.3
        eta.slope ~ 0.3
        eta.delay ~ 0.3
        add.err <- .1 # add. residual variability
      })
      model({
        Resp(0) <- 1 # default Resp(0) = 0
        emaxD <- exp(temaxD + eta.emaxD)
        ec50 <- exp(tec50)
        emaxT <- exp(temaxT + eta.emaxT)
        et50 <- exp(tet50)

        slope <- exp(tslope + TRX * TRXslope + eta.slope)

        delay <- exp(tdelay + eta.delay)
        kin <- exp(tkin)
        kout <- kin
        GAM <- exp(1) / 2 # Hill coef
        C2 <- centr / V
        CONC <- (C2)^2

        Stim1 <- emaxT * (time) / (time + et50)
        Stim2 <- emaxD * (CONC^GAM) / (CONC^GAM + ec50^GAM)
        Stim <- Stim1 * (1 + TRX * Stim2)

        Delta <- 1 / (1 + exp(-20 * (time - delay)))
        d / dt(depot) <- -KA * depot
        d / dt(centr) <- KA * depot - CL * C2
        d / dt(Resp) <- kin * (1 + Delta * slope) - kout * (1 + Stim) * Resp

        Resp ~ add(add.err)
      })
    }


    tmp <- nlmixr(Lesion7)
    expect_error(tmp$saem.model, NA)
  })

}, test="saem")
