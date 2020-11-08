nlmixrTest(
{
  test_that("focei 429", {

    pk.turnover.emax.lag <- function() {
      ini({
        talag <- c(0,0.1)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(10)
        ##
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        prop.err <- 0.1
        pkadd.err <- 0.1
        ##
        poplogit <- 2
        #temax <- 7.5
        tec50 <- log(0.5)
        tkout <- log(0.05)
        te0 <- log(100)
        ##
        eta.emax ~ .5
        eta.ec50  ~ .5
        eta.kout ~ .5
        eta.e0 ~ .5
        ##
        pdadd.err <- 10
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        alag1 <- talag
        ##
        #poplogit = log(temax/(1-temax))
        logit=exp(poplogit+eta.emax)
        #logit=temax+eta.emax
        emax = logit/(1+logit)
        ec50 =  exp(tec50 + eta.ec50)
        kout = exp(tkout + eta.kout)
        e0 = exp(te0 + eta.e0)
        ##
        DCP = center/v
        PD=1-emax*DCP/(ec50+DCP)
        ##
        effect(0) = e0
        kin = e0*kout
        ##
        d/dt(gut) =  -ka * gut
        d/dt(center) =  ka * gut - cl / v * center
        d/dt(effect) = kin*PD -kout*effect
        alag(gut) = alag1
        ##
        cp = center / v
        cp ~ prop(prop.err) + add(pkadd.err)
        pca = effect
        pca ~ add(pdadd.err) | pca
      })
    }

    fit.TOF <- expect_error(nlmixr(pk.turnover.emax.lag, warfarin, "focei", foceiControl(print=0, maxOuterIterations = 0)), NA)

    expect_true(inherits(fit.TOF, "nlmixrFitCore"))

  })
}, test="focei")
