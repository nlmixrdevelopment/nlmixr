ctnlmixrTest({

  context("Test SAEM functions")
  test_that("saem objective function tests", {

    mod <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45    # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        ipre = center / v
        ipre ~ add(add.sd)
      })
    }


    mod2 <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45    # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        lnorm.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        ipre = log(center / v)
        ipre ~ add(lnorm.sd)
      })
    }

    dat <- theo_sd
    dat <- dat[!(dat$TIME == 0 & dat$EVID == 0),]
    

    ctl1 <- saemControl(nEm=5, nBurn = 5, logLik=TRUE, print=0, addProp="combined1")
    ctl2 <- saemControl(nEm=5, nBurn = 5, logLik=TRUE, print=0, addProp="combined2")

    f <- suppressWarnings(nlmixr(mod, dat, "saem", control=ctl2))
    expect_equal(round(f$objective, 3), 114.124)

    mod %>%
      model(ipre ~ prop(prop.sd)) %>%
      ini(prop.sd=0.1) %>%
      nlmixr(dat, "saem", control=ctl2) -> f
    
    expect_equal(round(f$objective, 3), 127.313)

    mod %>%
      model(ipre ~ add(add.sd) + prop(prop.sd)) %>%
      ini(c(prop.sd=0.1, add.sd=0.1)) %>%
      nlmixr(dat, "saem", control=ctl1) -> f
    
    expect_equal(round(f$objective, 3), 473.996)

    mod %>%
      model(ipre ~ add(add.sd) + prop(prop.sd)) %>%
      ini(c(prop.sd=0.1, add.sd=0.1)) %>%
      nlmixr(dat, "saem", control=ctl2) -> f
    
    expect_equal(round(f$objective, 3), 473.276)

    mod %>%
      model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(c(prop.sd=0.1, add.sd=0.1, lambda=1)) %>%
      nlmixr(dat, "saem", control=ctl1) -> f
    
    expect_equal(round(f$objective, 3), 476.715)

    mod %>%
      model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(c(prop.sd=0.1, add.sd=0.1, lambda=1)) %>%
      nlmixr(dat, "saem", control=ctl2) -> f
    
    expect_equal(round(f$objective, 3), 490.495)
    
    mod %>%
      model(ipre ~ lnorm(lnorm.sd)) %>%
      ini(lnorm.sd=0.1) %>%
      nlmixr(dat, "saem", control=ctl2) -> f

    dat2 <- dat
    dat2$Y <- dat2$DV
    dat2$DV <- log(dat2$Y)

    nlmixr(mod2, dat2, est="saem", control=ctl2) -> f2

    expect_equal(f$theta, f2$theta)
    expect_equal(f$omega, f2$omega)

    expect_equal(f$PRED, exp(f2$PRED))
    expect_equal(f$IPRED, exp(f2$IPRED))
    expect_equal(f$IWRES, f2$IWRES)
    expect_equal(f$cov, f2$cov)

    expect_equal(f$objective , f2$objective - 2 * sum(dat2$DV[dat2$EVID == 0]))
    expect_equal(round(f$objective, 3), -606.185)

  })
}, test="saem")
