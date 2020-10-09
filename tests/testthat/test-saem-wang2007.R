ctnlmixrTest({

  context("Test SAEM functions")
  test_that("saem objective function tests", {

    mod <- function() {
      ini({
        tk <- log(0.5)
        eta.k ~ 0.04
        add.sd <- 0.1
      })
      model({
        ke <- exp(tk + eta.k)
        ipre <- 10 * exp(-ke * t)
        ipre ~ add(add.sd)
      })
    }

    mod2 <- function() {
      ini({
        tk <- log(0.5)
        eta.k ~ 0.04
        lnorm.sd <- 0.1
      })
      model({
        ke <- exp(tk + eta.k)
        ipre <- log(10 * exp(-ke * t))
        ipre ~ add(lnorm.sd)
      })
    }

    dat <- Wang2007
    dat$DV <- dat$Y

    ctl <- saemControl(nEm=5, nBurn = 5, logLik=TRUE, print=0)

    f <- suppressWarnings(nlmixr(mod, dat, "saem", control=ctl))
    expect_equal(round(f$objective, 3), -6.192)

    mod %>%
      model(ipre ~ prop(prop.sd)) %>%
      ini(prop.sd=0.1) %>%
      nlmixr(dat, "saem", control=ctl) -> f
    expect_equal(round(f$objective, 3), -3.862)

    mod %>%
      model(ipre ~ add(add.sd) + prop(prop.sd)) %>%
      ini(c(prop.sd=0.1, add.sd=0.1)) %>%
      nlmixr(dat, "saem", control=saemControl(nEm=5, nBurn = 5, logLik=TRUE, addProp="combined2")) -> f
    expect_equal(round(f$objective, 3), 87.124)

    mod %>%
      model(ipre ~ add(add.sd) + prop(prop.sd)) %>%
      ini(c(prop.sd=0.1, add.sd=0.1)) %>%
      nlmixr(dat, "saem", control=saemControl(nEm=5, nBurn = 5, logLik=TRUE, addProp="combined1")) -> f
    expect_equal(round(f$objective, 3), 82.95)

    mod %>%
      model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) %>%
      ini(c(prop.sd=0.1, add.sd=0.1, lambda=1)) %>%
      nlmixr(dat, "saem", control=saemControl(nEm=5, nBurn = 5, logLik=TRUE, addProp="combined2")) -> f
    expect_equal(round(f$objective, 3), 87.124)

    mod %>%
      model(ipre ~ lnorm(lnorm.sd)) %>%
      ini(lnorm.sd=0.1) %>%
      nlmixr(dat, "saem", control=ctl) -> f

    dat2 <- dat
    dat2$DV <- log(dat2$Y)

    nlmixr(mod2, dat2, est="saem", control=ctl) -> f2

    expect_equal(f$theta, f2$theta)
    expect_equal(f$omega, f2$omega)

    expect_equal(f$PRED, exp(f2$PRED))
    expect_equal(f$IPRED, exp(f2$IPRED))
    expect_equal(f$IWRES, f2$IWRES)
    expect_equal(f$cov, f2$cov)

    expect_equal(f$objective , f2$objective - 2 * sum(dat2$DV))

    expect_equal(round(f$objective, 3), -165.963)


  })
}, test="saem")
