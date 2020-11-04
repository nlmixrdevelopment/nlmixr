nlmixrTest(
  {
    pk.turnover.emax3 <- function() {
      ini({
        tktr <- log(1)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(10)
        ##
        eta.ktr ~ 1
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        prop.err <- 0.1
        pkadd.err <- 0.1
        ##
        poplogit <- 2
        # temax <- 7.5
        tec50 <- log(0.5)
        tkout <- log(0.05)
        te0 <- log(100)
        ##
        eta.emax ~ .5
        eta.ec50 ~ .5
        eta.kout ~ .5
        eta.e0 ~ .5
        ##
        pdadd.err <- 10
      })
      model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        ##
        # poplogit = log(temax/(1-temax))
        logit <- exp(poplogit + eta.emax)
        # logit=temax+eta.emax
        emax <- logit(poplogit + eta.emax)
        ec50 <- exp(tec50 + eta.ec50)
        kout <- exp(tkout + eta.kout)
        e0 <- exp(te0 + eta.e0)
        ##
        DCP <- center / v
        PD <- 1 - emax * DCP / (ec50 + DCP)
        ##
        effect(0) <- e0
        kin <- e0 * kout
        ##
        d / dt(depot) <- -ktr * depot
        d / dt(gut) <- ktr * depot - ka * gut
        d / dt(center) <- ka * gut - cl / v * center
        d / dt(effect) <- kin * PD - kout * effect
        ##
        cp <- center / v
        cp ~ prop(prop.err) + add(pkadd.err)
        effect ~ add(pdadd.err) | pca
      })
    }

    f <- nlmixr(pk.turnover.emax3)
    expect_equal(c(1L, 2L, 3L, 4L, NA, NA, 5L, 6L, 7L, 8L, NA), f$saem.theta.trans)
  },
  test = "cran"
)
