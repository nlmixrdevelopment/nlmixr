nlmixrTest({
  context("Multiple endpoint parsing")
  test_that("multiple-endpoint", {
    tmp <- theo_sd
    tmp$EVID <- ifelse(tmp$EVID == 0, 0, 1)
    tmp$CMT <- ifelse(tmp$CMT == 1, "depot", "central")
    tmp2 <- tmp[tmp$CMT == "central", ]
    tmp2$DV <- tmp2$DV + 20
    tmp2$CMT <- "c20"
    df <- rbind(tmp, tmp2)
    df$CMT <- factor(df$CMT)

    one.compartment.saem <- function() {
      ini({
        tka <- .5 # Log Ka
        tcl <- -3.2 # Log Cl
        tv <- -1 # Log V
        extra <- 20
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        add.err <- 0.1
        add.err2 <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp2 <- cp + extra
        cp ~ add(add.err) | center
        cp2 ~ add(add.err2) | c20
      })
    }

    tmp <- nlmixr(one.compartment.saem)

    expect_equal(tmp$predDf$cmt, c(2L, 3L))
    expect_equal(tmp$predDf$dvid, c(1L, 2L))
    expect_equal(tmp$predDf$cond, c("center", "c20"))
    expect_equal(tmp$predDf$var, c("cp", "cp2"))

    ## tmp2 <- RxODE::etTrans(df,tmp$rxode,TRUE)

    one.compartment.saem <- function() {
      ini({
        tka <- .5 # Log Ka
        tcl <- -3.2 # Log Cl
        tv <- -1 # Log V
        extra <- 20
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        add.err <- 0.1
        add.err2 <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp2 <- cp + extra
        cp ~ add(add.err)
        cp2 ~ add(add.err2)
      })
    }

    tmp <- nlmixr(one.compartment.saem)
    expect_true(inherits(tmp, "nlmixrUI"))



    one.compartment.saem <- function() {
      ini({
        tka <- .5 # Log Ka
        tcl <- -3.2 # Log Cl
        tv <- -1 # Log V
        extra <- 20
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        add.err <- 0.1
        add.err2 <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp2 <- cp + extra
        cp ~ add(add.err) | center
        cp2 ~ add(add.err2)
      })
    }

    d <- nlmixr(one.compartment.saem)
    expect_true(inherits(d, "nlmixrUI"))

    one.compartment.saem <- function() {
      ini({
        tka <- .5 # Log Ka
        tcl <- -3.2 # Log Cl
        tv <- -1 # Log V
        extra <- 20
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        add.err <- 0.1
        add.err2 <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp2 <- cp + extra
        cp ~ add(add.err)
        cp2 ~ add(add.err2)
      })
    }

    d <- nlmixr(one.compartment.saem)
    expect_true(inherits(d, "nlmixrUI"))

    pk.turnover.emax <- function() {
      ini({
        tktr <- log(0.00001)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(1)

        eta.ktr ~ 1
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        prop.err <- 1
        pkadd.err <- 0.00002

        poplogit <- 2
        tec50 <- log(0.5)
        tkout <- log(0.05)
        te0 <- log(100)

        eta.emax ~ .5
        eta.ec50 ~ .5
        eta.kout ~ .5
        eta.e0 ~ .5

        pdadd.err <- 4
      })
      model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)

        # poplogit = log(temax/(1-temax))
        logit <- exp(poplogit + eta.emax)
        # logit=temax+eta.emax
        emax <- logit / (1 + logit)
        ec50 <- exp(tec50 + eta.ec50)
        kout <- exp(tkout + eta.kout)
        e0 <- exp(te0 + eta.e0)

        DCP <- center / v
        PD <- 1 - emax * DCP / (ec50 + DCP)

        effect(0) <- e0
        kin <- e0 * kout

        d / dt(depot) <- -ktr * depot
        d / dt(gut) <- ktr * depot - ka * gut
        d / dt(center) <- ka * gut - cl / v * center
        d / dt(effect) <- kin * PD - kout * effect

        cp <- center / v
        cp ~ prop(prop.err) + add(pkadd.err)
        pca ~ add(pdadd.err)
      })
    }

    expect_error(nlmixr(pk.turnover.emax))

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
        emax <- logit / (1 + logit)
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

    pk.turnover.emax4 <- function() {
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
        emax <- logit / (1 + logit)
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
        ## Who would use this...
        log(cp) ~ prop(prop.err) + add(pkadd.err) | center
        effect * 1 / (1 + cp) ~ add(pdadd.err) | pca
      })
    }

    expect_true(inherits(nlmixr(pk.turnover.emax4), "nlmixrUI"))

    pk.turnover.emax4 <- function() {
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
        emax <- logit / (1 + logit)
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
        ## Who would use this...
        log(cp) ~ prop(prop.err) + add(pkadd.err)
        effect * 1 / (1 + cp) ~ add(pdadd.err) | pca
      })
    }

    expect_error(nlmixr(pk.turnover.emax4))

    pk.turnover.emax4 <- function() {
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
        emax <- logit / (1 + logit)
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
        ## Who would use this...
        log(cp) ~ prop(prop.err) + add(pkadd.err) | center
        log(pca) ~ add(pdadd.err) | cmt
      })
    }

    expect_error(nlmixr(pk.turnover.emax4))
  })
}, test="cran")
