## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  out.width = "100%"
)

## ------------------------------------------------------------------------
library(nlmixr)

## These options cache the models and the model simulations in R
## To run the actual models on your system, take the save options off.
options(nlmixr.save=TRUE,
        nlmixr.save.dir=system.file(package="nlmixr"));

invisible(memoise::cache_filesystem(file.path(system.file(package="nlmixr"),"memo")))

## ------------------------------------------------------------------------
pk.turnover.emax <- function() {
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
    ktr <- exp(tktr + eta.ktr)
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
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
    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    d/dt(effect) = kin*PD -kout*effect
    ##
    cp = center / v
    cp ~ prop(prop.err) + add(pkadd.err)
    effect ~ add(pdadd.err)
  })
}

## ------------------------------------------------------------------------
ui <- nlmixr(pk.turnover.emax)
ui

## ------------------------------------------------------------------------
ui$multipleEndpoint

## ------------------------------------------------------------------------
pk.turnover.emax2 <- function() {
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
    ktr <- exp(tktr + eta.ktr)
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
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
    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    d/dt(effect) = kin*PD -kout*effect
    ##
    cp = center / v
    cp ~ prop(prop.err) + add(pkadd.err) | center
    effect ~ add(pdadd.err)
  })
}
ui2 <- nlmixr(pk.turnover.emax2)
ui2$multipleEndpoint

## ------------------------------------------------------------------------
options(RxODE.combine.dvid=FALSE)
ui2$multipleEndpoint

## ------------------------------------------------------------------------
options(RxODE.combine.dvid=TRUE)
ui2$multipleEndpoint

## ------------------------------------------------------------------------
summary(warfarin)

