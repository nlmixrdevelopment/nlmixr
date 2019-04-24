## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
## Load libraries
library(ggplot2)
library(nlmixr)

## To allow nlmixr to reload runs without large run times
## To run the actual models on your system, take the save options off.
options(nlmixr.save=TRUE,
        nlmixr.save.dir=system.file(package="nlmixr"));


str(theo_sd)

ggplot(theo_sd, aes(TIME, DV)) + geom_line(aes(group=ID), col="red") + scale_x_continuous("Time (h)") + scale_y_continuous("Concentration") + labs(title="Theophylline single-dose", subtitle="Concentration vs. time by individual")


## ------------------------------------------------------------------------
one.cmt <- function() {
    ini({
        tka <- .5   # log Ka
        tcl <- -3.2 # log Cl
        tv <- -1    # log V
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        add.err <- 0.1
    })
    model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.err)
    })
}

## ------------------------------------------------------------------------
fit <- nlmixr(one.cmt, theo_sd, est="nlme")
print(fit)

## ------------------------------------------------------------------------
one.compartment <- function() {
    ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45    # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.err <- 0.7
    })
    model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.err)
    })
}

## ------------------------------------------------------------------------
fit <- nlmixr(one.compartment, theo_sd, est="saem")

## ------------------------------------------------------------------------
fitF <- nlmixr(one.compartment, theo_sd, est="focei")

## ------------------------------------------------------------------------
plot(fit)

## ------------------------------------------------------------------------
fit

## ------------------------------------------------------------------------
fit$eta

## ------------------------------------------------------------------------

iter <- fit$par.hist.stacked
iter$Parameter[iter$par=="add.err"] <- "Additive error"
iter$Parameter[iter$par=="eta.cl"]  <- "IIV CL/F"
iter$Parameter[iter$par=="eta.v"]   <- "IIV V/F"
iter$Parameter[iter$par=="eta.ka"]  <- "IIV ka"
iter$Parameter[iter$par=="tcl"]     <- "log(CL/F)"
iter$Parameter[iter$par=="tv"]      <- "log(V/F)"
iter$Parameter[iter$par=="tka"]     <- "log(ka)"
iter$Parameter <- ordered(iter$Parameter, c("log(CL/F)", "log(V/F)", "log(ka)",
                                            "IIV CL/F", "IIV V/F", "IIV ka",
                                            "Additive error"))

ggplot(iter, aes(iter, val)) +
  geom_line(col="red") + 
  scale_x_continuous("Iteration") +
  scale_y_continuous("Value") +
  facet_wrap(~ Parameter, scales="free_y") +
  labs(title="Theophylline single-dose", subtitle="Parameter estimation iterations")


## ------------------------------------------------------------------------

etas <- data.frame(eta = c(fit$eta$eta.ka, fit$eta$eta.cl, fit$eta$eta.v),
                   lab = rep(c("eta(ka)", "eta(CL/F)", "eta(V/F)"), each=nrow(fit$eta)))
etas$lab <- ordered(etas$lab, c("eta(CL/F)","eta(V/F)","eta(ka)"))

ggplot(etas, aes(eta)) +
  geom_histogram(fill="red", col="white") + 
  geom_vline(xintercept=0) +
  scale_x_continuous(expression(paste(eta))) +
  scale_y_continuous("Count") +
  facet_grid(~ lab) +
  coord_cartesian(xlim=c(-1.75,1.75)) +
  labs(title="Theophylline single-dose", subtitle="IIV distributions")


## ------------------------------------------------------------------------
## install.packages("xpose")
library(xpose)

## ------------------------------------------------------------------------
fit <- nlmixr(one.compartment) %>% saem.fit(data=theo_sd)

## ------------------------------------------------------------------------
fit2 <- nlmixr(one.compartment, data=theo_sd, est="saem")

## ------------------------------------------------------------------------
fit3 <- one.compartment %>% saem.fit(data=theo_sd)

## ------------------------------------------------------------------------
fit4 <- nlmixr(one.compartment, theo_sd,est="nlme",control = nlmeControl(pnlsTol = .5))

## ------------------------------------------------------------------------
fit5 <- nlmixr(one.compartment,theo_sd,est="saem",control=saemControl(n.burn=250,n.em=350,print=50))

## ------------------------------------------------------------------------
f <- function(){
    ini({
        lCl <- 1.6      #log Cl (L/hr)
        lVc <- log(90)  #log Vc (L)
        lKA <- 0.1      #log Ka (1/hr)
        prop.err <- c(0, 0.2, 1)
        eta.Cl ~ 0.1   # BSV Cl
        eta.Vc ~ 0.1   # BSV Vc
        eta.KA ~ 0.1   # BSV Ka
    })
    model({
        Cl <- exp(lCl + eta.Cl)
        Vc = exp(lVc + eta.Vc)
        KA <- exp(lKA + eta.KA)
        ## Instead of specifying the ODEs, you can use
        ## the linCmt() function to use the solved system.
        ##
        ## This function determines the type of PK solved system
        ## to use by the parameters that are defined.  In this case
        ## it knows that this is a one-compartment model with first-order
        ## absorption.
        linCmt() ~ prop(prop.err)
    })
}

## ------------------------------------------------------------------------
nlmixr(f)

