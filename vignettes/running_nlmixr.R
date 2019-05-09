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
print(fit)

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
xpdbLoc <- file.path(system.file(package="nlmixr"), "xpdb.rds");
if (file.exists(xpdbLoc)){
    load(xpdbLoc);
} else {
    devtools::install_github("nlmixrdevelopment/xpose.nlmixr")
    library(xpose.nlmixr)
    xp <- xpose_data_nlmixr(fit);
    save(xp, file=xpdbLoc)
}


## ------------------------------------------------------------------------
dv_vs_pred(xp)

## ------------------------------------------------------------------------
dv_vs_ipred(xp)

