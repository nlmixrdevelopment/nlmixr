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

