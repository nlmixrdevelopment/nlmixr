## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------

library(nlmixr)

## To allow nlmixr to reload runs without large run times
## To run the actual models on your system, take the save options off.
options(nlmixr.save=TRUE,
        nlmixr.save.dir=system.file(package="nlmixr"));

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
fit <- nlmixr(one.compartment, theo_sd, est="focei")

print(fit)

## ------------------------------------------------------------------------
## Example 2 -- Fix tka to 0.5 and re-estimate.
one.ka.0.5 <- fit %>%
    update(tka=fix(0.5)) %>%
    nlmixr(est="focei")

print(one.ka.0.5)

## ------------------------------------------------------------------------
## Example 3 -- Fix tka to model estimated value and re-estimate.
one.ka.0.5 <- fit %>%
    update(tka=fix) %>%
    nlmixr(est="focei")

print(one.ka.0.5)

## ------------------------------------------------------------------------
## Example 4 -- Change tka to 0.7 in orginal model function and then estimate
one.ka.0.7 <- one.compartment %>% 
    update(tka=0.7) %>% 
    nlmixr(data=theo_sd,est="focei")

print(one.ka.0.7)

## ------------------------------------------------------------------------
## Remove eta.ka on ka
noEta <- fit %>%
    update(ka <- exp(tka)) %>%
    nlmixr(est="focei")

print(noEta)

## ------------------------------------------------------------------------
addBackKa <- noEta %>%
    update({ka <- exp(tka + bsv.ka)}) %>%
    ini(bsv.ka=0.1) %>%
    nlmixr(est="focei")

print(addBackKa)

## ------------------------------------------------------------------------
addBackKa$omega

## ------------------------------------------------------------------------
## Note currently cov is needed as a prefix so nlmixr knows this is an
## estimated parameter not a parameter
wt70 <- fit %>% ## FIXME could be based on if it finds the covarite in the last used nlmixr data.
    update({cl <- exp(tcl + eta.cl)*(WT/70)^covWtPow}) %>%
    update(covWtPow=fix(0.75)) %>%
    nlmixr(est="focei")

print(wt70)

## ------------------------------------------------------------------------
## Since there are 0 observations in the data, these are changed to 0.0150 to show proportional error change.
d <- theo_sd
d$DV[d$EVID == 0 & d$DV == 0] <- 0.0150

addPropModel <- fit %>%
    update({cp ~ add(add.err)+prop(prop.err)}) %>%
    update(prop.err=0.1) %>%
    nlmixr(d,est="focei")

print(addPropModel)

