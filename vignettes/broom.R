## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(nlmixr)
## To allow nlmixr to reload runs without large run times
## To run the actual models on your system, take the save options off.
options(nlmixr.save=TRUE,
        nlmixr.save.dir=system.file(package="nlmixr"));

pheno <- function() {
    # Pheno with covariance
  ini({
    tcl <- log(0.008) # typical value of clearance
    tv <-  log(0.6)   # typical value of volume
    ## var(eta.cl)
    eta.cl + eta.v ~ c(1, 
                       0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
                      # interindividual variability on clearance and volume
    add.err <- 0.1    # residual variability
  })
  model({
    cl <- exp(tcl + eta.cl) # individual value of clearance
    v <- exp(tv + eta.v)    # individual value of volume
    ke <- cl / v            # elimination rate constant
    d/dt(A1) = - ke * A1    # model differential equation
    cp = A1 / v             # concentration in plasma
    cp ~ add(add.err)       # define error model
  })
}

## We will run it two ways to allow comparisons
fit.s <- nlmixr(pheno, pheno_sd, "saem", control=list(logLik=TRUE), table=list(cwres=TRUE))
fit.f <- nlmixr(pheno, pheno_sd, "focei")

## ------------------------------------------------------------------------
glance(fit.s)

## ------------------------------------------------------------------------
setOfv(fit.s,"FOCEi")

## ------------------------------------------------------------------------
glance(fit.s)

## ------------------------------------------------------------------------
setOfv(fit.s,"gauss3_1.6") # Setting objective function to gauss3_1.6

## ------------------------------------------------------------------------
glance(fit.s)

## ------------------------------------------------------------------------
glance(fit.s, type="FOCEi")

## ------------------------------------------------------------------------
tidy(fit.s)

## ------------------------------------------------------------------------
tidy(fit.s, exponentiate=FALSE) ## No transformation applied

## ------------------------------------------------------------------------
tidy(fit.s, exponentiate=TRUE) ## No transformation applied on every parameter

## ------------------------------------------------------------------------
options(broom.mixed.sep2="..")
tidy(fit.s)

## ------------------------------------------------------------------------
confint(fit.s)

## ------------------------------------------------------------------------
confint(fit.s, exponentiate=FALSE)

## ------------------------------------------------------------------------
tidy(fit.s, conf.level=0.9)

## ------------------------------------------------------------------------
tidy(fit.s, effects="fixed")

## ------------------------------------------------------------------------
tidy(fit.s, effects="ran_pars")

## ------------------------------------------------------------------------
tidy(fit.s, effects="ran_vals")

## ------------------------------------------------------------------------
tidy(fit.s, effects="ran_coef")

## ------------------------------------------------------------------------
tidy(fit.s, effects="ran_coef", exponentiate=FALSE)

## ------------------------------------------------------------------------
library(ggplot2)
library(dotwhisker)
library(dplyr)
fit.s %>%
    tidy() %>%
    filter(effect=="fixed") %>%
    dwplot()

## ------------------------------------------------------------------------
dwplot(list("SAEM"=fit.s, "FOCEi"=fit.f))

## ------------------------------------------------------------------------
library(huxtable)
tbl <- huxreg('Phenobarbitol'=fit.s)

tbl

## ------------------------------------------------------------------------
as_hux('Phenobarbitol'=fit.s)

## ------------------------------------------------------------------------
as_hux('SAEM'=fit.s, 'FOCEi'=fit.f)

## ------------------------------------------------------------------------
library(officer)
library(flextable)

ft  <- huxtable::as_flextable(tbl);
    
read_docx() %>%
    flextable::body_add_flextable(ft)  %>%
    print(target="pheno.docx")

