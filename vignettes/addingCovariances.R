## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  out.width = "100%"
)

## ------------------------------------------------------------------------
## Load Phenobarb data
library(nlmixr)

## These options cache the models and the model simulations in R
## To run the actual models on your system, take the save options off.
options(nlmixr.save=TRUE,
        nlmixr.save.dir=system.file(package="nlmixr"));

invisible(memoise::cache_filesystem(file.path(system.file(package="nlmixr"),"memo")))


## ------------------------------------------------------------------------

pheno <- function() {
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


## ------------------------------------------------------------------------
fit <- nlmixr(pheno, pheno_sd, "saem")

print(fit)

## ------------------------------------------------------------------------
plot(fit)

## ------------------------------------------------------------------------
plot(augPred(fit))

## ------------------------------------------------------------------------
library(ggplot2)
## A traditional VPC
p1 <- vpc(fit, show=list(obs_dv=TRUE)) + ylab("Concentrations")

## A prediction-corrected VPC
p2 <- vpc(fit, pred_corr = TRUE, show=list(obs_dv=TRUE)) + ylab("Prediction-Corrected Concentrations")

library(gridExtra)
grid.arrange(p1,p2)

