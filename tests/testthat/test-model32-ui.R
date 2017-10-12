library(testthat)
library(nlmixr)

context("UI-NLME32: two-compartment bolus, single-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {
    datr <- Bolus_2CPT

    datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
    datr <- datr[datr$EVID != 2,]

    runno <- "N032"

    dat <- datr[datr$SD == 1,]

    test_that("Closed-form", {

        uif <- function(){
            ini({
                lCL = 1.37
                lV = 4.19
                lCLD = 1.37
                lVT = 3.87
                prop.err = 1
                eta.Cl ~ 0.1
                eta.V ~ 0.1
                eta.Cld ~ 0.1
                eta.VT ~ 0.1
            })
            model({
                CL <- exp(lCL + eta.Cl)
                V  <- exp(lV + eta.V)
                CLD <- exp(lCLD + eta.Cld)
                VT <- exp(lVT + eta.VT)
                ## FIXME possibly include both?
                ## K10 <- CL / V
                ## K12 <- CLD / V
                ## K21 <- CLD / VT
                linCmt() ~ prop(prop.err)
            })
        }

        fit <- nlmixr(uif, dat, est="nlme",
                      control=nlmeControl(pnlsTol = .3, msVerbose = TRUE), focei.translate=FALSE)

    })

    test_that("ODE", {

        uif.ode <- function(){
            ini({
                lCL = 1.37
                lV = 4.19
                lCLD = 1.37
                lVT = 3.87
                prop.err = 1
                eta.Cl ~ 0.1
                eta.V ~ 0.1
                eta.Cld ~ 0.1
                eta.VT ~ 0.1
            })
            model({
                CL <- exp(lCL + eta.Cl)
                V  <- exp(lV + eta.V)
                CLD <- exp(lCLD + eta.Cld)
                VT <- exp(lVT + eta.VT)
                K10 <- CL / V
                K12 <- CLD / V
                K21 <- CLD / VT
                d/dt(centr)  = K21*periph-K12*centr-K10*centr;
                d/dt(periph) =-K21*periph+K12*centr;
                cp = centr / V
                cp ~ prop(prop.err)
            })
        }

        fit <- nlmixr(uif.ode, dat, est="nlme",
                      control=nlmeControl(pnlsTol = .3, msVerbose = TRUE), focei.translate=FALSE)


        uif.ode.add <- function(){
            ini({
                lCL = 1.37
                lV = 4.19
                lCLD = 1.37
                lVT = 3.87
                add.err = 1
                eta.Cl ~ 0.1
                eta.V ~ 0.1
                eta.Cld ~ 0.1
                eta.VT ~ 0.1
            })
            model({
                CL <- exp(lCL + eta.Cl)
                V  <- exp(lV + eta.V)
                CLD <- exp(lCLD + eta.Cld)
                VT <- exp(lVT + eta.VT)
                K10 <- CL / V
                K12 <- CLD / V
                K21 <- CLD / VT
                d/dt(centr)  = K21*periph-K12*centr-K10*centr;
                d/dt(periph) =-K21*periph+K12*centr;
                cp = centr / V
                cp ~ add(add.err)
            })
        }

        fit <- nlmixr(uif.ode.add, dat, est="nlme",
                      control=nlmeControl(pnlsTol = .3, msVerbose = TRUE), focei.translate=FALSE)


        uif.ode.add.prop <- function(){
            ini({
                lCL = 1.37
                lV = 4.19
                lCLD = 1.37
                lVT = 3.87
                err.add = 1
                err.prop = 1
                eta.Cl + eta.V ~ c(0.1,
                                   0, 0.01)
            })
            model({
                CL <- exp(lCL + eta.Cl)
                V  <- exp(lV + eta.V)
                CLD <- exp(lCLD)
                VT <- exp(lVT)
                K10 <- CL / V
                K12 <- CLD / V
                K21 <- CLD / VT
                d/dt(centr)  = K21*periph-K12*centr-K10*centr;
                d/dt(periph) =-K21*periph+K12*centr;
                cp = centr / V
                cp ~ add(err.add) + prop(err.prop)
            })
        }

        fit <- nlmixr(uif.ode.add.prop, dat, est="nlme",
                      control=nlmeControl(pnlsTol = .3, msVerbose = TRUE),
                      focei.translate=FALSE)

    })

}
