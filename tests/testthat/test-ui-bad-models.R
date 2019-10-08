rxPermissive({

    context("Bad UI models should raise errors")

    test_that("Duplicate parameters raise errors", {
        uif <- function(){
            ini({
                lCL = 1.37
                lV = 4.19
                lCLD = 1.37
                lVT = 3.87
                prop.err = 1
                eta.Cl ~ 0.1
                eta.V ~ 0.1
                ## Duplicate CLs
                eta.Cl ~ 0.1
                eta.VT ~ 0.1
            })
            model({
                CL <- exp(lCL + eta.Cl)
                V  <- exp(lV + eta.V)
                CLD <- exp(lCLD + eta.Cl)
                VT <- exp(lVT + eta.VT)
                ## FIXME possibly include both?
                ## K10 <- CL / V
                ## K12 <- CLD / V
                ## K21 <- CLD / VT
                linCmt() ~ prop(prop.err)
            })
        }
        expect_error(nlmixr(uif), rex::rex("The following parameter names were duplicated: eta.Cl."))
    })

    test_that("Un-estimated paramteres raise errors", {

        uif.ode <- function(){
            ini({
                lCL = 1.37
                lV = 4.19
                lCLD = 1.37
                lVT = 3.87
                ## Prop error isn't estimated
                prop.err = 1
                add.err = 0.1
                eta.Cl + eta.V~ c(0.1,
                                  0.01, 0.01)
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
                cp ~ add(add.err)
            })
        }

        expect_error(nlmixr(uif.ode), rex::rex("Model error: initial estimates provided without variables being used:r"))

        uif <- function(){
            ini({
                tka <- exp(0.5)
                tcl <- exp(-3.2)
                tv <- exp(1)
                eta.ka ~ 0.1
                ## Should be eta.cl
                eta.v ~ 0.2
                add.err ~ 0.1
            })
            model({
                ka <- tka + eta.ka
                cl <- tcl + eta.cl
                v <- tv
                d / dt(depot) = -ka * depot
                d / dt(center) = ka * depot - cl / v * center
                cp = center / v
                cp ~ add(add.err)
            })
        }

        expect_error(nlmixr(uif), rex::rex("The following parameter(s) were in the ini block but not in the model block: eta.v\nBad parsed model"))

    })


    test_that("Residuals are population parameters", {

        uif <- function(){
            ini({
                tka <- exp(0.5)
                tcl <- exp(-3.2)
                tv <- exp(1)
                eta.ka ~ 0.1
                eta.cl ~ 0.2
                ## Should be assign since it is a THETa, should I support it....?
                add.err ~ 0.1
            })
            model({
                ka <- tka + eta.ka
                cl <- tcl + eta.cl
                v <- tv
                d / dt(depot) = -ka * depot
                d / dt(center) = ka * depot - cl / v * center
                cp = center / v
                cp ~ add(add.err)
            })
        }

        expect_error(nlmixr(uif), rex::rex("Residual error component(s) need to be defined with assignment ('=' or '<-') in ini block (not '~'): add.err"))

    })

    test_that("Parameters need to be named", {

        uif <- function(){
            ini({
                tka <- exp(0.5)
                tcl <- exp(-3.2)
                tv <- exp(1)
                eta.ka ~ 0.1
                eta.cl ~ 0.2
                ## Should be assign since it is a THETa, should I support it....?
                0.1
            })
            model({
                ka <- tka + eta.ka
                cl <- tcl + eta.cl
                v <- tv
                d / dt(depot) = -ka * depot
                d / dt(center) = ka * depot - cl / v * center
                cp = center / v
                cp ~ add(add.err)
            })
        }
        ##, rex::rex("The following THETAs are unnamed: THETA[4]")

        expect_error(nlmixr(uif))

        uif <- function(){
            ini({
                tka <- exp(0.5)
                tcl <- exp(-3.2)
                tv <- exp(1)
                eta.ka ~ 0.1
                ~ 0.2
                ## Should be assign since it is a THETa, should I support it....?
                add.err = 0.1
            })
            model({
                ka <- tka + eta.ka
                cl <- tcl + eta.cl
                v <- tv
                d / dt(depot) = -ka * depot
                d / dt(center) = ka * depot - cl / v * center
                cp = center / v
                cp ~ add(add.err)
            })
        }

        ## rex::rex("The following ETAs are unnamed: ETA[2]")
        expect_error(nlmixr(uif))

    })

    test_that("Parameters cannot be missing or Infinite", {

        uif <- function(){
            ini({
                tka <- 1 / 0
                tcl <- exp(-3.2)
                tv <- exp(1)
                eta.ka ~ 0.1
                eta.cl ~ 0.2
                ## Should be assign since it is a THETa, should I support it....?
                add.err
            })
            model({
                ka <- tka + eta.ka
                cl <- tcl + eta.cl
                v <- tv
                d / dt(depot) = -ka * depot
                d / dt(center) = ka * depot - cl / v * center
                cp = center / v
                cp ~ add(add.err)
            })
        }
        expect_error(nlmixr(uif), rex::rex("The following parameters initial estimates are infinite: tka"))

        uif <- function(){
            ini({
                tka <- NA
                tcl <- exp(-3.2)
                tv <- exp(1)
                eta.ka ~ 0.1
                eta.cl ~ 0.2
                ## Should be assign since it is a THETa, should I support it....?
                add.err
            })
            model({
                ka <- tka + eta.ka
                cl <- tcl + eta.cl
                v <- tv
                d / dt(depot) = -ka * depot
                d / dt(center) = ka * depot - cl / v * center
                cp = center / v
                cp ~ add(add.err)
            })
        }
        expect_error(nlmixr(uif), rex::rex("The following parameters initial estimates are NA: tka"))

        two.cmt.pd <- function(){
            ini({
                tKa   <- log(0.64)
                tCl   <- log(5.22)
                tV2   <- log(41.3)
                tV3   <- log(115)
                tQ    <- log(11.96)
                BWef  <- log(1.87)
                tSlope     <- log(10) ## add for PD estimation
                tIntercept <- log(1)  ## add for PD estimation
                eta.Ka ~ 1.18
                eta.Cl ~ 0.09
                eta.V2 ~ 0.2
                eta.V3 ~ 0.12
                eta.Q  ~ 0.12
                eta.Slope     ~ 0.1 ## add for PD estimation
                eta.Intercept ~ 0.1 ## add for PD estimation

                prop.err1 <- 0.1  ## Cp
                prop.err2 <- 0.3  ## Ef

            })
            model({
                Ka <- exp(tKa + eta.Ka)
                Cl <- exp(tCl + BWef * log.BW.70 + eta.Cl)
                V2 <- exp(tV2 + eta.V2)
                V3 <- exp(tV3 + eta.V3)
                Q  <- exp(tQ + eta.Q)
                Slope     <- exp(tSlope + eta.Slope)         ## add for PD estimation
                Intercept <- exp(tIntercept + eta.Intercept) ## add for PD estimation

                d/dt(depot)  = -Ka * depot
                d/dt(center) = Ka * depot - Cl/V2 * center + Q/V3 * periph - Q/V2 * center
                d/dt(periph) = Q/V2 * center - Q/V3 * periph

                Cp = center / V2
                Ef = Cp * Slope + Intercept                  ## add for PD estimation

                Cp ~ prop(prop.err1) | Cp
                Ef ~ prop(prop.err2) | Ef                    ## add for PD estimation

            })
        }
        expect_error(nlmixr(two.cmt.pd),
                     rex::rex("The conditional statements (Cp, Ef) are not in terms of the RxODE states: depot, center, periph"))
    })

    test_that("Parameters cannot be missing or Infinite", {

        uif <- function(){
            ini({
                tka <- NA
                tcl <- exp(-3.2)
                tv <- exp(1)
                eta.ka ~ 0.1
                eta.cl ~ 0.2
                add.err = 1
            })
            model({
                ka <- tka + eta.ka
                cl <- tcl + eta.cl
                v <- tv
                d / dt(depot) = -ka * depot
                d / dt(center) = ka * depot - cl / v * center
                cp = center / v
                cp = add(add.err)
            })
        }
        expect_error(nlmixr(uif), "There must be at least one prediction")
    })

}, on.validate="NLMIXR_VALIDATION")
