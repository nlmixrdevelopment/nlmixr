rxPermissive({

    context("Good UI models should not raise errors")

    one.compartment.saem <- function() {
        ini({
            tka <- .5   # Log Ka
            tcl <- -3.2 # Log Cl
            tv <- -1    # Log V
            eta.ka ~ 1
            eta.cl ~ 2
            eta.v ~ 1
            add.err <- 0.1
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d/dt(depot) = -ka * depot + exp(-k0 * t)
            d/dt(center) = ka * depot - cl / v * center
            cp = center / v
            cp ~ add(add.err)
        })
    }

    expect_equal(nlmixr(one.compartment.saem)$all.covs, "k0")

}, on.validate="NLMIXR_VALIDATION")
