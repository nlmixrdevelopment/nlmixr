nlmixrTest({
  context("Back-transform vs scaling checks (#161)")
  test_that("Log-scaled vs Back-transformed parameters", {
    run7.mod <- function() {
      ini({
        tcl <- log(0.008) # typical value of clearance
        tv <- log(0.6) # typical value of volume
        all.cl <- 1 # allometric exponent on CL
        eta.cl ~ 1
        # interindividual variability on clearance and volume
        add.err <- 0.1 # residual variability
      })
      model({
        cl <- exp(tcl + all.cl * log_allo_wt + eta.cl) # individual value of clearance
        v <- exp(tv) # individual value of volume
        ke <- cl / v # elimination rate constant
        d / dt(A1) <- -ke * A1 # model differential equation
        cp <- A1 / v # concentration in plasma
        cp ~ add(add.err) # define error model
      })
    }

    f <- nlmixr(run7.mod)

    expect_equal(list(1:3, 1:2), f$logThetasList)
  })

}, test="cran")
