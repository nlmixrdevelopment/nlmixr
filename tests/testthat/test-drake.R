context("test drake parital compatibility fixes.");
oldOp <- options()
options(drake_make_menu = FALSE)
library(drake)
library(nlmixr)

.cacheDir <- file.path(tempdir(), ".drake")
if (file.exists(.cacheDir)) unlink(.cacheDir, recursive=TRUE)
.cache <- drake_cache(.cacheDir)

test_that("parsing functions fix drake compatibility", {
    ## It is best to define functions outside the plan.`
    ## Targets are for long-ish computations that produce data.
    model_fn <- function() {
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
            depot(0) = 3
            d/dt(center) = ka * depot - cl / v * center
            cp = center / v
            cp ~ add(add.err)
        })
    }

    ## By parsing the model first, it extracts interesting information
    ## about the model and makes it compatible with drake.  Note to track
    ## model changes you need at least drake > 7.4.0, right now it is in
    ## github at 7.4.0.9000
    model_p <- nlmixr(model_fn)

    ## This changes the model function but you cannot see it directly.
    ## It appears that it has not changed

    run_nlmixr <- function(...) {NULL}
    your_data <- NULL

    plan <- drake_plan(analysis = run_nlmixr(model_p, your_data),
                       try2=run_nlmixr(model_fn, your_data))

    expect_equal(make(plan, cache=.cache), NULL)
})

test_that("drake doesn't work without parsing the model", {
    model_fn <- function() {
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
            depot(0) = 3
            d/dt(center) = ka * depot - cl / v * center
            cp = center / v
            cp ~ add(add.err)
        })
    }

    run_nlmixr <- function(...) {NULL}
    your_data <- NULL

    plan <- drake_plan(analysis = run_nlmixr(model_fn, your_data))

    expect_error(make(plan, cache=.cache))
})

unlink(.cacheDir, recursive=TRUE)

options(oldOp)
