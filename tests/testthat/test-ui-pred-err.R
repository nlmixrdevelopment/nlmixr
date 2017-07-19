fn1 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ pois()
}


fn2 <- function() {
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop(Prop.Err, Prop.Err2)
}

fn3 <- function() {
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop()
}



fn4 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ pois(p1, p2)
}


fn5 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dpois();
}

fn6 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dpois(1, 2);
}

fn7 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbinom();
}

fn8 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbinom(13);
}

fn9 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbinom(13, 0.5, 32);
}

fn10 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbeta();
}

fn11 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbeta(0.1);
}

fn12 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbeta(0.1, 0.1, 0.3, 0.1);
}

fn13 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dt();
}

fn14 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dt(0.1, 0.1, 0.2);
}


fn15 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ binom();
}

fn16 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ binom(13);
}

fn17 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ binom(13, 0.5, 32);
}

fn18 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ beta();
}

fn19 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ beta(0.1);
}

fn20 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ beta(0.1, 0.1, 0.3, 0.1);
}

## Should ~ t() be supported?
## Currently t() is also transpose in R.

fn21 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ t();
}

fn22 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ t(0.1, 0.1, 0.2);
}


fn23 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ add();
}

fn24 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ add(0.1, 0.1);
}

fn25 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop();
}

fn26 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop(0.1, 0.1);
}

fn27 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ norm();
}

fn28 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ norm(0.1, 0.1);
}

fn29 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dnorm();
}

fn30 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dnorm(0.1, 0.1);
}

fn31 <- function(){
    KA = KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ nlmixrDist(0.1, 0.1);
}

context("Improperly specified residuals distributions throw errors")

test_that("Improper distribution functions throw errors", {
    expect_error(nlmixrUIModelPredErr(fn1), "The pois distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn2), "The prop distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn3), "The prop distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn4), "The pois distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn5), "The dpois distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn6), "The dpois distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn7), "The dbinom distribution requires 2 arguments.")
    expect_error(nlmixrUIModelPredErr(fn8), "The dbinom distribution requires 2 arguments.")
    expect_error(nlmixrUIModelPredErr(fn9), "The dbinom distribution requires 2 arguments.")
    expect_error(nlmixrUIModelPredErr(fn10), "The dbeta distribution requires 2-3 arguments.")
    expect_error(nlmixrUIModelPredErr(fn11), "The dbeta distribution requires 2-3 arguments.")
    expect_error(nlmixrUIModelPredErr(fn12), "The dbeta distribution requires 2-3 arguments.")
    expect_error(nlmixrUIModelPredErr(fn13), "The dt distribution requires 1-2 arguments.");
    expect_error(nlmixrUIModelPredErr(fn14), "The dt distribution requires 1-2 arguments.");
    expect_error(nlmixrUIModelPredErr(fn15), "The binom distribution requires 2 arguments.")
    expect_error(nlmixrUIModelPredErr(fn16), "The binom distribution requires 2 arguments.")
    expect_error(nlmixrUIModelPredErr(fn17), "The binom distribution requires 2 arguments.")
    expect_error(nlmixrUIModelPredErr(fn18), "The beta distribution requires 2-3 arguments.")
    expect_error(nlmixrUIModelPredErr(fn19), "The beta distribution requires 2-3 arguments.")
    expect_error(nlmixrUIModelPredErr(fn20), "The beta distribution requires 2-3 arguments.")
    expect_error(nlmixrUIModelPredErr(fn21), "The t distribution requires 1-2 arguments.");
    expect_error(nlmixrUIModelPredErr(fn22), "The t distribution requires 1-2 arguments.");
    expect_error(nlmixrUIModelPredErr(fn23), "The add distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn24), "The add distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn25), "The prop distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn26), "The prop distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn27), "The norm distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn28), "The norm distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn29), "The dnorm distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn30), "The dnorm distribution requires 1 arguments.")
    expect_error(nlmixrUIModelPredErr(fn31), "The nlmixrDist distribution is currently unsupported.")
})
