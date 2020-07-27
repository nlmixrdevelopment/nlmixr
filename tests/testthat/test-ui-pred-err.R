nlmixrTest({
  fn1 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ pois()
  }


  fn2 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop(Prop.Err, Prop.Err2)
  }

  fn3 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop()
  }

  fn4 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ pois(p1, p2)
  }


  fn5 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dpois()
  }

  fn6 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dpois(1, 2)
  }

  fn7 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbinom()
  }

  fn8 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbinom(13)
  }

  fn9 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbinom(13, 0.5, 32)
  }

  fn10 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbeta()
  }

  fn11 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbeta(par1)
  }

  fn12 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dbeta(par1, par2, par3, par4)
  }

  fn13 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dt()
  }

  fn14 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dt(par1, par2, par3)
  }


  fn15 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ binom()
  }

  fn16 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ binom(par1)
  }

  fn17 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ binom(par1, par2, par3)
  }

  fn18 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ beta()
  }

  fn19 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ beta(par1)
  }

  fn20 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ beta(par1, par2, par3, par4)
  }

  ## Should ~ t() be supported?
  ## Currently t() is also transpose in R.

  fn21 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ t()
  }

  fn22 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ t(par1, par2, par3)
  }


  fn23 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ add()
  }

  fn24 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ add(par1, par2)
  }

  fn25 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop()
  }

  fn26 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop(par1, par2)
  }

  fn27 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ norm()
  }

  fn28 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ norm(par1, par2)
  }

  fn29 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dnorm()
  }

  fn30 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ dnorm(par1, par2)
  }

  fn31 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ nlmixrDist(par1, par2)
  }

  fn32 <- function() {
    ini({
      KA <- c(0, 1)
      CL <- c(0, 0.5)
    })
    model({
      KA <- KA + eta.KA
      CL <- CL + eta.CL + add(par1)
      v1 <- 1
      linCmt() ~ pois()
    })
  }

  fn33 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ add(par1) + pois(par2)
  }

  fn36 <- function() {
    KA <- KA + eta.KA
    CL <- CL + eta.CL
    linCmt() ~ prop(par1) + add(par2) + pois(par3)
  }

  context("Improperly specified residuals distributions throw errors")

  test_that("Improper distribution functions throw errors", {
    ## expect_error(nlmixr:::nlmixrUIModel(fn1), "The pois distribution requires 1 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn2), "The prop distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn3), "The prop distribution requires 1 argument.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn4), "The pois distribution requires 1 argument.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn5), "The dpois distribution requires 1 arguments.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn6), "The dpois distribution requires 1 arguments.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn7), "The dbinom distribution requires 2 arguments.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn8), "The dbinom distribution requires 2 arguments.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn9), "The dbinom distribution requires 2 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn10), "The dbeta distribution requires 2-3 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn11), "The dbeta distribution requires 2-3 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn12), "The dbeta distribution requires 2-3 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn13), "The dt distribution requires 1-2 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn14), "The dt distribution requires 1-2 arguments.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn15), "The binom distribution requires 2 arguments.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn16), "The binom distribution requires 2 arguments.")
    ## expect_error(nlmixr:::nlmixrUIModel(fn17), "The binom distribution requires 2 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn18), "The beta distribution requires 2-3 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn19), "The beta distribution requires 2-3 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn20), "The beta distribution requires 2-3 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn21), "The t distribution requires 1-2 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn22), "The t distribution requires 1-2 arguments.")
    expect_error(nlmixr:::nlmixrUIModel(fn23), "The add distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn24), "The add distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn25), "The prop distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn26), "The prop distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn27), "The norm distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn28), "The norm distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn29), "The dnorm distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn30), "The dnorm distribution requires 1 argument.")
    expect_error(nlmixr:::nlmixrUIModel(fn31), "The nlmixrDist distribution is currently unsupported.")
    expect_error(nlmixr:::nlmixr(fn32), rex::rex("Distributions need to be on residual model lines (like f ~ add(add.err)).\nMisplaced Distribution(s): add"))
    expect_error(nlmixr:::nlmixrUIModel(fn33), rex::rex("The add and pois distributions cannot be combined\nCurrently can combine: add, prop"))
  })

  context("Proper Variances")

  fn1 <- function() {
    ini({
      KA <- c(0, 1)
      CL <- c(0, 0.5)
    })
    model({
      KA <- KA + eta.KA
      CL <- CL + eta.CL
      v1 <- 1
      linCmt() ~ pois()
    })
  }

  fn34 <- function() {
    ini({
      KA <- c(0, 1)
      CL <- c(0, 0.5)
      par1 <- 1
      par2 <- 2
    })
    model({
      KA <- KA + eta.KA
      CL <- CL + eta.CL
      v1 <- 1
      linCmt() ~ add(par1) + prop(par2)
    })
  }

  fn35 <- function() {
    ini({
      KA <- c(0, 1)
      CL <- c(0, 0.5)
      par1 <- 1
      par2 <- 2
    })
    model({
      KA <- KA + eta.KA
      CL <- CL + eta.CL
      v1 <- 1
      linCmt() ~ prop(par1) + add(par2)
    })
  }

  test_that("Good Parsing of proper variance specifications", {
    expect_equal(class(nlmixr:::nlmixr(fn1)), "nlmixrUI")
    expect_equal(class(nlmixr:::nlmixr(fn34)), "nlmixrUI")
    expect_equal(class(nlmixr:::nlmixr(fn35)), "nlmixrUI")
  })
},
test="cran")
