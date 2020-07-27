nlmixrTest({

  context("Back Transformation tests")

  test_that("PK106 mixed mu-referenced and non mu-referenced covariates", {
    PK106 <- function() {
      ini({
        ## Where initial conditions/variables are specified
        lktr <- log(1.15) # log k transit (/h)
        lcl <- log(0.135) # log Cl (L/hr)
        lv <- log(8) # log V (L)
        ALLC <- fix(0.75) # allometric exponent cl
        ALLV <- fix(1.00) # allometric exponent v
        Sex_V <- 0.1 # log Sex on v
        prop.err <- 0.15
        add.err <- 0.6
        eta.ktr ~ 0.5 # IIV ktr
        eta.cl ~ 0.1 # IIV Cl
        eta.v ~ 0.1 # IIV V
      })
      model({
        ## Allometric scaling on weight
        cl <- exp(lcl + eta.cl + ALLC * (log(WT / 70)))
        v <- exp(lv + eta.v + Sex_V * SEX + ALLV * (log(WT / 70)))
        ktr <- exp(lktr + eta.ktr)
        ## RxODE-style differential equations are supported
        d / dt(depot) <- -ktr * depot
        d / dt(gut) <- ktr * depot - ktr * gut
        d / dt(center) <- ktr * gut - (cl / v) * center
        ## Concentration is calculated
        cp <- center / v
        ## And is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }
    nm <- nlmixr(PK106)
    expect_equal(nm$logThetasList, list(1:6, c(1, 2, 3)))
  })

  test_that("PK106 with non mu-referenced exp(lktr)", {
    PK106 <- function() {
      ini({
        ## Where initial conditions/variables are specified
        lktr <- log(1.15) # log k transit (/h)
        lcl <- log(0.135) # log Cl (L/hr)
        lv <- log(8) # log V (L)
        ALLC <- fix(0.75) # allometric exponent cl
        ALLV <- fix(1.00) # allometric exponent v
        Sex_V <- 0.1 # log Sex on v
        prop.err <- 0.15
        add.err <- 0.6
        eta.cl ~ 0.1 # IIV Cl
        eta.v ~ 0.1 # IIV V
      })
      model({
        ## Allometric scaling on weight
        cl <- exp(lcl + eta.cl + ALLC * (log(WT / 70)))
        v <- exp(lv + eta.v + Sex_V * SEX + ALLV * (log(WT / 70)))
        ktr <- exp(lktr)
        ## RxODE-style differential equations are supported
        d / dt(depot) <- -ktr * depot
        d / dt(gut) <- ktr * depot - ktr * gut
        d / dt(center) <- ktr * gut - (cl / v) * center
        ## Concentration is calculated
        cp <- center / v
        ## And is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }
    nm <- nlmixr(PK106)
    expect_equal(nm$logThetasList, list(1:6, 1:3))
  })


  test_that("PK106 with non mu-referenced lktr", {
    PK106 <- function() {
      ini({
        # Where initial conditions/variables are specified
        lktr <- log(1.15) # log k transit (/h)
        lcl <- log(0.135) # log Cl (L/hr)
        lv <- log(8) # log V (L)
        ALLC <- fix(0.75) # allometric exponent cl
        ALLV <- fix(1.00) # allometric exponent v
        Sex_V <- 0.1 # log Sex on v
        prop.err <- 0.15
        add.err <- 0.6
        eta.cl ~ 0.1 # IIV Cl
        eta.v ~ 0.1 # IIV V
      })
      model({
        # Allometric scaling on weight
        cl <- exp(lcl + eta.cl + ALLC * (log(WT / 70)))
        v <- exp(lv + eta.v + Sex_V * SEX + ALLV * (log(WT / 70)))
        ktr <- lktr
        # RxODE-style differential equations are supported
        d / dt(depot) <- -ktr * depot
        d / dt(gut) <- ktr * depot - ktr * gut
        d / dt(center) <- ktr * gut - (cl / v) * center
        ## Concentration is calculated
        cp <- center / v
        # And is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }
    nm <- nlmixr(PK106)
    expect_equal(nm$logThetasList, list(2:6, 2:3))
  })

}, test="cran")
