nlmixrTest({

  one.compartment <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45 # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.err <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err)
    })
  }


  f <- nlmixr(one.compartment)

  testUi <- function(ui, has = NULL, exclude = NULL, values = NULL) {
    if (!is.null(has)) {
      expect_true(all(has %in% paste(ui$ini$name)))
    }
    if (!is.null(values) && !is.null(names(values))) {
      .vals <- setNames(ui$ini$est, paste(ui$ini$name))
      .vals <- .vals[names(values)]
      expect_equal(values, .vals)
    }
    if (!is.null(exclude)) {
      expect_false(any(exclude %in% paste(ui$ini$name)))
    }
  }

  test_that("UI updates work correctly", {
    context("update: Test Base model")
    testUi(
      f, c("tka", "tcl", "tv", "eta.ka", "eta.cl", "eta.v", "add.err"),
      "matt", c(tka = 0.45, tcl = 1, tv = 3.45, eta.ka = 0.6, eta.cl = 0.3, eta.v = 0.1, add.err = 0.7)
    )

    context("update: Multiple component change with c()")
    testUi(
      f %>% update(tka = 4, cl = exp(tcl), ka = exp(tka), c(tcl = 3, tv = 4)),
      c("tka", "tcl", "tv", "eta.v", "add.err"),
      c("eta.ka", "eta.cl"),
      c(tka = 4, tcl = 3, tv = 4, eta.v = 0.1, add.err = 0.7)
    )

    context("update: Multiple component change with list()")

    testUi(
      f %>% update(tka = 4, cl = exp(tcl), ka = exp(tka), list(tcl = 3, tv = 4)),
      c("tka", "tcl", "tv", "eta.v", "add.err"),
      c("eta.ka", "eta.cl"),
      c(tka = 4, tcl = 3, tv = 4, eta.v = 0.1, add.err = 0.7)
    )

    context("update: Multiple component change with assigned .tmp=list()")

    .tmp <- list(tcl = 3, tv = 4)
    .ui <- f %>% update(tka = 4, cl = exp(tcl), ka = exp(tka), .tmp)

    testUi(
      .ui,
      c("tka", "tcl", "tv", "eta.v", "add.err"),
      c("eta.ka", "eta.cl"),
      c(tka = 4, tcl = 3, tv = 4, eta.v = 0.1, add.err = 0.7)
    )

    context("update: Multiple component change with assigned .tmp=c()")

    .tmp <- c(tcl = 3, tv = 4)
    .ui <- f %>% update(tka = 4, cl = exp(tcl), ka = exp(tka), .tmp)

    testUi(
      .ui,
      c("tka", "tcl", "tv", "eta.v", "add.err"),
      c("eta.ka", "eta.cl"),
      c(tka = 4, tcl = 3, tv = 4, eta.v = 0.1, add.err = 0.7)
    )

    context("update: Multiple component change with assigned .tmp={}")

    .tmp <- quote({
      ka <- exp(tka)
    })
    .ui <- f %>% update(tka = 4, cl = exp(tcl), .tmp, c(tcl = 3, tv = 4))

    testUi(
      .ui,
      c("tka", "tcl", "tv", "eta.v", "add.err"),
      c("eta.ka", "eta.cl"),
      c(tka = 4, tcl = 3, tv = 4, eta.v = 0.1, add.err = 0.7)
    )


    testUi(
      f %>% update(
        tka = 4,
        cl = exp(tcl),
        {
          ka <- exp(tka)
        },
        c(tcl = 3, tv = 4)
      ),
      c("tka", "tcl", "tv", "eta.v", "add.err"),
      c("eta.ka", "eta.cl"),
      c(tka = 4, tcl = 3, tv = 4, eta.v = 0.1, add.err = 0.7)
    )

    testUi(
      f %>% update(ka = exp(tka)),
      c("tka", "tcl", "tv", "eta.cl", "eta.v", "add.err"),
      "eta.ka", c(tka = 0.45, tcl = 1, tv = 3.45, eta.cl = 0.3, eta.v = 0.1, add.err = 0.7)
    )

    ## Now test linCmt() issue #166
    one.cmt <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45 # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.err <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.err)
      })
    }

    .ui <- one.cmt %>% update({
      linCmt() ~ add(add.err) + prop(prop.err)
    })
    expect_true(inherits(.ui, "nlmixrUI"))
  })

  context("piping looks through parent environments")

  test_that("Looks through prior frames for the correct object", {
    fit <- nlmixr(one.compartment)
    fits <- lapply(seq(-1, -0.1, 0.1), function(kainit) {
      nlmixr(update(fit, tka = kainit))
    })

    expect_true(inherits(fits, "list"))

    expect_error(lapply(seq(-1, -0.1, 0.1), function(kainit) {
      nlmixr(update(fit, tka = matt))
    }), "object 'matt' not found")
  })

  f <- nlmixr(one.compartment)

  test_that("piping works for correlations #1", {
    testUi(f %>% ini(eta.ka + eta.cl ~ c(
      0.2,
      0.01, 0.2
    )),
    has = c("tka", "tcl", "tv", "eta.ka", "eta.cl", "eta.v", "add.err", "(eta.cl,eta.ka)"),
    exclude = "matt",
    values = c(
      tka = 0.45, tcl = 1, tv = 3.45, eta.ka = 0.2, eta.cl = 0.2, eta.v = 0.1, add.err = 0.7,
      `(eta.cl,eta.ka)` = 0.01
    )
    )
  })

  test_that("piping works for correlations #2", {
    expect_error(f %>% ini(eta.ka + eta.matt ~ c(
      0.2,
      0.01, 0.2
    )))
  })


  test_that("piping works for correlations #3", {
    testUi(
      f %>% update(eta.ka + eta.cl ~ c(
        0.2,
        0.01, 0.2
      )),
      c("tka", "tcl", "tv", "eta.ka", "eta.cl", "eta.v", "add.err", "(eta.cl,eta.ka)"),
      "matt", c(
        tka = 0.45, tcl = 1, tv = 3.45, eta.ka = 0.2, eta.cl = 0.2, eta.v = 0.1, add.err = 0.7,
        `(eta.cl,eta.ka)` = 0.01
      )
    )
  })

  test_that("piping works for correlations #4", {
    expect_error(f %>% update(eta.ka + eta.matt ~ c(
      0.2,
      0.01, 0.2
    )))
  })

}, test="cran")
