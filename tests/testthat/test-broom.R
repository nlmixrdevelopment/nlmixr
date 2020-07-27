nlmixrTest({
  tol <- 1e-5
  ## From https://raw.githubusercontent.com/bbolker/broom.mixed/master/tests/testthat/helper-checkers.R

  ##' test the basics of tidy/augment/glance output: is a data frame, no row names
  check_tidiness <- function(o) {
    testthat::expect_is(o, "tbl_df")
    testthat::expect_equal(rownames(o), as.character(seq_len(nrow(o))))
  }

  .nlmixr <- function(...) {
    suppressWarnings(nlmixr(...))
  }


  #' check the output of a tidy function
  check_tidy <- function(o, exp.row = NULL, exp.col = NULL, exp.names = NULL) {
    check_tidiness(o)

    if (!is.null(exp.row)) {
      testthat::expect_equal(nrow(o), exp.row)
    }
    if (!is.null(exp.col)) {
      testthat::expect_equal(ncol(o), exp.col)
    }
    if (!is.null(exp.names)) {
      testthat::expect_true(all(exp.names %in% colnames(o)))
    }
  }

  library(testthat)

  options(
    nlmixr.save = TRUE,
    nlmixr.save.dir = system.file(package = "nlmixr")
  )

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


  context("broom tidy nlmixr SAEM")

  fitS <- .nlmixr(one.compartment, theo_sd, est = "saem")
  test_that("tidy works on nlmixr fit SAEM fits", {
    td <- tidy(fitS, exponentiate = NA)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
      )
    )
    td <- tidy(fitS, conf.level = 0.9, exponentiate = NA)
    check_tidy(td, 7, 9, c(
      "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
      "conf.low", "conf.high"
    ))
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
      )
    )
    .est <- td$estimate
    .stdErr <- td$std.error
    .confLow <- td$conf.low

    td <- tidy(fitS, conf.level = 0.9, exponentiate = FALSE)
    check_tidy(td)
    expect_equal(td$estimate[1:3], log(.est)[1:3],
      tolerance = tol
    )
    expect_equal(td$estimate[-(1:3)], .est[-(1:3)],
      tolerance = tol
    )
    ## exp(.df$model.est[.exp])*.df$std.error[.exp]
    expect_equal(exp(td$estimate[1:3]) * td$std.error[1:3], .stdErr[1:3],
      tolerance = tol
    )
    expect_equal(exp(td$conf.low), .confLow,
      tolerance = tol
    )

    for (ef in c("ran_vals", "random")) {
      td <- tidy(fitS, effects = ef, exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- tidy(fitS, effects = ef, exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- tidy(fitS, effects = ef, exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      expect_equal(td1, td2)
      expect_equal(td2, td3)
    }

    td <- tidy(fitS, effects = "ran_coef", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    .td1 <- td1

    td <- tidy(fitS, effects = "ran_coef", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    td <- tidy(fitS, effects = "ran_coef", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(log(td1), td2)
    expect_equal(td2, log(td3))

    td <- tidy(fitS, effects = "ran_pars", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- tidy(fitS, effects = "ran_pars", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- tidy(fitS, effects = "ran_pars", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    expect_equal(td1, td2)
    expect_equal(td2, td3)
  })

  for (f in c("focei", "foce")) {
    context(sprintf("broom tidy nlmixr %s", f))
    fitF <- .nlmixr(one.compartment, theo_sd, est = f)
    test_that(sprintf("tidy works on nlmixr fit %s fits", f), {
      td <- tidy(fitF, exponentiate = NA)
      check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      td <- tidy(fitF, conf.level = 0.9, exponentiate = NA)
      check_tidy(td, 7, 9, c(
        "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
        "conf.low", "conf.high"
      ))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      .est <- td$estimate
      .stdErr <- td$std.error
      .confLow <- td$conf.low

      td <- tidy(fitF, conf.level = 0.9, exponentiate = FALSE)
      check_tidy(td)
      expect_equal(td$estimate[1:3], log(.est)[1:3],
        tolerance = tol
      )
      expect_equal(td$estimate[-(1:3)], .est[-(1:3)],
        tolerance = tol
      )
      expect_equal(exp(td$estimate[1:3]) * td$std.error[1:3], .stdErr[1:3],
        tolerance = tol
      )
      expect_equal(exp(td$conf.low), .confLow,
        tolerance = tol
      )

      for (ef in c("ran_vals", "random")) {
        td <- tidy(fitF, effects = ef, exponentiate = NA)
        td1 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- tidy(fitF, effects = ef, exponentiate = FALSE)
        td2 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- tidy(fitF, effects = ef, exponentiate = TRUE)
        td3 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        expect_equal(td1, td2)
        expect_equal(td2, td3)
      }
      ##
      td <- tidy(fitF, effects = "ran_coef", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_coef", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_coef", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      expect_equal(log(td1), td2)
      expect_equal(td2, log(td3))
      ##
      td <- tidy(fitF, effects = "ran_pars", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_pars", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_pars", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    })
  }

  for (f in c("foi", "fo")) {
    context(sprintf("broom tidy nlmixr %s", f))
    fitF <- .nlmixr(one.compartment, theo_sd, est = f)
    test_that(sprintf("tidy works on nlmixr fit %s fits", f), {
      td <- tidy(fitF, exponentiate = NA)
      check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      td <- tidy(fitF, conf.level = 0.9, exponentiate = NA)
      check_tidy(td, 7, 9, c(
        "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
        "conf.low", "conf.high"
      ))
      expect_equal(
        td$term,
        c(
          "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err"
        )
      )
      .est <- td$estimate
      .stdErr <- td$std.error
      .confLow <- td$conf.low

      td <- tidy(fitF, conf.level = 0.9, exponentiate = FALSE)
      check_tidy(td)
      expect_equal(td$estimate[1:3], log(.est)[1:3],
        tolerance = tol
      )
      expect_equal(td$estimate[-(1:3)], .est[-(1:3)],
        tolerance = tol
      )
      ## exp(.df$model.est[.exp])*.df$std.error[.exp]
      expect_equal(exp(td$estimate[1:3]) * td$std.error[1:3], .stdErr[1:3],
        tolerance = tol
      )
      expect_equal(exp(td$conf.low), .confLow,
        tolerance = tol
      )
      for (ef in c("ran_vals", "random")) {
        td <- tidy(fitF, effects = ef, exponentiate = NA)
        td1 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- tidy(fitF, effects = ef, exponentiate = FALSE)
        td2 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        td <- tidy(fitF, effects = ef, exponentiate = TRUE)
        td3 <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
        ##
        expect_equal(td1, td2, tolerance = tol)
        expect_equal(td2, td3, tolerance = tol)
      }
      ##
      td <- tidy(fitF, effects = "ran_coef", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_coef", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_coef", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
      ##
      expect_equal(log(td1), td2, tolerance = tol)
      expect_equal(td2, log(td3), tolerance = tol)
      ##
      td <- tidy(fitF, effects = "ran_pars", exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_pars", exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      td <- tidy(fitF, effects = "ran_pars", exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
      ##
      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    })
  }


  context("broom tidy nlmixr nlme")
  fitN <- .nlmixr(one.compartment, theo_sd, est = "nlme", control = nlmeControl(pnlsTol = 0.6))
  test_that("tidy works on nlmixr fit nlme fits", {
    td <- tidy(fitN, exponentiate = NA)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
      )
    )
    td <- tidy(fitN, conf.level = 0.9, exponentiate = NA)
    check_tidy(td, 7, 9, c(
      "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
      "conf.low", "conf.high"
    ))
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
      )
    )
    .est <- td$estimate
    .stdErr <- td$std.error
    .confLow <- td$conf.low

    td <- tidy(fitN, conf.level = 0.9, exponentiate = FALSE)
    check_tidy(td)
    expect_equal(td$estimate[1:3], log(.est)[1:3],
      tolerance = tol
    )
    expect_equal(td$estimate[-(1:3)], .est[-(1:3)],
      tolerance = tol
    )
    ## exp(.df$model.est[.exp])*.df$std.error[.exp]
    expect_equal(exp(td$estimate[1:3]) * td$std.error[1:3], .stdErr[1:3],
      tolerance = tol
    )
    expect_equal(exp(td$conf.low), .confLow,
      tolerance = tol
    )

    for (ef in c("ran_vals", "random")) {
      td <- tidy(fitN, effects = ef, exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- tidy(fitN, effects = ef, exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- tidy(fitN, effects = ef, exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      expect_equal(td1, td2)
      expect_equal(td2, td3)
    }

    td <- tidy(fitN, effects = "ran_coef", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    td <- tidy(fitN, effects = "ran_coef", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    td <- tidy(fitN, effects = "ran_coef", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(log(td1), td2, tolerance = tol)
    expect_equal(td2, log(td3), tolerance = tol)

    td <- tidy(fitN, effects = "ran_pars", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- tidy(fitN, effects = "ran_pars", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- tidy(fitN, effects = "ran_pars", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    expect_equal(td1, td2, tolerance = tol)
    expect_equal(td2, td3, tolerance = tol)
  })


  context("broom tidy nlmixr posthoc")
  fitP <- .nlmixr(one.compartment, theo_sd, est = "posthoc")

  test_that("tidy works on posthoc fit fits", {
    td <- tidy(fitP, exponentiate = NA)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error", "statistic", "p.value"))
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
      )
    )
    td <- tidy(fitP, conf.level = 0.9, exponentiate = NA)
    check_tidy(td, 7, 9, c(
      "effect", "group", "term", "estimate", "std.error", "statistic", "p.value",
      "conf.low", "conf.high"
    ))
    expect_equal(
      td$term,
      c(
        "tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
        "add.err"
      )
    )
    expect_equal(td$estimate, c(
      1.56831218549017, 2.71828182845905, 31.5003923087479, 0.774596669241483,
      0.547722557505166, 0.316227766016838, 0.7
    ),
    tolerance = tol
    )
    expect_equal(td$std.error, c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_
    ))
    expect_equal(td$conf.low, c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_
    ))

    td <- tidy(fitP, conf.level = 0.9, exponentiate = FALSE)
    check_tidy(td)
    expect_equal(td$estimate, c(
      0.45, 1, 3.45, 0.774596669241483, 0.547722557505166, 0.316227766016838,
      0.7
    ),
    tolerance = tol
    )
    expect_equal(td$std.error, c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_
    ),
    tolerance = tol
    )
    expect_equal(td$conf.low, c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_
    ),
    tolerance = tol
    )

    for (ef in c("ran_vals", "random")) {
      td <- tidy(fitP, effects = ef, exponentiate = NA)
      td1 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- tidy(fitP, effects = ef, exponentiate = FALSE)
      td2 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      td <- tidy(fitP, effects = ef, exponentiate = TRUE)
      td3 <- td$estimate
      check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

      expect_equal(td1, td2, tolerance = tol)
      expect_equal(td2, td3, tolerance = tol)
    }

    td <- tidy(fitP, effects = "ran_coef", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(td1, c(
      1.75611262206464, 1.92886734474446, 2.36872361414073, 1.18776870402381,
      1.48421770349779, 1.1603466406079, 0.728548560764397, 1.37496750947703,
      6.63504345280614, 0.728113370661147, 3.56192283422217, 0.87667797314966,
      1.62338691901516, 3.22874000491057, 2.81088025117325, 2.7068918775575,
      2.37701456492688, 4.04231470515951, 3.23451222726721, 3.26106788182658,
      2.88069603771061, 1.86803407932834, 3.73237924876953, 2.49388025078765,
      29.2305913871308, 31.8379070004274, 33.9077584685152, 31.2454399174589,
      27.0633487903421, 40.6859508443866, 33.6487106457298, 35.5118803877647,
      31.940988239621, 26.0728299589223, 37.2139995031868, 24.7515885649664
    ),
    tolerance = tol
    )

    td <- tidy(fitP, effects = "ran_coef", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    td <- tidy(fitP, effects = "ran_coef", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(log(td1), td2, tolerance = tol)
    expect_equal(td2, log(td3), tolerance = tol)

    td <- tidy(fitP, effects = "ran_pars", exponentiate = NA)
    td1 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
    expect_equal(td1, c(0.774596669241483, 0.547722557505166, 0.316227766016838, 0.7),
      tolerance = tol
    )

    td <- tidy(fitP, effects = "ran_pars", exponentiate = FALSE)
    td2 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td <- tidy(fitP, effects = "ran_pars", exponentiate = TRUE)
    td3 <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    expect_equal(td1, td2, tolerance = tol)
    expect_equal(td2, td3, tolerance = tol)
  })
}, test="broom")
