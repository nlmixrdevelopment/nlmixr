## From https://raw.githubusercontent.com/bbolker/broom.mixed/master/tests/testthat/helper-checkers.R

##' test the basics of tidy/augment/glance output: is a data frame, no row names
check_tidiness <- function(o) {
  testthat::expect_is(o, "tbl_df")
  testthat::expect_equal(rownames(o), as.character(seq_len(nrow(o))))
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


#' check the output of an augment function
check_augment <- function(au, original = NULL, exp.names = NULL,
                          same = NULL) {
  check_tidiness(au)

  if (!is.null(original)) {
    # check that all rows in original appear in output
    testthat::expect_equal(nrow(au), nrow(original))
    # check that columns are the same
    for (column in same) {
      testthat::expect_equal(au[[column]], original[[column]])
    }
  }

  if (!is.null(exp.names)) {
    testthat::expect_true(all(exp.names %in% colnames(au)))
  }
}


#' add NAs to a vector randomly
#'
#' @param v vector to add NAs to
#' @param number number of NAs to add
#'
#' @return vector with NAs added randomly
add_NAs <- function(v, number) {
  if (number >= length(v)) {
    stop("Would replace all or more values with NA")
  }
  v[sample(length(v), number)] <- NA
  v
}

#' check an augmentation function works as expected when given NAs
#'
#' @param func A modeling function that takes a dataset and additional
#' arguments, including na.action
#' @param .data dataset to test function on; must have at least 3 rows
#' and no NA values
#' @param column a column included in the model to be replaced with NULLs
#' @param column2 another column in the model; optional
#' @param ... extra arguments, not used
#'
#' @export


library(testthat)
library(nlmixr)

options(nlmixr.save=TRUE,
        nlmixr.save.dir=system.file(package="nlmixr"));

one.compartment <- function() {
    ini({
        tka <- .5   # Log Ka
        tcl <- -3.2 # Log Cl
        tv <- -0.6    # Log V
        eta.ka ~ 0.5
        eta.cl ~ 0.5
        eta.v ~ 0.5
        add.err <- 0.1
    })
    model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.err)
    })
}

fitS <- nlmixr(one.compartment, theo_sd, est="saem")
context("broom nlmixr SAEM")

test_that("tidy works on nlmixr fit SAEM fits", {
    td <- tidy(fitS)
    ## FIXME: fails if lmerTest has been loaded previously ...
    check_tidy(td,7,5,c("effect", "group", "term", "estimate", "std.error"));
    expect_equal(
        td$term,
        c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err")
    )
    td <- tidy(fitS, conf.level=0.9)
    check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error",
          "conf.low", "conf.high"))
    expect_equal(
        td$term,
        c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
          "add.err")
    )
    expect_equal(td$estimate,c(1.57250316684165, 2.75697373309734, 31.5119000868981, 0.641485375817318,
                               0.271453624934398, 0.133398409963236, 0.694922343196694))
    expect_equal(td$std.error,c(0.307124574041779, 0.235554957083684, 1.4343040419653, NA,
                                NA, NA, NA))
    expect_equal(td$conf.low, c(1.14043917332581, 2.39551403203811, 29.2388311287432, NA, NA,
                                NA, NA))

    td <- tidy(fitS, conf.level=0.9,exponentiate=FALSE)
    check_tidy(td)
    expect_equal(td$estimate,c(0.452668723480193, 1.01413360464674, 3.45036525502693, 0.641485375817318,
0.271453624934398, 0.133398409963236, 0.694922343196694))
    expect_equal(td$std.error, c(0.195309351687116, 0.0854396812910683, 0.0455162664901202,
NA, NA, NA, NA),tolerance=1e-5)
    expect_equal(td$conf.low, c(0.1314134279801, 0.873597834989551, 3.37549765900537, NA, NA,
                                NA, NA))

    for (ef in c("ran_vals", "random")){

        td  <- tidy(fitS, effects=ef)
        td1  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td  <- tidy(fitS, effects=ef, exponentiate=FALSE)
        td2  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td  <- tidy(fitS, effects=ef, exponentiate=TRUE)
        td3  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(td1, td2)
        expect_equal(td2, td3)
    }

    td  <- tidy(fitS, effects="ran_coef")
    td1  <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(td1, c(1.74706183526116, 1.95318511690162, 2.27639703775571, 1.19660268941229,
                        1.51057755792184, 1.0617571207384, 0.709150547688647, 1.30594405279382,
                        6.60135556415273, 0.763590560320776, 3.42355257477486, 0.936706165959318,
                        1.69243972945261, 3.17304918825295, 2.84360137604056, 2.71285514667348,
                        2.36741811727645, 4.02012422853366, 3.26362929310098, 3.26320683211502,
                        2.86177367480114, 1.8436536629697, 3.66043526991176, 2.42633680279661,
                        29.0361127660822, 32.0963402623569, 33.3988129997684, 31.3322798356558,
                        27.3827701153425, 38.2631456108251, 32.9216569074101, 34.5389598620118,
                        32.0015076568584, 26.8368511027911, 36.296341502626, 26.0163203497275
                        ))

    td  <- tidy(fitS, effects="ran_coef", exponentiate=FALSE)
    td2  <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    td  <- tidy(fitS, effects="ran_coef", exponentiate=TRUE)
    td3  <- td$estimate
    check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

    expect_equal(log(td1),td2)
    expect_equal(td2,log(td3))

    td  <- tidy(fitS, effects="ran_pars")
    td1  <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
    expect_equal(td1, c(0.641485375817318, 0.271453624934398, 0.133398409963236, 0.694922343196694))

    td  <- tidy(fitS, effects="ran_pars", exponentiate=FALSE)
    td2  <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    td  <- tidy(fitS, effects="ran_pars", exponentiate=TRUE)
    td3  <- td$estimate
    check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

    expect_equal(td1, td2)
    expect_equal(td2, td3)
  })



## fitF <- nlmixr(one.compartment, theo_sd, est="focei")
## fitN <- nlmixr(one.compartment, theo_sd, est="nlme")
