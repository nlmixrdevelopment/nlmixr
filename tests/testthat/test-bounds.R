context("Test bounds extraction")
rxPermissive({

    ref <- structure(list(ntheta = c(1, 2, 3, 4, 5, 6, 7, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9, 10, 11, 12, 13, 14), neta1 = c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12, 12, NA, NA, NA, NA, NA, NA), neta2 = c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 2, 3, 4, 5, 5, 6, 7, 7, 8, 7, 8, 9, 10, 10, 11, 10, 11, 12, NA, NA, NA, NA, NA, NA), name = structure(c(1L, 2L, 3L, 4L, NA, NA, NA, NA, 5L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 6L, 7L, 8L, 9L, 10L, 11L, NA, NA, NA, 12L, 13L, 14L), .Label = c("a", "b", "c", "d", "et1", "et2", "(et3,et2)", "et3", "(et4,et2)", "(et4,et3)", "et4", "a5", "a6", "a7"), class = "factor"), lower = c(0, 0, -Inf, -Inf, 0, 0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 9, 11, -Inf, 9, 11), est = c(1, 3, 4, 4, 1, 1, 1, 1, 10, 20, 30, 40, 40, 0.1, 20, 40, 0.1, 20, 0.1, 0.1, 30, 40, 0.1, 20, 0.1, 0.1, 30, 8, 10, 12, 8, 10, 12), upper = c(2, Inf, Inf, Inf, 2, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 13, Inf, Inf, 13), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), err = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), label = c("A", NA, NA, NA, "e", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "labels", NA, NA, NA), condition = c(NA, NA, NA, NA, NA, NA, NA, NA, "ID", "ID", "ID", "ID", "STUD", "STUD", "STUD", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", NA, NA, NA, NA, NA, NA)), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -33L), class = c("nlmixrBounds", "data.frame"));

    testbounds <- function(){
        a = c(0, 1, 2) # A
        b = c(0, 3);
        c <- 4
        d <- c(4);
        c(0, 1, 2) # e
        c(0, 1)
        c(1)
        1
        et1 ~ 10
        ~ 20
        ~ 30
        ~ c(40)
        ~ c(40,
            0.1, 20) | STUD;
        ~ c(40,
            0.1, 20,
            0.1, 0.1, 30);
        et2 + et3 + et4 ~ c(40,
                            0.1, 20,
                            0.1, 0.1, 30)

        ## new test fixed parameters...
        c(8, fixed)
        c(9, 10, fixed)
        c(11, 12, 13, fixed) # labels

        a5 = c(8, fixed)
        a6 = c(9, 10, fixed)
        a7 = c(11, 12, 13, fixed)
    }

    test_that("bounds are extracted correctly", {
        expect_equal(nlmixrBounds(testbounds), ref);
    })

    bnd <- function(){
        c(1, 2, 3, 4, 5)
    }

    bnda <- function(){
        a = c(1, 2, 3, 4, 5)
    }
    bndb <- function(){
        a <- c(1, 2, 3, 4, 5)
    }

    test_that("Theta Bounds above 5 don't work", {
        expect_error(nlmixrBounds(bnd), rex::rex("c(1, 2, 3, 4, 5) syntax is not supported for thetas") )
        expect_error(nlmixrBounds(bnda), rex::rex("a = c(1, 2, 3, 4, 5) syntax is not supported for thetas"))
        expect_error(nlmixrBounds(bndb), rex::rex("a <- c(1, 2, 3, 4, 5) syntax is not supported for thetas"))
    })

    bnd <- function(){
        c(1, 2, 3, 4)
    }
    bnda <- function(){
        a = c(1, 2, 3, 4)
    }
    bndb <- function(){
        a <- c(1, 2, 3, 4)
    }
    test_that("Theta Bounds above 4 don't work", {
        expect_error(nlmixrBounds(bnd), rex::rex("c(1, 2, 3, 4) syntax is not supported for thetas") )
        expect_error(nlmixrBounds(bnda), rex::rex("a = c(1, 2, 3, 4) syntax is not supported for thetas"))
        expect_error(nlmixrBounds(bndb), rex::rex("a <- c(1, 2, 3, 4) syntax is not supported for thetas"))
    })

    bnd1 <- function(){
        ~ c(1)
    }
    ref1 <- structure(list(ntheta = NA, neta1 = 1, neta2 = 1, name = NA,     lower = -Inf, est = 1, upper = Inf, fix = FALSE, err = NA,     label = NA, condition = structure(1L, .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -1L), class = c("nlmixrBounds", "data.frame"))

    bnd2 <- function(){
        ~ c(1, 2)
    }

    bnd3 <- function(){
        ~ c(1,
            2, 3)
    }

    ref3 <- structure(list(ntheta = c(NA, NA, NA), neta1 = c(1, 2, 2), neta2 = c(1, 1, 2), name = c(NA, NA, NA), lower = c(-Inf, -Inf, -Inf), est = c(1, 2, 3), upper = c(Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE),     err = c(NA, NA, NA), label = c(NA, NA, NA), condition = structure(c(1L,     1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -3L), class = c("nlmixrBounds", "data.frame"))

    bnd4  <- function(){
        ~ c(1,
            2, 3,
            4)
    }

    bnd5  <- function(){
        ~ c(1,
            2, 3,
            4, 5)
    }

    bnd6  <- function(){
        ~ c(1,
            2, 3,
            4, 5, 6)
    }

    ref6 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(1, 2, 3, 4, 5, 6), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))


    bnd7  <- function(){
        ~ c(1,
            2, 3,
            4, 5, 6,
            7)
    }
    bnd8  <- function(){
        ~ c(1,
            2, 3,
            4, 5, 6,
            7, 8)
    }
    bnd9  <- function(){
        ~ c(1,
            2, 3,
            4, 5, 6,
            7, 8, 9)
    }
    bnd10  <- function(){
        ~ c(1,
            2, 3,
            4, 5, 6,
            7, 8, 9, 10)
    }

    ref10 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4), neta2 = c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4), name = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -10L), class = c("nlmixrBounds", "data.frame"))


    test_that("Bad Lower trianglar matrices throw errors.", {
        expect_equal(nlmixrBounds(bnd1), ref1);
        expect_error(nlmixrBounds(bnd2), rex::rex("~c(1, 2) does not have the right dimensions for a lower triangular matrix."))
        expect_equal(nlmixrBounds(bnd3), ref3);
        expect_error(nlmixrBounds(bnd4), rex::rex("~c(1, 2, 3, 4) does not have the right dimensions for a lower triangular matrix."))
        expect_error(nlmixrBounds(bnd5), rex::rex("~c(1, 2, 3, 4, 5) does not have the right dimensions for a lower triangular matrix."))
        expect_equal(nlmixrBounds(bnd6), ref6);
        expect_error(nlmixrBounds(bnd7), rex::rex("~c(1, 2, 3, 4, 5, 6, 7) does not have the right dimensions for a lower triangular matrix."))
        expect_error(nlmixrBounds(bnd8), rex::rex("~c(1, 2, 3, 4, 5, 6, 7, 8) does not have the right dimensions for a lower triangular matrix."))
        expect_error(nlmixrBounds(bnd9), rex::rex("~c(1, 2, 3, 4, 5, 6, 7, 8, 9) does not have the right dimensions for a lower triangular matrix."))
        expect_equal(nlmixrBounds(bnd10), ref10)
    })


    bnd1 <- function(){
        eta1 ~ c(1)
    }

    ref1 <- structure(list(ntheta = NA, neta1 = 1, neta2 = 1, name = structure(1L, .Label = "eta1", class = "factor"),     lower = -Inf, est = 1, upper = Inf, fix = FALSE, err = NA,     label = NA, condition = structure(1L, .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -1L), class = c("nlmixrBounds", "data.frame"))

    bnd2 <- function(){
        eta1 ~ c(1, 2)
    }

    bnd3 <- function(){
        eta1 + eta2~ c(1,
                       2, 3)
    }

    ref3 <- structure(list(ntheta = c(NA, NA, NA), neta1 = c(1, 2, 2), neta2 = c(1, 1, 2), name = structure(1:3, .Label = c("eta1", "(eta2,eta1)", "eta2"), class = "factor"), lower = c(-Inf, -Inf, -Inf), est = c(1, 2, 3), upper = c(Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE),     err = c(NA, NA, NA), label = c(NA, NA, NA), condition = structure(c(1L,     1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -3L), class = c("nlmixrBounds", "data.frame"));

    bnd4  <- function(){
        eta1 + eta2~ c(1,
                       2, 3,
                       4)
    }

    bnd5  <- function(){
        eta1 + eta2 ~ c(1,
                        2, 3,
                        4, 5)
    }

    bnd6  <- function(){
        eta1 + eta2 + eta3 ~ c(1,
                               2, 3,
                               4, 5, 6)
    }

    ref6 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(1,     2, 3, 4, 5, 6), upper = c(Inf, Inf, Inf, Inf, Inf, Inf),     fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd7  <- function(){
        eta1 + eta2 + eta3 ~ c(1,
                               2, 3,
                               4, 5, 6,
                               7)
    }
    bnd8  <- function(){
        eta1 + eta2 + eta3 ~ c(1,
                               2, 3,
                               4, 5, 6,
                               7, 8)
    }
    bnd9  <- function(){
        eta1 + eta2 + eta3 ~ c(1,
                               2, 3,
                               4, 5, 6,
                               7, 8, 9)
    }
    bnd10  <- function(){
        eta1 + eta2 + eta3 + eta4 ~ c(1,
                                      2, 3,
                                      4, 5, 6,
                                      7, 8, 9, 10)
    }

    ref10 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4), neta2 = c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4), name = structure(1:10, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3", "(eta4,eta1)", "(eta4,eta2)", "(eta4,eta3)", "eta4"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,     -Inf, -Inf), est = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), upper = c(Inf,     Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE,     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE    ), err = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), label = c(NA,     NA, NA, NA, NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -10L), class = c("nlmixrBounds", "data.frame"))


    test_that("Bad Lower trianglar matrices (with labels) throw errors.", {
        expect_equal(nlmixrBounds(bnd1), ref1);
        expect_error(nlmixrBounds(bnd2), rex::rex("eta1 ~ c(1, 2) does not have the right dimensions for a lower triangular matrix."))
        expect_equal(nlmixrBounds(bnd3), ref3);
        expect_error(nlmixrBounds(bnd4), rex::rex("eta1 + eta2 ~ c(1, 2, 3, 4) does not have the right dimensions for a lower triangular matrix."))
        expect_error(nlmixrBounds(bnd5), rex::rex("eta1 + eta2 ~ c(1, 2, 3, 4, 5) does not have the right dimensions for a lower triangular matrix."))
        expect_equal(nlmixrBounds(bnd6), ref6);
        expect_error(nlmixrBounds(bnd7), rex::rex("eta1 + eta2 + eta3 ~ c(1, 2, 3, 4, 5, 6, 7) does not have the right dimensions for a lower triangular matrix."))
        expect_error(nlmixrBounds(bnd8), rex::rex("eta1 + eta2 + eta3 ~ c(1, 2, 3, 4, 5, 6, 7, 8) does not have the right dimensions for a lower triangular matrix."))
        expect_error(nlmixrBounds(bnd9), rex::rex("eta1 + eta2 + eta3 ~ c(1, 2, 3, 4, 5, 6, 7, 8, 9) does not have the right dimensions for a lower triangular matrix."))
        expect_equal(nlmixrBounds(bnd10), ref10)
    })

    bnd10.a  <- function(){
        eta1 + eta2 + eta3 + eta4 + eta5~ c(1,
                                            2, 3,
                                            4, 5, 6,
                                            7, 8, 9, 10)
    }

    bnd10.b  <- function(){
        eta1 + eta2 + eta3~ c(1,
                              2, 3,
                              4, 5, 6,
                              7, 8, 9, 10)
    }

    test_that("Number of eta variables must match", {
        expect_error(nlmixrBounds(bnd10.a), rex::rex("The left handed side of the expression must match the number of ETAs in the lower triangular matrix."))
        expect_error(nlmixrBounds(bnd10.b), rex::rex("The left handed side of the expression must match the number of ETAs in the lower triangular matrix."))
    })

    bnd3 <- function(){
        eta1 + eta2~ c(1, # Comment here.
                       2, 3)
    }

    test_that("Comments inside bounds are not supported!", {
        expect_error(nlmixrBounds(bnd3), rex::rex("Error parsing bounds; Perhaps there is an (unsupported) comment/condition inside the bounds themselves."))
    })

    bnd1 <- function(){
        eta0 ~ 0.3
        eta1 + eta2~ c(1,
                       2, 3) | STUD
        ~ c(1,
            2, 3)
    }

    ref <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 3, 3, 4, 5, 5), neta2 = c(1, 2, 2, 3, 4, 4, 5), name = structure(c(1L, 2L, 3L, 4L, NA, NA, NA), .Label = c("eta0", "eta1", "(eta2,eta1)", "eta2"), class = "factor"), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(0.3, 1, 2, 3, 1, 2, 3), upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA, NA), condition = c("ID", "STUD", "STUD", "STUD", "ID", "ID", "ID")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -7L), class = c("nlmixrBounds", "data.frame"))

    test_that("Conditional statments are supported correctly.", {
        expect_equal(nlmixrBounds(bnd1), ref)
    })

    ## Now test err/pred type parsing.

    bnd1 <- function(){
        err ~ add(0.1)
    }

    ref1 <- structure(list(ntheta = 1, neta1 = NA, neta2 = NA, name = structure(1L, .Label = "err", class = "factor"),     lower = 0, est = 0.1, upper = Inf, fix = TRUE, err = structure(1L, .Label = "add", class = "factor"),     label = NA, condition = NA), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -1L), class = c("nlmixrBounds", "data.frame"))

    bnd2 <- function(){
        err ~ add(0.1) + prop(0.1)
    }

    ref2 <- structure(list(ntheta = c(1, 2), neta1 = c(NA, NA), neta2 = c(NA, NA), name = structure(c(1L, 1L), .Label = "err", class = "factor"),     lower = c(0, 0), est = c(0.1, 0.1), upper = c(Inf, Inf),     fix = c(TRUE, TRUE), err = structure(1:2, .Label = c("add",     "prop"), class = "factor"), label = c(NA, NA), condition = c(NA,     NA)), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -2L), class = c("nlmixrBounds", "data.frame"))

    bnd3 <- function(){
        err ~ prop(0.1) + add(0.1)
    }

    ref3 <- structure(list(ntheta = c(1, 2), neta1 = c(NA, NA), neta2 = c(NA, NA), name = structure(c(1L, 1L), .Label = "err", class = "factor"),     lower = c(0, 0), est = c(0.1, 0.1), upper = c(Inf, Inf),     fix = c(TRUE, TRUE), err = structure(1:2, .Label = c("prop",     "add"), class = "factor"), label = c(NA, NA), condition = c(NA,     NA)), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -2L), class = c("nlmixrBounds", "data.frame"));

    bnd4 <- function(){
        err ~ prop(0.1)
    }

    ref4 <-structure(list(ntheta = 1, neta1 = NA, neta2 = NA, name = structure(1L, .Label = "err", class = "factor"),     lower = 0, est = 0.1, upper = Inf, fix = TRUE, err = structure(1L, .Label = "prop", class = "factor"),     label = NA, condition = NA), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -1L), class = c("nlmixrBounds", "data.frame"))


    bnd5 <- function(){
        err ~ prop(0.1) + add(0.2) | NVS == 1
        err ~ add(0.3) + prop(0.4) | NVS == 2
        err ~ add(0.5)   | NVS == 3
    }

    ref5 <- structure(list(ntheta = c(1, 2, 3, 4, 5), neta1 = c(NA, NA, NA, NA, NA), neta2 = c(NA, NA, NA, NA, NA), name = structure(c(1L, 1L, 1L, 1L, 1L), .Label = "err", class = "factor"), lower = c(0, 0, 0, 0, 0), est = c(0.1, 0.2, 0.3, 0.4, 0.5), upper = c(Inf, Inf, Inf, Inf, Inf), fix = c(TRUE, TRUE, TRUE, TRUE, TRUE), err = structure(c(1L, 2L, 2L, 1L, 2L), .Label = c("prop", "add"), class = "factor"),     label = c(NA, NA, NA, NA, NA), condition = c("NVS == 1",     "NVS == 1", "NVS == 2", "NVS == 2", "NVS == 3")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -5L), class = c("nlmixrBounds", "data.frame"));

    test_that("Error parsing is reasonable", {
        expect_equal(nlmixrBounds(bnd1), ref1)
        ## expect_equal(nlmixrBounds(bnd2), ref2) ## add+prop is taken care of in ui now
        ## expect_equal(nlmixrBounds(bnd3), ref3) ## add prop is taken care of in ui now
        expect_equal(nlmixrBounds(bnd4), ref4)
        ## expect_equal(nlmixrBounds(bnd5), ref5)
    })


    ## Now try a fixed parameter block

    ref1 <- structure(list(ntheta = c(1, 2, 3, 4), neta1 = c(NA, NA, NA, NA), neta2 = c(NA, NA, NA, NA), name = structure(1:4, .Label = c("a", "b", "c", "d"), class = "factor"), lower = c(0, 0, -Inf, -Inf), est = c(1, 3, 4, 4), upper = c(2, Inf, Inf, Inf), fix = c(TRUE, TRUE, FALSE, TRUE), err = c(NA, NA, NA, NA), label = c("A", NA, NA, NA), condition = c(NA, NA, NA, NA)), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -4L), class = c("nlmixrBounds", "data.frame"));

    bnd1 <- function(){
        a = fix(0, 1, 2) # A
        b = fix(0, 3);
        c <- 4
        d <- fix(4);
    }

    bnd2 <- function(){
        a = FIX(0, 1, 2) # A
        b = FIX(0, 3);
        c <- 4
        d <- FIX(4);
    }

    bnd3 <- function(){
        a = fixed(0, 1, 2) # A
        b = fixed(0, 3);
        c <- 4
        d <- fixed(4);
    }

    bnd4 <- function(){
        a = FIXED(0, 1, 2) # A
        b = FIXED(0, 3);
        c <- 4
        d <- FIXED(4);
    }

    bnd5 <- function(){
        a = c(0, fix(1), 2) # A
        b = c(0, fix(3));
        c <- 4
        d <- fix(4);
    }

    bnd6 <- function(){
        a = c(0, FIX(1), 2) # A
        b = c(0, FIX(3));
        c <- 4
        d <- FIX(4);
    }

    bnd7 <- function(){
        a = c(0, fixed(1), 2) # A
        b = c(0, fixed(3));
        c <- 4
        d <- fixed(4);
    }

    bnd8 <- function(){
        a = c(0, FIXED(1), 2) # A
        b = c(0, FIXED(3));
        c <- 4
        d <- FIXED(4);
    }
    test_that("Theta fix fixed are reasonable", {
        expect_equal(nlmixrBounds(bnd1), ref1)
        expect_equal(nlmixrBounds(bnd2), ref1)
        expect_equal(nlmixrBounds(bnd3), ref1)
        expect_equal(nlmixrBounds(bnd4), ref1)
        expect_equal(nlmixrBounds(bnd5), ref1)
        expect_equal(nlmixrBounds(bnd6), ref1)
        expect_equal(nlmixrBounds(bnd7), ref1)
        expect_equal(nlmixrBounds(bnd8), ref1)
    })

    ref1 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd1 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, 30)
    }

    ref2 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),     err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA,     NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd2 <- function(){
        ~ fix(40,
              0.1, 20,
              0.1, 0.1, 30)
    }

    bnd3 <- function(){
        ~ fixed(40,
                0.1, 20,
                0.1, 0.1, 30)
    }

    bnd4 <- function(){
        ~ FIX(40,
              0.1, 20,
              0.1, 0.1, 30)
    }

    bnd5 <- function(){
        ~ FIXED(40,
                0.1, 20,
                0.1, 0.1, 30)
    }

    test_that("Total ETA fixed (unnamed)", {
        expect_equal(nlmixrBounds(bnd1), ref1)
        expect_equal(nlmixrBounds(bnd2), ref2)
        expect_equal(nlmixrBounds(bnd3), ref2)
        expect_equal(nlmixrBounds(bnd4), ref2)
        expect_equal(nlmixrBounds(bnd5), ref2)
    })

    ref6 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"));

    bnd6 <- function(){
        ~ c(fix(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    ref7 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd7 <- function(){
        ~ c(40,
            fix(0.1), 20,
            0.1, 0.1, 30)
    }

    ref8 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd8 <- function(){
        ~ c(40,
            0.1, fix(20),
            0.1, 0.1, 30)
    }

    ref9 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd9 <- function(){
        ~ c(40,
            0.1, 20,
            fix(0.1), 0.1, 30)
    }

    ref10 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))


    bnd10 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, fix(0.1), 30)
    }

    ref11 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"));

    bnd11 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, fix(30))
    }


    bnd6 <- function(){
        ~ c(fix(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    ref7 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd7 <- function(){
        ~ c(40,
            fix(0.1), 20,
            0.1, 0.1, 30)
    }

    ref8 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd8 <- function(){
        ~ c(40,
            0.1, fix(20),
            0.1, 0.1, 30)
    }

    ref9 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"));

    bnd9 <- function(){
        ~ c(40,
            0.1, 20,
            fix(0.1), 0.1, 30)
    }

    ref10 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"));

    bnd10 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, fix(0.1), 30)
    }

    ref11 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c(NA, NA, NA, NA, NA, NA), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd11 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, fix(30))
    }


    test_that("Total ETA fixed (unnamed)", {
        expect_equal(nlmixrBounds(bnd6), ref6)
        expect_equal(nlmixrBounds(bnd7), ref7)
        expect_equal(nlmixrBounds(bnd8), ref8)
        expect_equal(nlmixrBounds(bnd9), ref9)
        expect_equal(nlmixrBounds(bnd10), ref10)
        expect_equal(nlmixrBounds(bnd11), ref11)
    })

    bnd6a <- function(){
        ~ c(fixed(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    bnd7a <- function(){
        ~ c(40,
            fixed(0.1), 20,
            0.1, 0.1, 30)
    }

    bnd8a <- function(){
        ~ c(40,
            0.1, fixed(20),
            0.1, 0.1, 30)
    }

    bnd9a <- function(){
        ~ c(40,
            0.1, 20,
            fixed(0.1), 0.1, 30)
    }

    bnd10a <- function(){
        ~ c(40,
            0.1, 20,
            0.1, fixed(0.1), 30)
    }

    bnd11a <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, fixed(30))
    }

    bnd6a <- function(){
        ~ c(fixed(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    bnd7a <- function(){
        ~ c(40,
            fixed(0.1), 20,
            0.1, 0.1, 30)
    }

    bnd8a <- function(){
        ~ c(40,
            0.1, fixed(20),
            0.1, 0.1, 30)
    }

    bnd9a <- function(){
        ~ c(40,
            0.1, 20,
            fixed(0.1), 0.1, 30)
    }

    bnd10a <- function(){
        ~ c(40,
            0.1, 20,
            0.1, fixed(0.1), 30)
    }

    bnd11a <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, fixed(30))
    }

    test_that("Total ETA fixed (unnamed) #a", {
        expect_equal(nlmixrBounds(bnd6a), ref6)
        expect_equal(nlmixrBounds(bnd7a), ref7)
        expect_equal(nlmixrBounds(bnd8a), ref8)
        expect_equal(nlmixrBounds(bnd9a), ref9)
        expect_equal(nlmixrBounds(bnd10a), ref10)
        expect_equal(nlmixrBounds(bnd11a), ref11)
    })

    bnd6b <- function(){
        ~ c(FIX(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    bnd7b <- function(){
        ~ c(40,
            FIX(0.1), 20,
            0.1, 0.1, 30)
    }

    bnd8b <- function(){
        ~ c(40,
            0.1, FIX(20),
            0.1, 0.1, 30)
    }

    bnd9b <- function(){
        ~ c(40,
            0.1, 20,
            FIX(0.1), 0.1, 30)
    }

    bnd10b <- function(){
        ~ c(40,
            0.1, 20,
            0.1, FIX(0.1), 30)
    }

    bnd11b <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, FIX(30))
    }

    bnd6b <- function(){
        ~ c(FIX(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    bnd7b <- function(){
        ~ c(40,
            FIX(0.1), 20,
            0.1, 0.1, 30)
    }

    bnd8b <- function(){
        ~ c(40,
            0.1, FIX(20),
            0.1, 0.1, 30)
    }

    bnd9b <- function(){
        ~ c(40,
            0.1, 20,
            FIX(0.1), 0.1, 30)
    }

    bnd10b <- function(){
        ~ c(40,
            0.1, 20,
            0.1, FIX(0.1), 30)
    }

    bnd11b <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, FIX(30))
    }

    ##

    test_that("Total ETA fixed (unnamed) #b", {
        expect_equal(nlmixrBounds(bnd6b), ref6)
        expect_equal(nlmixrBounds(bnd7b), ref7)
        expect_equal(nlmixrBounds(bnd8b), ref8)
        expect_equal(nlmixrBounds(bnd9b), ref9)
        expect_equal(nlmixrBounds(bnd10b), ref10)
        expect_equal(nlmixrBounds(bnd11b), ref11)
    })

    bnd6c <- function(){
        ~ c(FIXED(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    bnd7c <- function(){
        ~ c(40,
            FIXED(0.1), 20,
            0.1, 0.1, 30)
    }

    bnd8c <- function(){
        ~ c(40,
            0.1, FIXED(20),
            0.1, 0.1, 30)
    }

    bnd9c <- function(){
        ~ c(40,
            0.1, 20,
            FIXED(0.1), 0.1, 30)
    }

    bnd10c <- function(){
        ~ c(40,
            0.1, 20,
            0.1, FIXED(0.1), 30)
    }

    bnd11c <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, FIXED(30))
    }

    bnd6c <- function(){
        ~ c(FIXED(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    bnd7c <- function(){
        ~ c(40,
            FIXED(0.1), 20,
            0.1, 0.1, 30)
    }

    bnd8c <- function(){
        ~ c(40,
            0.1, FIXED(20),
            0.1, 0.1, 30)
    }

    bnd9c <- function(){
        ~ c(40,
            0.1, 20,
            FIXED(0.1), 0.1, 30)
    }

    bnd10c <- function(){
        ~ c(40,
            0.1, 20,
            0.1, FIXED(0.1), 30)
    }

    bnd11c <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, FIXED(30))
    }

    test_that("Total ETA fixed (unnamed) #c", {
        expect_equal(nlmixrBounds(bnd6c), ref6)
        expect_equal(nlmixrBounds(bnd7c), ref7)
        expect_equal(nlmixrBounds(bnd8c), ref8)
        expect_equal(nlmixrBounds(bnd9c), ref9)
        expect_equal(nlmixrBounds(bnd10c), ref10)
        expect_equal(nlmixrBounds(bnd11c), ref11)
    })

    ref12 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),     err = c(NA, NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA,     NA, NA), condition = structure(c(1L, 1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd12 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, 0.1, 30)
    }

    ref13 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))


    bnd13 <- function(){
        eta1 + eta2 + eta3 ~ fix(40,
                                 0.1, 20,
                                 0.1, 0.1, 30)
    }

    bnd14 <- function(){
        eta1 + eta2 + eta3 ~ fixed(40,
                                   0.1, 20,
                                   0.1, 0.1, 30)
    }

    bnd15 <- function(){
        eta1 + eta2 + eta3 ~ FIX(40,
                                 0.1, 20,
                                 0.1, 0.1, 30)
    }

    bnd16 <- function(){
        eta1 + eta2 + eta3 ~ FIXED(40,
                                   0.1, 20,
                                   0.1, 0.1, 30)
    }

    test_that("Total ETA fixed (named)", {
        expect_equal(nlmixrBounds(bnd12), ref12)
        expect_equal(nlmixrBounds(bnd13), ref13)
        expect_equal(nlmixrBounds(bnd14), ref13)
        expect_equal(nlmixrBounds(bnd15), ref13)
        expect_equal(nlmixrBounds(bnd16), ref13)
    })

    ref17 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd17 <- function(){
        eta1 + eta2 + eta3 ~ c(fix(40),
                               0.1, 20,
                               0.1, 0.1, 30)
    }


    ref18 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"));

    bnd18 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               fix(0.1), 20,
                               0.1, 0.1, 30)
    }

    ref19 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd19 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, fix(20),
                               0.1, 0.1, 30)
    }

    ref20 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd20 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               fix(0.1), 0.1, 30)
    }

    ref21 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd21 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, fix(0.1), 30)
    }

    ref22 <- structure(list(ntheta = c(NA, NA, NA, NA, NA, NA), neta1 = c(1, 2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = structure(1:6, .Label = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"), class = "factor"),     lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,     0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,     Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), err = c(NA,     NA, NA, NA, NA, NA), label = c(NA, NA, NA, NA, NA, NA), condition = structure(c(1L,     1L, 1L, 1L, 1L, 1L), .Label = "ID", class = "factor")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"))

    bnd22 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, 0.1, fix(30))
    }
    ##

    bnd17a <- function(){
        eta1 + eta2 + eta3 ~ c(FIX(40),
                               0.1, 20,
                               0.1, 0.1, 30)
    }

    bnd18a <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               FIX(0.1), 20,
                               0.1, 0.1, 30)
    }

    bnd19a <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, FIX(20),
                               0.1, 0.1, 30)
    }

    bnd20a <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               FIX(0.1), 0.1, 30)
    }

    bnd21a <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, FIX(0.1), 30)
    }

    bnd22a <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, 0.1, FIX(30))
    }
    ##

    bnd17b <- function(){
        eta1 + eta2 + eta3 ~ c(fixed(40),
                               0.1, 20,
                               0.1, 0.1, 30)
    }

    bnd18b <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               fixed(0.1), 20,
                               0.1, 0.1, 30)
    }

    bnd19b <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, fixed(20),
                               0.1, 0.1, 30)
    }

    bnd20b <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               fixed(0.1), 0.1, 30)
    }

    bnd21b <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, fixed(0.1), 30)
    }

    bnd22b <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, 0.1, fixed(30))
    }

    ##

    bnd17c <- function(){
        eta1 + eta2 + eta3 ~ c(FIX(40),
                               0.1, 20,
                               0.1, 0.1, 30)
    }

    bnd18c <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               FIX(0.1), 20,
                               0.1, 0.1, 30)
    }

    bnd19c <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, FIX(20),
                               0.1, 0.1, 30)
    }

    bnd20c <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               FIX(0.1), 0.1, 30)
    }

    bnd21c <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, FIX(0.1), 30)
    }

    bnd22c <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, 0.1, FIX(30))
    }

    ##
    bnd17d <- function(){
        eta1 + eta2 + eta3 ~ c(FIXED(40),
                               0.1, 20,
                               0.1, 0.1, 30)
    }

    bnd18d <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               FIXED(0.1), 20,
                               0.1, 0.1, 30)
    }

    bnd19d <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, FIXED(20),
                               0.1, 0.1, 30)
    }

    bnd20d <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               FIXED(0.1), 0.1, 30)
    }

    bnd21d <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, FIXED(0.1), 30)
    }

    bnd22d <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, 0.1, FIXED(30))
    }


    test_that("Total ETA FIXED (named)", {
        expect_equal(nlmixrBounds(bnd17), ref17)
        expect_equal(nlmixrBounds(bnd18), ref18)
        expect_equal(nlmixrBounds(bnd19), ref19)
        expect_equal(nlmixrBounds(bnd20), ref20)
        expect_equal(nlmixrBounds(bnd21), ref21)
        expect_equal(nlmixrBounds(bnd22), ref22)
        expect_equal(nlmixrBounds(bnd17a), ref17)
        expect_equal(nlmixrBounds(bnd18a), ref18)
        expect_equal(nlmixrBounds(bnd19a), ref19)
        expect_equal(nlmixrBounds(bnd20a), ref20)
        expect_equal(nlmixrBounds(bnd21a), ref21)
        expect_equal(nlmixrBounds(bnd22a), ref22)
        expect_equal(nlmixrBounds(bnd17b), ref17)
        expect_equal(nlmixrBounds(bnd18b), ref18)
        expect_equal(nlmixrBounds(bnd19b), ref19)
        expect_equal(nlmixrBounds(bnd20b), ref20)
        expect_equal(nlmixrBounds(bnd21b), ref21)
        expect_equal(nlmixrBounds(bnd22b), ref22)
        expect_equal(nlmixrBounds(bnd17c), ref17)
        expect_equal(nlmixrBounds(bnd18c), ref18)
        expect_equal(nlmixrBounds(bnd19c), ref19)
        expect_equal(nlmixrBounds(bnd20c), ref20)
        expect_equal(nlmixrBounds(bnd21c), ref21)
        expect_equal(nlmixrBounds(bnd22c), ref22)
        expect_equal(nlmixrBounds(bnd17d), ref17)
        expect_equal(nlmixrBounds(bnd18d), ref18)
        expect_equal(nlmixrBounds(bnd19d), ref19)
        expect_equal(nlmixrBounds(bnd20d), ref20)
        expect_equal(nlmixrBounds(bnd21d), ref21)
        expect_equal(nlmixrBounds(bnd22d), ref22)

    })

    f1 <- function(){
        lCl = c(5, 5, 5) # A
    }

    f2 <- function(){
        lCl = c(0, -1.3) # lCl
    }

    f3 <- function(){
        lCl = c(0, -1.3, -10) # lCl
    }

    f4 <- function(){
        lCl = c(0, 5, 5) # A
    }

    f5 <- function(){
        lCl = c(0, 0, 5) # A
    }

    f6 <- function(){
        lCl = c(5, 5) # A
    }


    test_that("Invalid bounds raise errors",{
        expect_error(nlmixrBounds(f1), rex::rex("The estimate, and upper and lower bounds are the same for the following parameters: lCl\nTo fix parameters use lCl=fix(5) instead."))
        expect_error(nlmixrBounds(f2), rex::rex("The lower bound is higher than the estimate for these parameters: lCl.\nYou can adjust by lCl=c(-1.3, 0) # c(lower, est)"))
        expect_error(nlmixrBounds(f3), rex::rex("The bounds make no sense for these parameters: lCl.\nThey should be ordered as follows: lCl=c(-10, -1.3, 0) # c(lower, est, upper)"))
        expect_error(nlmixrBounds(f4), rex::rex("The estimate is the same as a boundary for the following parameters: lCl\nInstead use lCl=c(0, 5) # c(lower, est)"))
        expect_error(nlmixrBounds(f5), rex::rex("The estimate is the same as a boundary for the following parameters: lCl\nInstead use lCl=c(0, 5) # c(lower, est)"))
        expect_error(nlmixrBounds(f6), rex::rex("The estimate is the same as a boundary for the following parameter: lCl\nInstead use lCl=5 # est"))
    })
}, cran=TRUE)
