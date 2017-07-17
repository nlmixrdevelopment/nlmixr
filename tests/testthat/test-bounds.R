context("Test bounds extraction")

ref <-structure(list(ntheta = c(1, 2, 3, 4, 5, 6, 7, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19), neta1 = c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12, 12, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), neta2 = c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 2, 3, 4, 5, 5, 6, 7, 7, 8, 7, 8, 9, 10, 10, 11, 10, 11, 12, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), name = structure(c(1L, 2L, 3L, 4L, NA, NA, NA, NA, 5L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 6L, 7L, 8L, 9L, 10L, 11L, NA, NA, NA, 12L, 13L, 14L, 15L, 15L, 15L, 15L, 15L), .Label = c("a", "b", "c", "d", "et1", "et2", "(et3,et2)", "et3", "(et4,et2)", "(et4,et3)", "et4", "a5", "a6", "a7", "err"), class = "factor"), lower = c(0, 0, -Inf, -Inf, 0, 0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 9, 11, -Inf, 9, 11, 0, 0, 0, 0, 0), est = c(1, 3, 4, 4, 1, 1, 1, 1, 10, 20, 30, 40, 40, 0.1, 20, 40, 0.1, 20, 0.1, 0.1, 30, 40, 0.1, 20, 0.1, 0.1, 30, 8, 10, 12, 8, 10, 12, 0.1, 0.2, 0.3, 0.4, 0.5), upper = c(2, Inf, Inf, Inf, 2, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 13, Inf, Inf, 13, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), err = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "prop", "add", "add", "prop", "add"), label = c("A", NA, NA, NA, "e", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "labels", NA, NA, NA, NA, NA, NA, NA, NA), condition = c(NA, NA, NA, NA, NA, NA, NA, NA, "ID", "ID", "ID", "ID", "STUD", "STUD", "STUD", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", NA, NA, NA, NA, NA, NA, "NVS == 1", "NVS == 1", "NVS == 2", "NVS == 2", "NVS == 3")), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -38L), class = c("nlmixrBounds", "data.frame"))

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

    err ~ prop(0.1) + add(0.2) | NVS == 1
    err ~ add(0.3) + prop(0.4) | NVS == 2
    err ~ add(0.5)   | NVS == 3
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
    expect_equal(nlmixrBounds(bnd2), ref2)
    expect_equal(nlmixrBounds(bnd3), ref3)
    expect_equal(nlmixrBounds(bnd4), ref4)
    expect_equal(nlmixrBounds(bnd5), ref5)
})


## Now try a fixed parameter block

bnd1 <- function(){
    a = fix(0, 1, 2) # A
    b = fix(0, 3);
    c <- 4
    d <- fix(4);
}

ref1 <- structure(list(ntheta = c(1, 2, 3, 4), neta1 = c(NA, NA, NA, NA), neta2 = c(NA, NA, NA, NA), name = structure(1:4, .Label = c("a", "b", "c", "d"), class = "factor"), lower = c(0, 0, -Inf, -Inf), est = c(1, 3, 4, 4), upper = c(2, Inf, Inf, Inf), fix = c(TRUE, TRUE, FALSE, TRUE), err = c(NA, NA, NA, NA), label = c("A", NA, NA, NA), condition = c(NA, NA, NA, NA)), .Names = c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper", "fix", "err", "label", "condition"), row.names = c(NA, -4L), class = c("nlmixrBounds", "data.frame"))


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
test_that("Theta fix fixed are reasonable", {
    expect_equal(nlmixrBounds(bnd1), ref1)
    expect_equal(nlmixrBounds(bnd2), ref1)
    expect_equal(nlmixrBounds(bnd3), ref1)
    expect_equal(nlmixrBounds(bnd4), ref1)
})

nlmixrBounds(bnd1)
