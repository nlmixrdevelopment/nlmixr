context("Test bounds extraction")

test_that("as.nlmixrBounds, data.frame to bounds creation works", {
  expect_error(
    as.nlmixrBounds(data.frame()),
    regexp="no parameter information"
  )
  expect_error(
    as.nlmixrBounds(data.frame(ntheta=1)),
    regexp=
      paste(
        "columns missing:",
        paste0("'", setdiff(names(nlmixr:::nlmixrBoundsTemplate), "ntheta"), "'", collapse=", ")
      )
  )
  {
    zero_bound <- nlmixrBoundsTemplate[1:2,]
    zero_bound$ntheta <- 1:2
    zero_bound$lower <- c(-Inf, 0)
    zero_bound$est <- c(-5, 5)
    zero_bound$upper <- c(0, Inf)
    expect_equal(
      as.data.frame(as.nlmixrBounds(zero_bound)[, c("lower", "upper")]),
      data.frame(
        lower=c(-Inf, sqrt(.Machine$double.eps)),
        upper=c(-sqrt(.Machine$double.eps), Inf)
      ),
      # row.names will not be equal
      check.attributes=FALSE
    )
  }
})

rxPermissive({
    ref <- structure(list(ntheta = c(1, 2, 3, 4, 5, 6, 7, 8, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
9, 10, 11, 12, 13, 14), neta1 = c(NA, NA, NA, NA, NA, NA, NA,
NA, 1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12,
12, NA, NA, NA, NA, NA, NA), neta2 = c(NA, NA, NA, NA, NA, NA,
NA, NA, 1, 2, 3, 4, 5, 5, 6, 7, 7, 8, 7, 8, 9, 10, 10, 11, 10,
11, 12, NA, NA, NA, NA, NA, NA), name = c("a", "b", "c", "d",
NA, NA, NA, NA, "et1", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, "et2", "(et3,et2)", "et3", "(et4,et2)", "(et4,et3)",
"et4", NA, NA, NA, "a5", "a6", "a7"), lower = c(1.49011611938477e-08,
1.49011611938477e-08, -Inf, -Inf, 1.49011611938477e-08, 1.49011611938477e-08,
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,
9, 11, -Inf, 9, 11), est = c(1, 3, 4, 4, 1, 1, 1, 1, 10, 20,
30, 40, 40, 0.1, 20, 40, 0.1, 20, 0.1, 0.1, 30, 40, 0.1, 20,
0.1, 0.1, 30, 8, 10, 12, 8, 10, 12), upper = c(2, Inf, Inf, Inf,
2, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf,
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 13,
Inf, Inf, 13), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), err = as.character(c(NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
)), label = c("A", NA, NA, NA, "e", NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, "labels", NA, NA, NA), condition = c(NA, NA, NA, NA, NA,
NA, NA, NA, "ID", "ID", "ID", "ID", "STUD", "STUD", "STUD", "ID",
"ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID",
NA, NA, NA, NA, NA, NA)), row.names = c(NA, -33L), class = c("nlmixrBounds",
"data.frame"))

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
        expect_error(nlmixrBounds(bnd), regexp="Syntax is not supported for thetas: c(1, 2, 3, 4, 5)", fixed=TRUE)
        expect_error(nlmixrBounds(bnda), regexp="Syntax is not supported for thetas: c(1, 2, 3, 4, 5)", fixed=TRUE)
        expect_error(nlmixrBounds(bndb), regexp="Syntax is not supported for thetas: c(1, 2, 3, 4, 5)", fixed=TRUE)
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
        expect_error(nlmixrBounds(bnd), regexp="Syntax is not supported for thetas: c(1, 2, 3, 4)", fixed=TRUE)
        expect_error(nlmixrBounds(bnda), regexp="Syntax is not supported for thetas: c(1, 2, 3, 4)", fixed=TRUE)
        expect_error(nlmixrBounds(bndb), regexp="Syntax is not supported for thetas: c(1, 2, 3, 4)", fixed=TRUE)
    })

    bnd1 <- function(){
        ~ c(1)
    }
    ref1 <- structure(list(ntheta = NA_real_, neta1 = 1, neta2 = 1, name = NA_character_,
    lower = -Inf, est = 1, upper = Inf, fix = FALSE, err = NA_character_,
    label = NA_character_, condition = "ID"), row.names = c(NA, -1L), class = c("nlmixrBounds",
"data.frame"))

    bnd2 <- function(){
        ~ c(1, 2)
    }

    bnd3 <- function(){
        ~ c(1,
            2, 3)
    }

    ref3 <- structure(list(ntheta = as.numeric(c(NA, NA, NA)), neta1 = c(1, 2, 2), neta2 = c(1,
1, 2), name = as.character(c(NA, NA, NA)), lower = c(-Inf, -Inf, -Inf), est = c(1,
2, 3), upper = c(Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE),
    err = as.character(c(NA, NA, NA)), label = as.character(c(NA, NA, NA)), condition = c("ID",
    "ID", "ID")), row.names = c(NA, -3L), class = c("nlmixrBounds",
"data.frame"))

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

    ref6 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(1, 2, 3, 4, 5, 6), upper = c(Inf, Inf, Inf, Inf, Inf,
Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = as.character(c(NA,
NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
"ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))


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

    ref10 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA)), neta1 = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4), neta2 = c(1, 1,
2, 1, 2, 3, 1, 2, 3, 4), name = as.character(c(NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,
-Inf, -Inf, -Inf), est = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), upper = c(Inf,
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE,
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)), label = as.character(c(NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA)), condition = c("ID", "ID",
"ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-10L), class = c("nlmixrBounds", "data.frame"))


    test_that("Bad Lower trianglar matrices throw errors.", {
        expect_equal(nlmixrBounds(bnd1), ref1);
        expect_error(nlmixrBounds(bnd2), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2)", fixed=TRUE)
        expect_equal(nlmixrBounds(bnd3), ref3);
        expect_error(nlmixrBounds(bnd4), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4)", fixed=TRUE)
        expect_error(nlmixrBounds(bnd5), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5)", fixed=TRUE)
        expect_equal(nlmixrBounds(bnd6), ref6);
        expect_error(nlmixrBounds(bnd7), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7)", fixed=TRUE)
        expect_error(nlmixrBounds(bnd8), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8)", fixed=TRUE)
        expect_error(nlmixrBounds(bnd9), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8, 9)", fixed=TRUE)
        expect_equal(nlmixrBounds(bnd10), ref10)
    })


    bnd1 <- function(){
        eta1 ~ c(1)
    }

    ref1 <- structure(list(ntheta = NA_real_, neta1 = 1, neta2 = 1, name = "eta1",
    lower = -Inf, est = 1, upper = Inf, fix = FALSE, err = NA_character_,
    label = NA_character_, condition = "ID"), row.names = c(NA, -1L), class = c("nlmixrBounds",
"data.frame"))

    bnd2 <- function(){
        eta1 ~ c(1, 2)
    }

    bnd3 <- function(){
        eta1 + eta2~ c(1,
                       2, 3)
    }

    ref3 <- structure(list(ntheta = as.numeric(c(NA, NA, NA)), neta1 = c(1, 2, 2), neta2 = c(1,
1, 2), name = c("eta1", "(eta2,eta1)", "eta2"), lower = c(-Inf,
-Inf, -Inf), est = c(1, 2, 3), upper = c(Inf, Inf, Inf), fix = c(FALSE,
FALSE, FALSE), err = as.character(c(NA, NA, NA)), label = as.character(c(NA, NA, NA)), condition = c("ID",
"ID", "ID")), row.names = c(NA, -3L), class = c("nlmixrBounds",
"data.frame"))

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

    ref6 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(1,
    2, 3, 4, 5, 6), upper = c(Inf, Inf, Inf, Inf, Inf, Inf),
    fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))

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

    ref10 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA)), neta1 = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4), neta2 = c(1, 1,
2, 1, 2, 3, 1, 2, 3, 4), name = c("eta1", "(eta2,eta1)", "eta2",
"(eta3,eta1)", "(eta3,eta2)", "eta3", "(eta4,eta1)", "(eta4,eta2)",
"(eta4,eta3)", "eta4"), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf,
-Inf, -Inf, -Inf, -Inf, -Inf), est = c(1, 2, 3, 4, 5, 6, 7, 8,
9, 10), upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf,
Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
FALSE, FALSE, FALSE), err = as.character(c(NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)), condition = c("ID",
"ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-10L), class = c("nlmixrBounds", "data.frame"))


    test_that("Bad Lower trianglar matrices (with labels) throw errors.", {
        expect_equal(nlmixrBounds(bnd1), ref1);
        expect_error(nlmixrBounds(bnd2), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2)", fixed=TRUE)
        expect_equal(nlmixrBounds(bnd3), ref3);
        expect_error(nlmixrBounds(bnd4), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4)", fixed=TRUE)
        expect_error(nlmixrBounds(bnd5), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5)", fixed=TRUE)
        expect_equal(nlmixrBounds(bnd6), ref6);
        expect_error(nlmixrBounds(bnd7), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7)", fixed=TRUE)
        expect_error(nlmixrBounds(bnd8), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8)", fixed=TRUE)
        expect_error(nlmixrBounds(bnd9), regexp="incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8, 9)", fixed=TRUE)
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
        expect_error(nlmixrBounds(bnd10.a), regexp="omega assignment left handed side must match lower triangular matrix size", fixed=TRUE)
        expect_error(nlmixrBounds(bnd10.b), regexp="omega assignment left handed side must match lower triangular matrix size", fixed=TRUE)
    })

    bnd3 <- function(){
        eta1 + eta2~ c(1, # Comment here.
                       2, 3)
    }

    test_that("Comments inside bounds are not supported!", {
        expect_error(nlmixrBounds(bnd3), regexp="error parsing bounds: possible (unsupported) comment/condition inside bounds", fixed=TRUE)
    })

    bnd1 <- function(){
        eta0 ~ 0.3
        eta1 + eta2~ c(1,
                       2, 3) | STUD
        ~ c(1,
            2, 3)
    }

    ref <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 3, 3, 4, 5, 5), neta2 = c(1, 2, 2, 3, 4, 4, 5), name = c("eta0",
"eta1", "(eta2,eta1)", "eta2", NA, NA, NA), lower = c(-Inf, -Inf,
-Inf, -Inf, -Inf, -Inf, -Inf), est = c(0.3, 1, 2, 3, 1, 2, 3),
    upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf), fix = c(FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), err = as.character(c(NA, NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA, NA
    )), condition = c("ID", "STUD", "STUD", "STUD", "ID", "ID",
    "ID")), row.names = c(NA, -7L), class = c("nlmixrBounds",
"data.frame"))

    test_that("Conditional statments are supported correctly.", {
        expect_equal(nlmixrBounds(bnd1), ref)
    })

    ## Now test err/pred type parsing.

    bnd1 <- function(){
        err ~ add(0.1)
    }

    ref1 <- structure(list(ntheta = 1, neta1 = NA, neta2 = NA, name = "err",
    lower = 1.49011611938477e-08, est = 0.1, upper = Inf, fix = TRUE,
    err = "add", label = NA, condition = NA), row.names = c(NA,
-1L), class = c("nlmixrBounds", "data.frame"))


    bnd4 <- function(){
        err ~ prop(0.1)
    }

    ref4 <- structure(list(ntheta = 1, neta1 = NA, neta2 = NA, name = "err",
    lower = 1.49011611938477e-08, est = 0.1, upper = Inf, fix = TRUE,
    err = "prop", label = NA, condition = NA), row.names = c(NA,
-1L), class = c("nlmixrBounds", "data.frame"))

    test_that("Error parsing is reasonable", {
        expect_equal(nlmixrBounds(bnd1), ref1)
        expect_equal(nlmixrBounds(bnd4), ref4)
    })


    ## Now try a fixed parameter block

    ref1 <- structure(list(ntheta = c(1, 2, 3, 4), neta1 = as.numeric(c(NA, NA, NA,
NA)), neta2 = as.numeric(c(NA, NA, NA, NA)), name = c("a", "b", "c", "d"),
    lower = c(1.49011611938477e-08, 1.49011611938477e-08, -Inf,
    -Inf), est = c(1, 3, 4, 4), upper = c(2, Inf, Inf, Inf),
    fix = c(TRUE, TRUE, FALSE, TRUE), err = as.character(c(NA, NA, NA, NA)),
    label = c("A", NA, NA, NA), condition = as.character(c(NA, NA, NA, NA))), row.names = c(NA,
-4L), class = c("nlmixrBounds", "data.frame"))

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

    ref1 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd1 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, 0.1, 30)
    }

    ref2 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
    NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID"
    )), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"
))

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

    ref6 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd6 <- function(){
        ~ c(fix(40),
            0.1, 20,
            0.1, 0.1, 30)
    }

    ref7 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd7 <- function(){
        ~ c(40,
            fix(0.1), 20,
            0.1, 0.1, 30)
    }

    ref8 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd8 <- function(){
        ~ c(40,
            0.1, fix(20),
            0.1, 0.1, 30)
    }

    ref9 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd9 <- function(){
        ~ c(40,
            0.1, 20,
            fix(0.1), 0.1, 30)
    }

    ref10 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd10 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, fix(0.1), 30)
    }

    ref11 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

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

    ref7 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd7 <- function(){
        ~ c(40,
            fix(0.1), 20,
            0.1, 0.1, 30)
    }

    ref8 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd8 <- function(){
        ~ c(40,
            0.1, fix(20),
            0.1, 0.1, 30)
    }

    ref9 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd9 <- function(){
        ~ c(40,
            0.1, 20,
            fix(0.1), 0.1, 30)
    }

    ref10 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

    bnd10 <- function(){
        ~ c(40,
            0.1, 20,
            0.1, fix(0.1), 30)
    }

    ref11 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = as.character(c(NA, NA,
NA, NA, NA, NA)), lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf
), est = c(40, 0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf,
Inf, Inf, Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
), err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID")), row.names = c(NA,
-6L), class = c("nlmixrBounds", "data.frame"))

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

    ref12 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    err = as.character(c(NA, NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA,
    NA, NA)), condition = c("ID", "ID", "ID", "ID", "ID", "ID"
    )), row.names = c(NA, -6L), class = c("nlmixrBounds", "data.frame"
))

    bnd12 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, 0.1, 30)
    }

    ref13 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))


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

    ref17 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))

    bnd17 <- function(){
        eta1 + eta2 + eta3 ~ c(fix(40),
                               0.1, 20,
                               0.1, 0.1, 30)
    }


    ref18 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))

    bnd18 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               fix(0.1), 20,
                               0.1, 0.1, 30)
    }

    ref19 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))

    bnd19 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, fix(20),
                               0.1, 0.1, 30)
    }

    ref20 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))

    bnd20 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               fix(0.1), 0.1, 30)
    }

    ref21 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))

    bnd21 <- function(){
        eta1 + eta2 + eta3 ~ c(40,
                               0.1, 20,
                               0.1, fix(0.1), 30)
    }

    ref22 <- structure(list(ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA)), neta1 = c(1,
2, 2, 3, 3, 3), neta2 = c(1, 1, 2, 1, 2, 3), name = c("eta1",
"(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
    lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), est = c(40,
    0.1, 20, 0.1, 0.1, 30), upper = c(Inf, Inf, Inf, Inf, Inf,
    Inf), fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), err = as.character(c(NA,
    NA, NA, NA, NA, NA)), label = as.character(c(NA, NA, NA, NA, NA, NA)), condition = c("ID",
    "ID", "ID", "ID", "ID", "ID")), row.names = c(NA, -6L), class = c("nlmixrBounds",
"data.frame"))

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

    f7 <- function(){
        lCl < 3
    }

    test_that("Invalid bounds raise errors",{
        expect_error(nlmixrBounds(f1), regexp="consider fixing these:\n     lCl = fixed(5)", fixed=TRUE)
        expect_error(nlmixrBounds(f2), regexp="reorder bounds:\n     lCl = c(-1.3, 0)", fixed=TRUE)
        expect_error(nlmixrBounds(f3), regexp="reorder bounds:\n     lCl = c(-10, -1.3, 0)", fixed=TRUE)
        expect_error(nlmixrBounds(f4), regexp="consider fixing these:\n     lCl = fixed(5)", fixed=TRUE)
        expect_error(nlmixrBounds(f5), regexp="consider fixing these:\n     lCl = fixed(0)", fixed=TRUE)
        expect_error(nlmixrBounds(f6), regexp="consider fixing these:\n     lCl = fixed(5)", fixed=TRUE)
        expect_error(nlmixrBounds(f7), regexp="invalid call in initial conditions: lCl < 3", fixed=TRUE)
    })
}, test="cran")


# nlmixrBounds ####

test_that("nlmixrBounds", {
  expect_equal(
    nlmixrBoundsParser(
      function() {
        ({({a = 1;b=2})})
        {c = 3}
        d <- 4
      }
    ),
    list(
      list(
        operation=c("assign", "theta"),
        varname="a",
        value=1
      ),
      list(
        operation=c("assign", "theta"),
        varname="b",
        value=2
      ),
      list(
        operation=c("assign", "theta"),
        varname="c",
        value=3
      ),
      list(
        operation=c("assign", "theta"),
        varname="d",
        value=4
      )
    ),
    info="Nested assignments are unnested"
  )
})

# nlmixrBoundsValueFixed ####

test_that("nlmixrBoundsValueFixed", {
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~1)[[2]]),
    list(value=1, fixed=FALSE)
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1))[[2]]),
    list(value=1, fixed=FALSE)
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1, 2))[[2]]),
    list(value=c(1, 2), fixed=rep(FALSE, 2))
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1, 2, 3))[[2]]),
    list(value=c(1, 2, 3), fixed=rep(FALSE, 3))
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1, fixed))[[2]]),
    list(value=1, fixed=TRUE)
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1, 2, fixed))[[2]]),
    list(value=c(1, 2), fixed=rep(TRUE, 2))
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1, 2, 3, fixed))[[2]]),
    list(value=c(1, 2, 3), fixed=rep(TRUE, 3))
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1, fixed(2), 3, fixed))[[2]]),
    list(value=c(1, 2, 3), fixed=rep(TRUE, 3)),
    info="Fixed is specified two ways, but they are not in conflict"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(1, fixed(2), 3))[[2]]),
    list(value=c(1, 2, 3), fixed=c(FALSE, TRUE, FALSE))
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(fixed(1), 2, 3))[[2]]),
    list(value=c(1, 2, 3), fixed=c(TRUE, FALSE, FALSE))
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(fixed(1), 2, 3))[[2]]),
    list(value=c(1, 2, 3), fixed=c(TRUE, FALSE, FALSE))
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(fixed(1), log(2), 3))[[2]]),
    list(value=c(1, log(2), 3), fixed=c(TRUE, FALSE, FALSE)),
    info="Function evaluation works (though it may have issues related to environment precedence). (Fix #253)"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(log(0), log(1.5), log(20)))[[2]]),
    list(value=log(c(0, 1.5, 20)), fixed=rep(FALSE, 3)),
    info="Function evaluation works (though it may have issues related to environment precedence). (Fix #253)"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~FIX(log(0), log(1.5), log(20)))[[2]]),
    list(value=log(c(0, 1.5, 20)), fixed=rep(TRUE, 3)),
    info="Function evaluation works (though it may have issues related to environment precedence). (Fix #253)"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsValueFixed((~c(FIX(log(0), 1/log(1.5)), 1/log(20)))[[2]]),
    list(value=c(log(0), 1/log(1.5), 1/log(20)), fixed=c(TRUE, TRUE, FALSE)),
    info="Function evaluation works and arbitrary complexity may be within the fixed() call (or outside of it). (Fix #253)"
  )
  expect_error(
    expect_warning(
      nlmixr:::nlmixrBoundsValueFixed((~FIX(sqrt(-1)))[[2]]),
      regexp="NaNs produced"
    ),
    regexp="NaN values in initial condition: FIX(sqrt(-1))",
    fixed=TRUE,
    info="Invalid math stops execution"
  )
  expect_error(
    nlmixr:::nlmixrBoundsValueFixed((~FIX("A"))[[2]]),
    regexp='non-numeric values in initial condition: FIX("A")',
    fixed=TRUE,
    info="Values must be numbers"
  )
  expect_error(
    nlmixr:::nlmixrBoundsValueFixed((~a)[[2]]),
    regexp="error parsing initial condition 'a': object 'a' not found",
    fixed=TRUE,
    info="No variable substitutions are performed for parsing."
  )
})

# nlmixrBoundsReplaceFixed ####

test_that("nlmixrBoundsReplaceFixed, testing replacement of fixed names within calls", {
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~a)[[2]]),
    list(
      call=(~a)[[2]],
      fixed=FALSE
    )
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~c(1, fixed))[[2]]),
    list(
      call=(~c(1))[[2]],
      fixed=TRUE
    )
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~c(1, c(1, fixed)))[[2]]),
    list(
      call=(~c(1, c(1)))[[2]],
      fixed=TRUE
    )
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~1)[[2]]),
    list(
      call=1,
      fixed=FALSE
    )
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~fixed(1))[[2]], replacementFun="c"),
    list(
      call=(~c(1))[[2]],
      fixed=FALSE
    )
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~fixed(1))[[2]], replacementFun="c"),
    list(
      call=(~c(1))[[2]],
      fixed=FALSE
    )
  )
  # This is weird syntax to use, but it is logically okay
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~c(1, 1, c(fixed)))[[2]]),
    list(
      call=(~c(1, 1, c()))[[2]],
      fixed=TRUE
    ),
    info="Fixed can only be at the end of a vector of values, and detection of that works even when it is in a sub-expression."
  )
  expect_error(
    nlmixr:::nlmixrBoundsReplaceFixed((~c(fixed, 1))[[2]]),
    regexp="'fixed' may only be the last item in a list: c(fixed, 1)",
    fixed=TRUE,
    info="Fixed can only be at the end of a vector of values"
  )
  expect_error(
    nlmixr:::nlmixrBoundsReplaceFixed((~c(1, c(fixed), 1))[[2]]),
    regexp="'fixed' may only be the last item in a list: c(1, c(fixed), 1)",
    fixed=TRUE,
    info="Fixed can only be at the end of a vector of values, and detection of that works even when it is at the end of its sub-expression, but it is not at the overall-end of the expression."
  )
})

test_that("nlmixrBoundsReplaceFixed, testing replacement of fixed function calls within calls", {
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~fixed(a))[[2]]),
    list(
      call=(~fixed(a))[[2]],
      fixed=FALSE
    ),
    info="`fixed()` is returned unchanged"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~fix(a))[[2]]),
    list(
      call=(~fixed(a))[[2]],
      fixed=FALSE
    ),
    info="`fix()` is changed to `fixed()`"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~FIX(a))[[2]]),
    list(
      call=(~fixed(a))[[2]],
      fixed=FALSE
    ),
    info="`FIX()` is changed to `fixed()`"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~FIXED(a))[[2]]),
    list(
      call=(~fixed(a))[[2]],
      fixed=FALSE
    ),
    info="`FIXED()` is changed to `fixed()`"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~1/FIX(a))[[2]]),
    list(
      call=(~1/fixed(a))[[2]],
      fixed=FALSE
    ),
    info="`FIX()` is changed to `fixed()` inside another expression"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsReplaceFixed((~1/FIX(a))[[2]], replacementFun="c"),
    list(
      call=(~1/c(a))[[2]],
      fixed=FALSE
    ),
    info="`FIX()` is changed to `c()` inside another expression to allow for arbitrary calculations as part of initial conditions."
  )
})

# nlmixrBoundsPrepareFun ####

test_that("preparation of the function for bound extraction", {
  expect_equal(
    nlmixr:::nlmixrBoundsPrepareFun(function() {1}),
    function() {1}
  )
  expect_equal(
    expect_message(
      nlmixr:::nlmixrBoundsPrepareFun(
        function() {
          1 # foo
        }
      ),
      regexp="parameter labels from comments will be replaced by 'label()'",
      fixed=TRUE
    ),
    function() {
      1
      label("foo")
    },
    # Env and srcref attributes will not be equal
    check.attributes=FALSE,
    info="comment lines are converted to labels"
  )
})

# nlmixrBoundsPrepareFunComments ####

test_that("Extraction of comments to labels with nlmixrBoundsPrepareFunComments", {
  nlmixrTestFunToChar <- function(x) {
    as.character(attr(x, "srcref"), useSource = TRUE)
  }
  expect_equal(
    nlmixr:::nlmixrBoundsPrepareFunComments(nlmixrTestFunToChar(
      function() {
        # hello
      }
    )),
    function() {},
    # Env and srcref attributes will not be equal
    check.attributes=FALSE,
    info="comment lines without other information are dropped"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsPrepareFunComments(nlmixrTestFunToChar(
      function() {
        1 # hello
      }
    )),
    function() {
      1
      label("hello")
    },
    # Env and srcref attributes will not be equal
    check.attributes=FALSE,
    info="comment lines with other information are converted to label()"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsPrepareFunComments(nlmixrTestFunToChar(
      function() {
        1|STUDY # hello
      }
    )),
    function() {
      1|STUDY
      label("hello")
    },
    # Env and srcref attributes will not be equal
    check.attributes=FALSE,
    info="comment lines with other information are converted to label() (even if they are on a line with a condition)"
  )
  expect_equal(
    nlmixr:::nlmixrBoundsPrepareFunComments(nlmixrTestFunToChar(
      function() {
        1 # label 1
        label("# hash in a quote may try to be detected as a label, but that is wrong")
      }
    )),
    function() {
      1
      label("label 1")
      label("# hash in a quote may try to be detected as a label, but that is wrong")
    },
    # Env and srcref attributes will not be equal
    check.attributes=FALSE,
    info="This is challenging to parse, and it was formerly a bug.  It is the reason that we are moving to parsing and not string extraction."
  )
})
