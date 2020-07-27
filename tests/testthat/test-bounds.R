nlmixrTest({

  context("Test bounds extraction")

  test_that("as.nlmixrBounds, data.frame to bounds creation works", {
    expect_error(
      nlmixr:::as.nlmixrBounds(data.frame()),
      regexp = "no parameter information"
    )
    expect_error(
      nlmixr:::as.nlmixrBounds(data.frame(ntheta = 1)),
      regexp =
        paste(
          "columns missing:",
          paste0("'", setdiff(names(nlmixr:::nlmixrBoundsTemplate), "ntheta"), "'", collapse = ", ")
        )
    )
    ref <- nlmixr:::nlmixrBoundsTemplate
    ref$ntheta <- 1
    ref$lower <- 0
    ref$est <- 1
    ref$upper <- 2
    expect_equal(
      as.data.frame(nlmixr:::as.nlmixrBounds(
        data.frame(ntheta = 1, est = 1, lower = 0, upper = 2),
        addMissingCols = TRUE
      )),
      ref,
      info = "Missing column addition works"
    )
    {
      zero_bound <- nlmixr:::nlmixrBoundsTemplate[1:2, ]
      zero_bound$ntheta <- 1:2
      zero_bound$lower <- c(-Inf, 0)
      zero_bound$est <- c(-5, 5)
      zero_bound$upper <- c(0, Inf)
      expect_equal(
        as.data.frame(nlmixr:::as.nlmixrBounds(zero_bound)[, c("lower", "upper")]),
        data.frame(
          lower = c(-Inf, sqrt(.Machine$double.eps)),
          upper = c(-sqrt(.Machine$double.eps), Inf)
        ),
        # row.names will not be equal
        check.attributes = FALSE
      )
    }
  })

  test_that("bounds are extracted correctly", {
    ref <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = c(1, 2, 3, 4, 5, 6, 7, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9, 10, 11, 12, 13, 14),
          neta1 = c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12, 12, NA, NA, NA, NA, NA, NA),
          neta2 = c(NA, NA, NA, NA, NA, NA, NA, NA, 1, 2, 3, 4, 5, 5, 6, 7, 7, 8, 7, 8, 9, 10, 10, 11, 10, 11, 12, NA, NA, NA, NA, NA, NA),
          name = c("a", "b", "c", "d", NA, NA, NA, NA, "et1", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "et2", "(et3,et2)", "et3", "(et4,et2)", "(et4,et3)", "et4", NA, NA, NA, "a5", "a6", "a7"),
          lower = c(1.49011611938477e-08, 1.49011611938477e-08, -Inf, -Inf, 1.49011611938477e-08, 1.49011611938477e-08, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 9, 11, -Inf, 9, 11),
          est = c(1, 3, 4, 4, 1, 1, 1, 1, 10, 20, 30, 40, 40, 0.1, 20, 40, 0.1, 20, 0.1, 0.1, 30, 40, 0.1, 20, 0.1, 0.1, 30, 8, 10, 12, 8, 10, 12),
          upper = c(2, Inf, Inf, Inf, 2, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 13, Inf, Inf, 13),
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
          err = NA_character_,
          label = c("A", NA, NA, NA, "e", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "labels", NA, NA, NA),
          condition = c(NA, NA, NA, NA, NA, NA, NA, NA, "ID", "ID", "ID", "ID", "STUD", "STUD", "STUD", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", "ID", NA, NA, NA, NA, NA, NA),
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    testbounds <- function() {
      a <- c(0, 1, 2) # A
      b <- c(0, 3)
      c <- 4
      d <- c(4)
      c(0, 1, 2) # e
      c(0, 1)
      c(1)
      1
      et1 ~ 10
      ~20
      ~30
      ~ c(40)
      ~ c(
        40,
        0.1, 20
      ) | STUD
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
      et2 + et3 + et4 ~ c(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )

      ## new test fixed parameters...
      c(8, fixed)
      c(9, 10, fixed)
      c(11, 12, 13, fixed) # labels

      a5 <- c(8, fixed)
      a6 <- c(9, 10, fixed)
      a7 <- c(11, 12, 13, fixed)
    }

    expect_equal(nlmixrBounds(testbounds), ref)
  })

  test_that("Theta Bounds above 5 don't work", {
    bnd <- function() {
      c(1, 2, 3, 4, 5)
    }

    bnda <- function() {
      a <- c(1, 2, 3, 4, 5)
    }
    bndb <- function() {
      a <- c(1, 2, 3, 4, 5)
    }

    expect_error(nlmixrBounds(bnd), regexp = "Syntax is not supported for thetas: c(1, 2, 3, 4, 5)", fixed = TRUE)
    expect_error(nlmixrBounds(bnda), regexp = "Syntax is not supported for thetas: c(1, 2, 3, 4, 5)", fixed = TRUE)
    expect_error(nlmixrBounds(bndb), regexp = "Syntax is not supported for thetas: c(1, 2, 3, 4, 5)", fixed = TRUE)
  })

  test_that("Theta Bounds above 4 don't work", {
    bnd <- function() {
      c(1, 2, 3, 4)
    }
    bnda <- function() {
      a <- c(1, 2, 3, 4)
    }
    bndb <- function() {
      a <- c(1, 2, 3, 4)
    }
    expect_error(nlmixrBounds(bnd), regexp = "Syntax is not supported for thetas: c(1, 2, 3, 4)", fixed = TRUE)
    expect_error(nlmixrBounds(bnda), regexp = "Syntax is not supported for thetas: c(1, 2, 3, 4)", fixed = TRUE)
    expect_error(nlmixrBounds(bndb), regexp = "Syntax is not supported for thetas: c(1, 2, 3, 4)", fixed = TRUE)
  })

  test_that("Bad Lower trianglar matrices throw errors.", {
    bnd1 <- function() {
      ~ c(1)
    }
    ref1 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = 1,
          neta2 = 1,
          name = NA_character_,
          lower = -Inf,
          est = 1,
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd2 <- function() {
      ~ c(1, 2)
    }

    bnd3 <- function() {
      ~ c(
        1,
        2, 3
      )
    }

    ref3 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2),
          neta2 = c(1, 1, 2),
          name = NA_character_,
          lower = -Inf,
          est = c(1, 2, 3),
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd4 <- function() {
      ~ c(
        1,
        2, 3,
        4
      )
    }

    bnd5 <- function() {
      ~ c(
        1,
        2, 3,
        4, 5
      )
    }

    bnd6 <- function() {
      ~ c(
        1,
        2, 3,
        4, 5, 6
      )
    }

    ref6 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = 1:6,
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd7 <- function() {
      ~ c(
        1,
        2, 3,
        4, 5, 6,
        7
      )
    }
    bnd8 <- function() {
      ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8
      )
    }
    bnd9 <- function() {
      ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8, 9
      )
    }
    bnd10 <- function() {
      ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8, 9, 10
      )
    }

    ref10 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4),
          neta2 = c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4),
          name = NA_character_,
          lower = -Inf,
          est = 1:10,
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    expect_equal(nlmixrBounds(bnd1), ref1)
    expect_error(nlmixrBounds(bnd2), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2)", fixed = TRUE)
    expect_equal(nlmixrBounds(bnd3), ref3)
    expect_error(nlmixrBounds(bnd4), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4)", fixed = TRUE)
    expect_error(nlmixrBounds(bnd5), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5)", fixed = TRUE)
    expect_equal(nlmixrBounds(bnd6), ref6)
    expect_error(nlmixrBounds(bnd7), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7)", fixed = TRUE)
    expect_error(nlmixrBounds(bnd8), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8)", fixed = TRUE)
    expect_error(nlmixrBounds(bnd9), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8, 9)", fixed = TRUE)
    expect_equal(nlmixrBounds(bnd10), ref10)
  })

  test_that("Bad Lower trianglar matrices (with labels) throw errors.", {
    bnd1 <- function() {
      eta1 ~ c(1)
    }

    ref1 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = 1,
          neta2 = 1,
          name = "eta1",
          lower = -Inf,
          est = 1,
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd2 <- function() {
      eta1 ~ c(1, 2)
    }

    bnd3 <- function() {
      eta1 + eta2 ~ c(
        1,
        2, 3
      )
    }

    ref3 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2),
          neta2 = c(1, 1, 2),
          name = c("eta1", "(eta2,eta1)", "eta2"),
          lower = -Inf,
          est = 1:3,
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd4 <- function() {
      eta1 + eta2 ~ c(
        1,
        2, 3,
        4
      )
    }

    bnd5 <- function() {
      eta1 + eta2 ~ c(
        1,
        2, 3,
        4, 5
      )
    }

    bnd6 <- function() {
      eta1 + eta2 + eta3 ~ c(
        1,
        2, 3,
        4, 5, 6
      )
    }

    ref6 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = 1:6,
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd7 <- function() {
      eta1 + eta2 + eta3 ~ c(
        1,
        2, 3,
        4, 5, 6,
        7
      )
    }
    bnd8 <- function() {
      eta1 + eta2 + eta3 ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8
      )
    }
    bnd9 <- function() {
      eta1 + eta2 + eta3 ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8, 9
      )
    }
    bnd10 <- function() {
      eta1 + eta2 + eta3 + eta4 ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8, 9, 10
      )
    }

    ref10 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4),
          neta2 = c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3", "(eta4,eta1)", "(eta4,eta2)", "(eta4,eta3)", "eta4"),
          lower = -Inf,
          est = 1:10,
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    expect_equal(nlmixrBounds(bnd1), ref1)
    expect_error(nlmixrBounds(bnd2), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2)", fixed = TRUE)
    expect_equal(nlmixrBounds(bnd3), ref3)
    expect_error(nlmixrBounds(bnd4), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4)", fixed = TRUE)
    expect_error(nlmixrBounds(bnd5), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5)", fixed = TRUE)
    expect_equal(nlmixrBounds(bnd6), ref6)
    expect_error(nlmixrBounds(bnd7), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7)", fixed = TRUE)
    expect_error(nlmixrBounds(bnd8), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8)", fixed = TRUE)
    expect_error(nlmixrBounds(bnd9), regexp = "incorrect lower triangular matrix dimensions: ~c(1, 2, 3, 4, 5, 6, 7, 8, 9)", fixed = TRUE)
    expect_equal(nlmixrBounds(bnd10), ref10)
  })

  test_that("Number of eta variables must match", {
    bnd10.a <- function() {
      eta1 + eta2 + eta3 + eta4 + eta5 ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8, 9, 10
      )
    }

    bnd10.b <- function() {
      eta1 + eta2 + eta3 ~ c(
        1,
        2, 3,
        4, 5, 6,
        7, 8, 9, 10
      )
    }
    expect_error(nlmixrBounds(bnd10.a), regexp = "omega assignment left handed side must match lower triangular matrix size", fixed = TRUE)
    expect_error(nlmixrBounds(bnd10.b), regexp = "omega assignment left handed side must match lower triangular matrix size", fixed = TRUE)
  })

  test_that("Comments inside bounds are not supported!", {
    bnd3 <- function() {
      eta1 + eta2 ~ c(
        1, # Comment here.
        2, 3
      )
    }

    expect_error(nlmixrBounds(bnd3), regexp = "error parsing bounds: possible (unsupported) comment/condition inside bounds", fixed = TRUE)
  })

  test_that("Conditional statments are supported correctly.", {
    bnd1 <- function() {
      eta0 ~ 0.3
      eta1 + eta2 ~ c(
        1,
        2, 3
      ) | STUD
      ~ c(
        1,
        2, 3
      )
    }

    ref <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = as.numeric(c(NA, NA, NA, NA, NA, NA, NA)),
          neta1 = c(1, 2, 3, 3, 4, 5, 5),
          neta2 = c(1, 2, 2, 3, 4, 4, 5),
          name = c("eta0", "eta1", "(eta2,eta1)", "eta2", NA, NA, NA),
          lower = -Inf,
          est = c(0.3, 1, 2, 3, 1, 2, 3),
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = c("ID", "STUD", "STUD", "STUD", "ID", "ID", "ID"),
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    expect_equal(nlmixrBounds(bnd1), ref)
  })

  test_that("Theta fix fixed are reasonable", {
    ref1 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = 1:4,
          neta1 = NA_real_,
          neta2 = NA_real_,
          name = c("a", "b", "c", "d"),
          lower = c(1.49011611938477e-08, 1.49011611938477e-08, -Inf, -Inf),
          est = c(1, 3, 4, 4),
          upper = c(2, Inf, Inf, Inf),
          fix = c(TRUE, TRUE, FALSE, TRUE),
          err = NA_character_,
          label = c("A", NA, NA, NA),
          condition = NA_character_,
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd1 <- function() {
      a <- fix(0, 1, 2) # A
      b <- fix(0, 3)
      c <- 4
      d <- fix(4)
    }

    bnd2 <- function() {
      a <- FIX(0, 1, 2) # A
      b <- FIX(0, 3)
      c <- 4
      d <- FIX(4)
    }

    bnd3 <- function() {
      a <- fixed(0, 1, 2) # A
      b <- fixed(0, 3)
      c <- 4
      d <- fixed(4)
    }

    bnd4 <- function() {
      a <- FIXED(0, 1, 2) # A
      b <- FIXED(0, 3)
      c <- 4
      d <- FIXED(4)
    }

    bnd5 <- function() {
      a <- c(0, fix(1), 2) # A
      b <- c(0, fix(3))
      c <- 4
      d <- fix(4)
    }

    bnd6 <- function() {
      a <- c(0, FIX(1), 2) # A
      b <- c(0, FIX(3))
      c <- 4
      d <- FIX(4)
    }

    bnd7 <- function() {
      a <- c(0, fixed(1), 2) # A
      b <- c(0, fixed(3))
      c <- 4
      d <- fixed(4)
    }

    bnd8 <- function() {
      a <- c(0, FIXED(1), 2) # A
      b <- c(0, FIXED(3))
      c <- 4
      d <- FIXED(4)
    }
    expect_equal(nlmixrBounds(bnd1), ref1)
    expect_equal(nlmixrBounds(bnd2), ref1)
    expect_equal(nlmixrBounds(bnd3), ref1)
    expect_equal(nlmixrBounds(bnd4), ref1)
    expect_equal(nlmixrBounds(bnd5), ref1)
    expect_equal(nlmixrBounds(bnd6), ref1)
    expect_equal(nlmixrBounds(bnd7), ref1)
    expect_equal(nlmixrBounds(bnd8), ref1)

    ref2 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = 1,
          neta1 = NA_real_,
          neta2 = NA_real_,
          name = "a",
          lower = 1.49011611938477e-08,
          est = 2,
          upper = 3,
          fix = TRUE,
          err = NA_character_,
          label = NA_character_,
          condition = NA_character_,
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    bnd1 <- function() {
      a <- fixed(0, 2, 3)
    }
    bnd2 <- function() {
      a <- c(0, fixed(2), 3)
    }
    expect_equal(nlmixrBounds(bnd1), ref2)
    expect_equal(nlmixrBounds(bnd2), ref2)
  })

  test_that("Total ETA fixed (unnamed)", {
    ref1 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = FALSE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd1 <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    ref2 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = TRUE,
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd2 <- function() {
      ~ fix(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd3 <- function() {
      ~ fixed(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd4 <- function() {
      ~ FIX(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd5 <- function() {
      ~ FIXED(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }
    expect_equal(nlmixrBounds(bnd1), ref1)
    expect_equal(nlmixrBounds(bnd2), ref2)
    expect_equal(nlmixrBounds(bnd3), ref2)
    expect_equal(nlmixrBounds(bnd4), ref2)
    expect_equal(nlmixrBounds(bnd5), ref2)
  })

  test_that("Total ETA fixed (unnamed)", {
    ref6 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd6 <- function() {
      ~ c(
        fix(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    ref7 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd7 <- function() {
      ~ c(
        40,
        fix(0.1), 20,
        0.1, 0.1, 30
      )
    }

    ref8 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd8 <- function() {
      ~ c(
        40,
        0.1, fix(20),
        0.1, 0.1, 30
      )
    }

    ref9 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd9 <- function() {
      ~ c(
        40,
        0.1, 20,
        fix(0.1), 0.1, 30
      )
    }

    ref10 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd10 <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, fix(0.1), 30
      )
    }

    ref11 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd11 <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, fix(30)
      )
    }


    bnd6 <- function() {
      ~ c(
        fix(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    ref7 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd7 <- function() {
      ~ c(
        40,
        fix(0.1), 20,
        0.1, 0.1, 30
      )
    }

    ref8 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd8 <- function() {
      ~ c(
        40,
        0.1, fix(20),
        0.1, 0.1, 30
      )
    }

    ref9 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd9 <- function() {
      ~ c(
        40,
        0.1, 20,
        fix(0.1), 0.1, 30
      )
    }

    ref10 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd10 <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, fix(0.1), 30
      )
    }

    ref11 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd11 <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, fix(30)
      )
    }

    expect_equal(nlmixrBounds(bnd6), ref6)
    expect_equal(nlmixrBounds(bnd7), ref7)
    expect_equal(nlmixrBounds(bnd8), ref8)
    expect_equal(nlmixrBounds(bnd9), ref9)
    expect_equal(nlmixrBounds(bnd10), ref10)
    expect_equal(nlmixrBounds(bnd11), ref11)
  })

  test_that("Total ETA fixed (unnamed) #a", {
    ref6 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref7 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref8 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref9 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref10 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref11 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd6a <- function() {
      ~ c(
        fixed(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd7a <- function() {
      ~ c(
        40,
        fixed(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd8a <- function() {
      ~ c(
        40,
        0.1, fixed(20),
        0.1, 0.1, 30
      )
    }

    bnd9a <- function() {
      ~ c(
        40,
        0.1, 20,
        fixed(0.1), 0.1, 30
      )
    }

    bnd10a <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, fixed(0.1), 30
      )
    }

    bnd11a <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, fixed(30)
      )
    }

    bnd6a <- function() {
      ~ c(
        fixed(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd7a <- function() {
      ~ c(
        40,
        fixed(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd8a <- function() {
      ~ c(
        40,
        0.1, fixed(20),
        0.1, 0.1, 30
      )
    }

    bnd9a <- function() {
      ~ c(
        40,
        0.1, 20,
        fixed(0.1), 0.1, 30
      )
    }

    bnd10a <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, fixed(0.1), 30
      )
    }

    bnd11a <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, fixed(30)
      )
    }
    expect_equal(nlmixrBounds(bnd6a), ref6)
    expect_equal(nlmixrBounds(bnd7a), ref7)
    expect_equal(nlmixrBounds(bnd8a), ref8)
    expect_equal(nlmixrBounds(bnd9a), ref9)
    expect_equal(nlmixrBounds(bnd10a), ref10)
    expect_equal(nlmixrBounds(bnd11a), ref11)
  })

  test_that("Total ETA fixed (unnamed) #b", {
    ref6 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref7 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref8 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref9 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref10 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref11 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd6b <- function() {
      ~ c(
        FIX(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd7b <- function() {
      ~ c(
        40,
        FIX(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd8b <- function() {
      ~ c(
        40,
        0.1, FIX(20),
        0.1, 0.1, 30
      )
    }

    bnd9b <- function() {
      ~ c(
        40,
        0.1, 20,
        FIX(0.1), 0.1, 30
      )
    }

    bnd10b <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, FIX(0.1), 30
      )
    }

    bnd11b <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, FIX(30)
      )
    }

    bnd6b <- function() {
      ~ c(
        FIX(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd7b <- function() {
      ~ c(
        40,
        FIX(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd8b <- function() {
      ~ c(
        40,
        0.1, FIX(20),
        0.1, 0.1, 30
      )
    }

    bnd9b <- function() {
      ~ c(
        40,
        0.1, 20,
        FIX(0.1), 0.1, 30
      )
    }

    bnd10b <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, FIX(0.1), 30
      )
    }

    bnd11b <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, FIX(30)
      )
    }
    expect_equal(nlmixrBounds(bnd6b), ref6)
    expect_equal(nlmixrBounds(bnd7b), ref7)
    expect_equal(nlmixrBounds(bnd8b), ref8)
    expect_equal(nlmixrBounds(bnd9b), ref9)
    expect_equal(nlmixrBounds(bnd10b), ref10)
    expect_equal(nlmixrBounds(bnd11b), ref11)
  })

  test_that("Total ETA fixed (unnamed) #c", {
    ref6 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref7 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref8 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref9 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref10 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )
    ref11 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = NA_character_,
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd6c <- function() {
      ~ c(
        FIXED(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd7c <- function() {
      ~ c(
        40,
        FIXED(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd8c <- function() {
      ~ c(
        40,
        0.1, FIXED(20),
        0.1, 0.1, 30
      )
    }

    bnd9c <- function() {
      ~ c(
        40,
        0.1, 20,
        FIXED(0.1), 0.1, 30
      )
    }

    bnd10c <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, FIXED(0.1), 30
      )
    }

    bnd11c <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, FIXED(30)
      )
    }

    bnd6c <- function() {
      ~ c(
        FIXED(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd7c <- function() {
      ~ c(
        40,
        FIXED(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd8c <- function() {
      ~ c(
        40,
        0.1, FIXED(20),
        0.1, 0.1, 30
      )
    }

    bnd9c <- function() {
      ~ c(
        40,
        0.1, 20,
        FIXED(0.1), 0.1, 30
      )
    }

    bnd10c <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, FIXED(0.1), 30
      )
    }

    bnd11c <- function() {
      ~ c(
        40,
        0.1, 20,
        0.1, 0.1, FIXED(30)
      )
    }

    expect_equal(nlmixrBounds(bnd6c), ref6)
    expect_equal(nlmixrBounds(bnd7c), ref7)
    expect_equal(nlmixrBounds(bnd8c), ref8)
    expect_equal(nlmixrBounds(bnd9c), ref9)
    expect_equal(nlmixrBounds(bnd10c), ref10)
    expect_equal(nlmixrBounds(bnd11c), ref11)
  })

  test_that("Total ETA FIXED (named)", {
    ref12 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd12 <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    ref13 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd13 <- function() {
      eta1 + eta2 + eta3 ~ fix(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd14 <- function() {
      eta1 + eta2 + eta3 ~ fixed(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd15 <- function() {
      eta1 + eta2 + eta3 ~ FIX(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd16 <- function() {
      eta1 + eta2 + eta3 ~ FIXED(
        40,
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    test_that("Total ETA fixed (named)", {
      expect_equal(nlmixrBounds(bnd12), ref12)
      expect_equal(nlmixrBounds(bnd13), ref13)
      expect_equal(nlmixrBounds(bnd14), ref13)
      expect_equal(nlmixrBounds(bnd15), ref13)
      expect_equal(nlmixrBounds(bnd16), ref13)
    })

    ref17 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd17 <- function() {
      eta1 + eta2 + eta3 ~ c(
        fix(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }


    ref18 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd18 <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        fix(0.1), 20,
        0.1, 0.1, 30
      )
    }

    ref19 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd19 <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, fix(20),
        0.1, 0.1, 30
      )
    }

    ref20 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd20 <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        fix(0.1), 0.1, 30
      )
    }

    ref21 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd21 <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, fix(0.1), 30
      )
    }

    ref22 <-
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = NA_real_,
          neta1 = c(1, 2, 2, 3, 3, 3),
          neta2 = c(1, 1, 2, 1, 2, 3),
          name = c("eta1", "(eta2,eta1)", "eta2", "(eta3,eta1)", "(eta3,eta2)", "eta3"),
          lower = -Inf,
          est = c(40, 0.1, 20, 0.1, 0.1, 30),
          upper = Inf,
          fix = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          err = NA_character_,
          label = NA_character_,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      )

    bnd22 <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, 0.1, fix(30)
      )
    }

    bnd17a <- function() {
      eta1 + eta2 + eta3 ~ c(
        FIX(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd18a <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        FIX(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd19a <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, FIX(20),
        0.1, 0.1, 30
      )
    }

    bnd20a <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        FIX(0.1), 0.1, 30
      )
    }

    bnd21a <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, FIX(0.1), 30
      )
    }

    bnd22a <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, 0.1, FIX(30)
      )
    }

    bnd17b <- function() {
      eta1 + eta2 + eta3 ~ c(
        fixed(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd18b <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        fixed(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd19b <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, fixed(20),
        0.1, 0.1, 30
      )
    }

    bnd20b <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        fixed(0.1), 0.1, 30
      )
    }

    bnd21b <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, fixed(0.1), 30
      )
    }

    bnd22b <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, 0.1, fixed(30)
      )
    }

    bnd17c <- function() {
      eta1 + eta2 + eta3 ~ c(
        FIX(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd18c <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        FIX(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd19c <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, FIX(20),
        0.1, 0.1, 30
      )
    }

    bnd20c <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        FIX(0.1), 0.1, 30
      )
    }

    bnd21c <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, FIX(0.1), 30
      )
    }

    bnd22c <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, 0.1, FIX(30)
      )
    }

    bnd17d <- function() {
      eta1 + eta2 + eta3 ~ c(
        FIXED(40),
        0.1, 20,
        0.1, 0.1, 30
      )
    }

    bnd18d <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        FIXED(0.1), 20,
        0.1, 0.1, 30
      )
    }

    bnd19d <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, FIXED(20),
        0.1, 0.1, 30
      )
    }

    bnd20d <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        FIXED(0.1), 0.1, 30
      )
    }

    bnd21d <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, FIXED(0.1), 30
      )
    }

    bnd22d <- function() {
      eta1 + eta2 + eta3 ~ c(
        40,
        0.1, 20,
        0.1, 0.1, FIXED(30)
      )
    }
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

  test_that("Invalid bounds raise errors", {
    f1 <- function() {
      lCl <- c(5, 5, 5) # A
    }

    f2 <- function() {
      lCl <- c(0, -1.3) # lCl
    }

    f3 <- function() {
      lCl <- c(0, -1.3, -10) # lCl
    }

    f4 <- function() {
      lCl <- c(0, 5, 5) # A
    }

    f5 <- function() {
      lCl <- c(0, 0, 5) # A
    }

    f6 <- function() {
      lCl <- c(5, 5) # A
    }

    f7 <- function() {
      lCl < 3
    }

    expect_error(nlmixrBounds(f1), regexp = "consider fixing these:\n     lCl = fixed(5)", fixed = TRUE)
    expect_error(nlmixrBounds(f2), regexp = "reorder bounds:\n     lCl = c(-1.3, 0)", fixed = TRUE)
    expect_error(nlmixrBounds(f3), regexp = "reorder bounds:\n     lCl = c(-10, -1.3, 0)", fixed = TRUE)
    expect_error(nlmixrBounds(f4), regexp = "consider fixing these:\n     lCl = fixed(5)", fixed = TRUE)
    expect_error(nlmixrBounds(f5), regexp = "consider fixing these:\n     lCl = fixed(0)", fixed = TRUE)
    expect_error(nlmixrBounds(f6), regexp = "consider fixing these:\n     lCl = fixed(5)", fixed = TRUE)
    expect_error(nlmixrBounds(f7), regexp = "invalid call in initial conditions: lCl < 3", fixed = TRUE)
  })

  # nlmixrBoundsParser ####

  test_that("nlmixrBoundsParser", {
    expect_equal(
      nlmixrBoundsParser(
        function() {
          ({
            ({
              a <- 1
              b <- 2
            })
          })
          {
            c <- 3
          }
          d <- 4
        }
      ),
      list(
        list(
          operation = c("assign", "theta"),
          varname = "a",
          value = 1
        ),
        list(
          operation = c("assign", "theta"),
          varname = "b",
          value = 2
        ),
        list(
          operation = c("assign", "theta"),
          varname = "c",
          value = 3
        ),
        list(
          operation = c("assign", "theta"),
          varname = "d",
          value = 4
        )
      ),
      info = "Nested assignments are unnested"
    )
  })

  # nlmixrBoundsParserOmega ####

  test_that("Implicitly test nlmixrBoundsParserOmega", {
    expect_equal(
      nlmixrBounds(function() {
        ~1
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          neta1 = 1,
          neta2 = 1,
          name = NA_character_,
          lower = -Inf,
          est = 1,
          upper = Inf,
          fix = FALSE,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Unnamed omega scalar assignment"
    )
    expect_equal(
      nlmixrBounds(function() {
        a ~ 1
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          neta1 = 1,
          neta2 = 1,
          name = "a",
          lower = -Inf,
          est = 1,
          upper = Inf,
          fix = FALSE,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Named omega scalar assignment"
    )
    expect_equal(
      expect_warning(
        nlmixrBounds(function() {
          a ~ cor(1)
        }),
        regexp = "'cor(...)' with a single value is ignored: ~cor(1)",
        fixed = TRUE
      ),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          neta1 = 1,
          neta2 = 1,
          name = "a",
          lower = -Inf,
          est = 1,
          upper = Inf,
          fix = FALSE,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Named scalar omega assignment with correlation"
    )
    expect_equal(
      nlmixrBounds(function() {
        a + b ~ cor(2, -0.5, 3)
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          neta1 = c(1, 2, 2),
          neta2 = c(1, 1, 2),
          name = c("a", "(b,a)", "b"),
          lower = -Inf,
          est = c(4, -3, 9),
          upper = Inf,
          fix = FALSE,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Named vector omega assignment with correlation"
    )
    expect_equal(
      nlmixrBounds(function() {
        a + b ~ fixed(cor(2, -0.5, 3))
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          neta1 = c(1, 2, 2),
          neta2 = c(1, 1, 2),
          name = c("a", "(b,a)", "b"),
          lower = -Inf,
          est = c(4, -3, 9),
          upper = Inf,
          fix = TRUE,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Named vector omega assignment with correlation, all fixed with an outer function"
    )
    expect_equal(
      nlmixrBounds(function() {
        a + b ~ cor(fixed(2, -0.5, 3))
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          neta1 = c(1, 2, 2),
          neta2 = c(1, 1, 2),
          name = c("a", "(b,a)", "b"),
          lower = -Inf,
          est = c(4, -3, 9),
          upper = Inf,
          fix = TRUE,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Named vector omega assignment with correlation, all fixed with an inner function"
    )
    expect_equal(
      nlmixrBounds(function() {
        a + b ~ cor(fixed(2), fixed(-0.5), fixed(3))
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          neta1 = c(1, 2, 2),
          neta2 = c(1, 1, 2),
          name = c("a", "(b,a)", "b"),
          lower = -Inf,
          est = c(4, -3, 9),
          upper = Inf,
          fix = TRUE,
          condition = "ID",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Named vector omega assignment with correlation, all fixed with individual functions (unusual, but acceptable)"
    )
    expect_error(
      nlmixrBounds(function() {
        a + b ~ cor(2, fixed(-0.5), fixed(3))
      }),
      regexp = "either all or none of the elements may be fixed with cor(...): ~cor(2, fixed(-0.5), fixed(3))",
      fixed = TRUE,
      info = "Named vector omega assignment with correlation, some fixed"
    )
  })

  # nlmixrBoundsParserAttribute ####

  test_that("nlmixrBoundsParserAttribute backTransform", {
    expect_equal(
      nlmixrBounds(function() {
        0.1
        backTransform(exp())
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = 1,
          lower = -Inf,
          est = 0.1,
          upper = Inf,
          fix = FALSE,
          backTransform = "exp()",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Simple back-transform"
    )
    expect_equal(
      nlmixrBounds(function() {
        0.1
        backTransform(exp)
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = 1,
          lower = -Inf,
          est = 0.1,
          upper = Inf,
          fix = FALSE,
          backTransform = "exp()",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Simple back-transform as a name is replaced by the function name"
    )
    expect_equal(
      nlmixrBounds(function() {
        0.1
        backTransform("exp")
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = 1,
          lower = -Inf,
          est = 0.1,
          upper = Inf,
          fix = FALSE,
          backTransform = "exp()",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Simple back-transform as a character string is replaced by the function name"
    )
    expect_equal(
      expect_warning(
        nlmixrBounds(function() {
          0.1
          backTransform("exp")
          backTransform("log")
        }),
        regexp = 'only last backTransform used: backTransform("log")',
        fixed = TRUE
      ),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = 1,
          lower = -Inf,
          est = 0.1,
          upper = Inf,
          fix = FALSE,
          backTransform = "log()",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "Multiple backTransform()s keep the last one"
    )
    expect_equal(
      nlmixrBounds(function() {
        0.1
        backTransform(function(x) x^2)
      }),
      nlmixr:::as.nlmixrBounds(
        data.frame(
          ntheta = 1,
          lower = -Inf,
          est = 0.1,
          upper = Inf,
          fix = FALSE,
          backTransform = "function(x) x^2",
          stringsAsFactors = FALSE
        ),
        addMissingCols = TRUE
      ),
      info = "A function may be defined within the backTransform"
    )
    expect_error(
      nlmixrBounds(function() {
        0.1
        backTransform(exp(), foo())
      }),
      regexp = "'backTransform()' must have zero or one arguments: backTransform(exp(), foo())",
      fixed = TRUE,
      info = "Multiple arguments to backTransform() are not allowed."
    )
  })

  # nlmixrBoundsValueFixed ####

  test_that("nlmixrBoundsValueFixed", {
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~1)[[2]]),
      list(value = 1, fixed = FALSE)
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1))[[2]]),
      list(value = 1, fixed = FALSE)
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1, 2))[[2]]),
      list(value = c(1, 2), fixed = rep(FALSE, 2))
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1, 2, 3))[[2]]),
      list(value = c(1, 2, 3), fixed = rep(FALSE, 3))
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1, fixed))[[2]]),
      list(value = 1, fixed = TRUE)
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1, 2, fixed))[[2]]),
      list(value = c(1, 2), fixed = rep(TRUE, 2))
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1, 2, 3, fixed))[[2]]),
      list(value = c(1, 2, 3), fixed = rep(TRUE, 3))
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1, fixed(2), 3, fixed))[[2]]),
      list(value = c(1, 2, 3), fixed = rep(TRUE, 3)),
      info = "Fixed is specified two ways, but they are not in conflict"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(1, fixed(2), 3))[[2]]),
      list(value = c(1, 2, 3), fixed = c(FALSE, TRUE, FALSE))
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(fixed(1), 2, 3))[[2]]),
      list(value = c(1, 2, 3), fixed = c(TRUE, FALSE, FALSE))
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(fixed(1), 2, 3))[[2]]),
      list(value = c(1, 2, 3), fixed = c(TRUE, FALSE, FALSE))
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(fixed(1), log(2), 3))[[2]]),
      list(value = c(1, log(2), 3), fixed = c(TRUE, FALSE, FALSE)),
      info = "Function evaluation works (though it may have issues related to environment precedence). (Fix #253)"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(log(0), log(1.5), log(20)))[[2]]),
      list(value = log(c(0, 1.5, 20)), fixed = rep(FALSE, 3)),
      info = "Function evaluation works (though it may have issues related to environment precedence). (Fix #253)"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ FIX(log(0), log(1.5), log(20)))[[2]]),
      list(value = log(c(0, 1.5, 20)), fixed = rep(TRUE, 3)),
      info = "Function evaluation works (though it may have issues related to environment precedence). (Fix #253)"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueFixed((~ c(FIX(log(0), 1 / log(1.5)), 1 / log(20)))[[2]]),
      list(value = c(log(0), 1 / log(1.5), 1 / log(20)), fixed = c(TRUE, TRUE, FALSE)),
      info = "Function evaluation works and arbitrary complexity may be within the fixed() call (or outside of it). (Fix #253)"
    )
    expect_error(
      expect_warning(
        nlmixr:::nlmixrBoundsValueFixed((~ FIX(sqrt(-1)))[[2]]),
        regexp = "NaNs produced"
      ),
      regexp = "NaN values in initial condition: FIX(sqrt(-1))",
      fixed = TRUE,
      info = "Invalid math stops execution"
    )
    expect_error(
      nlmixr:::nlmixrBoundsValueFixed((~ FIX("A"))[[2]]),
      regexp = 'non-numeric values in initial condition: FIX("A")',
      fixed = TRUE,
      info = "Values must be numbers"
    )
    expect_error(
      nlmixr:::nlmixrBoundsValueFixed((~a)[[2]]),
      regexp = "error parsing initial condition 'a': object 'a' not found",
      fixed = TRUE,
      info = "No variable substitutions are performed for parsing."
    )
  })

  # nlmixrBoundsReplaceFixed ####

  test_that("nlmixrBoundsReplaceFixed, testing replacement of fixed names within calls", {
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~a)[[2]]),
      list(
        call = (~a)[[2]],
        fixed = FALSE
      )
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ c(1, fixed))[[2]]),
      list(
        call = (~ c(1))[[2]],
        fixed = TRUE
      )
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ c(1, c(1, fixed)))[[2]]),
      list(
        call = (~ c(1, c(1)))[[2]],
        fixed = TRUE
      )
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~1)[[2]]),
      list(
        call = 1,
        fixed = FALSE
      )
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ fixed(1))[[2]], replacementFun = "c"),
      list(
        call = (~ c(1))[[2]],
        fixed = FALSE
      )
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ fixed(1))[[2]], replacementFun = "c"),
      list(
        call = (~ c(1))[[2]],
        fixed = FALSE
      )
    )
    # This is weird syntax to use, but it is logically okay
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ c(1, 1, c(fixed)))[[2]]),
      list(
        call = (~ c(1, 1, c()))[[2]],
        fixed = TRUE
      ),
      info = "Fixed can only be at the end of a vector of values, and detection of that works even when it is in a sub-expression."
    )
    expect_error(
      nlmixr:::nlmixrBoundsReplaceFixed((~ c(fixed, 1))[[2]]),
      regexp = "'fixed' may only be the last item in a list: c(fixed, 1)",
      fixed = TRUE,
      info = "Fixed can only be at the end of a vector of values"
    )
    expect_error(
      nlmixr:::nlmixrBoundsReplaceFixed((~ c(1, c(fixed), 1))[[2]]),
      regexp = "'fixed' may only be the last item in a list: c(1, c(fixed), 1)",
      fixed = TRUE,
      info = "Fixed can only be at the end of a vector of values, and detection of that works even when it is at the end of its sub-expression, but it is not at the overall-end of the expression."
    )
  })

  test_that("nlmixrBoundsReplaceFixed, testing replacement of fixed function calls within calls", {
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ fixed(a))[[2]]),
      list(
        call = (~ fixed(a))[[2]],
        fixed = FALSE
      ),
      info = "`fixed()` is returned unchanged"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ fix(a))[[2]]),
      list(
        call = (~ fixed(a))[[2]],
        fixed = FALSE
      ),
      info = "`fix()` is changed to `fixed()`"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ FIX(a))[[2]]),
      list(
        call = (~ fixed(a))[[2]],
        fixed = FALSE
      ),
      info = "`FIX()` is changed to `fixed()`"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ FIXED(a))[[2]]),
      list(
        call = (~ fixed(a))[[2]],
        fixed = FALSE
      ),
      info = "`FIXED()` is changed to `fixed()`"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ 1 / FIX(a))[[2]]),
      list(
        call = (~ 1 / fixed(a))[[2]],
        fixed = FALSE
      ),
      info = "`FIX()` is changed to `fixed()` inside another expression"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsReplaceFixed((~ 1 / FIX(a))[[2]], replacementFun = "c"),
      list(
        call = (~ 1 / c(a))[[2]],
        fixed = FALSE
      ),
      info = "`FIX()` is changed to `c()` inside another expression to allow for arbitrary calculations as part of initial conditions."
    )
  })

  # nlmixrBoundsReplaceCor ####

  test_that("", {
    expect_equal(
      nlmixr:::nlmixrBoundsValueCor(x = (~ cor(1, 2, 3))[[2]]),
      list(
        value = 1:3,
        fixed = rep(FALSE, 3),
        cor = rep(TRUE, 3)
      )
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueCor(x = (~ cor(1, fixed(2), 3))[[2]]),
      list(
        value = 1:3,
        fixed = c(FALSE, TRUE, FALSE),
        cor = rep(TRUE, 3)
      )
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueCor(x = (~ fixed(cor(1, fixed(2), 3)))[[2]]),
      list(
        value = 1:3,
        fixed = rep(TRUE, 3),
        cor = rep(TRUE, 3)
      ),
      info = "Unusual syntax works"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueCor(x = (~ c(1, fixed(2), cor(3)))[[2]]),
      list(
        value = 1:3,
        fixed = c(FALSE, TRUE, FALSE),
        cor = c(FALSE, FALSE, TRUE)
      ),
      info = "Legitimacy of syntax checking will be confirmed elsewhere"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsValueCor(x = (~ c(1, 2, cor(3), fixed))[[2]]),
      list(
        value = 1:3,
        fixed = c(TRUE, TRUE, TRUE),
        cor = c(FALSE, FALSE, TRUE)
      ),
      info = "Trailing 'fixed' is handled correctly"
    )
  })

  # nlmixrBoundsPrepareFun ####

  test_that("preparation of the function for bound extraction", {
    expect_equal(
      nlmixr:::nlmixrBoundsPrepareFun(function() {
        1
      }),
      function() {
        1
      }
    )
    expect_equal(
      expect_message(
        nlmixr:::nlmixrBoundsPrepareFun(
          function() {
            1 # foo
          }
        ),
        regexp = "parameter labels from comments will be replaced by 'label()'",
        fixed = TRUE
      ),
      function() {
        1
        label("foo")
      },
      # Env and srcref attributes will not be equal
      check.attributes = FALSE,
      info = "comment lines are converted to labels"
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
      check.attributes = FALSE,
      info = "comment lines without other information are dropped"
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
      check.attributes = FALSE,
      info = "comment lines with other information are converted to label()"
    )
    expect_equal(
      nlmixr:::nlmixrBoundsPrepareFunComments(nlmixrTestFunToChar(
        function() {
          1 | STUDY # hello
        }
      )),
      function() {
        1 | STUDY
        label("hello")
      },
      # Env and srcref attributes will not be equal
      check.attributes = FALSE,
      info = "comment lines with other information are converted to label() (even if they are on a line with a condition)"
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
      check.attributes = FALSE,
      info = "This is challenging to parse, and it was formerly a bug.  It is the reason that we are moving to parsing and not string extraction."
    )
  })

  # Test call and name replacement ####
  test_that("call replacement", {
    expect_equal(
      nlmixr:::replaceCallName(x = a ~ b(), replacementFun = "c", sourceNames = "b"),
      a ~ c(),
      check.attributes = FALSE,
      info = "Simple replacement works"
    )
    expect_equal(
      nlmixr:::replaceCallName(x = a ~ b(c + d * b(e)), replacementFun = "c", sourceNames = "b"),
      a ~ c(c + d * c(e)),
      check.attributes = FALSE,
      info = "Nested replacement works"
    )
    expect_equal(
      nlmixr:::replaceCallName(x = a ~ b, replacementFun = "c", sourceNames = "b"),
      a ~ b,
      check.attributes = FALSE,
      info = "Names that are not calls are not replaced"
    )
    expect_equal(
      nlmixr:::replaceCallName(x = a ~ b(b), replacementFun = "c", sourceNames = "b"),
      a ~ c(b),
      check.attributes = FALSE,
      info = "Nested names that are not calls are not replaced"
    )
    expect_equal(
      nlmixr:::replaceCallName(x = a ~ b(1 + "A"), replacementFun = "c", sourceNames = "b"),
      a ~ c(1 + "A"),
      check.attributes = FALSE,
      info = "Non-name values are permitted"
    )
  })

  test_that("name replacement", {
    expect_equal(
      nlmixr:::replaceNameName(x = a ~ b, replacementName = "c", sourceNames = "b"),
      a ~ c,
      check.attributes = FALSE,
      info = "Names that are not calls are replaced"
    )
    expect_equal(
      nlmixr:::replaceNameName(x = a ~ b(), replacementName = "c", sourceNames = "b"),
      a ~ b(),
      check.attributes = FALSE,
      info = "Function calls are skipped"
    )
    expect_equal(
      nlmixr:::replaceNameName(x = a ~ b(c + d * b(e)), replacementName = "c", sourceNames = "b"),
      a ~ b(c + d * b(e)),
      check.attributes = FALSE,
      info = "Nested replacement still ignores function calls"
    )
    expect_equal(
      nlmixr:::replaceNameName(x = a ~ b(b), replacementName = "c", sourceNames = "b"),
      a ~ b(c),
      check.attributes = FALSE,
      info = "Nested names that are not calls are replaced"
    )
    expect_equal(
      nlmixr:::replaceNameName(x = a ~ b(1 + "A"), replacementName = "c", sourceNames = "b"),
      a ~ b(1 + "A"),
      check.attributes = FALSE,
      info = "Non-name values are permitted"
    )
    expect_equal(
      nlmixr:::replaceNameName(x = a ~ b(b), replacementName = NULL, sourceNames = "b"),
      a ~ b(),
      check.attributes = FALSE,
      info = "Null removes the name"
    )
    expect_equal(
      nlmixr:::replaceNameName(x = a ~ b * c, replacementName = NULL, sourceNames = "b"),
      {
        comparison <- a ~ b * c
        comparison[[3]][[2]] <- NULL
        comparison
      },
      check.attributes = FALSE,
      info = "This is weird but consistent (and no longer a valid formula), you can remove anything"
    )
  })

}, test="cran")
