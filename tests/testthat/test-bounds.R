context("Test bounds extraction")

ref <- structure(list(theta = c(1, 2, 3, 4, 5, 6, 7, 8, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
9, 10, 11, 12, 13, 14), eta1 = c(NA, NA, NA, NA, NA, NA, NA,
NA, 1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12,
12, NA, NA, NA, NA, NA, NA), eta2 = c(NA, NA, NA, NA, NA, NA,
NA, NA, 1, 2, 3, 4, 5, 5, 6, 7, 7, 8, 7, 8, 9, 10, 10, 11, 10,
11, 12, NA, NA, NA, NA, NA, NA), name = structure(c(1L, 2L, 3L,
4L, NA, NA, NA, NA, 5L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, 6L, 7L, 8L, 9L, 10L, 11L, NA, NA, NA, 12L, 13L, 14L), .Label = c("a",
"b", "c", "d", "et1", "et2", "(et3,et2)", "et3", "(et4,et2)",
"(et4,et3)", "et4", "a5", "a6", "a7"), class = "factor"), lower = c(0,
0, -Inf, -Inf, 0, 0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,
-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,
-Inf, -Inf, -Inf, -Inf, 9, 11, -Inf, 9, 11), est = c(1, 3, 4,
4, 1, 1, 1, 1, 10, 20, 30, 40, 40, 0.1, 20, 40, 0.1, 20, 0.1,
0.1, 30, 40, 0.1, 20, 0.1, 0.1, 30, 8, 10, 12, 8, 10, 12), upper = c(2,
Inf, Inf, Inf, 2, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf,
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf,
Inf, Inf, 13, Inf, Inf, 13), fix = c(FALSE, FALSE, FALSE, FALSE,
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
TRUE), label = c("A", NA, NA, NA, "e", NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, "labels", NA, NA, NA)), .Names = c("theta", "eta1", "eta2",
"name", "lower", "est", "upper", "fix", "label"), row.names = c(NA,
-33L), class = "data.frame");


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
        0.1, 20);
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
