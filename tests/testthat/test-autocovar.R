context("autocovar")

# Test cases
# 1. addCovariate -  categorical, HS, non-categorical, log/not-log -- done
# 2. removecovariate - categorical, HS, non-categorical, log/not-log -- done
# 3. addCovVar --> performNorm() -- done
# 4. performNorm - norm-median/mean, norm type - mul, div, add, sub, initial Est  - done
# 5. addCovMultiple - function not complete (yet!)

# === addCovariate
testthat::test_that("adding non-categorical covariate to funstring with log-transformation", {
  funstring <- "ka <- exp(tka + eta.ka)"
  varName <- "ka"
  covariate <- "cov_WT_ka*WT"
  theta <- list("tka")
  isLog <- TRUE

  funstring1 <- addCovariate(funstring, varName, covariate, theta, isLog)
  funstring2 <- "ka<-exp(tka+cov_WT_ka*WT+eta.ka)"

  testthat::expect_equal(funstring1, funstring2)
})

testthat::test_that("adding non-categorical covariate to funstring without log-transformation", {
  funstring <- "ka <- tka + eta.ka"
  varName <- "ka"
  covariate <- "cov_WT_ka*WT"
  theta <- list("tka")
  isLog <- FALSE

  funstring1 <- addCovariate(funstring, varName, covariate, theta, isLog)
  funstring2 <- "ka<-(tka+eta.ka)*(cov_WT_ka*WT)"

  testthat::expect_equal(funstring1, funstring2)
})


testthat::test_that("adding categorical covariate to funstring with log-transformation", {
  funstring <- "ka <- exp(tka + eta.ka)"
  varName <- "ka"
  covariate <- "cov_factor_1_ka*factor_1+cov_factor_2_ka*factor_2+cov_factor_3_ka*factor_3"
  theta <- list("tka")
  isLog <- TRUE

  funstring1 <- addCovariate(funstring, varName, covariate, theta, isLog)
  funstring2 <- "ka<-exp(tka+cov_factor_1_ka*factor_1+cov_factor_2_ka*factor_2+cov_factor_3_ka*factor_3+eta.ka)"

  testthat::expect_equal(funstring1, funstring2)
})

testthat::test_that("adding categorical covariate to funstring without log-transformation", {
  funstring <- "ka <- tka + eta.ka"
  varName <- "ka"
  covariate <- "cov_factor_1_ka*factor_1+cov_factor_2_ka*factor_2+cov_factor_3_ka*factor_3"
  theta <- list("tka")
  isLog <- FALSE

  funstring1 <- addCovariate(funstring, varName, covariate, theta, isLog)
  funstring2 <- "ka<-(tka+eta.ka)*(cov_factor_1_ka*factor_1+cov_factor_2_ka*factor_2+cov_factor_3_ka*factor_3)"

  testthat::expect_equal(funstring1, funstring2)
})


# ==== removeCovariate

testthat::test_that("removing non-categorical covariate from funstring with log-transformation", {
  funstring <- "ka<-exp(tka+cov_WT_ka*WT+eta.ka)"
  varName <- "ka"
  covariate <- "cov_WT_ka*WT"
  theta <- list("tka")

  funstring1 <- removeCovariate(funstring, varName, covariate, theta)
  funstring2 <- "ka<-exp(tka+eta.ka)"

  testthat::expect_equal(funstring1, funstring2)
})

testthat::test_that("removing non-categorical covariate from funstring without log-transformation", {
  funstring <- "ka<-(tka+eta.ka)*(cov_WT_ka*WT)"
  varName <- "ka"
  covariate <- "cov_WT_ka*WT"
  theta <- list("tka")

  funstring1 <- removeCovariate(funstring, varName, covariate, theta)
  funstring2 <- "ka<-(tka+eta.ka)"

  testthat::expect_equal(funstring1, funstring2)
})


testthat::test_that("removing categorical covariate from funstring with log-transformation", {
  funstring <- "ka<-exp(tka+cov_factor_1_ka*factor_1+cov_factor_2_ka*factor_2+cov_factor_3_ka*factor_3+eta.ka)"
  varName <- "ka"
  covariate <- "cov_factor_1_ka*factor_1 + cov_factor_2_ka*factor_2 + cov_factor_3_ka*factor_3"
  theta <- list("tka")

  funstring1 <- removeCovariate(funstring, varName, covariate, theta)
  funstring2 <- "ka<-exp(tka+eta.ka)"

  testthat::expect_equal(funstring1, funstring2)
})

testthat::test_that("removing categorical covariate from funstring without log-transformation", {
  funstring <- "ka<-(tka+eta.ka)*(cov_factor_1_ka*factor_1+cov_factor_2_ka*factor_2+cov_factor_3_ka*factor_3)"
  varName <- "ka"
  covariate <- "cov_factor_1_ka*factor_1 + cov_factor_2_ka*factor_2 + cov_factor_3_ka*factor_3"
  theta <- list("tka")

  funstring1 <- removeCovariate(funstring, varName, covariate, theta)
  funstring2 <- "ka<-(tka+eta.ka)"

  testthat::expect_equal(funstring1, funstring2)
})

# ==== performNorm

testthat::test_that("normalize non-categorical covariates using mean and norm type mul", {
  data <- theo_sd
  covariate <- "WT"
  normOp <- `*`
  normValVec <- mean(data[, "WT"])
  varName <- "ka"
  isLog <- FALSE
  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- theo_sd[, "WT"] * normValVec
  covNameMod2 <- "centered_WT*cov_WT_ka"
  covNames2 <- "cov_WT_ka"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})

testthat::test_that("normalize non-categorical covariates using mean and norm type div", {
  data <- theo_sd
  covariate <- "WT"
  normOp <- `/`
  normValVec <- mean(data[, "WT"])
  varName <- "ka"
  isLog <- FALSE
  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- theo_sd[, "WT"] / normValVec
  covNameMod2 <- "centered_WT*cov_WT_ka"
  covNames2 <- "cov_WT_ka"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})

testthat::test_that("normalize non-categorical covariates using mean and norm type sub", {
  data <- theo_sd
  covariate <- "WT"
  normOp <- `-`
  varName <- "ka"
  normValVec <- mean(data[, "WT"])
  isLog <- FALSE
  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- theo_sd[, "WT"] - normValVec
  covNameMod2 <- "centered_WT*cov_WT_ka"
  covNames2 <- "cov_WT_ka"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})

testthat::test_that("normalize non-categorical covariates using mean and norm type add", {
  data <- theo_sd
  covariate <- "WT"
  varName <- "ka"
  normOp <- `+`
  normValVec <- mean(data[, "WT"])
  isLog <- FALSE
  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- theo_sd[, "WT"] + normValVec
  covNameMod2 <- "centered_WT*cov_WT_ka"
  covNames2 <- "cov_WT_ka"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})


testthat::test_that("normalize non-categorical covariates with log-transformation using mean and norm type mul", {
  data <- theo_sd
  covariate <- "WT"
  varName <- "ka"
  normOp <- `*`
  normValVec <- mean(data[, "WT"])
  isLog <- TRUE
  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- log(theo_sd[, "WT"] * normValVec)
  covNameMod2 <- "centered_WT*cov_WT_ka"
  covNames2 <- "cov_WT_ka"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})

testthat::test_that("normalize non-categorical covariates with log-transformation using mean and norm type div", {
  data <- theo_sd
  covariate <- "WT"
  varName <- "ka"
  normOp <- `/`
  normValVec <- mean(data[, "WT"])
  isLog <- TRUE

  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- log(theo_sd[, "WT"] / normValVec)
  covNameMod2 <- "centered_WT*cov_WT_ka"
  covNames2 <- "cov_WT_ka"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})


testthat::test_that("normalize non-categorical covariates with log-transformation using mean and norm type add", {
  data <- theo_sd
  covariate <- "WT"
  varName <- "ka"
  normOp <- `+`
  normValVec <- mean(data[, "WT"])
  isLog <- TRUE

  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- log(theo_sd[, "WT"] + normValVec)
  covNameMod2 <- "centered_WT*cov_WT_ka"
  covNames2 <- "cov_WT_ka"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})

testthat::test_that("normalize (with prefactor for Cl) for non-categorical covariates with log-transformation using mean and norm type add", {
  data <- theo_sd
  covariate <- "WT"
  varName <- "cl"
  normOp <- `+`
  normValVec <- mean(data[, "WT"])
  isLog <- TRUE

  isCat <- FALSE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1 <- res1[[1]][, "centered_WT"]
  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2 <- 0.75 * log(theo_sd[, "WT"] + normValVec)
  covNameMod2 <- "centered_WT*cov_WT_cl"
  covNames2 <- "cov_WT_cl"

  testthat::expect_equal(dat1, dat2)
  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})

testthat::test_that("categorical covariates without normalization", {
  data <- theo_sd
  data[, "factor"] <- rep(c(1, 2, 3, 4), nrow(data) / 4)
  covariate <- "factor"
  varName <- "ka"
  normOp <- NULL
  normValVec <- NULL
  isLog <- FALSE
  isCat <- TRUE
  isHS <- FALSE

  res1 <- performNorm(data, covariate, varName, normOp, normValVec, isLog, isCat, isHS)

  dat1a <- res1[[1]][, "categorical_factor_1"]
  dat1b <- res1[[1]][, "categorical_factor_2"]
  dat1c <- res1[[1]][, "categorical_factor_3"]

  covNameMod1 <- res1[[2]]
  covNames1 <- res1[[3]]

  dat2a <- 1L * (data[, "factor"] == 1)
  dat2b <- 1L * (data[, "factor"] == 2)
  dat2c <- 1L * (data[, "factor"] == 3)

  covNameMod2 <- "categorical_factor_1*cov_factor_1_ka+categorical_factor_2*cov_factor_2_ka+categorical_factor_3*cov_factor_3_ka"
  covNames2 <- c("cov_factor_1_ka", "cov_factor_2_ka", "cov_factor_3_ka")

  testthat::expect_equal(dat1a, dat2a)
  testthat::expect_equal(dat1b, dat2b)
  testthat::expect_equal(dat1c, dat2c)

  testthat::expect_equal(covNameMod1, covNameMod2)
  testthat::expect_equal(covNames1, covNames2)
})

# ==== makeHockeyStick
testthat::test_that("creating hockey stick variables", {
  data <- theo_sd
  covariate <- "WT"
  varName = 'ka'
  med <- median(data[, covariate])

  res1 <- makeHockeyStick(data, covariate, varName)
  dat1a <- res1[[1]][, "centered_WT_lower"]
  dat1b <- res1[[1]][, "centered_WT_upper"]

  covModExpr1 <- res1[[2]]
  covNames1 <- res1[[3]]
  colNames <- res1[[4]]


  dat2a <- 1L * (data[, covariate] < med) * data[, covariate]
  dat2b <- 1L * (data[, covariate] >= med) * data[, covariate]

  covModExpr2 <- c("centered_WT_lower*cov_WT_lower_ka", "centered_WT_upper*cov_WT_upper_ka")
  covNames2 <- c("cov_WT_lower_ka", "cov_WT_upper_ka")

  testthat::expect_equal(dat1a, dat2a)
  testthat::expect_equal(dat1b, dat2b)
  testthat::expect_equal(covModExpr1, covModExpr2)
  testthat::expect_equal(covNames1, covNames2)
})
