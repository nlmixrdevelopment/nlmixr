nlmixrTest({
  context("M2 -- all observations")

  .nlmixr <- function(...) {
    suppressWarnings(nlmixr(...))
  }

  dat <- Wang2007
  dat$DV <- dat$Y # Add the required DV data item

  f <- function() {
    ini({
      tvK <- 0.5 # Typical Value of K
      bsvK ~ 0.04 # Between Subject Variance of K
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tvK * exp(bsvK)
      v <- 1
      ipre <- 10 * exp(-ke * t)
      ipre ~ prop(prop.sd)
    })
  }

  dat2 <- dat
  dat2$limit <- 0

  dat3 <- dat
  dat3$limit <- 3

  dat4 <- dat
  dat4$limit <- 12

  f.focei <- .nlmixr(f, dat, "posthoc")

  f.focei2 <- .nlmixr(f, dat2, "posthoc")

  expect_false(isTRUE(all.equal(f.focei$objf, f.focei2$objf)))

  ## Limit affects values

  f.focei3 <- .nlmixr(f, dat3, "posthoc")

  expect_false(isTRUE(all.equal(f.focei2$objf, f.focei3$objf)))


  f.focei4 <- .nlmixr(f, dat4, "posthoc")

  f.foce <- .nlmixr(f, dat, "posthoc", control = list(interaction = FALSE))

  f.foce2 <- .nlmixr(f, dat2, "posthoc", control = list(interaction = FALSE))

  expect_false(isTRUE(all.equal(f.foce$objf, f.foce2$objf)))

  f.foce3 <- .nlmixr(f, dat3, "posthoc", control = list(interaction = FALSE))

  expect_false(isTRUE(all.equal(f.foce2$objf, f.foce3$objf)))

  context("M3/M4 -- Missing, assume LLOQ=3 at t=1.5")

  datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
  datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
  datL <- datL[order(datL$ID, datL$Time), ]

  datL4 <- datL
  datL4$limit <- 0

  f.foceiL <- .nlmixr(f, datL, "posthoc")
  expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL$objf)))

  f.foceiL4 <- .nlmixr(f, datL4, "posthoc")
  expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL4$objf)))
  expect_false(isTRUE(all.equal(f.foceiL$objf, f.foceiL4$objf)))

  datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
  datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
  datL <- datL[order(datL$ID, datL$Time), ]

  datL4 <- datL
  datL4$limit <- 0

  f.foceiL <- .nlmixr(f, datL, "posthoc")
  expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL$objf)))

  f.foceiL4 <- .nlmixr(f, datL4, "posthoc")
  expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL4$objf)))
  expect_false(isTRUE(all.equal(f.foceiL$objf, f.foceiL4$objf)))

  ## foce

  datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
  datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
  datL <- datL[order(datL$ID, datL$Time), ]

  datL4 <- datL
  datL4$limit <- 0

  f.foceL <- .nlmixr(f, datL, "posthoc", control = list(interaction = FALSE))
  expect_false(isTRUE(all.equal(f.foce$objf, f.foceL$objf)))

  f.foceL4 <- .nlmixr(f, datL4, "posthoc", control = list(interaction = FALSE))
  expect_false(isTRUE(all.equal(f.foce$objf, f.foceL4$objf)))
  expect_false(isTRUE(all.equal(f.foceL$objf, f.foceL4$objf)))
},
test="focei")
