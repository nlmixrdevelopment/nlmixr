nlmixrTest({

  f <- function() {
    ini({
      popKe <- 0.5
      etaKe ~ 0.04
      etaF ~ 0.1
      fp <- 2
      prop.sd <- 0.1
    })
    model({
      ke <- popKe * exp(etaKe)
      d / dt(ipre) <- -ke * ipre
      f(ipre) <- fp * exp(etaF)
      ipre ~ prop(prop.sd)
    })
  }

  set.seed(42)
  dat <- Wang2007
  dat$DV <- dat$Y
  dat2 <- dat[dat$Time == 0, ]
  dat2$EVID <- 101
  dat2$AMT <- 10
  dat2 <- rbind(dat2, data.frame(dat, EVID = 0, AMT = 0))
  dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$Time)), ]
  dat2 <- dat2[, names(dat2) != "Y"]
  dat3 <- data.frame(ID = 1:10, f0 = 2 * exp(rnorm(10, sd = 0.1)))
  dat2 <- merge(dat2, dat3)
  dat2$DV <- dat2$DV * dat2$f0

  testIt <- function(meth) {
    context(sprintf("test finite difference sensitivities in %s", meth))
    test_that(sprintf("finite difference %s", meth), {
      suppressWarnings({
        fit <- nlmixr(f, dat2, "focei",
          control = foceiControl(maxOuterIterations = 0, covMethod = "")
        )
        expect_false(all(fit$eta$etaF == 0))

        fit <- nlmixr(f, dat2, "focei",
          control = foceiControl(maxOuterIterations = 0, covMethod = "", eventCentral = FALSE)
        )
        expect_false(all(fit$eta$etaF == 0))
      })
    })
  }

  lapply(c("focei", "foce", "foi", "fo"), testIt)
}, test="focei")
