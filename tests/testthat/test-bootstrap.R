nlmixrTest(
  {
    context("bootstrap")

    samp_dat <- theo_sd

    test_that("sampling should return different datasets at each call", {
      a <- digest::digest(nlmixr:::sampling(samp_dat))
      b <- digest::digest(nlmixr:::sampling(samp_dat))
      testthat::expect_false(isTRUE(all.equal(a, b)))
    })


    test_that("resuming the fit should not return the same datasets as before", {
      one.cmt <- function() {
        ini({
          ## You may label each parameter with a comment
          tka <- 0.45 # Log Ka
          tcl <- 1 # Log Cl
          ## This works with interactive models
          ## You may also label the preceding line with label("label text")
          tv <- 3.45
          label("log V")
          ## the label("Label name") works with all models
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd)
        })
      }

      fit <- suppressWarnings(nlmixr(
        one.cmt,
        samp_dat,
        est = "focei",
        control = list(print = 0),
        table = list(npde = TRUE, cwres = TRUE)
      ))


      fit1 <- nlmixr:::bootstrapFit(fit, nboot = 2, restart = TRUE)
      fit2 <- nlmixr:::bootstrapFit(fit, nboot = 4, restart = FALSE)

      output_dir <-
        paste0("nlmixrBootstrapCache_", "fit", "_", fit$bootstrapMd5)

      fnameBootDataPattern <- paste0("boot_data",
        "_", "[0-9]+", ".RData",
        sep = ""
      )
      fileExists <- list.files(paste0("./", output_dir), pattern = fnameBootDataPattern)

      ## print(output_dir)

      fitdata <- lapply(fileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })

      a <- digest::digest(fitdata[[1]])
      b <- digest::digest(fitdata[[3]])
      expect_false(isTRUE(all.equal(a, b)))

      a <- digest::digest(fitdata[[2]])
      b <- digest::digest(fitdata[[4]])
      expect_false(isTRUE(all.equal(a, b)))

      unlink(output_dir, recursive = TRUE, force = TRUE)
    })

    testthat::test_that("different confidence levels should result in different bands", {
      one.cmt <- function() {
        ini({
          ## You may label each parameter with a comment
          tka <- 0.45 # Log Ka
          tcl <- 1 # Log Cl
          ## This works with interactive models
          ## You may also label the preceding line with label("label text")
          tv <- 3.45
          label("log V")
          ## the label("Label name") works with all models
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd)
        })
      }

      fit <- suppressWarnings(nlmixr(
        one.cmt,
        samp_dat,
        est = "focei",
        control = list(print = 0),
        table = list(npde = TRUE, cwres = TRUE)
      ))

      fitlist <- nlmixr:::modelBootstrap(fit, nboot = 4, restart = TRUE)
      bootSummary1 <- nlmixr:::getBootstrapSummary(fitlist, ci = 0.95)
      bootSummary2 <- nlmixr:::getBootstrapSummary(fitlist, ci = 0.75)

      a <- digest::digest(bootSummary1$parFixedDf$confLower)
      b <- digest::digest(bootSummary2$parFixedDf$confLower)
      testthat::expect_false(isTRUE(all.equal(a, b)))

      a <- digest::digest(bootSummary1$parFixedDf$confUpper)
      b <- digest::digest(bootSummary2$parFixedDf$confUpper)
      testthat::expect_false(isTRUE(all.equal(a, b)))
    })

    testthat::test_that("expected columns in fit$parFixedDf object should match", {
      one.cmt <- function() {
        ini({
          ## You may label each parameter with a comment
          tka <- 0.45 # Log Ka
          tcl <- 1 # Log Cl
          ## This works with interactive models
          ## You may also label the preceding line with label("label text")
          tv <- 3.45
          label("log V")
          ## the label("Label name") works with all models
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd)
        })
      }

      fit <- suppressWarnings(nlmixr(
        one.cmt,
        samp_dat,
        est = "focei",
        control = list(print = 0),
        table = list(npde = TRUE, cwres = TRUE)
      ))

      colsBefore <- colnames(fit$parFixedDf)
      fitlist <- nlmixr:::modelBootstrap(fit, nboot = 4, restart = TRUE)

      bootSummary <- nlmixr:::getBootstrapSummary(fitlist, ci = 0.95)

      colsAfter <- colnames(fit$parFixedDf)

      testthat::expect_equal(colsAfter, colsBefore)
    })
  },
  silent = TRUE,
  test = "bootstrap"
)
