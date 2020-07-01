context('bootstrap')
samp_dat = theo_sd

testthat::test_that("sampling should return different datasets at each call", {
  a <- digest::digest(nlmixr:::sampling(samp_dat))
  b <- digest::digest(nlmixr:::sampling(samp_dat))
  testthat::expect_false(isTRUE(all.equal(a, b)))
})


testthat::test_that("resuming the fit should return not return the same datasets as before",
                    {
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
                        est = "saem",
                        control = list(print = 0),
                        table = list(npde = TRUE, cwres = TRUE)
                      ))
                      
                      fit1 <- nlmixr:::bootstrapFit(fit, numModels = 2, resume = FALSE)
                      fit2 <- nlmixr:::bootstrapFit(fit, numModels = 4, resume = TRUE)
                      
                      fnamebootdata <- paste0(getwd(),
                                              "/",
                                              "nlmixrBootstrapCache_fit/",
                                              "boot_data.Rdata")
                      fitdata <- readRDS(fnamebootdata)
                      
                      a <- digest::digest(fitdata[[1]])
                      b <- digest::digest(fitdata[[3]])
                      testthat::expect_false(isTRUE(all.equal(a, b)))
                      
                      a <- digest::digest(fitdata[[2]])
                      b <- digest::digest(fitdata[[4]])
                      testthat::expect_false(isTRUE(all.equal(a, b)))
                      
                      unlink(paste0(getwd(),
                             "/",
                             "nlmixrBootstrapCache_fit/"), recursive = TRUE, force=TRUE)
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
  
  fitlist <- modelBootstrap(fit, numModels = 4, resume = FALSE)
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
  
  modelBootstrap(fit, numModels = 4, resume = FALSE)
  
  
  bootSummary1 <- nlmixr:::getBootstrapSummary(fitlist, ci = 0.95)
  bootSummary2 <- nlmixr:::getBootstrapSummary(fitlist, ci = 0.75)
  
  a <- digest::digest(bootSummary1$parFixedDf$confLower)
  b <- digest::digest(bootSummary2$parFixedDf$confLower)
  testthat::expect_false(isTRUE(all.equal(a, b)))
  
  a <- digest::digest(bootSummary1$parFixedDf$confUpper)
  b <- digest::digest(bootSummary2$parFixedDf$confUpper)
  testthat::expect_false(isTRUE(all.equal(a, b)))
})
