nlmixrTest({
  context("testing saem without table can add focei objf")

  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45; label("log V")
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

  fit <- nlmixr(one.cmt, theo_sd, est="saem", control=list(calcTables=FALSE, print=0))

  expect_true(inherits(fit, "nlmixrFitCore"))
  expect_false(inherits(fit, "data.frame"))
  expect_false(inherits(fit, "nlmixrFitData"))
  expect_error(setOfv(fit, "focei"), NA)
}, test="saem")
