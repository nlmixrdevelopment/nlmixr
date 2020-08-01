nlmixrTest({
  context("precondition tests")

  one.compartment <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
  }

  fit2 <- nlmixr(one.compartment, theo_sd, est="focei",control=list(print=0))

  df1 <- fit2$parFixedDf
  cov1 <- fit2$cov

  ## Simply re-evaluate with no estimation (including inner estimation)
  preconditionFit(fit2, estType = "none")

  df2 <- fit2$parFixedDf
  cov2 <- fit2$cov

  ## In this case there isn't a theta/omega estimate so these should be the same
  expect_equal(df1$Estimate, df2$Estimate)
  expect_equal(df1$`Back-transformed`, df2$`Back-transformed`)
  expect_equal(df1$`BSV(CV%)`, df2$`BSV(CV%)`)
  expect_equal(df1$`Shrink(SD)%`, df2$`Shrink(SD)%`)

  expect_false(isTRUE(all.equal(df1$SE, df2$SE)))
  expect_false(isTRUE(all.equal(df1$`%RSE`, df2$`%RSE`)))
  expect_false(isTRUE(all.equal(df1$`CI Lower`, df2$`CI Lower`)))
  expect_false(isTRUE(all.equal(df1$`%RSE`, df2$`%RSE`)))
  expect_false(isTRUE(all.equal(cov1, cov2)))

}, test="precondition")
