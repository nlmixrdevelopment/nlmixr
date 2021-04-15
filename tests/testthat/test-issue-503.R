nlmixrTest({

  PKdata <- warfarin[warfarin$dvid == "cp", ]

  One.comp.KA.solved <- function() {
    ini({
      # Where initial conditions/variables are specified
      lka  <- log(1.15)  #log ka (1/h)
      lcl  <- log(0.135) #log Cl (L/h)
      lv   <- log(8)     #log V (L)
      prop.err <- 0.15   #proportional error (SD/mean)
      add.err  <- 0.6    #additive error (mg/L)
      eta.ka ~ 0.5   #IIV ka
      eta.cl ~ 0.1   #IIV cl
      eta.v  ~ 0.1   #IIV v
    })
    model({
      # Where the model is specified
      cl <- exp(lcl + eta.cl)
      v  <- exp(lv + eta.v)
      ka <- exp(lka + eta.ka)
      ## solved system example
      ## where residual error is assumed to follow proportional and additive error
      linCmt() ~ propT(prop.err) + add(add.err)
    })
  }


  fitOne.comp.KA.solved_S2 <-
    nlmixr(
      One.comp.KA.solved,    #the model definition
      PKdata,                #the data set
      est = "saem",          #the estimation algorithm (SAEM)
      saemControl(nBurn = 200, #200 SAEM burn-in iterations (the default)
                  nEm   = 300, #300 EM iterations (the default)
                  print = 50,
                  #type="newuoa",
                  addProp="combined1"),
      tableControl(npde=TRUE, cwres=TRUE)
    )

  expect_true(fitOne.comp.KA.solved_S2$theta["add.err"] > 0.4)
},
test = "lvl2"
)
