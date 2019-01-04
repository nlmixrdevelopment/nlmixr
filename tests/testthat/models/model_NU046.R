source("helper-prep_fit.R")
context("NLME46: two-compartment infusion Michaelis-Menten, multiple-dose")
runno <- "N046"

datr <-
  read.csv("../Infusion_2CPT.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- datr[datr$EVID != 2,]

datIV <- datr[datr$AMT > 0,]
datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
datIV$AMT  <- -1*datIV$AMT

datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME),]

two.compartment.IV.model <- function(){
    ini({ # Where initial conditions/variables are specified
          # '<-' or '=' defines population parameters
          # Simple numeric expressions are supported
        lCl <- 1.6      #log Cl (L/hr)
        lVc <- 4.5      #log Vc (L)
        lQ  <- 1.6      #log Q (L/hr)
        lVp <- 4        #log Vp (L)

        # Bounds may be specified by c(lower, est, upper), like NONMEM:
        # Residuals errors are assumed to be population parameters
        prop.err <- c(0, 0.3, 1)
        # Between subject variability estimates are specified by '~'
        # Semicolons are optional
        eta.Vc ~ 0.15
        eta.Cl ~ 0.15
        eta.Vp ~ 0.15
        eta.Q  ~ 0.15
    })
    model({ # Where the model is specified
        # The model uses the ini-defined variable names
        Vc <- exp(lVc + eta.Vc)
        Cl <- exp(lCl + eta.Cl)
        Vp <- exp(lVp + eta.Vp)
        Q <- exp(lQ + eta.Q)
        # And is assumed to follow proportional error estimated by prop.err
        linCmt() ~ prop(prop.err)
    })
}


dat <- datr

mod <- nlmixr(two.compartment.IV.model);

opts <- c("nlme", "saem", "fo", "foi", "foce", "focei")
for (opt in opts){
    context(sprintf("%s-UI-046: two-compartment infusion, single-dose", opt))
    runno <- paste0(opt, "N046")
    fit[[runno]] <- nlmixr(mod, datr, opt, control=defaultControl(opt), table=tableControl(cwres=TRUE))
    source(genIfNeeded())
}
