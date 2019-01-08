source("helper-prep_fit.R")

datr <-
  read.csv("../Infusion_2CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- datr[datr$EVID != 2,]
datIV <- datr[datr$AMT > 0,]
datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
datIV$AMT <- -1*datIV$AMT
datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME),]

dat <- datr[datr$SD == 1,]



two.compartment.IV.MM.model2 <- function(){
    ini({ # Where initial conditions/variables are specified
          # '<-' or '=' defines population parameters
          # Simple numeric expressions are supported
        lVM <- 7      #log Vmax (mg/hr)
        lKM <- 5.7    #log KM (mg/L)
        lVc <- 4.5    #log Vc (L)
        lQ  <- 1.5    #log Q (L/hr)
        lVp <- 4      #log Vp (L)
        # Bounds may be specified by c(lower, est, upper), like NONMEM:
        # Residuals errors are assumed to be population parameters
        prop.err <- c(0, 0.3, 1)
        # Between subject variability estimates are specified by '~'
        # Semicolons are optional
        eta.VM ~ 0.15
        eta.KM ~ 0.15
        eta.Vc ~ 0.15
        eta.Q  ~ 0.15
        eta.Vp ~ 0.15
    })
    model({ # Where the model is specified
        # The model uses the ini-defined variable names
        VM <- exp(lVM + eta.VM)
        KM <- exp(lKM + eta.KM)
        Vc <- exp(lVc + eta.Vc)
        Q  <- exp(lQ  + eta.Q)
        Vp <- exp(lVp + eta.Vp)
        # RxODE-style differential equations are supported
        K12<- Q/Vc
        K21<- Q/Vp
        d/dt(centr)  = K21*periph-K12*centr-(VM*centr/Vc)/(KM+centr/Vc);
        d/dt(periph) =-K21*periph+K12*centr;
        ## Concentration is calculated
        cp = centr / Vc;
        # And is assumed to follow proportional error estimated by prop.err
        cp ~ prop(prop.err)
    })
}

mod <- nlmixr(two.compartment.IV.MM.model2);

opts <- c("nlme", "saem", "fo", "foi", "foce", "focei")
for (opt in opts){
    context(sprintf("%s-UI-054: two-compartment infusion Michaelis-Menten, single-dose", opt))
    runno <- paste0(opt, "N054_ode")
    fit[[runno]] <- nlmixr(mod, dat, opt, control=defaultControl(opt), table=tableControl(cwres=TRUE))
    source(genIfNeeded())
}
