source("helper-prep_fit.R")
context("NLME02: one-compartment bolus, multiple-dose")
datr <- Bolus_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2, ]
dat <- datr[datr$SD == 0, ]

one.compartment.IV.model <- function(){
    ini({ # Where initial conditions/variables are specified
          # '<-' or '=' defines population parameters
          # Simple numeric expressions are supported
        lCl <- 1.6      #log Cl (L/hr)
        lVc <- 4.5      #log V (L)
        # Bounds may be specified by c(lower, est, upper), like NONMEM:
        # Residuals errors are assumed to be population parameters
        prop.err <- c(0, 0.3, 1)
        # Between subject variability estimates are specified by '~'
        # Semicolons are optional
        eta.Vc ~ 0.1   #IIV V
        eta.Cl ~ 0.1   #IIV Cl
    })
    model({ # Where the model is specified
        # The model uses the ini-defined variable names
        Vc <- exp(lVc + eta.Vc)
        Cl <- exp(lCl + eta.Cl)
        # RxODE-style differential equations are supported
        d / dt(centr) = -(Cl / Vc) * centr;
        ## Concentration is calculated
        cp = centr / Vc;
        # And is assumed to follow proportional error estimated by prop.err
        cp ~ prop(prop.err)
    })
}


mod <- nlmixr(one.compartment.IV.model);

opts <- c("nlme", "saem", "fo", "foi", "foce", "focei")
for (opt in opts){
    context(sprintf("%s-UI-002: one-compartment bolus, single-dose", opt))
    runno <- paste0(opt, "U002ode")
    fit[[runno]] <- nlmixr(mod, dat, opt, control=defaultControl(opt), table=tableControl(cwres=TRUE))
    source(genIfNeeded())
}
