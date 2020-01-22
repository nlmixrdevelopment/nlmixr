rxPermissive({
    context("M4 -- all observations")

    .nlmixr <- function(...){suppressWarnings(nlmixr(...))}
    
    dat <- Wang2007
    dat$DV <- dat$Y # Add the required DV data item
    
    f <- function(){
        ini({
            tvK <- 0.5 # Typical Value of K
            bsvK ~ 0.04 # Between Subject Variance of K
            prop.sd <- sqrt(0.1)
        })
        model({
            ke <- tvK * exp(bsvK);
            v <- 1
            ipre <- 10*exp(-ke*t)
            ipre ~ prop(prop.sd)
        })
    }

    dat2 <- dat;
    dat2$limit <- 0

    dat3 <- dat;
    dat3$limit <- 3

    f.focei <- .nlmixr(f, dat, "posthoc")
    
    f.focei2 <- .nlmixr(f, dat2, "posthoc")

    expect_false(isTRUE(all.equal(f.focei$objf, f.focei2$objf)))

    ## Limit affects values

    f.focei3 <- .nlmixr(f, dat3, "posthoc")

    expect_false(isTRUE(all.equal(f.focei2$objf, f.focei3$objf)))
    
    f.foce <- .nlmixr(f, dat, "posthoc", control=list(interaction=FALSE))
    
    f.foce2 <- .nlmixr(f, dat2, "posthoc", control=list(interaction=FALSE))
    
    expect_false(isTRUE(all.equal(f.foce$objf, f.foce2$objf)))

    f.foce3 <- .nlmixr(f, dat3, "posthoc", control=list(interaction=FALSE))

    expect_false(isTRUE(all.equal(f.foce2$objf, f.foce3$objf)))
    
}, on.validate="NLMIXR_VALIDATION", silent=TRUE)
