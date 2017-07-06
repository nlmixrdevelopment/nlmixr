context("Test FOCEI with NONMEM")
rxPermissive({

    mypar1 = function ()
    {
        ke = theta[1] * exp(eta[1]);
    }

    mypar2 = function ()
    {
        k = theta[1] * exp(eta[1]);
        v = 1
    }

    mod <- RxODE({
        ipre = 10 * exp(-ke * t)
    })

    pred = function() ipre

    inits = list(THTA=c(0.5))
    inits$OMGA=list(ETA[1]~.04)
    inits$ERROR=list(~prop(.1))

    dat <- Wang2007
    dat$DV <- dat$Y

    fit.prop <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))

    test_that("Matches NONMEM objective proportional function; (Based on Wang2007)", {
        expect_equal(round(fit.prop$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
    })

    m1 <- RxODE({
        d / dt(ipre) = -ke * ipre
    })

    ## Enhance data frame to include dosing records.
    dat2 <- dat[dat$Time == 0, ]
    dat2$EVID <- 101
    dat2$AMT <- 10;
    dat2 <- rbind(dat2, data.frame(dat, EVID=0, AMT=0))
    dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$Time)), ]

    fit.prop2 <- focei.fit(data=dat2, inits, mypar1, model=m1,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))

    test_that("Matches NONMEM objective proportional function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.prop2$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
    })

    fit.prop3 <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))

    test_that("Matches NONMEM objective proportional function; pred GRAD (Based on Wang2007)", {
        expect_equal(round(fit.prop3$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
    })

    fit.prop4 <- focei.fit(data=dat2, inits, mypar1, model=m1,pred=pred, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))

    test_that("Matches NONMEM objective proportional function; ODE GRAD (Based on Wang2007)", {
        expect_equal(round(fit.prop4$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
    })

    fit.prop5 <- focei.fit(data=dat2, inits, mypar2, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Solved GRAD (Based on Wang2007)", {
        expect_equal(round(fit.prop5$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
    })

    fit.prop6 <- focei.fit(data=dat2, inits, mypar2, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Solved (Based on Wang2007)", {
        expect_equal(round(fit.prop6$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
    })

    inits = list(THTA=c(0.5))
    inits$OMGA=list(ETA[1]~.04)
    inits$ERROR=list(~add(.1))

    fit.add <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCE=FOCEI for NONMEM
    })

    fit.add <- focei.fit(data=dat2, inits, mypar1, model=m1,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCE=FOCEI for NONMEM
    })

    fit.add <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Grad (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCE=FOCEI for NONMEM
    })

    fit.add <- focei.fit(data=dat2, inits, mypar1, model=m1,pred=pred, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Grad ODE (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCE=FOCEI for NONMEM
    })

    fit.add <- focei.fit(data=dat2, inits, mypar2, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Solved Grad (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCE=FOCEI for NONMEM
    })

    fit.add <- focei.fit(data=dat2, inits, mypar2, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Solved (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCE=FOCEI for NONMEM
    })

    inits = list(THTA=c(0.5))
    inits$OMGA=list(ETA[1]~.04)
    inits$ERROR=list(~add(.1) + prop(0.2))

    fit.add.prop <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))

    test_that("Matches NONMEM objective proportional function; (Based on Wang2007)", {
        expect_equal(round(fit.add.prop$objective, 3), 51.882); # New entry that does not match Wang2007
    })

    fit.add.prop2 <- focei.fit(data=dat2, inits, mypar1, model=m1,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))

    test_that("Matches NONMEM objective proportional function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.add.prop2$objective, 3), 51.882); # New entry that does not match Wang2007
    })

    fit.add.prop3 <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))

    test_that("Matches NONMEM objective proportional function; Grad (Based on Wang2007)", {
        expect_equal(round(fit.add.prop$objective, 3), 51.882); # New entry that does not match Wang2007
    })

    fit.add.prop3 <- focei.fit(data=dat2, inits, mypar1, model=m1,pred=pred, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Grad ODE (Based on Wang2007)", {
        expect_equal(round(fit.add.prop3$objective, 3), 51.882); # New entry that does not match Wang2007
    })

    fit.add.prop4 <- focei.fit(data=dat2, inits, mypar2, optim="nlminb", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Grad Solved (Based on Wang2007)", {
        expect_equal(round(fit.add.prop4$objective, 3), 51.882); # New entry that does not match Wang2007
    })

    fit.add.prop5 <- focei.fit(data=dat2, inits, mypar2, optim="bobyqa", control=list(NOTRUN=TRUE, cores=1))
    test_that("Matches NONMEM objective proportional function; Grad Solved (Based on Wang2007)", {
        expect_equal(round(fit.add.prop5$objective, 3), 51.882); # New entry that does not match Wang2007
    })

}, silent=TRUE)
