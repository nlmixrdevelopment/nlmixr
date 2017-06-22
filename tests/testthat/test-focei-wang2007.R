context("Test FOCEI with NONMEM")
rxPermissive({

    mypar1 = function ()
    {
        ke = theta[1] * exp(eta[1]);
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

    fit.prop <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE))

    inits = list(THTA=c(0.5))
    inits$OMGA=list(ETA[1]~.04)
    inits$ERROR=list(~add(.1))

    fit.add <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE))

    inits = list(THTA=c(0.5))
    inits$OMGA=list(ETA[1]~.04)
    inits$ERROR=list(~add(.1) + prop(0.2))

    fit.add.prop <- focei.fit(data=dat, inits, mypar1, model=mod,pred=pred, optim="bobyqa", control=list(NOTRUN=TRUE))

    test_that("Matches NONMEM objective function; (Based on Wang2007)", {
        expect_equal(round(fit.prop$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCE=FOCEI for NONMEM
        expect_equal(round(fit.add.prop$objective, 3), 51.882); # New entry that does not match Wang2007
    })

}, silent=TRUE)
