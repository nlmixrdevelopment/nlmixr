rxPermissive({
    context("Test FO, FOCE, and FOCEi objective functions")

    ## For some reason the ODE and solved FOCE proportional models
    ## give quite different results.  However, FOCE doesn't work as
    ## well with prop models.

    .foceiFit <- function(...) {
        suppressWarnings(foceiFit(...));
    }

    mypar1 = function () {
        ke = theta[1] * exp(eta[1]);
    }

    mypar2 = function () {
        k = theta[1] * exp(eta[1]);
        v = 1
    }

    mod <- RxODE({
        ipre = 10 * exp(-ke * t)
    })

    out.focei.prop <- structure(list(IPRE = c(10, 5.8445, 10, 6.048, 10, 5.9892, 10, 5.9676, 10, 6.0798, 10, 6.0532, 10, 6.101, 10, 6.0474, 10, 6.0748, 10, 6.0673), ID = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10), NPRED = c(10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653), NRES = c(0.68, -2.3816, 0.402, 0.38869, -0.1186, -0.20881, -0.6592, -0.44441, 0.082, 0.69299, -0.1062, 0.43959, -0.1092, 0.89039, 0.234, 0.38349, -0.0118, 0.64589, -0.3264, 0.57489), NWRES = c(0.21503, -1.1839, 0.12712, 0.19322, -0.037505, -0.1038, -0.20846, -0.22092, 0.025931, 0.34449, -0.033583, 0.21853, -0.034532, 0.44262, 0.073997, 0.19064, -0.0037315, 0.32108, -0.10322, 0.28578), PREDI = c(10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653), RESI = c(0.68, -2.3816, 0.402, 0.38869, -0.1186, -0.20881, -0.6592, -0.44441, 0.082, 0.69299, -0.1062, 0.43959, -0.1092, 0.89039, 0.234, 0.38349, -0.0118, 0.64589, -0.3264, 0.57489), WRESI = c(0.21503, -1.1786, 0.12712, 0.19235, -0.037505, -0.10333, -0.20846, -0.21992, 0.025931, 0.34294, -0.033583, 0.21754, -0.034532, 0.44062, 0.073997, 0.18978, -0.0037315, 0.31963, -0.10322, 0.28449), CPRED = c(10, 6.0691, 10, 6.0653, 10, 6.0658, 10, 6.0661, 10, 6.0653, 10, 6.0653, 10, 6.0654, 10, 6.0653, 10, 6.0653, 10, 6.0653), CRES = c(0.68, -2.3854, 0.402, 0.38867, -0.1186, -0.20927, -0.6592, -0.44517, 0.082, 0.69298, -0.1062, 0.43958, -0.1092, 0.89029, 0.234, 0.38347, -0.0118, 0.64589, -0.3264, 0.57489), CWRES = c(0.21503, -1.2221, 0.12712, 0.19366, -0.037505, -0.10511, -0.20846, -0.22425, 0.025931, 0.34381, -0.033583, 0.21888, -0.034532, 0.44045, 0.073997, 0.19108, -0.0037315, 0.32066, -0.10322, 0.28571), CPREDI = c(10, 6.0691, 10, 6.0653, 10, 6.0658, 10, 6.0661, 10, 6.0653, 10, 6.0653, 10, 6.0654, 10, 6.0653, 10, 6.0653, 10, 6.0653), CRESI = c(0.68, -2.3854, 0.402, 0.38867, -0.1186, -0.20927, -0.6592, -0.44517, 0.082, 0.69298, -0.1062, 0.43958, -0.1092, 0.89029, 0.234, 0.38347, -0.0118, 0.64589, -0.3264, 0.57489), CWRESI = c(0.21503, -1.1756, 0.12712, 0.19228, -0.037505, -0.10343, -0.20846, -0.21993, 0.025931, 0.34301, -0.033583, 0.21749, -0.034532, 0.44082, 0.073997, 0.18971, -0.0037315, 0.31967, -0.10322, 0.2845), EPRED = c(10, 5.9937, 10, 6.0137, 10, 5.9991, 10, 6.072, 10, 5.975, 10, 6.0976, 10, 6.0263, 10, 6.0775, 10, 6.0322, 10, 6.0261), ERES = c(0.68, -2.31, 0.402, 0.44029, -0.1186, -0.14263, -0.6592, -0.45105, 0.082, 0.78326, -0.1062, 0.40725, -0.1092, 0.92936, 0.234, 0.37133, -0.0118, 0.67904, -0.3264, 0.6141), EWRES = c(0.21503, -1.1547, 0.12712, 0.21859, -0.037505, -0.071498, -0.20846, -0.22253, 0.025931, 0.39309, -0.033583, 0.19972, -0.034532, 0.46276, 0.073997, 0.18277, -0.0037315, 0.33789, -0.10322, 0.30445), ECWRES = c(0.21503, -1.1868, 0.12712, 0.21855, -0.037505, -0.071919, -0.20846, -0.2271, 0.025931, 0.3887, -0.033583, 0.20201, -0.034532, 0.45969, 0.073997, 0.18447, -0.0037315, 0.33725, -0.10322, 0.30406), NPDE = c(0.15943, -1.2628, 0.14252, 0.23613, 0.041789, -0.15097, -0.23613, -0.30548, -0.14252, 0.4677, 0.041789, 0.20189, -0.17637, 0.59277, 0.14252, 0.23613, 0.0083555, 0.4677, -0.14252, 0.34069), NPD = c(0.21503, -1.1678, 0.12712, 0.24692, -0.037505, -0.045664, -0.20846, -0.19727, 0.025931, 0.41814, -0.033584, 0.22571, -0.034532, 0.48374, 0.073997, 0.21048, -0.0037315, 0.3651, -0.10322, 0.32815), TIME = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1), ETA1 = c(0.071545, 0.071545, 0.0057045, 0.0057045, 0.024947, 0.024947, 0.03195, 0.03195, -0.0047874, -0.0047874, 0.0039789, 0.0039789, -0.011804, -0.011804, 0.0058801, 0.0058801, -0.0031367, -0.0031367, -0.00066641, -0.00066641),     IRES = c(0.68, -2.1608, 0.402, 0.40602, -0.1186, -0.13268,     -0.6592, -0.34674, 0.082, 0.67849, -0.1062, 0.45167, -0.1092,     0.8547, 0.234, 0.40135, -0.0118, 0.63639, -0.3264, 0.57287    ), IWRES = c(0.21503, -1.1691, 0.12712, 0.21229, -0.037505,     -0.070055, -0.20846, -0.18374, 0.025931, 0.3529, -0.033583,     0.23596, -0.034532, 0.44301, 0.073997, 0.20987, -0.0037315,     0.33128, -0.10322, 0.29858), OBJI = c(5.2018, 5.2018, 3.776,     3.776, 3.7206, 3.7206, 3.7958, 3.7958, 3.848, 3.848, 3.7726,     3.7726, 3.9287, 3.9287, 3.7643, 3.7643, 3.831, 3.831, 3.8188,     3.8188), DV = c(10.68, 3.6837, 10.402, 6.454, 9.8814, 5.8565,     9.3408, 5.6209, 10.082, 6.7583, 9.8938, 6.5049, 9.8908, 6.9557,     10.234, 6.4488, 9.9882, 6.7112, 9.6736, 6.6402), PRED = c(10,     6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10,     6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653, 10, 6.0653),     RES = c(0.68, -2.3816, 0.402, 0.38869, -0.1186, -0.20881,     -0.6592, -0.44441, 0.082, 0.69299, -0.1062, 0.43959, -0.1092,     0.89039, 0.234, 0.38349, -0.0118, 0.64589, -0.3264, 0.57489    ), WRES = c(0.21503, -1.1786, 0.12712, 0.19235, -0.037505,     -0.10333, -0.20846, -0.21992, 0.025931, 0.34294, -0.033583,     0.21754, -0.034532, 0.44062, 0.073997, 0.18978, -0.0037315,     0.31963, -0.10322, 0.28449)), class = "data.frame", .Names = c("IPRE", "ID", "NPRED", "NRES", "NWRES", "PREDI", "RESI", "WRESI", "CPRED", "CRES", "CWRES", "CPREDI", "CRESI", "CWRESI", "EPRED", "ERES", "EWRES", "ECWRES", "NPDE", "NPD", "TIME", "ETA1", "IRES", "IWRES", "OBJI", "DV", "PRED", "RES", "WRES"), row.names = c(NA, -20L));

    pred = function() ipre

    inits = list(THTA=c(0.5))
    inits$OMGA=list(ETA[1]~.04)

    dat <- Wang2007
    dat$DV <- dat$Y

    fit.prop <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod=""))

    fit.prop2 <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1))},
                           control=foceiControl(maxOuterIterations=0,covMethod="", optExpression = FALSE))


    test_that("Matches NONMEM objective proportional function; (Based on Wang2007)", {
        expect_equal(round(fit.prop$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
        expect_equal(round(fit.prop$`ETA[1]`, 4), round(out.focei.prop$ETA1, 4)) # match NONMEM output
        ## Individual properties
        expect_equal(round(fit.prop$IPRED, 4), round(out.focei.prop$IPRE, 4))
        expect_equal(round(fit.prop$IRES, 4), round(out.focei.prop$IRES, 4))
        expect_equal(round(fit.prop$IWRES, 4), round(out.focei.prop$IWRES, 4))
        ## WRES variants
        expect_equal(round(fit.prop$PRED, 4), round(out.focei.prop$NPRED, 4)) # matches output of PRED from NONMEM
        expect_equal(round(fit.prop$PRED, 4), round(out.focei.prop$PRED, 4)) # matches output of PRED from NONMEM
        expect_equal(round(fit.prop$RES, 4), round(out.focei.prop$RES, 4)) # match NONMEM output
        expect_equal(round(fit.prop$RES, 4), round(out.focei.prop$NRES, 4)) # match NONMEM output
        ## FOI equivalents
        expect_equal(round(fit.prop$PRED, 4), round(out.focei.prop$PREDI, 4)) # matches output of PRED from NONMEM
        ## CWRES variants
        expect_equal(round(fit.prop$CRES, 4), round(out.focei.prop$CRES, 4)) # match NONMEM output
        expect_equal(round(fit.prop$CPRED, 4), round(out.focei.prop$CPRED, 4)) # match NONMEM output
        expect_equal(round(fit.prop$CWRES, 4), round(out.focei.prop$CWRES, 4)) # match NONMEM output
        ## Note that E[x] for CPRED and CPREDI are equal
        expect_equal(round(fit.prop$CRES, 4), round(out.focei.prop$CRESI, 4)) # match NONMEM output
        expect_equal(round(fit.prop$CPRED, 4), round(out.focei.prop$CPREDI, 4)) # match NONMEM output
    })
    test_that("Matches NONMEM objective proportional function; (Based on Wang2007; unoptimized)", {
        # Check unoptimized expression
        expect_equal(round(fit.prop2$objective, 3), 39.458); # Matches Table 2 Prop FOCEI for NONMEM
        expect_equal(round(fit.prop2$`ETA[1]`, 4), round(out.focei.prop$ETA1, 4)) # match NONMEM output
        ## Individual properties
        expect_equal(round(fit.prop2$IPRED, 4), round(out.focei.prop$IPRE, 4))
        expect_equal(round(fit.prop2$IRES, 4), round(out.focei.prop$IRES, 4))
        expect_equal(round(fit.prop2$IWRES, 4), round(out.focei.prop$IWRES, 4))
        ## WRES variants
        expect_equal(round(fit.prop2$PRED, 4), round(out.focei.prop$NPRED, 4)) # matches output of PRED from NONMEM
        expect_equal(round(fit.prop2$PRED, 4), round(out.focei.prop$PRED, 4)) # matches output of PRED from NONMEM
        expect_equal(round(fit.prop2$RES, 4), round(out.focei.prop$RES, 4)) # match NONMEM output
        expect_equal(round(fit.prop2$RES, 4), round(out.focei.prop$NRES, 4)) # match NONMEM output
        ## FOI equivalents
        expect_equal(round(fit.prop2$PRED, 4), round(out.focei.prop$PREDI, 4)) # matches output of PRED from NONMEM
        ## CWRES variants
        expect_equal(round(fit.prop2$CRES, 4), round(out.focei.prop$CRES, 4)) # match NONMEM output
        expect_equal(round(fit.prop2$CPRED, 4), round(out.focei.prop$CPRED, 4)) # match NONMEM output
        expect_equal(round(fit.prop2$CWRES, 4), round(out.focei.prop$CWRES, 4)) # match NONMEM output
        ## Note that E[x] for CPRED and CPREDI are equal
        expect_equal(round(fit.prop2$CRES, 4), round(out.focei.prop$CRESI, 4)) # match NONMEM output
        expect_equal(round(fit.prop2$CPRED, 4), round(out.focei.prop$CPREDI, 4)) # match NONMEM output
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

    fit.prop2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1))},
                      control=foceiControl(maxOuterIterations=0,covMethod=""))


    test_that("Matches NONMEM objective proportional function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.prop2$objective, 3), 39.458);
    })


    fit.prop <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective proportional error FOCE (Based on Wang2007)", {
        expect_equal(round(fit.prop$objective, 3), 39.207);
    })

    etaMat <- as.matrix(ranef(fit.prop)[, -1, drop = FALSE])

    fit.prop2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective proportional error FOCE; ODE (Based on Wang2007)", {
        expect_equal(round(fit.prop2$objective, 3), 39.279);
    })


    fit.prop <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective proportional error FO (Based on Wang2007)", {
        expect_equal(round(fit.prop$objective, 3), 39.213);
    })

    fit.prop2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective proportional error FO; ODE (Based on Wang2007)", {
        expect_equal(round(fit.prop2$objective, 3), 39.213); # Matches Table 2 Prop FOCEI for NONMEM
    })


#### Add

    fit.add <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(add(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective additive error; (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059); # Matches Table 2 Add FOCEI for NONMEM
    })

    fit.add2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(add(.1))},
                      control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective additive error; ODE (Based on Wang2007)", {
        expect_equal(round(fit.add2$objective, 3), -2.059);
    })

    fit.add <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(add(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective additive error FOCE (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), -2.059);
    })


    fit.add2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(add(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective additive error FOCE; ODE (Based on Wang2007)", {
        expect_equal(round(fit.add2$objective, 3), -2.059);
    })


    fit.add <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(add(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective additive error FO (Based on Wang2007)", {
        expect_equal(round(fit.add$objective, 3), 0.026);
    })

    fit.add2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(add(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective additive error FO; ODE (Based on Wang2007)", {
        expect_equal(round(fit.add2$objective, 3), 0.026);
    })

### Add+Prop

    fit.addprop <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(add(.1)+prop(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective additive+proportional function; (Based on Wang2007)", {
        expect_equal(round(fit.addprop$objective, 3), 39.735);
    })

    fit.addprop2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(add(.1)+prop(.1))},
                      control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective additive+proportional function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.addprop2$objective, 3), 39.735);
    })

    fit.addprop <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(add(.1)+prop(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective additive+proportional error FOCE (Based on Wang2007)", {
        expect_equal(round(fit.addprop$objective, 3), 39.499);
    })

    fit.addprop2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(add(.1)+prop(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective additive+proportional error FOCE; ODE (Based on Wang2007)", {
        expect_equal(round(fit.addprop2$objective, 3), 39.563);
    })

    fit.addprop <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(add(.1)+prop(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective additive+proportional error FO (Based on Wang2007)", {
        expect_equal(round(fit.addprop$objective, 3), 39.505);
    })

    fit.addprop2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(add(.1)+prop(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective additive+proportional error FO; ODE (Based on Wang2007)", {
        expect_equal(round(fit.addprop2$objective, 3), 39.505);
    })

    ## lognormal -- equivalent to add on log-space and back-transformed.

    ## Next run on the log-transformed space
    datl <- dat; datl$DV <- log(datl$DV)
    datl2 <- dat2; datl2$DV <- log(datl2$DV)
    predl = function() log(ipre)


    fit.lnorm <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(lnorm(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective lognormal function; (Based on Wang2007)", {
        expect_equal(round(fit.lnorm$objective, 3), 40.039);
    })

    fit.lnorm0 <- .foceiFit(datl, inits, mypar1,mod,predl,function(){return(add(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective lognormal function; (Based on Wang2007)", {
        expect_equal(round(fit.lnorm$objective, 3), 40.039);
        expect_equal(fit.lnorm0$objective + 2*sum(datl$DV), fit.lnorm$objective);
        expect_equal(round(fit.lnorm0$objective, 3), -42.106);
    })

    fit.lnorm2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(lnorm(.1))},
                           control=foceiControl(maxOuterIterations=0,covMethod=""))

    fit.lnorm20 <- .foceiFit(datl2, inits, mypar1, m1, predl,function(){return(add(.1))},
                      control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective lognormal function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.lnorm2$objective, 3), 40.039);
        expect_equal(fit.lnorm20$objective + 2*sum(datl$DV), fit.lnorm2$objective);
        expect_equal(round(fit.lnorm20$objective, 3), -42.106);
    })

    fit.lnorm <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(lnorm(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    fit.lnorm0 <- .foceiFit(datl, inits, mypar1,mod,predl,function(){return(add(.1))},
                           control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective lognormal error FOCE (Based on Wang2007)", {
        expect_equal(round(fit.lnorm$objective, 3), 40.039);
        expect_equal(fit.lnorm0$objective + 2*sum(datl$DV), fit.lnorm$objective);
        expect_equal(round(fit.lnorm0$objective, 3), -42.106);
    })

    fit.lnorm2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(lnorm(.1))},
                           control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    fit.lnorm20 <- .foceiFit(datl2, inits, mypar1,m1,predl,function(){return(add(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective lognormal error FOCE; ODE (Based on Wang2007)", {
        expect_equal(round(fit.lnorm2$objective, 3), 40.039);
        expect_equal(fit.lnorm20$objective + 2*sum(datl$DV), fit.lnorm2$objective);
        expect_equal(round(fit.lnorm20$objective, 3), -42.106);
    })

    fit.lnorm <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(lnorm(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    fit.lnorm0 <- .foceiFit(datl, inits, mypar1,mod,predl,function(){return(add(.1))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective lognormal error FO (Based on Wang2007)", {
        expect_equal(round(fit.lnorm$objective, 3), 40.055);
        expect_equal(fit.lnorm0$objective + 2*sum(datl$DV), fit.lnorm$objective);
        expect_equal(round(fit.lnorm0$objective, 3), -42.09);
    })

    fit.lnorm2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(lnorm(.1))},
                           control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    fit.lnorm20 <- .foceiFit(datl2, inits, mypar1,m1,predl,function(){return(add(.1))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective lognormal error FO; ODE (Based on Wang2007)", {
        expect_equal(round(fit.lnorm2$objective, 3), 40.055);
        expect_equal(fit.lnorm20$objective + 2*sum(datl$DV), fit.lnorm2$objective);
        expect_equal(round(fit.lnorm20$objective, 3), -42.09);
    })

### Now try  TBS boxCox

    fit.boxCox <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1)+boxCox(.5))},
                         control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective Box-Cox function; (Based on Wang2007)", {
        expect_equal(round(fit.boxCox$objective, 3), 61.473);
    })

    fit.boxCox2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1)+boxCox(.5))},
                      control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective Box-Cox function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.boxCox2$objective, 3), 61.473);
    })

    fit.boxCox <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1)+boxCox(.5))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective Box-Cox error FOCE (Based on Wang2007)", {
        expect_equal(round(fit.boxCox$objective, 3), 61.324);
    })

    fit.boxCox2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1)+boxCox(.5))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective Box-Cox error FOCE; ODE (Based on Wang2007)", {
        expect_equal(round(fit.boxCox2$objective, 3), 61.28);
    })

    fit.boxCox <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1)+boxCox(.5))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective Box-Cox error FO (Based on Wang2007)", {
        expect_equal(round(fit.boxCox$objective, 3), 61.325);
    })

    fit.boxCox2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1)+boxCox(.5))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective Box-Cox error FO; ODE (Based on Wang2007)", {
        expect_equal(round(fit.boxCox2$objective, 3), 61.325);
    })

    ### Now try  TBS

    fit.yeoJohnson <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1)+yeoJohnson(.5))},
                         control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective Yeo-Johnson function; (Based on Wang2007)", {
        expect_equal(round(fit.yeoJohnson$objective, 3), 62.821);
    })

    fit.yeoJohnson2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1)+yeoJohnson(.5))},
                      control=foceiControl(maxOuterIterations=0,covMethod=""))

    test_that("Matches NONMEM objective Yeo-Johnson function; ODE (Based on Wang2007)", {
        expect_equal(round(fit.yeoJohnson2$objective, 3), 62.821);
    })

    fit.yeoJohnson <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1)+yeoJohnson(.5))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective Yeo-Johnson error FOCE (Based on Wang2007)", {
        expect_equal(round(fit.yeoJohnson$objective, 3), 62.676);
    })

    fit.yeoJohnson2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1)+yeoJohnson(.5))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE))

    test_that("Matches NONMEM objective Yeo-Johnson error FOCE; ODE (Based on Wang2007)", {
        expect_equal(round(fit.yeoJohnson2$objective, 3), 62.627);
    })

    fit.yeoJohnson <- .foceiFit(dat, inits, mypar1,mod,pred,function(){return(prop(.1)+yeoJohnson(.5))},
                         control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective Yeo-Johnson error FO (Based on Wang2007)", {
        expect_equal(round(fit.yeoJohnson$objective, 3), 62.677);
    })

    fit.yeoJohnson2 <- .foceiFit(dat2, inits, mypar1,m1,pred,function(){return(prop(.1)+yeoJohnson(.5))},
                          control=foceiControl(maxOuterIterations=0,covMethod="", interaction=FALSE, fo=TRUE))

    test_that("Matches NONMEM objective Yeo-Johnson error FO; ODE (Based on Wang2007)", {
        expect_equal(round(fit.yeoJohnson2$objective, 3), 62.677);
    })

}, on.validate="NLMIXR_VALIDATION", silent=TRUE)
