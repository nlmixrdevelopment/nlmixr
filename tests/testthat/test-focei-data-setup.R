rxPermissive({
    ## Test the behavior
    context("Test FOCEi Data Setup (for RcppParallel-style for loop); 0 cov")
    dat <- suppressWarnings({
        nmDataConvert(Bolus_1CPT[, names(Bolus_1CPT) != "SS"])
    });
    test_that("conversion without covariates", {
        convert1 <- foceiDataSetup(dat);
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        for (i in unique(dat$ID)){
            w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert1$obs[w, ])))
            w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert1$dose[w, ])))
        }
    })

    cn <- c("V")
    convert2 <- foceiDataSetup(dat, cn);
    context("Test FOCEi Data Setup (for RcppParallel-style for loop); 1 cov")
    test_that("conversion with 1 covariate", {
        expect_equal(length(convert2$cov.names), 1);
        for (i in unique(dat$ID)){
            w <- seq(convert2$ids$posObs[i] + 1,
                     convert2$ids$posObs[i] + convert2$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert2$obs[w, ])))
            w <- seq(convert2$ids$posDose[i] + 1,
                     convert2$ids$posDose[i] + convert2$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert2$dose[w, ])));
            w <- seq(convert2$ids$posCov[i] + 1,
                     convert2$ids$posCov[i] + convert2$ids$nCov[i])
            expect_equal(as.double(convert2$cov[w]),
                         as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,cn])));
        }
    })

    cn <- c("V", "CL")
    convert2 <- foceiDataSetup(dat, cn);
    context("Test FOCEi Data Setup (for RcppParallel-style for loop); 2 cov")
    test_that("conversion with 2 covariate", {
        expect_equal(length(convert2$cov.names), 2);
        for (i in unique(dat$ID)){
            w <- seq(convert2$ids$posObs[i] + 1,
                     convert2$ids$posObs[i] + convert2$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert2$obs[w, ])))
            w <- seq(convert2$ids$posDose[i] + 1,
                     convert2$ids$posDose[i] + convert2$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert2$dose[w, ])));
            w <- seq(convert2$ids$posCov[i] + 1,
                     convert2$ids$posCov[i] + convert2$ids$nCov[i])
            expect_equal(as.double(convert2$cov[w]),
                         as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,cn])));
        }
    })

    cn <- c("V", "CL", "DOSE")
    convert2 <- foceiDataSetup(dat, cn);
    context("Test FOCEi Data Setup (for RcppParallel-style for loop); 3 cov")
    test_that("conversion with 3 covariate", {
        expect_equal(length(convert2$cov.names), 3);
        for (i in unique(dat$ID)){
            w <- seq(convert2$ids$posObs[i] + 1,
                     convert2$ids$posObs[i] + convert2$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert2$obs[w, ])))
            w <- seq(convert2$ids$posDose[i] + 1,
                     convert2$ids$posDose[i] + convert2$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert2$dose[w, ])));
            w <- seq(convert2$ids$posCov[i] + 1,
                     convert2$ids$posCov[i] + convert2$ids$nCov[i])
            expect_equal(as.double(convert2$cov[w]),
                         as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,cn])));
        }
    })
})
