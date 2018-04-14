rxPermissive({

    context("Test NONMEM conversion");

    d1 <- data.frame(DATE=c("10-1-86", "10-1-86", "10-2-86"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)

    test_that("Throw error with missing DV", {
        expect_error(nmDataConvert(d1))
    })

    d1$DV <- 0;

    d2 <- rbind(data.frame(ID=1, d1, stringsAsFactors=F), data.frame(ID=2, d1, stringsAsFactors=F))
    d2[d2$ID == 2, "DATE"] <- gsub("^10", "11", d2[d2$ID == 2, "DATE"]);

    d3 <- d1;
    d3$DATE <- c("10-1-1986", "10-1-1986", "10-2-1986");

    d4 <- d1;
    d4$DATE <- c("10 1 1986", "10/1/86", "10-2-1986");

    test_that("DATE conversion works correctly", {
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d1))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25, 0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d2))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d3))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d4))$TIME)
    })


    ## Dat1= day month year
    d1 <- data.frame(DV=0, DAT1=c("1-10-86", "1-10-86", "2-10-86"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)

    d2 <- rbind(data.frame(ID=1, d1, stringsAsFactors=F), data.frame(ID=2, d1, stringsAsFactors=F))
    d2[d2$ID == 2, "DAT1"] <- gsub("-10-", "-11-", d2[d2$ID == 2, "DAT1"]);

    d3 <- d1;
    d3$DAT1 <- c("1-10-1986", "1-10-1986", "2-10-1986");

    d4 <- d1;
    d4$DAT1 <- c("1-10-1986", "1-10-86", "2-10-1986");

    test_that("DAT1 conversion works correctly", {
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d1))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25, 0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d2))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d3))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d4))$TIME)
    })

    ## Dat2 = year month day

    d1 <- data.frame(DAT2=c("86-10-1", "86-10-1", "86-10-2"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)
    d1$DV <- 0;

    d2 <- rbind(data.frame(ID=1, d1, stringsAsFactors=F), data.frame(ID=2, d1, stringsAsFactors=F))
    d2[d2$ID == 2, "DAT2"] <- gsub("-10-", "-11-", d2[d2$ID == 2, "DAT2"]);

    d3 <- d1;
    d3$DAT2 <- c("1986-10-1", "1986-10-1", "1986-10-2");

    d4 <- d1;
    d4$DAT2 <- c("1986-10-1", "86-10-1", "1986-10-2");

    test_that("DAT2 conversion works correctly", {
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d1))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25, 0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d2))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d3))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d4))$TIME)
    })

    ## DAT3 conversion
    d1 <- data.frame(DAT3=c("86-1-10", "86-1-10", "86-2-10"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)
    d1$DV <- 0;

    d2 <- rbind(data.frame(ID=1, d1, stringsAsFactors=F), data.frame(ID=2, d1, stringsAsFactors=F))
    d2[d2$ID == 2, "DAT3"] <- gsub("-10$", "-11", d2[d2$ID == 2, "DAT3"]);

    d3 <- d1;
    d3$DAT3 <- c("1986-1-10", "1986-1-10", "1986-2-10");

    d4 <- d1;
    d4$DAT3 <- c("1986-1-10", "86-1-10", "1986-2-10");

    test_that("DAT3 conversion works correctly", {
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d1))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25, 0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d2))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d3))$TIME)
        expect_equal(c(0, 5.41666666666667, 23.25), suppressWarnings(nmDataConvert(d4))$TIME)
    })

    d1 <- data.frame(DV=0, DATE=c("10-1-86", "10-1-86", "10-2-86"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)

    d2 <- d1;
    d2$DAT1 <- d2$DATE

    d3 <- d1;
    d3$DAT2 <- d3$DATE

    d4 <- d1;
    d4$DAT3 <- d4$DATE

    test_that("Multiple DATE errors", {
        expect_error(nmDataConvert(d2), rex::rex("Dates can only be specified by one of DATE, DAT1, DAT2, DAT3 / TIME;  This dataset has multiple DATE columns."))
        expect_error(nmDataConvert(d3), rex::rex("Dates can only be specified by one of DATE, DAT1, DAT2, DAT3 / TIME;  This dataset has multiple DATE columns."))
        expect_error(nmDataConvert(d4), rex::rex("Dates can only be specified by one of DATE, DAT1, DAT2, DAT3 / TIME;  This dataset has multiple DATE columns."))
    })


    d1 <- data.frame(DV=0, DATE=c("10-1-86", "10-1-86", "10-2-86"), TIME=c("9.15", "14:40", "8:30"), stringsAsFactors=F)

    test_that("Bad Date/Time combination", {
        expect_error(nmDataConvert(d1), rex::rex("The date time format was not correctly specified."))
    })

    d1 <- data.frame(DATE=c("10-1-86", "10-1-86", "10-2-86"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)
    d1$DV <- 0;

    d2 <- rbind(data.frame(ID=1, d1, stringsAsFactors=F), data.frame(ID=2, d1, stringsAsFactors=F))
    d2[d2$ID == 2, "DATE"] <- gsub("^10", "11", d2[d2$ID == 2, "DATE"]);


    d3 <- d2
    d3$AMT <- 0;

    test_that("EVID warning", {
        expect_warning(nmDataConvert(d3), rex::rex("Assumed all DV values are observations. (EVID=0)"));
    })

    d3 <- d2
    d3$EVID <- 0;
    test_that("AMT warning", {
        expect_warning(nmDataConvert(d3), rex::rex("Assuming AMT=0."));
    })

    test_that("data frame error", {
        expect_error(nmDataConvert(as.matrix(d3)), rex::rex("Input must be a data frame."));
    })

    d3 <- d2
    d3$EVID <- 5;
    test_that("Bad EVID", {
        expect_error(nmDataConvert(d3), "EVIDs do not match NONMEM-style EVIDs.");
    })


    d3 <- d2
    d3$EVID <- 3;
    d4 <- d2
    d4$EVID <- 4;

    test_that("Reset EVID", {
        expect_error(nmDataConvert(d3), "A reset of the ODE system is not currently supported.");
        expect_error(nmDataConvert(d4), "A reset of the ODE system is not currently supported.");
    })

    d3 <- d2;
    d3$EVID <- 0
    d3$EVID[3] <- 2;
    d3$AMT <- 0;

    test_that("Ignore EVID=2", {
        expect_warning(nmDataConvert(d3), "EVID=2 is dropped from RxODE/nlmixr datasets.");
        expect_equal(suppressWarnings(nmDataConvert(d3))$TIME, c(0, 5.41666666666667, 0, 5.41666666666667, 23.25));
    })

    d3 <- d2;
    d3$MDV <- 0;
    d3$MDV[c(1, 4)] <- 1
    d3$AMT <- 0;
    d3$AMT[c(1, 4)] <- 1

    test_that("Add CMT=1", {
        expect_equal(suppressWarnings(nmDataConvert(d3))$EVID, c(101, 0, 0, 101, 0, 0));
    })


    d1 <- data.frame(DV=0, DATE=c("10-1-86", "10-1-86", "10-2-86"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)
    d1$SS <- 0;

    d2 <- d1;
    d2$SS <- 0;
    d2$SS[1] <- 1;

    test_that("Steady State Checks", {
        expect_error(suppressWarnings(nmDataConvert(d1)), NA)
        expect_error(nmDataConvert(d2), rex::rex("Steady State dosing events are not currently supported."))
    })


    d1 <- data.frame(DV=0, DATE=c("10-1-86", "10-1-86", "10-2-86"), TIME=c("9:15", "14:40", "8:30"), stringsAsFactors=F)
    d1$RATE <- 0;

    d2 <- d1;
    d2$RATE <- 0;
    d2$RATE[1] <- -1;


    test_that("Calculated Rate checks", {
        expect_error(suppressWarnings(nmDataConvert(d1)), NA)
        expect_error(nmDataConvert(d2), rex::rex("RxODE does not currently estimating infusion rates or durations. (Found RATE < 0)"))
    })


    test_that("Rik's datasets", {
        expect_error(nmDataConvert(Bolus_1CPT), "Steady State dosing events are not currently supported.");
        expect_equal(suppressWarnings(unique(nmDataConvert(Bolus_1CPT[,names(Bolus_1CPT) != "SS"])$EVID)), c(101, 0))
        expect_equal(suppressWarnings(unique(nmDataConvert(Infusion_1CPT[,names(Infusion_1CPT) != "SS"])$EVID)), c(10101, 0))
        expect_warning(nmDataConvert(Bolus_1CPT[,names(Bolus_1CPT) != "SS"]), "EVID=2 is dropped from RxODE/nlmixr datasets.")
        expect_warning(nmDataConvert(Infusion_1CPT[,names(Infusion_1CPT) != "SS"]), "EVID=2 is dropped from RxODE/nlmixr datasets.")
    })

    test_that("100+ cmt", {
        tmp <- Infusion_1CPT[,names(Infusion_1CPT) != "SS"]
        tmp$CMT[tmp$CMT == 1] <- 123;
        expect_equal(suppressWarnings(unique(nmDataConvert(tmp)$EVID)), c(112301L, 0L))
        tmp <- Bolus_1CPT[,names(Bolus_1CPT) != "SS"]
        tmp$CMT[tmp$CMT == 1] <- 123;
        expect_equal(suppressWarnings(unique(nmDataConvert(tmp)$EVID)), c(102301L, 0L))
    })

}, cran=TRUE)
