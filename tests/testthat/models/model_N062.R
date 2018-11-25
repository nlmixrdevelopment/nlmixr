source("helper-prep_fit.R")
context("NLME62: two-compartment oral, steady-state, multiple-dose")
runno <- "N062"

datr <- read.csv("../Oral_2CPT.csv",
                 header = TRUE,
                 stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2, ]

datSS <- datr[datr$SS == 0, ]
datSD <- datr[datr$SS == 1, ]

#general solution to allow different times of SS dose and different II values per subject:
datSSD <- datr[datr$SS == 1, c("ID", "TIME", "II")]

datSSD$V0 <- datSSD$TIME
datSSD$V1 <- datSSD$TIME - datSSD$II
datSSD$V2 <- datSSD$TIME - 2 * datSSD$II
datSSD$V3 <- datSSD$TIME - 3 * datSSD$II
datSSD$V4 <- datSSD$TIME - 4 * datSSD$II
datSSD$V5 <- datSSD$TIME - 5 * datSSD$II
datSSD$V6 <- datSSD$TIME - 6 * datSSD$II
datSSD$TIME <- NULL
datSSD$II <- NULL

index <- reshape2::melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
index$variable <- NULL
index <- index[index$TIMED > 0, ]
index <- index[order(index$ID, index$TIMED), ]

#much easier solution if you know the time of SS dose and the II and if it is the same for all
#index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))

datSD2 <- merge(datSD, index, by = c("ID"), all = T)
datSD2$TIME <- datSD2$TIMED
datSD2$TIMED <- NULL

datSS <- rbind(datSS, datSD2)
datSS <- datSS[order(datSS$ID, datSS$TIME), ]
dat <- datSS

specs8i <-
  list(
    fixed = lCL + lV + lCLD + lVT + lKA ~ 1,
    random = pdDiag(value = diag(c(6, 6, 6, 6, 6)), form = lCL + lV + lCLD +
                      lVT + lKA ~ 1),
    start = c(
      lCL = 1.4,
      lV = 4.2,
      lCLD = 1.3,
      lVT = 3.9,
      lKA = 0.1
    )
  )

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs8i,
    ncmt = 2,
    oral = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-13243.36, 26508.73, 26571.78),
    param=c(1.3489, 4.2246, 1.4426, 3.8546, -0.051931),
    stdev_param=c(0.32225, 0.30297, 0.0021809, 0.30263, 0.0010639),
    sigma=0.20137
  )

