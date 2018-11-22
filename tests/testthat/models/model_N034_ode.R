source("tests/testthat/models/helper-prep_fit.R")
context("NLME34: two-compartment bolus, steady-state")
runno <- "N034_ode"

datr <- Bolus_2CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2, ]

ode2 <- "
d/dt(centr)  = K21*periph-K12*centr-K10*centr;
d/dt(periph) =-K21*periph+K12*centr;
"

mypar6 <- function(lCL, lV, lCLD, lVT)
{
  CL <- exp(lCL)
  V  <- exp(lV)
  CLD <- exp(lCLD)
  VT <- exp(lVT)
  K10 <- CL / V
  K12 <- CLD / V
  K21 <- CLD / VT
}

specs6 <-
  list(
    fixed = lCL + lV + lCLD + lVT ~ 1,
    random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
    start = c(
      lCL = 1.3,
      lV = 4.19,
      lCLD = 1.5,
      lVT = 3.9
    )
  )
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

index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
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

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode2,
    par_model = specs6,
    par_trans = mypar6,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values <-
  list(
    lik=c(-13285.23, 26580.47, 26609.12),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84730, 0),
    sigma=0.43887
  )
