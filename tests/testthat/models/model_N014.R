source("helper-prep_fit.R")
context("NLME14: one-compartment infusion, steady state")
runno <- "N014"

datr <-
  read.csv("../Infusion_1CPT.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)
datr <- datr[datr$EVID != 2,]

datSSobs <- datr[datr$SS == 0,]
datSD <- datr[datr$SS == 1,]

#general solution to allow different times of SS dose and different II values per subject:
datSSD <- datSD[, c("ID","TIME","II")]
specs1 <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.5, lV = 4)
  )

#updates datSSD with 7 columns to account for the new dosing times
datSSD$V0<-datSSD$TIME
datSSD$V1<-datSSD$TIME-datSSD$II
datSSD$V2<-datSSD$TIME-2*datSSD$II
datSSD$V3<-datSSD$TIME-3*datSSD$II
datSSD$V4<-datSSD$TIME-4*datSSD$II
datSSD$V5<-datSSD$TIME-5*datSSD$II
datSSD$V6<-datSSD$TIME-6*datSSD$II
datSSD$TIME<-NULL
datSSD$II<-NULL

index <- reshape2::melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
index$variable <- NULL
index <- index[index$TIMED > 0,]
index<-index[order(index$ID,index$TIMED),]

datSD2 <- merge(datSD, index, by = c("ID"), all = TRUE)
datSD2$TIME <- datSD2$TIMED
datSD2$TIMED <- NULL

datSDoff <- datSD2
datSDoff$TIME <- datSDoff$TIME + datSDoff$AMT / datSDoff$RATE
datSDoff$AMT <- -1 * datSDoff$AMT

datSD2 <- rbind(datSD2, datSDoff)
datSS <- rbind(datSSobs, datSD2)
datSS <- datSS[order(datSS$ID,datSS$TIME),]

dat <- datSS

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs1,
    ncmt = 1,
    oral = FALSE,
    infusion = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-14247.56, 28505.12, 28533.77),
    param=c(1.2821, 4.155),
    stdev_param=c(5.5601e-05, 0.49595),
    sigma=0.41285
  )
