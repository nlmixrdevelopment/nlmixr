source("helper-prep_fit.R")
context("NLME14: one-compartment infusion, steady state")
runno <- "N014_ode"

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
    start = c(lCL = 1.4, lV = 4.3)
  )

ode1 <- "
d/dt(centr)  = -(CL/V)*centr;
"

mypar1 <- function(lCL, lV)
{
  CL <- exp(lCL)
  V <- exp(lV)
}

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
  nlme_ode(
    dat,
    model = ode1,
    par_model = specs1,
    par_trans = mypar1,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-12668.12, 25346.24, 25374.89),
    param=c(1.3862, 4.2647),
    stdev_param=c(1.4077, 1.518),
    sigma=c(0.1963)
  )
