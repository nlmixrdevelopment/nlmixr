source("tests/testthat/models/helper-prep_fit.R")
context("NLME03: one-compartment bolus, steady-state")
runno <- "N003"

datr <- Bolus_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

specs1 <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.6, lV = 4.5)
  )

datSS <- datr[datr$SS == 0,]
datSD <- datr[datr$SS == 1,]

# general solution to allow different times of SS dose and different II values per subject:
datSSD <- datr[datr$SS == 1, c("ID", "TIME", "II")]

nvar <- function(TIME, II) {
  for (i in seq(1, 7)) {
    r = TIME - (i - 1) * II
    assign(paste("ret", i, sep = ""), r)
  }
  return(list(
    r1 = ret1,
    r2 = ret2,
    r3 = ret3,
    r4 = ret4,
    r5 = ret5,
    r6 = ret6,
    r7 = ret7
  ))
}

#updates datSSD with 7 columns to account for the new dosing times
datSSD[,(paste("V",seq(1,7)-1,sep=""))] <- nvar(datSSD$TIME, datSSD$II)
datSSD <- datSSD[,c(1,4:10)]

index <- reshape2::melt(datSSD, id.vars="ID",value.name="TIMED")
index <- index[,-2]
index <- index[index$TIMED>0,]

#much easier solution if you know the time of SS dose and the II and if it is the same for all
#index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))

datSD2 <- merge(datSD,index,by=c("ID"),all=TRUE)
datSD2$TIME <- datSD2$TIMED
datSD2 <- datSD2[,-15]

datSS <- rbind(datSS,datSD2)
datSS <- datSS[order(datSS$ID, datSS$TIME),]

dat <- datSS

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs1,
    ncmt = 1,
    oral = FALSE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values <-
  list(
    lik=c(-12854.18, 25718.37, 25747.02),
    param=c(1.3585, 4.1934),
    stdev_param=c(1.3479, 1.5912),
    sigma=c(0.20027)
  )
