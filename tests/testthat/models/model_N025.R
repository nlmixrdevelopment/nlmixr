source("tests/testthat/models/helper-prep_fit.R")
context("NLME25: one-compartment oral, multiple-dose")
runno<-"N025"

datr <- Oral_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)

datr <- datr[datr$EVID!=2,]

specs4 <-
  list(
    fixed = lCL + lV + lKA ~ 1,
    random = pdDiag(lCL + lV + lKA ~ 1),
    start = c(lCL = 1, lV = 4, lKA = 0)
  )

datSS <- datr[datr$SS==0,]
datSD <- datr[datr$SS==1,]

#general solution to allow different times of SS dose and different II values per subject:
datSSD <- datr[datr$SS==1, c("ID","TIME","II")]

datSSD$V0<-datSSD$TIME
datSSD$V1<-datSSD$TIME-datSSD$II
datSSD$V2<-datSSD$TIME-2*datSSD$II
datSSD$V3<-datSSD$TIME-3*datSSD$II
datSSD$V4<-datSSD$TIME-4*datSSD$II
datSSD$V5<-datSSD$TIME-5*datSSD$II
datSSD$V6<-datSSD$TIME-6*datSSD$II
datSSD$TIME<-NULL
datSSD$II<-NULL

index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
index$variable <- NULL
index <- index[index$TIMED > 0,]
index<-index[order(index$ID,index$TIMED),]

#much easier solution if you know the time of SS dose and the II and if it is the same for all
#index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))

datSD2 <- merge(datSD, index, by = c("ID"), all=T)
datSD2$TIME <- datSD2$TIMED
datSD2$TIMED <- NULL

datSS <- rbind(datSS, datSD2)
datSS <- datSS[order(datSS$ID, datSS$TIME),]

dat <- datSS

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs4,
    ncmt = 1,
    oral = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values <-
  list(
    lik=c(-12633.52, 25281.05, 25321.15),
    param=c(1.3964, 4.215, 0.0013108),
    stdev_param=c(0.2706, 0.28435, 0.00041952),
    sigma=0.20936
  )
