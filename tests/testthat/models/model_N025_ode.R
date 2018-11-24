source("helper-prep_fit.R")
context("NLME25: one-compartment oral, multiple-dose")
runno<-"N025_ode"

datr <- Oral_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)

datr <- datr[datr$EVID!=2,]

ode1KA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs-(CL/V)*centr;
    "

mypar4 <- function(lCL, lV, lKA)
{
  CL <- exp(lCL)
  V <- exp(lV)
  KA <- exp(lKA)
}

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
  nlme_ode(
    dat,
    model = ode1KA,
    par_model = specs4,
    par_trans = mypar4,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-13285.23, 26580.47, 26609.12),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84730, 0),
    sigma=0.43887
  )
