library(nlmixr)
library(data.table)

################################################################################################
################################################################################################
##
## 1 compartment IV Bolus
##
## For a bolus into compartment 1 (for ODEs) EVID needs to be 101
##
################################################################################################
################################################################################################

datr <- Bolus_1CPT;
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

specs1<-list(fixed=lCL+lV~1, random = pdDiag(lCL+lV~1), start=c(lCL=1.6,lV=4.5))

################################################################################################
## 1 compartment IV Bolus: single dose, closed form solution
################################################################################################

runno<-"N001"
dat<-datr[SD==1]
fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
intervals(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus: single dose, ODE solution
################################################################################################

ode1 <- "
d/dt(centr)  = -(CL/V)*centr;
"

mypar1 <- function(lCL, lV )
{
    CL = exp(lCL)
    V  = exp(lV)
}

fitODE <- nlme_ode(dat, model=ode1, par_model=specs1, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
intervals(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus: multiple dose, closed form solution
################################################################################################

runno<-"N002"
dat<-datr[SD==0]
fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
intervals(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus: multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1, par_model=specs1, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
intervals(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus: single and multiple dose, closed form solution
################################################################################################

runno<-"N004"
dat<-datr
fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
intervals(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus: single and multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1, par_model=specs1, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
intervals(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus: steady state data implemented using 7 doses, closed form solution
################################################################################################

runno<-"N003"
datSS<-datr[SS==0]
datSD<-datr[SS==1]
#general solution to allow different times of SS dose and different II values per subject:
datSSD<-datr[SS==1,.(ID,TIME,II)]

nvar<-function(TIME,II){
 for(i in seq(1,7)){
   r=TIME-(i-1)*II
   assign(paste("ret",i,sep=""),r)
 }
 return(list(r1 = ret1, r2 = ret2,r3=ret3,r4=ret4,r5=ret5,r6=ret6,r7=ret7))
}

#updates datSSD with 7 columns to account for the new dosing times
datSSD[,(paste("V",seq(1,7)-1,sep="")):=nvar(TIME,II)][,TIME:=NULL][,II:=NULL]
index<-melt(datSSD,id.vars="ID",value.name="TIMED")
index[,variable:=NULL]
index<-index[TIMED>0]

datSD2<-merge(datSD,index,by=c("ID"),all=TRUE)
datSD2[,TIME:=TIMED][,TIMED:=NULL]
datSS<-rbind(datSS,datSD2)
setkey(datSS,ID,TIME)
dat<-datSS

fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
intervals(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus: steady state data implemented using 7 doses, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1, par_model=specs1, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
intervals(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 1 compartment IV Infusion
##
## Note that infusions need to be coded both with a record to start infusions
## and another one to stop infusions!
## For infusions into compartment 1 (for ODEs) EVID needs to be 10101
##
################################################################################################
################################################################################################

datr <- read.csv("Infusion_1CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,10101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]
datIV<-datr[AMT>0][,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
datr<-rbind(datr,datIV)
setkey(datr,ID,TIME)

specs1<-list(fixed=lCL+lV~1, random = pdDiag(lCL+lV~1), start=c(lCL=1.5,lV=4))


################################################################################################
## 1 compartment IV Infusion: single dose, closed form solution
################################################################################################

runno<-"N012"
dat<-datr[SD==1]
fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
intervals(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion: single dose, ODE solution
################################################################################################

#slightly tweaked starting values to ensure convergence
specs1m<-list(fixed=lCL+lV~1, random = pdDiag(lCL+lV~1), start=c(lCL=1.3,lV=4))
fitODE <- nlme_ode(dat, model=ode1, par_model=specs1m, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
intervals(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion: multiple dose, closed form solution
################################################################################################

runno<-"N013"
dat<-datr[SD==0]
fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion: multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1, par_model=specs1, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion: single and multiple dose, closed form solution
################################################################################################

runno<-"N015"
dat<-datr
fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion: single and multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1, par_model=specs1, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion: steady state data implemented using 7 doses, closed form solution
################################################################################################

runno<-"N014"
datr <- read.csv("Infusion_1CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,10101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]
datSSobs<-datr[SS==0]
datSD<-datr[SS==1]
setkey(datSD,ID,TIME)
#general solution to allow different times of SS dose and different II values per subject:
datSSD<-datSD[,.(ID,TIME,II)]

nvar<-function(TIME,II){
 for(i in seq(1,7)){
   r=TIME-(i-1)*II
   assign(paste("ret",i,sep=""),r)
 }
 return(list(r1 = ret1, r2 = ret2,r3=ret3,r4=ret4,r5=ret5,r6=ret6,r7=ret7))
}

#updates datSSD with 7 columns to account for the new dosing times
datSSD[,(paste("V",seq(1,7)-1,sep="")):=nvar(TIME,II)][,TIME:=NULL][,II:=NULL]

index<-melt(datSSD,id.vars=c("ID"),value.name="TIMED")
index[,variable:=NULL]
index<-index[TIMED>0]
setkey(index,ID,TIMED)

datSD2<-merge(datSD,index,by=c("ID"),all=TRUE)
datSD2[,TIME:=TIMED][,TIMED:=NULL]

datSDoff<-copy(datSD2)
datSDoff[,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
datSD2<-rbind(datSD2,datSDoff)

datSS<-rbind(datSSobs,datSD2)
setkey(datSS,ID,TIME)

dat<-datSS
fit <- nlme_lin_cmpt(dat, par_model=specs1, ncmt=1, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion: steady state data implemented using 7 doses, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1, par_model=specs1, par_trans=mypar1, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 1 compartment IV Bolus, Michaelis-Menten elimination
##
## For a bolus into compartment 1 EVID needs to be 101
##
################################################################################################
################################################################################################

ode1MM <- "
d/dt(centr)  = -(VM*centr/V)/(KM+centr/V);
"

mypar3 <- function(lVM, lKM, lV )
{
    VM = exp(lVM)
    KM = exp(lKM)
    V  = exp(lV)
}
specs3 <- list(fixed=lVM+lKM+lV~1, random = pdDiag(lVM+lKM+lV~1), start=c(lVM=7, lKM=6, lV=4))

datr <- read.csv("BOLUS_1CPTMM.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]


################################################################################################
## 1 compartment IV Bolus, Michaelis-Menten elimination: single dose
################################################################################################

runno<-"N009"
dat<-datr[SD==1]
fitODE <- nlme_ode(dat, model=ode1MM, par_model=specs3, par_trans=mypar3, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus, Michaelis-Menten elimination: multiple dose
################################################################################################

runno<-"N010"
dat<-datr[SD==0]
fitODE <- nlme_ode(dat, model=ode1MM, par_model=specs3, par_trans=mypar3, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Bolus, Michaelis-Menten elimination: single and multiple dose
################################################################################################

runno<-"N011"
dat<-datr
fitODE <- nlme_ode(dat, model=ode1MM, par_model=specs3, par_trans=mypar3, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 1 compartment IV Infusion, Michaelis-Menten elimination
##
## Note that infusions need to be coded both with a record to start infusions
## and another one to stop infusions!
## For infusions into compartment 1 (for ODEs) EVID needs to be 10101
##
################################################################################################
################################################################################################

datr <- read.csv("Infusion_1CPTMM.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,10101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]
datIV<-datr[AMT>0][,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
datr<-rbind(datr,datIV)
setkey(datr,ID,TIME)


################################################################################################
## 1 compartment IV Infusion, Michaelis-Menten elimination: single dose
################################################################################################

runno<-"N020"
dat<-datr[SD==1]
fitODE <- nlme_ode(dat, model=ode1MM, par_model=specs3, par_trans=mypar3, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion, Michaelis-Menten elimination: multiple dose
################################################################################################

runno<-"N021"
dat<-datr[SD==0]
fitODE <- nlme_ode(dat, model=ode1MM, par_model=specs3, par_trans=mypar3, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment IV Infusion, Michaelis-Menten elimination: single and multiple dose
################################################################################################

runno<-"N022"
dat<-datr
fitODE <- nlme_ode(dat, model=ode1MM, par_model=specs3, par_trans=mypar3, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))
summary(fit)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 1 compartment elimination, first order absorption
##
## For a bolus into compartment 1 (for ODEs) EVID needs to be 101
##
################################################################################################
################################################################################################

datr <- read.csv("ORAL_1CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

specs4 <- list(fixed=lCL+lV+lKA~1, random = pdDiag(lCL+lV+lKA~1), start=c(lCL=1,lV=4,lKA=0))

################################################################################################
## 1 compartment oral administration: single dose, closed form solution
################################################################################################

runno<-"N023"
dat<-datr[SD==1]
fit <- nlme_lin_cmpt(dat, par_model=specs4, ncmt=1, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment oral administration: single dose, ODE solution
################################################################################################

ode1KA <- "
d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs-(CL/V)*centr;
"

mypar4 <- function(lCL, lV, lKA )
{
    CL = exp(lCL)
    V  = exp(lV)
    KA = exp(lKA)
}

fitODE <- nlme_ode(dat, model=ode1KA, par_model=specs4, par_trans=mypar4, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3,tolerance=1e-3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment oral administration: multiple dose, closed form solution
################################################################################################

runno<-"N024"
dat<-datr[SD==0]

#slightly tweaked initial estimates:
specs4i <- list(fixed=lCL+lV+lKA~1, random = pdDiag(value=diag(c(3,3,3)), form = lCL+lV+lKA~1), start=c(lCL=1.6,lV=4.5,lKA=0.2))

fit <- nlme_lin_cmpt(dat, par_model=specs4i, ncmt=1, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment oral administration: multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1KA, par_model=specs4i, par_trans=mypar4, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment oral administration: single and multiple dose, closed form solution
################################################################################################

runno<-"N026"
dat<-datr
fit <- nlme_lin_cmpt(dat, par_model=specs4, ncmt=1, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment oral administration: single and multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1KA, par_model=specs4, par_trans=mypar4, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment oral administration:
## steady state data implemented using 7 doses, closed form solution
################################################################################################

runno<-"N025"
datSS<-datr[SS==0]
datSD<-datr[SS==1]

#general solution to allow different times of SS dose and different II values per subject:
datSSD<-datr[SS==1,.(ID,TIME,II)]

nvar<-function(TIME,II){
 for(i in seq(1,7)){
   r=TIME-(i-1)*II
   assign(paste("ret",i,sep=""),r)
 }
 return(list(r1 = ret1, r2 = ret2,r3=ret3,r4=ret4,r5=ret5,r6=ret6,r7=ret7))
}

#updates datSSD with 7 columns to account for the new dosing times
datSSD[,(paste("V",seq(1,7)-1,sep="")):=nvar(TIME,II)][,TIME:=NULL][,II:=NULL]
index<-melt(datSSD,id.vars="ID",value.name="TIMED")
index[,variable:=NULL]
index<-index[TIMED>0]

datSD2<-merge(datSD,index,by=c("ID"),all=TRUE)
datSD2[,TIME:=TIMED][,TIMED:=NULL]
datSS<-rbind(datSS,datSD2)
setkey(datSS,ID,TIME)
dat<-datSS

fit <- nlme_lin_cmpt(dat, par_model=specs4, ncmt=1, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 1 compartment oral administration: steady state data implemented using 7 doses, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode1KA, par_model=specs4, par_trans=mypar4, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## first order absorption, 1 compartment distribution, Michaelis-Menten elimination
##
## For a bolus into compartment 1 (for ODEs) EVID needs to be 101
##
################################################################################################
################################################################################################


ode1MMKA <- "
d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs-(VM*centr/V)/(KM+centr/V);
"

mypar5 <- function(lVM, lKM, lV, lKA )
{
    VM = exp(lVM)
    KM = exp(lKM)
    V  = exp(lV)
    KA = exp(lKA)

}
specs5 <- list(fixed=lVM+lKM+lV+lKA~1, random = pdDiag(lVM+lKM+lV+lKA~1), start=c(lVM=7, lKM=6, lV=4,lKA=0))

datr <- read.csv("ORAL_1CPTMM.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]



################################################################################################
## first order absorption, 1 compartment distribution, Michaelis-Menten elimination: single dose
################################################################################################

runno<-"N029"
dat<-datr[SD==1]
fitODE <- nlme_ode(dat, model=ode1MMKA, par_model=specs5, par_trans=mypar5, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 1 compartment distribution, Michaelis-Menten elimination:
## multiple dose
################################################################################################

#slightly tweaked initial estimates, including starting estimates for IIV.
#The value is a scaled variance and is derived as: 30%IIV, 20% residual error
#Initial value=(0.3/0.2)^2=2.25:
specs5i <- list(fixed=lVM+lKM+lV+lKA~1, random = pdDiag(value=diag(c(2,2,2,2)), form = lVM+lKM+lV+lKA~1), start=c(lVM=7, lKM=6.2, lV=4.5,lKA=-0.2))

runno<-"N030"
dat<-datr[SD==0]
fitODE <- nlme_ode(dat, model=ode1MMKA, par_model=specs5i, par_trans=mypar5, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1,msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 1 compartment distribution, Michaelis-Menten elimination:
## single and multiple dose
################################################################################################

runno<-"N031"
dat<-datr
fitODE <- nlme_ode(dat, model=ode1MMKA, par_model=specs5, par_trans=mypar5, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 2 compartment IV Bolus
##
## For a bolus into compartment 1 (for ODEs) EVID needs to be 101
##
################################################################################################
################################################################################################

datr <- read.csv("BOLUS_2CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

specs6<-list(fixed=lCL+lV+lCLD+lVT~1, random = pdDiag(lCL+lV+lCLD+lVT~1), start=c(lCL=1.6,lV=4.5,lCLD=1.5,lVT=3.9))


################################################################################################
## 2 compartment IV Bolus: single dose, closed form solution
################################################################################################

runno<-"N032"
dat<-datr[SD==1]
fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus: single dose, ODE solution
################################################################################################

ode2 <- "
d/dt(centr)  = K21*periph-K12*centr-K10*centr;
d/dt(periph) =-K21*periph+K12*centr;
"

mypar6 <- function(lCL, lV, lCLD, lVT)
{
    CL = exp(lCL)
    V  = exp(lV)
    CLD= exp(lCLD)
    VT = exp(lVT)
    K10= CL/V
    K12= CLD/V
    K21= CLD/VT
}

fitODE <- nlme_ode(dat, model=ode2, par_model=specs6, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus: multiple dose, closed form solution
################################################################################################

runno<-"N033"
dat<-datr[SD==0]
fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus: multiple dose, ODE solution
################################################################################################

specs6i<-list(fixed=lCL+lV+lCLD+lVT~1, random = pdDiag(value=diag(c(3,3,3,3)), form = lCL+lV+lCLD+lVT~1), start=c(lCL=1.6,lV=4.5,lCLD=1.6,lVT=4))
fitODE <- nlme_ode(dat, model=ode2, par_model=specs6i, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1,msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus: single and multiple dose, closed form solution
################################################################################################

runno<-"N035"
dat<-datr
fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus: single and multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2, par_model=specs6, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus: steady state data implemented using 7 doses, closed form solution
################################################################################################

runno<-"N034"
datSS<-datr[SS==0]
datSD<-datr[SS==1]

#general solution to allow different times of SS dose and different II values per subject:
datSSD<-datr[SS==1,.(ID,TIME,II)]

nvar<-function(TIME,II){
 for(i in seq(1,7)){
   r=TIME-(i-1)*II
   assign(paste("ret",i,sep=""),r)
 }
 return(list(r1 = ret1, r2 = ret2,r3=ret3,r4=ret4,r5=ret5,r6=ret6,r7=ret7))
}

#updates datSSD with 7 columns to account for the new dosing times
datSSD[,(paste("V",seq(1,7)-1,sep="")):=nvar(TIME,II)][,TIME:=NULL][,II:=NULL]
index<-melt(datSSD,id.vars="ID",value.name="TIMED")
index[,variable:=NULL]
index<-index[TIMED>0]

datSD2<-merge(datSD,index,by=c("ID"),all=TRUE)
datSD2[,TIME:=TIMED][,TIMED:=NULL]
datSS<-rbind(datSS,datSD2)
setkey(datSS,ID,TIME)
dat<-datSS

fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus: steady state data implemented using 7 doses, ODE solution
################################################################################################

specs6<-list(fixed=lCL+lV+lCLD+lVT~1, random = pdDiag(lCL+lV+lCLD+lVT~1), start=c(lCL=1.3,lV=4.2,lCLD=1.6,lVT=4))

fitODE <- nlme_ode(dat, model=ode2, par_model=specs6, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 2 compartment IV Infusion
##
## Note that infusions need to be coded both with a record to start infusions
## and another one to stop infusions!
## For infusions into compartment 1 (for ODEs) EVID needs to be 10101
##
################################################################################################
################################################################################################

datr <- read.csv("INFUSION_2CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,10101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]
datIV<-datr[AMT>0][,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
datr<-rbind(datr,datIV)
setkey(datr,ID,TIME)

specs6<-list(fixed=lCL+lV+lCLD+lVT~1, random = pdDiag(lCL+lV+lCLD+lVT~1), start=c(lCL=1.6,lV=4.5,lCLD=1.5,lVT=3.9))


################################################################################################
## 2 compartment IV Infusion: single dose, closed form solution
################################################################################################

runno<-"N046"
dat<-datr[SD==1]
fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)),
       control = nlmeControl(pnlsTol = .1, msVerbose = TRUE, maxIter=200))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion: single dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2, par_model=specs6, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion: multiple dose, closed form solution
################################################################################################

runno<-"N047"
dat<-datr[SD==0]
fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion: multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2, par_model=specs6, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion: single and multiple dose, closed form solution
################################################################################################

runno<-"N049"
dat<-datr
fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion: single and multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2, par_model=specs6, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion: steady state data implemented using 7 doses, closed form solution
################################################################################################

runno<-"N048"
datr <- read.csv("INFUSION_2CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,10101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]
datSSobs<-datr[SS==0]
datSD<-datr[SS==1]
setkey(datSD,ID,TIME)
#general solution to allow different times of SS dose and different II values per subject:
datSSD<-datSD[,.(ID,TIME,II)]

nvar<-function(TIME,II){
 for(i in seq(1,7)){
   r=TIME-(i-1)*II
   assign(paste("ret",i,sep=""),r)
 }
 return(list(r1 = ret1, r2 = ret2,r3=ret3,r4=ret4,r5=ret5,r6=ret6,r7=ret7))
}

#updates datSSD with 7 columns to account for the new dosing times
datSSD[,(paste("V",seq(1,7)-1,sep="")):=nvar(TIME,II)][,TIME:=NULL][,II:=NULL]

index<-melt(datSSD,id.vars=c("ID"),value.name="TIMED")
index[,variable:=NULL]
index<-index[TIMED>0]
setkey(index,ID,TIMED)

datSD2<-merge(datSD,index,by=c("ID"),all=TRUE)
datSD2[,TIME:=TIMED][,TIMED:=NULL]

datSDoff<-copy(datSD2)
datSDoff[,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
datSD2<-rbind(datSD2,datSDoff)

datSS<-rbind(datSSobs,datSD2)
setkey(datSS,ID,TIME)

dat<-datSS
fit <- nlme_lin_cmpt(dat, par_model=specs6, ncmt=2, verbose=TRUE,oral=FALSE,infusion=TRUE,weight=varPower(fixed=c(1)),
       control = nlmeControl(pnlsTol = .3, msVerbose = TRUE, maxIter=200))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion: steady state data implemented using 7 doses, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2, par_model=specs6, par_trans=mypar6, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 2 compartment IV Bolus, Michaelis-Menten elimination
##
## For a bolus into compartment 1 (for ODEs) EVID needs to be 101
##
################################################################################################
################################################################################################

datr <- read.csv("BOLUS_2CPTMM.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

ode2MM <- "
d/dt(centr)  = K21*periph-K12*centr-(VM*centr/V)/(KM+centr/V);
d/dt(periph) =-K21*periph+K12*centr;
"

mypar7 <- function(lVM, lKM, lV, lCLD, lVT )
{
    VM = exp(lVM)
    KM = exp(lKM)
    V  = exp(lV)
    CLD= exp(lCLD)
    VT = exp(lVT)
    K12= CLD/V
    K21= CLD/VT
}
specs7 <- list(fixed=lVM+lKM+lV+lCLD+lVT~1, random = pdDiag(lVM+lKM+lV+lCLD+lVT~1),
          start=c(lVM=7, lKM=6, lV=4, lCLD=1.5, lVT=4))

################################################################################################
## 2 compartment IV Bolus, Michaelis-Menten elimination: single dose
################################################################################################

runno<-"N040"
dat<-datr[SD==1]
fit <- nlme_ode(dat, model=ode2MM, par_model=specs7, par_trans=mypar7, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus, Michaelis-Menten elimination: multiple dose
################################################################################################

runno<-"N041"
dat<-datr[SD==0]
fitODE <- nlme_ode(dat, model=ode2MM, par_model=specs7, par_trans=mypar7, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Bolus, Michaelis-Menten elimination: single and multiple dose
################################################################################################

runno<-"N042"
dat<-datr
fitODE <- nlme_ode(dat, model=ode2MM, par_model=specs7, par_trans=mypar7, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## 2 compartment IV Infusion, Michaelis-Menten elimination
##
## Note that infusions need to be coded both with a record to start infusions
## and another one to stop infusions!
## For infusions into compartment 1 (for ODEs) EVID needs to be 10101
##
################################################################################################
################################################################################################

datr <- read.csv("INFUSION_2CPTMM.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,10101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]
datIV<-datr[AMT>0][,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
datr<-rbind(datr,datIV)
setkey(datr,ID,TIME)


################################################################################################
## 2 compartment IV Infusion, Michaelis-Menten elimination: single dose
################################################################################################

runno<-"N054"
dat<-datr[SD==1]
fitODE <- nlme_ode(dat, model=ode2MM, par_model=specs7, par_trans=mypar7, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion, Michaelis-Menten elimination: multiple dose
################################################################################################

specs7i <- list(fixed=lVM+lKM+lV+lCLD+lVT~1, random = pdDiag(value=diag(c(3,3,3,3,3)), form = lVM+lKM+lV+lCLD+lVT~1),
           start=c(lVM=7, lKM=5, lV=4, lCLD=1.2, lVT=4))

runno<-"N055"
dat<-datr[SD==0]
fitODE <- nlme_ode(dat, model=ode2MM, par_model=specs7i, par_trans=mypar7, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## 2 compartment IV Infusion, Michaelis-Menten elimination: single and multiple dose
################################################################################################

runno<-"N056"
dat<-datr
fitODE <- nlme_ode(dat, model=ode2MM, par_model=specs7, par_trans=mypar7, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## first order absorption, 2 compartment distribution, linear elimination
##
## For a bolus into compartment 1 (for ODEs) EVID needs to be 101
##
################################################################################################
################################################################################################

datr <- read.csv("ORAL_2CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

specs8<-list(fixed=lCL+lV+lCLD+lVT+lKA~1, random = pdDiag(lCL+lV+lCLD+lVT+lKA~1), start=c(lCL=1.6,lV=4.5,lCLD=1.5,lVT=3.9,lKA=0.1))

################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## single dose, closed form solution
################################################################################################

runno<-"N060"
dat<-datr[SD==1]
fit <- nlme_lin_cmpt(dat, par_model=specs8, ncmt=2, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)),
       control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## single dose, ODE solution
################################################################################################

ode2KA <- "
d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs+K21*periph-K12*centr-K10*centr;
d/dt(periph) =        -K21*periph+K12*centr;
"

mypar8 <- function(lCL, lV, lCLD, lVT,lKA)
{
    CL = exp(lCL)
    V  = exp(lV)
    CLD= exp(lCLD)
    VT = exp(lVT)
    KA = exp(lKA)
    K10= CL/V
    K12= CLD/V
    K21= CLD/VT
}

fitODE <- nlme_ode(dat, model=ode2KA, par_model=specs8, par_trans=mypar8, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .3, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## multiple dose, closed form solution
################################################################################################

runno<-"N061"
dat<-datr[SD==0]
fit <- nlme_lin_cmpt(dat, par_model=specs8, ncmt=2, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)),
       control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2KA, par_model=specs8, par_trans=mypar8, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## single and multiple dose, closed form solution
################################################################################################

runno<-"N063"
dat<-datr
fit <- nlme_lin_cmpt(dat, par_model=specs8, ncmt=2, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)),
       control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## single and multiple dose, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2KA, par_model=specs8, par_trans=mypar8, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## steady state data implemented using 7 doses, closed form solution
################################################################################################

runno<-"N062"
datr <- read.csv("ORAL_2CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

datSS<-datr[SS==0]
datSD<-datr[SS==1]

#general solution to allow different times of SS dose and different II values per subject:
datSSD<-datr[SS==1,.(ID,TIME,II)]

nvar<-function(TIME,II){
 for(i in seq(1,7)){
   r=TIME-(i-1)*II
   assign(paste("ret",i,sep=""),r)
 }
 return(list(r1 = ret1, r2 = ret2,r3=ret3,r4=ret4,r5=ret5,r6=ret6,r7=ret7))
}

#updates datSSD with 7 columns to account for the new dosing times
datSSD[,(paste("V",seq(1,7)-1,sep="")):=nvar(TIME,II)][,TIME:=NULL][,II:=NULL]
index<-melt(datSSD,id.vars="ID",value.name="TIMED")
index[,variable:=NULL]
index<-index[TIMED>0]

datSD2<-merge(datSD,index,by=c("ID"),all=TRUE)
datSD2[,TIME:=TIMED][,TIMED:=NULL]
datSS<-rbind(datSS,datSD2)
setkey(datSS,ID,TIME)
dat<-datSS

specs8i<-list(fixed=lCL+lV+lCLD+lVT+lKA~1, random = pdDiag(value=diag(c(6,6,6,6,6)),
         form = lCL+lV+lCLD+lVT+lKA~1), start=c(lCL=1.4,lV=4.2,lCLD=1.3,lVT=3.9,lKA=0.1))

fit <- nlme_lin_cmpt(dat, par_model=specs8i, ncmt=2, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)),
       control = nlmeControl(pnlsTol = .15, msVerbose = TRUE))
summary(fit)
assign(paste("fit",runno,sep=""),fit, pos = .GlobalEnv)
save(fit,file=paste("fit",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, linear elimination:
## steady state data implemented using 7 doses, ODE solution
################################################################################################

fitODE <- nlme_ode(dat, model=ode2KA, par_model=specs8i, par_trans=mypar8, response="centr", response.scaler="V",
          verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .15, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
################################################################################################
##
## first order absorption, 2 compartment distribution, Michaelis-Menten elimination
##
## For a bolus into compartment 1 (for ODEs) EVID needs to be 101
##
################################################################################################
################################################################################################

datr <- read.csv("ORAL_2CPTMM.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

ode2MMKA <- "
d/dt(abs)    =-KA*abs;
d/dt(centr)  = KA*abs+K21*periph-K12*centr-(VM*centr/V)/(KM+centr/V);
d/dt(periph) =-K21*periph+K12*centr;
"

mypar9 <- function(lVM, lKM, lV, lCLD, lVT, lKA )
{
    VM = exp(lVM)
    KM = exp(lKM)
    V  = exp(lV)
    CLD= exp(lCLD)
    VT = exp(lVT)
    KA = exp(lKA)
    K12= CLD/V
    K21= CLD/VT
}
specs9 <- list(fixed=lVM+lKM+lV+lCLD+lVT+lKA~1, random = pdDiag(lVM+lKM+lV+lCLD+lVT+lKA~1),
          start=c(lVM=7, lKM=6, lV=4, lCLD=1.5, lVT=4,lKA=0.1))


################################################################################################
## first order absorption, 2 compartment distribution, Michaelis-Menten elimination: single dose
################################################################################################

runno<-"N068"
dat<-datr[SD==1]
fitODE <- nlme_ode(dat, model=ode2MMKA, par_model=specs9, par_trans=mypar9, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, Michaelis-Menten elimination:
## multiple dose
################################################################################################

runno<-"N069"
dat<-datr[SD==0]
fitODE <- nlme_ode(dat, model=ode2MMKA, par_model=specs9, par_trans=mypar9, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))


################################################################################################
## first order absorption, 2 compartment distribution, Michaelis-Menten elimination:
## single and multiple dose
################################################################################################

runno<-"N070"
dat<-datr
fitODE <- nlme_ode(dat, model=ode2MMKA, par_model=specs9, par_trans=mypar9, response="centr", response.scaler="V",
       verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))
summary(fitODE)
assign(paste("fitODE",runno,sep=""),fitODE, pos = .GlobalEnv)
save(fitODE,file=paste("fitODE",runno,".Rdata",sep=""))

