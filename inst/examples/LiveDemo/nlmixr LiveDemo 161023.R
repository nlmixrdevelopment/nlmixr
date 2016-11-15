#Install Rtools first: https://cran.r-project.org/bin/windows/Rtools/Rtools34.exe
#go to your r installation directory (e.g. C:\Program Files\R\R-3.3.1) and change 
#the properties/security/change permissions/ and give full control to 
#"Users (YourComputerName/Users)"
 
#then run:
#install.packages("devtools")
#library(devtools)
#install_github('nlmixrdevelopment/RxODE')
#install_github('nlmixrdevelopment/nlmixr')
#install.packages("data.table")

### remove all pre-existing information
rm(list=ls())

library(nlmixr, quietly = TRUE)
library(data.table)
source("print.summary.lme.R")

datr <- read.csv("ORAL_1CPT.csv", header=TRUE,stringsAsFactors=F)
datr$EVID<-ifelse(datr$EVID==1,101,datr$EVID)
datr<-data.table(datr)
datr<-datr[EVID!=2]

mypar4 <- function(lCL, lV, lKA )
{   CL = exp(lCL) 
    V  = exp(lV)
    KA = exp(lKA)
}

specs4 <- list(fixed=lCL+lV+lKA~1, random = pdDiag(lCL+lV+lKA~1), start=c(lCL=1,lV=4,lKA=0))

fit <- nlme_lin_cmpt(datr, par_model=specs4, ncmt=1, verbose=TRUE,oral=TRUE,weight=varPower(fixed=c(1)))

summary(fit)
intervals(fit)
#coef(fit)



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


datrMM <- read.csv("Infusion_1CPTMM.csv", header=TRUE,stringsAsFactors=F)
datrMM$EVID<-ifelse(datrMM$EVID==1,10101,datrMM$EVID)
datrMM<-data.table(datrMM)
datrMM<-datrMM[EVID!=2]
datIV<-datrMM[AMT>0][,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
datrMM<-rbind(datrMM,datIV)
setkey(datrMM,ID,TIME)

fit <- nlme_ode(datrMM, model=ode1MM, par_model=specs3, par_trans=mypar3, response="centr", response.scaler="V", 
verbose=TRUE,weight=varPower(fixed=c(1)),control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))

summary(fit)
intervals(fit)




#SAEM instead of nlme... 

#use when you have installed the MKL library with Microsoft R Open:
#setMKLthreads(1)

source("summary.saemFit.R")

#if you have run the saem_fit command during an R-session you will need to 
#unload the saem dll prior to a new compile:

#dyn.unload("saem_main.dll")

saem_fit <- gen_saem_user_fn(model=lincmt(ncmt=1, oral=TRUE))

#temporary work around for specifying covariates
datr$WT<-1
model = list(saem_mod=saem_fit, res.mod=2,covars="WT")
inits = list(theta=c(5,90,1),omega=c(0.1,0.1,0.1),bres=0.2)
cfg   = configsaem(model, datr, inits)
cfg$print = 50
fit = saem_fit(cfg)

summary(fit)





m3 = RxODE(ode1MM, modName="m3")

PRED = function() centr / V

dyn.unload("saem_main.dll")
saem_fit <- gen_saem_user_fn(model=m3, PKpars=mypar3, pred=PRED)


datrMM$WT<-1
model = list(saem_mod=saem_fit, res.mod=2,covars="WT")
inits = list(theta=c(1000,250,70),omega=c(0.1,0.1,0.1),bres=0.2)
cfg   = configsaem(model, datrMM, inits, mcmc = list(niter = c(50, 100), nmc = 3, nu =  c(2, 2, 2)),seed=17624764)
cfg$print = 10

fit = saem_fit(cfg)

summary(fit)

