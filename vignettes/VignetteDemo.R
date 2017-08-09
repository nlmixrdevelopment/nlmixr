library(nlmixr)

## -- nlme examples
dat <- read.table(system.file("examples/theo_md.txt", package = "nlmixr"), head=TRUE)

specs <- list(fixed=lKA+lCL+lV~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lCL=-3.2, lV=-1))

fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1)
summary(fit)
plot(augPred(fit,level=0:1))

## ---------------------
specs <- list(
	fixed=list(lKA~1, lCL+lV~WT),
	random = pdDiag(lKA+lCL~1),
    start=c(0.5, -3.2, 0, -1, 0))
fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1)
#plot(augPred(fit,level=0:1))
#fit


## ---------------------
mypar <- function(lKA, lKE, lCL)
{
    KA <- exp(lKA)
    KE <- exp(lKE)
    CL <- exp(lCL)
    V  <- CL/KE
}
specs <- list(
    fixed=lKA+lCL+lKE~1,
	random = pdDiag(lKA+lCL~1),
	start=c(0.5, -2.5, -3.2)
)
fit <- nlme_lin_cmpt(
	dat, par_model=specs,
	ncmt=1, parameterization=2, par_trans=mypar)
#plot(augPred(fit,level=0:1))
fit


## ------------------------------------------------------------------------
ode <- "
d/dt(depot) =-KA*depot;
d/dt(centr) = KA*depot - KE*centr;
"
m1 = RxODE(ode)

dat <- read.table("theo_md.txt", head=TRUE)
dat$WG <- dat$WT>70
mypar <- function(lKA, lKE, lCL)
{
    KA <- exp(lKA)
    KE <- exp(lKE)
    CL <- exp(lCL)
    V  <- CL/KE
}
specs <- list(fixed=lKA+lKE+lCL~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lKE=-2.5, lCL=-3.2))
fit <- nlme_ode(dat, model=m1, par_model=specs, par_trans=mypar, response="centr", response.scaler="V")
nlme_gof(fit)
fit


## ------------------------------------------------------------------------
## -- vpc examples
vpc(fit, 100)
par(mfrow=c(1,2))
vpc(fit, 100, condition="WG")


## ------------------------------------------------------------------------
## -- boostrap example
dat <- read.table("theo_md.txt", head=TRUE)
specs <- list(fixed=lKA+lCL+lV~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lCL=-3.2, lV=-1))
set.seed(99); nboot = 20;

cat("generating", nboot, "bootstrap samples...\n")
cmat <- matrix(NA, nboot, 3)
for (i in 1:nboot)
{
	#print(i)
	bd <- bootdata(dat)
	fit <- nlme_lin_cmpt(bd, par_model=specs, ncmt=1)
	cmat[i,] = fit$coefficients$fixed
}
dimnames(cmat)[[2]] <- names(fit$coefficients$fixed)
print(head(cmat))

require(lattice)
df <- do.call("make.groups", split(cmat, col(cmat)))
df$grp <- dimnames(cmat)[[2]][df$which]
print(bwplot(grp~exp(data), df))


## ------------------------------------------------------------------------
## -- covariate search examples
dat <- read.table("theo_md.txt", head=TRUE)
dat$LOGWT <- log(dat$WT)
dat$TG <- (dat$ID < 6) + 0    #dummy covariate

specs <- list(
	fixed=list(lKA=lKA~1, lCL=lCL~1, lV=lV~1),
	random = pdDiag(lKA+lCL~1),
	start=c(0.5, -3.2, -1))
fit0 <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1)
cv <- list(lCL=c("WT", "TG", "LOGWT"), lV=c("WT", "TG", "LOGWT"))
fit <- frwd_selection(fit0, cv, dat)
print(summary(fit))



## ------------------------------------------------------------------------
## -- SAEM

uif <- function(){
    ini({
        lCl <- 0.05
        lV <- 0.2
        lKA = 0.02
    })
    model({
        CL = exp(lCl)
        V = exp(lV)
        KA = exp(lKA)
        KE = CL / V
    })
}

fit <- nlmixr(uif, theo_sd, est="saem")

PKpars = function()
{
  CL = exp(lCL)
  V = exp(lV)
  KA = exp(lKA)
  KE = CL / V
  xxx = 0;
                                        #initCondition = c(0,xxx)
}

#closed form solution:
saem_fit <- gen_saem_user_fn(model=lincmt(ncmt=1, oral=T))

##--- saem cfg
nmdat = read.table("theo_sd.dat",  head=T)
model = list(saem_mod=saem_fit, covars="WT")
inits = list(theta=c(.05, .5, 2))
cfg   = configsaem(model, nmdat, inits)

#cfg$Rfn = nlmixr:::Ruser_function_cmt
#dyn.load("m1.d/m1.dll");cfg$Rfn = nlmixr:::Ruser_function_ode
fit = saem_fit(cfg)
fit

df = plot(fit)


#ODE implementation:
ode <-
    "d/dt(depot) =-KA*depot;
     d/dt(centr) = KA*depot - KE*centr;"
PRED = function() centr / V
m1 = RxODE(ode, modName="m1")

saem_fit <- gen_saem_user_fn(model=m1, PKpars, pred=PRED)

#--- saem cfg
nmdat = read.table("theo_sd.dat",  head=T)
model = list(saem_mod=saem_fit, covars="WT")
inits = list(theta=c(.05, .5, 2))
#change additional settings for the SAEM algorithm:
cfg   = configsaem(model, nmdat, inits, mcmc = list(niter = c(50, 100), nmc = 3, nu =  c(2, 2, 2)),seed=17624764)


#cfg$Rfn = nlmixr:::Ruser_function_cmt
#dyn.load("m1.d/m1.dll");cfg$Rfn = nlmixr:::Ruser_function_ode
fit = saem_fit(cfg)
fit
df = plot(fit)


#Calculate C2 inside the ODE:
ode <-
  "C2 = centr/V;
  d/dt(depot) =-KA*depot;
  d/dt(centr) = KA*depot - KE*centr;"
m2 = RxODE(ode, modName="m2")
PRED2 = function() C2

saem_fit <- gen_saem_user_fn(model=m2, PKpars, pred=PRED2)


#--- saem cfg
nmdat = read.table("theo_sd.dat",  head=T)
model = list(saem_mod=saem_fit, covars="WT")
inits = list(theta=c(.05, .5, 2))
cfg   = configsaem(model, nmdat, inits)
#only print every 50th iteration:
cfg$print = 50

#cfg$Rfn = nlmixr:::Ruser_function_cmt
#dyn.load("m1.d/m1.dll");cfg$Rfn = nlmixr:::Ruser_function_ode
fit = saem_fit(cfg)
fit
df = plot(fit)




## ------------------------------------------------------------------------
## gnlmm examples
llik <- function()
{
	if (group==1) lp = THETA[1]+THETA[2]*logtstd+ETA[1]
	else          lp = THETA[3]+THETA[4]*logtstd+ETA[1]
	lam = exp(lp)
	dpois(y, lam, log=TRUE)
}
inits = list(THTA=c(1,1,1,1), OMGA=list(ETA[1]~1))

fit = gnlmm(llik, pump, inits,
	control=list(
	    reltol.outer=1e-4,
		optim.outer="nmsimplex",
            nAQD=5 ))

cv = calcCov(fit)
Rinv = attr(cv,"RinvS")$Rinv
S    = attr(cv,"RinvS")$S
Rinv*2					#inverse hessian matrix
solve(S)*4			    #inverse of score function product sum
Rinv %*% S %*% Rinv		#sandwich estimate


## ------------------------------------------------------------------------
llik <- function()
{
	lp = THETA[1]*x1+THETA[2]*x2+(x1+x2*THETA[3])*ETA[1]
	p = pnorm(lp)
	dbinom(x, m, p, log=TRUE)
}
inits = list(THTA=c(1,1,1), OMGA=list(ETA[1]~1))

gnlmm(llik, rats, inits, control=list(nAQD=7))


## ------------------------------------------------------------------------
ode <- "
d/dt(depot) =-KA*depot;
d/dt(centr) = KA*depot - KE*centr;
"

#same model as the one run above for SAEM generated using
# m1 = RxODE(ode, modName="m1")
#and so if run above, will not need to be run again,
# but if not, do run the line below:
#m1 = RxODE(ode)

pars <- function()
{
	CL = exp(THETA[1] + ETA[1])#; if (CL>100) CL=100
	KA = exp(THETA[2] + ETA[2])#; if (KA>20) KA=20
	KE = exp(THETA[3])
	V  = CL/KE
	sig2 = exp(THETA[4])
}
llik <- function() {
	pred = centr/V
	dnorm(DV, pred, sd=sqrt(sig2), log=TRUE)
}
inits = list(THTA=c(-3.22, 0.47, -2.45, 0))
inits$OMGA=list(ETA[1]~.027, ETA[2]~.37)
#inits$OMGA=list(ETA[1]+ETA[2]~c(.027, .01, .37))
theo <- read.table("theo_md.txt", head=TRUE)

fit = gnlmm(llik, theo, inits, pars, m1,
	control=list(trace=TRUE, nAQD=5))

pred = function() centr/V
s = prediction(fit, pred)
plot(s$p, s$dv); abline(0,1,col="red")



## ------------------------------------------------------------------------
#-- dynmodel examples
ode <- "
   dose=200;
   pi = 3.1415926535897931;

   if (t<=0) {
      fI = 0;
   } else {
      fI = F*dose*sqrt(MIT/(2.0*pi*CVI2*t^3))*exp(-(t-MIT)^2/(2.0*CVI2*MIT*t));
   }

   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(centr) = fI - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =              Q*C2 - Q*C3;
"
sys1 <- RxODE(model = ode)


## ------------------------------------------------------------------------
dat <- read.table("invgaussian.txt", header=TRUE)
mod <- cp ~ C2 + prop(.1)
inits <- c(MIT=190, CVI2=.65, F=.92)
fixPars <- c(CL=.0793, V2=.64, Q=.292, V3=9.63)
ev <- eventTable()
ev$add.sampling(c(0, dat$time))
(fit <- dynmodel(sys1, mod, ev, inits, dat, fixPars))
summary(fit)
par(mar=c(4,4,1,1), mfrow=c(1,3))
plot(fit, cex=2)


## ------------------------------------------------------------------------
ode <- "
Cp = centr/Vp;
Cm = meta/Vm;
d/dt(centr) = -k12*centr + k21*peri -kp*centr;
d/dt(peri)  =  k12*centr - k21*peri;
d/dt(meta)  =                        kp*centr - km*meta;
"
sys2 <- RxODE(model = ode)

dat <- read.table("metabolite.txt", header=TRUE)
mod <- list(y1 ~ Cp+prop(.1), y2 ~ Cm+prop(.15))
inits <- c(kp=0.4, Vp=10., k12=0.2, k21=0.1, km=0.2, Vm=30.)
ev <- eventTable()
ev$add.dosing(100, rate=100)
ev$add.sampling(c(0, dat$time))
(fit <- dynmodel(sys2, mod, ev, inits, dat))

## ----------------
mod <- list(y1 ~ Cp+add(.2)+prop(.1), y2 ~ Cm+prop(.15))
(fit <- dynmodel(sys2, mod, ev, inits, dat))

## ----------------
mod <- list(y1 ~ Cp+prop(.1), y2 ~ Cm+prop(.15))
(fit <- dynmodel.mcmc(sys2, mod, ev, inits, dat))

par(mfrow=c(4,2), mar=c(2,4,1,1))
s <- lapply(1:dim(fit)[2], function(k)
     plot(fit[,k], type="l", col="red", ylab=dimnames(fit)[[2]][k]))

s
