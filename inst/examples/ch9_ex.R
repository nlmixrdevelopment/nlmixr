library(nlmixr)
library(RxODE)

ode <- "
Cp = centr/Vp;
Cm = meta/Vm;
d/dt(centr) = -k12*centr + k21*peri -kp*centr;
d/dt(peri)  =  k12*centr - k21*peri;
d/dt(meta)  =                        kp*centr - km*meta;
"
sys1 <- RxODE(model = ode, modName = "s1")

dat = read.table("metabolite.txt", header=TRUE)
mod = list(y1 ~ Cp+prop(.1), y2 ~ Cm+prop(.15))
#mod = list(y1 ~ Cp+add(.2)+prop(.1), y2 ~ Cm+prop(.1))
#mod = list(y1 ~ Cp+prop(.1)+add(.2), y2 ~ Cm+prop(.1))
#mod = list(y1 ~ Cp+add(.1), y2 ~ Cm+prop(.1))
inits = c(kp=0.4, Vp=10., k12=0.2, k21=0.1, km=0.2, Vm=30.)
ev = eventTable()
ev$add.dosing(100, rate=100)
ev$add.sampling(c(0, dat$time))
#source("utils.R")
#(fit = dynmodel(sys1, mod, ev, inits, dat, method="PORT"))
(fit = dynmodel(sys1, mod, ev, inits, dat))
plot(fit)


#----------------
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
sys2 <- RxODE(model = ode, modName = "s2")

dat = read.table("invgaussian.txt", header=TRUE)
mod = cp ~ C2 + prop(.1)
inits = c(MIT=190, CVI2=.65, F=.92)
fixPars = c(CL=.0793, V2=.64, Q=.292, V3=9.63)
ev = eventTable()
ev$add.sampling(c(0, dat$time))
#source("utils.R")
#(fit = dynmodel(sys2, mod, ev, inits, dat, fixPars, method="PORT"))
(fit = dynmodel(sys2, mod, ev, inits, dat, fixPars))
plot(fit)

