require(nlmixr)
ode = "d/dt(depot) =-KA*depot;\nd/dt(centr) = KA*depot - KE*centr;"
PKpars = function(lCL, lV, lKA)
{
  CL = exp(lCL)
  V = exp(lV)
  KA = exp(lKA)
  KE = CL / V
}

#--- gen user fn
 saem_fit <- gen_saem_user_fn(model=lincmt(ncmt=1, oral=T))  
#saem_fit <- gen_saem_user_fn(model=ode, PKpars, depvar="centr", scaler="V")

#--- saem cfg
nmdat = read.table("theo_sd.dat",  head=T)
model = list(N.eta=3, omega=diag(3)*6, covars="WT", res.mod=1)
inits = list(theta=c(.05, .5, 2), ares=4, bres=1)
cfg   = configsaem(model, nmdat, inits)

fit = saem_fit(cfg)
df = simple.gof(fit)
xyplot(DV~TIME|ID, df, type=c("p", "l"), lwd=c(NA,1), pch = c(1, NA), groups=grp)
diag(solve(fit$Ha))
