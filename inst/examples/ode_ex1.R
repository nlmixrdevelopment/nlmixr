require(nlmixr)
require(nlme)

ode <- "
d/dt(depot) =-KA*depot;
d/dt(centr) = KA*depot - KE*centr;
"
dat <- read.table("theo_md.txt", head=TRUE)
mypar <- function(lKA, lKE, lCL)
{
    KA <- exp(lKA) 
    KE <- exp(lKE) 
    CL <- exp(lCL)
    V  <- CL/KE
}
specs <- list(fixed=lKA+lKE+lCL~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lKE=-2.5, lCL=-3.2))
fit <- nlme_ode(dat, model=ode, par_model=specs, par_trans=mypar, response="centr", response.scaler="V")
nlme_gof(fit)
