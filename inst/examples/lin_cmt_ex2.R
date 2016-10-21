require(nlmixr)
require(nlme)

dat <- read.table("theo_inf.txt", head=TRUE)
mypar <- function(lKE, lCL)
{
     KE=exp(lKE)
     CL=exp(lCL)
     V = CL/KE
}

specs <- list(fixed=lKE+lCL~1, random = pdDiag(lCL~1), start=c(lKE=-2.5, lCL=-3.2))
fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1, oral=FALSE, infusion=TRUE, parameterization=2, par_trans=mypar)
plot(augPred(fit,level=0:1))
summary(fit)

