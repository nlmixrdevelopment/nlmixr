require(nlmixr)
require(nlme)

dat <- read.table("theo_md.txt", head=TRUE)


#-- clearance parameterization: 1 cmt
specs <- list(fixed=lKA+lCL+lV~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lCL=-3.2, lV=-1))
fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1, verbose=TRUE)
plot(augPred(fit,level=0:1))
summary(fit)


#-- micro parameterization: 1 cmt
specs <- list(fixed=lKA+lKE+lV~1, random = pdDiag(lKA+lKE~1), start=c(lKA=0.5, lKE=-3.2, lV=-1))
fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1, parameterization=2)
plot(augPred(fit,level=0:1))
summary(fit)


#-- clearance parameterization: 2 cmt
specs <- list(fixed=lKA+lCL+lV+lCLD+lVT~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lCL=-3.2, lV=-1, lCLD=-1, lVT=2))
fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=2)
plot(augPred(fit,level=0:1))
summary(fit)


#-- micro parameterization: 2 cmt
specs <- list(fixed=lKA+lKE+lV+lK12+lK21~1, random = pdDiag(lKA+lKE~1), start=c(lKA=0.5, lKE=-2.2, lV=-1, lK12=.1, lK21=.1))
fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=2, parameterization=2)
plot(augPred(fit,level=0:1))
summary(fit)


#-- clearance parameterization: 1 cmt; covariate analysis
specs <- list(
	fixed=list(lKA~1, lCL+lV~WT), 
	random = pdDiag(lKA+lCL~1), 
	start=c(0.5, -3.2, 0, -1, 0))
fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1)
plot(augPred(fit,level=0:1))
summary(fit)


#-- non-standard parameterization & covariate analysis
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
plot(augPred(fit,level=0:1))
summary(fit)

