require(nlmixr)
require(nlme)

ode <- "
d/dt(centr) = -KE*centr;
"
dat <- read.table("theo_inf.txt", head=TRUE)
mypar <- function(lKE, lCL)
{
    KE=exp(lKE) 
    CL=exp(lCL)
    V = CL/KE
}

specs <- list(fixed=lKE+lCL~1, random = pdDiag(lCL~1), start=c(lKE=-2.5, lCL=-3.2))
fit <- nlme_ode(dat, model=ode, par_model=specs, par_trans=mypar, response="centr", response.scaler="V")

require(lattice)
df <- getData(fit)
df <- rbind(
	cbind(df[,c("ID", "TIME", "DV")], grp=0),
	cbind(df[,c("ID", "TIME")], DV=fit$fitted[,1], grp=1),
	cbind(df[,c("ID", "TIME")], DV=fit$fitted[,2], grp=2)
)
xyplot(DV~TIME|ID, df, group=grp, type="b", lwd=c(NA, 1, 1), pch=c(1,NA,NA))

