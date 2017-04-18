#' GOF plots for nlme::nlme-based mixed-effect models
#'
#' Generates basic goodness-of-fit plots for nlme::nlme-based mixed-effect models
#'
#' @param fit nlme::nlme fit object
#' @param ... optional additional arguments
#' @return NULL
#' @useDynLib nlmixr
#' @importFrom graphics abline lines matplot plot points title
#' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var
#' @export
nlme_gof <- function(fit, ...){
	df <- nlme::getData(fit)
	df <- rbind(
		cbind(df[,c("ID", "TIME", "DV")], grp=0),
		cbind(df[,c("ID", "TIME")], DV=fit$fitted[,1], grp=1),
		cbind(df[,c("ID", "TIME")], DV=fit$fitted[,2], grp=2)
	)

	p = lattice::xyplot(DV~TIME|ID, df, group=grp, type="b", lwd=c(NA, 1, 1), pch=c(1,NA,NA), col=lattice::trellis.par.get("superpose.line")$col[c(1,1,2)], ...)
	p
}
