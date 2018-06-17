## vpc.R: population PK/PD modeling library
##
## Copyright (C) 2014 - 2016  Wenping Wang
##
## This file is part of nlmixr.
##
## nlmixr is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## nlmixr is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with nlmixr.  If not, see <http:##www.gnu.org/licenses/>.


sim.one = function(zz, x) {
    ..ModList <- nlmeModList();
       nsub=length(unique(x$groups[[1]]))
	om = as.matrix(x$modelStruct$reStruct$ID) * x$sigma^2
	eta = multi2(rep(0, dim(om)[1]), om, nsub)
	th = x$coefficients$fixed

	m = sapply(seq(length(th)), function(k) {
		ix = match(names(th)[k], dimnames(eta)[[1]], nomatch=0)
		if (ix) th[k]+eta[ix,]
		else rep(th[k], nsub)
	})
	dimnames(m)[[2]] = names(th)

	r = as.integer(as.character(x$groups[[1]]))
	m = m[r,]
	m = cbind(m, TIME=..ModList$dat.g$TIME, ID=..ModList$dat.g$ID)
	m = as.data.frame(m)
    res = do.call(..ModList$user_fn, m)
	res+rnorm(res, 0, x$sigma)
}


##' Vpc function for nlmixr
##'
##' @param sim Observed data frame or fit object
##' @param ... Other parameters
##'
##' @export
vpc <- function (sim, ...)
{
    UseMethod("vpc")
}
##' @rdname vpc
##' @export
vpc.default <- function(sim, ...){
    ns <- loadNamespace("vpc");
    if (exists("vpc_vpc",ns)){
        vpcn <- "vpc_vpc"
    } else {
        vpcn <- "vpc"
    }
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    call <- call[names(call) %in% methods::formalArgs(getFromNamespace(vpcn,"vpc"))]
    p = do.call(getFromNamespace(vpcn,"vpc"), call, envir = parent.frame(1))
}

#' Visual predictive check (VPC) for nlmixr nlme objects
#'
#' Do visual predictive check (VPC) plots for nlme-based non-linear mixed effect models
#'
#' @param fit nlme fit object
#' @param nsim number of simulations
#' @param condition conditional variable
#' @param ... Additional arguments
#' @inheritParams vpc::vpc
#' @return NULL
#' @examples
#' specs <- list(fixed=lKA+lCL+lV~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lCL=-3.2, lV=-1))
#' fit <- nlme_lin_cmpt(theo_md, par_model=specs, ncmt=1, verbose=TRUE)
#' vpc(fit, nsim = 100, condition = NULL)
#' @export
vpc_nlmixr_nlme = function(fit, nsim=100, condition=NULL, ...)
{
    nlmeModList(fit$env);
    on.exit({nlmeModList(new.env(parent=emptyenv()))})
    ..ModList <- nlmeModList();

	options(warn=-1)

    s = sapply(1:nsim, sim.one, x=fit)

    cond.var = if(is.null(condition)) rep(1, dim(..ModList$dat.g)[1]) else ..ModList$dat.g[, condition]
	levels = sort(unique(cond.var))
	for (k in 1:length(levels)) {
		sel = cond.var == levels[k]
		xs = s[sel, ]
		xd = ..ModList$dat.g[sel, ]
		matplot(xd$TIME, xs, col="#33FF66", pch=19, xlab="TIME", ylab="DV")
		points(xd$TIME, xd$DV, col="#000066")

		if(!is.null(condition)) {
			title(paste0(condition, ": ", levels[k]))
		}
	}

	options(warn=0)
	invisible(NULL)
}

##' @rdname vpc_nlmixr_nlme
##' @export
vpcNlmixrNlme <- vpc_nlmixr_nlme

#' @rdname vpc_nlmixr_nlme
#' @export
vpc.nlmixrNlme <- function(sim, ...){
    vpc_nlmixr_nlme(sim, ...);
}

#vpc(fit, 100)

multi2 <- function (mu, vmat, n)
{
    eta <- matrix(rnorm(length(mu) * n), ncol = n, nrow = length(mu))
    Q <- chol(vmat, pivot = TRUE)
    pivot <- attr(Q, "pivot")
    oo <- order(pivot)
    para <- t(Q[, oo]) %*% eta
    sweep(para, 1, mu, "+")
}


#' Bootstrap data
#'
#' Bootstrap data by sampling the same number of subjects from the original dataset by sampling with replacement.
#'
#' @param dat model data to be bootstrapped
#' @return Bootstrapped data
#' @examples
#'
#' specs <- list(fixed=lKA+lCL+lV~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lCL=-3.2, lV=-1))
#' set.seed(99); nboot = 5;
#'
#' cat("generating", nboot, "bootstrap samples...\n")
#' cmat <- matrix(NA, nboot, 3)
#' for (i in 1:nboot)
#' {
#' 	#print(i)
#' 	bd <- bootdata(theo_md)
#' 	fit <- nlme_lin_cmpt(bd, par_model=specs, ncmt=1)
#' 	cmat[i,] = fit$coefficients$fixed
#' }
#' dimnames(cmat)[[2]] <- names(fit$coefficients$fixed)
#' print(head(cmat))
#'
#' @export
bootdata = function(dat)
{
    id = unique(dat$ID)
	nsub = length(id)
	s = sample(id, nsub, replace=TRUE)

	do.call("rbind",
		lapply(1:nsub, function(ix)
		{
			k = s[ix]
			d = dat[dat$ID == k, ];
			d$ID = ix
			d
		})
	)
}


#' Forward covariate selection for nlme-base non-linear mixed effect models
#'
#' Implements forward covariate selection for nlme-based non-linear mixed effect models
#'
#' @param base base model
#' @param cv a list of candidate covariate to model parameters
#' @param dat model data
#' @param cutoff significance level
#' @return an nlme object of the final model
#' @examples
#' dat <- theo_md
#' dat$LOGWT <- log(dat$WT)
#' dat$TG <- (dat$ID < 6) + 0    #dummy covariate
#'
#' specs <- list(
#' 	fixed=list(lKA=lKA~1, lCL=lCL~1, lV=lV~1),
#' 	random = pdDiag(lKA+lCL~1),
#' 	start=c(0.5, -3.2, -1))
#' fit0 <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1)
#' cv <- list(lCL=c("WT", "TG"), lV=c("WT"))
#' fit <- frwd_selection(fit0, cv, dat)
#' print(summary(fit))
#' @export
frwd_selection = function(base, cv, dat, cutoff=.05)
{
	#dat = getData(base)
	fixed.save = base$call$fixed
	start.save = as.list(base$call$start)
	names(start.save) = names(fixed.save)


	cat("covariate selection process:\n")
	while(1)
	{
		rl = NULL; pval=NULL
		for(par in names(cv))
		{
			fixed = fixed.save
			start = start.save
			for (wh in cv[[par]])
			{
				fixed[[par]] = as.formula(sprintf("%s+%s", deparse(fixed.save[[par]]), wh))
				start[[par]] = c(start.save[[par]], 0)
				specs <- list(
					fixed=fixed,
					random = pdDiag(lKA+lCL~1),
					start=unlist(start))
				fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1)
				aov = anova(base, fit)[2, "p-value"]
				cat("\nadding", wh, "to", par, ": p-val =", aov)
				pval = c(pval, aov)
				rl = c(rl, list(list(par, wh, fixed[[par]], start, fit)))
			}
		}

		if (min(pval)>cutoff) {
			cat("\n\ncovariate selection finished.\n\n\n")
			break
		}

		wh = match(min(pval), pval)
		s = rl[[wh]]
		ix = match(s[[2]], cv[[s[[1]]]])

		cat("\n", cv[[s[[1]]]][ix], "added to", s[[1]], "\n")

		cv[[s[[1]]]] = cv[[s[[1]]]][-ix]
		if (length(cv[[s[[1]]]])==0) cv[[s[[1]]]] = NULL
		fixed.save[[s[[1]]]] = s[[3]]
		start.save[[s[[1]]]] = s[[4]][[s[[1]]]]
		base = s[[5]]
	}

	base
}


sim.one = function (zz, x)
{
    ..ModList = nlmeModList();
    nsub = length(unique(x$groups[[1]]))
    om = as.matrix(x$modelStruct$reStruct$ID) * x$sigma^2
    eta = t(multi2(rep(0, dim(om)[1]), om, nsub))
    dimnames(eta)[[1]] = dimnames(x$coefficients$random$ID)[[1]]
    x$coefficients$random$ID = eta
    pred = predict(x, getData(x))
    if (is.null(x$call$weights)){
        sd = 1
    }
    else if (is(x$call$weights, "varPower")) {
        sd = abs(pred)^as.double(coef(x$modelStruct$varStruct, allCoef=TRUE))
    }
    else if (class(x$call$weights)[1] == "varConstPower") {
        sd = exp(x$modelStruct$varStruct$const) + abs(pred)^x$modelStruct$varStruct$power
    } else {
        stop("residual model not implemented")
    }
    ##pred + rnorm(pred, 0, x$sigma)
    pred + rnorm(pred, 0, sd*x$sigma)
}


