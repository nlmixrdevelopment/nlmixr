## utils.R: population PK/PD modeling library
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


thresh = function(x, cut=.Machine$double.xmin)
{
	ifelse(x>cut, x, cut)
}

gof = function(x, ...)
{
    gof_str = "
{
    theta = th[1:npar]
    names(theta) = names(inits)[1:npar]
    theta = c(theta, fixPars)
    system$solve(theta, ev, atol = 1e-08, rtol = 1e-08)
}
"
	f = x$obj
	dat = x$data
	body(f) = parse(text=gof_str)
	x = f(x$par)

	rows = environment(f)$rows
	m = environment(f)$model
	for(s in m) {
		time = x[,"time"]
		yo = dat[, s["dv"]]		#FIXME
		yp = x[, s["pred"]]

		#dv vs pred
		plot(time, yp, type="n", xlab="time", ylab=s["dv"])
		points(dat[,"time"], yo, ...)
		lines(time, yp)

		#pred vs res
		yp = yp[rows]
		res = yo - yp
		plot(yp, yp, xlab="pred", ylab=s["dv"], ...)
		abline(0, 1, col="red", lty=2)

		#pred vs res
		plot(yp, res, xlab="pred", ylab="res", ...)
		abline(h=0, col="red", lty=2)
	}
}


#' Plot of a non-population dynamic model fit
#'
#' Plot of a non-population dynamic model fit
#'
#' @param x a dynamodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
plot.dyn.ID = gof

#' Print a non-population dynamic model fit object
#'
#' Print a non-population dynamic model fit object
#'
#' @param x a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
print.dyn.ID = function(x, ...)
{
	print(x$res[, c(1,3)])
	cat("\n")
	aic = 2*x$value + 2*x$npar
	bic = 2*x$value + log(x$nobs)*x$npar
	print(c("-loglik"=x$value, "AIC"=aic, "BIC"=bic))
	cat("\n")
}

#' Summary of a non-population dynamic model fit
#'
#' Summary of a non-population dynamic model fit
#'
#' @param object a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
summary.dyn.ID = function(object, ...)
{
	print(object$res)
	cat("\n")
	aic = 2*object$value + 2*object$npar
	bic = 2*object$value + log(object$nobs)*object$npar
	print(c("-loglik"=object$value, "AIC"=aic, "BIC"=bic))
	cat("\niter:", object$iter, "\n")
	cat(object$message, "\n")
}



err.msg = function(x, pre="", post="")
{
	msg = paste0(x, collapse=", ")
	paste0(pre, msg, post)
}


mymin = function(start, fr, rho=NULL, control=list())
{
	if (is.null(rho)) rho = environment(fr)
	step = -.2*start

    con <- list(maxeval=999, ftol_rel=1e-6, rcoeff=1., ecoeff=2., ccoeff=.5, trace=F)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

	.Call(neldermead_wrap, fr, rho, length(start), start, step,
		  as.integer(con$maxeval), con$ftol_rel, con$rcoeff, con$ecoeff, con$ccoeff,
		  as.integer(con$trace),
		  PACKAGE = 'nlmixr')
}

#' Nelder-Mead of simplex search
#'
#' Nelder-Mead of simplex search
#'
#' @param start initials
#' @param fr objective function
#' @param rho evaluation environment
#' @param control additional optimization options
#' @return a list of ...
#' @export
nmsimplex = function(start, fr, rho=NULL, control=list())
{
	if (is.null(rho)) rho = environment(fr)
	step = -.2*start

    con <- list(maxeval=999, reltol=1e-6, rcoeff=1., ecoeff=2., ccoeff=.5, trace=F)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

	.Call(neldermead_wrap, fr, rho, length(start), start, step,
		  as.integer(con$maxeval), con$reltol, con$rcoeff, con$ecoeff, con$ccoeff,
		  as.integer(con$trace),
		  PACKAGE = 'nlmixr')
}


#' Fit a non-population dynamic model
#'
#' Fit a non-population dynamic model
#'
#' @param system an RxODE object
#' @param model a list of statistical meaurement models
#' @param evTable an Event Table object
#' @param inits initial values of system parameters
#' @param data input data
#' @param fixPars fixed system parameters
#' @param method estimation method: choice of Nelder-Mead, L-BFGS-B, and PORT.
#' @param control optional minimization control parameters
#' @param squared if parameters be squared during estimation
#' @return NULL
#' @author Wenping Wang
#' @examples
#' ode <- "
#'    dose=200;
#'    pi = 3.1415926535897931;
#'
#'    if (t<=0) {
#'       fI = 0;
#'    } else {
#'       fI = F*dose*sqrt(MIT/(2.0*pi*CVI2*t^3))*exp(-(t-MIT)^2/(2.0*CVI2*MIT*t));
#'    }
#'
#'    C2 = centr/V2;
#'    C3 = peri/V3;
#'    d/dt(centr) = fI - CL*C2 - Q*C2 + Q*C3;
#'    d/dt(peri)  =              Q*C2 - Q*C3;
#' "
#' sys1 <- RxODE(model = ode)
#'
#'
#' ## ------------------------------------------------------------------------
#' dat <- invgaussian
#' mod <- cp ~ C2 + prop(.1)
#' inits <- c(MIT=190, CVI2=.65, F=.92)
#' fixPars <- c(CL=.0793, V2=.64, Q=.292, V3=9.63)
#' ev <- eventTable()
#' ev$add.sampling(c(0, dat$time))
#' (fit <- dynmodel(sys1, mod, ev, inits, dat, fixPars))
#'
#' @export
dynmodel = function(system, model, evTable, inits, data, fixPars=NULL,
	method=c("Nelder-Mead", "L-BFGS-B", "PORT"),
	control=list(ftol_rel=1e-6, maxeval=999), squared=T)
{
	if (class(model)=="formula") {
		model = list(model)
	}
	inits.err = NULL
	model = lapply(model, function(f) {
		s = unlist(lapply(attr(terms(f),"variables"), as.list))
		s = sapply(s, deparse)

		ix.add = match("add",  s, nomatch=0)
		ix.pro = match("prop", s, nomatch=0)
		err.type = c("add", "prop", "combo")[(ix.add>0)+2*(ix.pro>0)]

		sig.add = if (ix.add>0) as.numeric(s[ix.add+1]) else NULL
		sig.pro = if (ix.pro>0) as.numeric(s[ix.pro+1]) else NULL

		inits.err <<- c(inits.err, sig.add, sig.pro)

		if (any(is.na(inits.err) | inits.err<=0)) stop("error model misspecification")

		s = c(s[2:3], err.type)
		names(s) = c("dv", "pred", "err")
		s
	})
	names(inits.err) = rep("err", length(inits.err))
	inits = c(inits, inits.err)
	#print(model)

	vars = names(data)
	nodef = setdiff(sapply(model, function(x) x["dv"]), vars)
	if (length(nodef)) {
		msg = err.msg(nodef, pre="var(s) not found in data: ")
		stop(msg)
	}

	modelVars = system$cmpMgr$get.modelVars()
	vars = c(modelVars$state, modelVars$lhs)
	nodef = setdiff(sapply(model, function(x) x["pred"]), vars)
	if (length(nodef)) {
		msg = err.msg(nodef, pre="var(s) not found in model: ")
		stop(msg)
	}

	pars = modelVars$params
	nodef = setdiff(pars, c(names(inits), names(fixPars)))
	if (length(nodef)) {
		msg = err.msg(nodef, pre="par(s) not found: ")
		stop(msg)
	}

	npar = length(pars) - length(fixPars)
	have_zero = min(data$time) <= 0
	rows = if(have_zero) T else -1

	if (squared) inits = sqrt(inits)

	obj = function(th)
	{
		#squared = get("squared", envir=sys.parent(n = 1))
		if (squared) th = th^2
		.ixpar = npar
		theta = th[1:npar]
		names(theta) = names(inits)[1:npar]
            theta = c(theta, fixPars)

            s = system$solve(theta, evTable, atol=1e-06, rtol=1e-06)

		l = lapply(model, function(x) {
			err.combo = (x["err"]=="combo")+0
			.ixpar <<- .ixpar+1
			sig = th[.ixpar:(.ixpar+err.combo)]
			sig = if (x["err"]=="add") {
				c(sig, 0)
			} else if (x["err"]=="prop") {
				c(0, sig)
			} else {
				.ixpar <<- .ixpar+1
				sig
			}
			#print(sig)

			yp = s[rows,x["pred"]]
			sgy = thresh(sig[1]+yp*sig[2])
			yo = data[, x["dv"]]
			ll = .5*((yo - yp)^2/sgy^2 + 2*log(sgy) + log(2*pi))
			sum(ll)
		})

		do.call("sum", l)
	}

	method <- match.arg(method)
	if (method=="Nelder-Mead") {
		fit = mymin(as.vector(inits), obj, control=control)
		fit$message=c("NON-CONVERGENCE", "NELDER_FTOL_REACHED")[1+fit$convergence]
	} else if (method=="L-BFGS-B") {
		fit = lbfgs(as.vector(inits), obj, control=control)
	} else {
		if ("ftol_rel" %in% names(control)) {
			control$rel.tol = control$ftol_rel
			control$ftol_rel = NULL
		}
		if ("maxeval" %in% names(control)) {
			control$eval.max = control$maxeval
			control$maxeval = NULL
		}

		fit = nlminb(as.vector(inits), obj, control=control)
	}

	if (squared) fit$par = fit$par^2
	squared = F
	fit$hessian = optimHess(fit$par, obj)
	se = sqrt(diag(solve(fit$hessian)))
	res = cbind(fit$par, se, se/fit$par*100)
	dimnames(res) = list(names(inits), c("est", "se", "%cv"))

	nobs = 0
	l = lapply(model, function(x) {
		yo = data[, x["dv"]]
		nobs <<- nobs + length(yo)
	})

	if (!is.null(fit$objective)) fit$value = fit$objective

	res = c(list(res=res, obj=obj, npar=length(fit$par), nobs=nobs, data=data), fit)
	class(res) = "dyn.ID"

	res
}


#---------------------
uni_slice = function(x0, fr, rho=NULL, w=1, m=1000, lower=-1.0e20, upper=1.0e20)
{
	if (is.null(rho)) rho = environment(fr)
	.Call(slice_wrap, fr, rho, x0, w, as.integer(m), lower, upper, PACKAGE = 'nlmixr')$x1
}

genobj = function(system, model, evTable, inits, data, fixPars=NULL,
	squared=T)
{
	if (class(model)=="formula") {
		model = list(model)
	}
	inits.err = NULL
	model = lapply(model, function(f) {
		s = unlist(lapply(attr(terms(f),"variables"), as.list))
		s = sapply(s, deparse)

		ix.add = match("add",  s, nomatch=0)
		ix.pro = match("prop", s, nomatch=0)
		err.type = c("add", "prop", "combo")[(ix.add>0)+2*(ix.pro>0)]

		sig.add = if (ix.add>0) as.numeric(s[ix.add+1]) else NULL
		sig.pro = if (ix.pro>0) as.numeric(s[ix.pro+1]) else NULL

		inits.err <<- c(inits.err, sig.add, sig.pro)

		if (any(is.na(inits.err) | inits.err<=0)) stop("error model misspecification")

		s = c(s[2:3], err.type)
		names(s) = c("dv", "pred", "err")
		s
	})
	names(inits.err) = paste0("err", 1:length(inits.err))
	inits = c(inits, inits.err)
	#print(model)

	vars = names(data)
	nodef = setdiff(sapply(model, function(x) x["dv"]), vars)
	if (length(nodef)) {
		msg = err.msg(nodef, pre="var(s) not found in data: ")
		stop(msg)
	}

	modelVars = system$cmpMgr$get.modelVars()
	vars = c(modelVars$state, modelVars$lhs)
	nodef = setdiff(sapply(model, function(x) x["pred"]), vars)
	if (length(nodef)) {
		msg = err.msg(nodef, pre="var(s) not found in model: ")
		stop(msg)
	}

	pars = modelVars$params
	nodef = setdiff(pars, c(names(inits), names(fixPars)))
	if (length(nodef)) {
		msg = err.msg(nodef, pre="par(s) not found: ")
		stop(msg)
	}

	npar = length(pars) - length(fixPars)
	have_zero = min(data$time) <= 0
	rows = if(have_zero) T else -1

	if (squared) inits = sqrt(inits)
	s.save = NULL

	obj = function(th, do.ode.solving=T, negation=F)
	{
		#squared = get("squared", envir=sys.parent(n = 1))
		if (squared) th = th^2
		.ixpar = npar
		theta = th[1:npar]
		names(theta) = names(inits)[1:npar]
		theta = c(theta, fixPars)
		if (do.ode.solving) {
			s = system$solve(theta, evTable, atol=1e-06, rtol=1e-06)
			s.save <<- s
		} else {
			s = s.save
		}

		l = lapply(model, function(x) {
			err.combo = (x["err"]=="combo")+0
			.ixpar <<- .ixpar+1
			sig = th[.ixpar:(.ixpar+err.combo)]
			sig = if (x["err"]=="add") {
				c(sig, 0)
			} else if (x["err"]=="prop") {
				c(0, sig)
			} else {
				.ixpar <<- .ixpar+1
				sig
			}
			#print(sig)

			yp = s[rows,x["pred"]]
			sgy = thresh(sig[1]+yp*sig[2])
			yo = data[, x["dv"]]
			ll = .5*((yo - yp)^2/sgy^2 + 2*log(sgy) + log(2*pi))
			sum(ll)
		})

		res = do.call("sum", l)
		if (negation) -res else res
	}
	list(obj=obj, inits=inits)
}


#-- mcmc
error.terms = paste0("err", 1:40)

do.slice = function(pars, fr0)
{
	rho = environment()
	lapply(names(pars), function(wh, fr0) {
		do.ode.solving = match(wh, error.terms, nomatch=0) == 0
		pars.cp = get("pars", rho)
		x0 = pars.cp[wh]
		fr = function(x) {
			pars.cp[wh] = x
			fr0(pars.cp, do.ode.solving=do.ode.solving, negation=T)
		}
		pars.cp[wh] = uni_slice(x0, fr, lower=0)
		assign("pars", pars.cp, rho)
		NULL
	}, fr0=fr0)

	pars
}


#' Fit a non-population dynamic model using mcmc
#'
#' Fit a non-population dynamic model using mcmc
#'
#' @param system an RxODE object
#' @param model a list of statistical meaurement models
#' @param evTable an Event Table object
#' @param inits initial values of system parameters
#' @param data input data
#' @param fixPars fixed system paraameters
#' @param nsim number of mcmc iteractions
#' @param squared if parameters be squared during estimation
#' @param seed random number seed
#' @author Wenping Wang
#' @return NULL
#' @examples
#' \dontrun{
#'
#' ode <- "
#'    dose=200;
#'    pi = 3.1415926535897931;
#'
#'    if (t<=0) {
#'       fI = 0;
#'    } else {
#'       fI = F*dose*sqrt(MIT/(2.0*pi*CVI2*t^3))*exp(-(t-MIT)^2/(2.0*CVI2*MIT*t));
#'    }
#'
#'    C2 = centr/V2;
#'    C3 = peri/V3;
#'    d/dt(centr) = fI - CL*C2 - Q*C2 + Q*C3;
#'    d/dt(peri)  =              Q*C2 - Q*C3;
#' "
#' sys1 <- RxODE(model = ode)
#'
#'
#' ## ------------------------------------------------------------------------
#' dat <- read.table("invgaussian.txt", header=TRUE)
#' mod <- cp ~ C2 + prop(.1)
#' inits <- c(MIT=190, CVI2=.65, F=.92)
#' fixPars <- c(CL=.0793, V2=.64, Q=.292, V3=9.63)
#' ev <- eventTable()
#' ev$add.sampling(c(0, dat$time))
#' (fit <- dynmodel.mcmc(sys1, mod, ev, inits, dat, fixPars))
#'
#' }
#' @export
dynmodel.mcmc = function(system, model, evTable, inits, data,
	fixPars=NULL, nsim = 500, squared=T, seed=NULL)
{
	calls = match.call()

	l = genobj(system, model, evTable, inits, data, fixPars, squared)
	rho = environment()
	pars = l$inits
	fr0 = l$obj

	if (is.null(seed)) seed=99
	set.seed(seed)
	s = t(sapply(1:nsim, function(k,rho) {
		pars = do.slice(get("pars", rho), fr0)
		assign("pars", pars, rho)
	}, rho=rho))

	if (squared) s = s*s
	attr(s, "calls") <- calls
	attr(s, "obj") <- fr0
	attr(s, "class") <- "dyn.mcmc"
	s
}

#' Summary of a non-population dynamic model fit using mcmc
#'
#' Summary of a non-population dynamic model fit using mcmc
#'
#' @param object a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
summary.dyn.mcmc = function(object, ...)
{
	s <- t(apply(object, 2, function(x) c(mean(x), sd(x), sd(x)/mean(x)*100)))
	dimnames(s)[[2]] <- c("mean", "sd", "cv%")
	print(s)
	cat("\n# samples:", dim(object)[1], "\n")
}

#' Summary of a non-population dynamic model fit using mcmc
#'
#' Summary of a non-population dynamic model fit using mcmc
#'
#' @param x a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
print.dyn.mcmc = function(x, ...)
{
	s <- t(apply(x, 2, function(x) c(mean(x), sd(x), sd(x)/mean(x)*100)))
	dimnames(s)[[2]] <- c("mean", "sd", "cv%")
	print(s)
	cat("\n# samples:", dim(x)[1], "\n")
}

#' Plot of a non-population dynamic model fit using mcmc
#'
#' Plot of a non-population dynamic model fit using mcmc
#'
#' @param x a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
plot.dyn.mcmc = function(x, ...)
{
	fit = list(obj=attr(x, "obj"), par=apply(x, 2, mean))
	gof(fit)
}



## Utilities for building nlmixr

refresh <- function(){
    ## nocov start
    source(devtools::package_file("build/refresh.R"))
    ## nocov end
}

nsis <- function(){ ## build installer...
    ## nocov start
    source(devtools::package_file("build/nsis.R"))
    ## nocov end
}
##' Collect warnings and just warn once.
##'
##' @param expr R expression
##' @return The value of the expression
##' @author Matthew L. Fidler
.collectWarnings <- function(expr){
    ws <- c();
    this.env <- environment()
    ret <- suppressWarnings(withCallingHandlers(expr,warning=function(w){assign("ws", unique(c(w$message, ws)), this.env)}))
    for (w in ws){
        warning(w)
    }
    return(ret);
}

##' Print x using the message facility
##'
##' This allows the suppressMessages to work on print functions.  This
##' captures the output via R.Util's captureOutput function and then
##' sends it through the message routine.
##'
##' catpureOutput was used since it is much faster than the internal
##' capture.output see https://www.r-bloggers.com/performance-captureoutput-is-much-faster-than-capture-output/
##' @param x object to print
##' @param ... Other things output
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
nlmixrPrint <- function(x, ...){
    this.env <- environment();
    message(invisible(paste(R.utils::captureOutput(assign("x", print(x, ...), this.env)), collapse="\n")), appendLF=TRUE);
    invisible(x)
}

.dontRun <- function(...){
    ## This is for r checks, though they need to be loaded...
    vpc::vpc(...)
    dparser::dparse(...)
}
