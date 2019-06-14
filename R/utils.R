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
    #if (length(noNms <- namc[!namc %in% nmsC]))
    #    warning("unknown names in control: ", paste(noNms, collapse = ", "))

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
dynmodel = function(system, model, evTable, inits, data, fixPars=NULL, lower = -Inf, upper = Inf,
  method=c("bobyqa", "Nelder-Mead", "lbfgsb3c", "PORT"),
	control=list(ftol_rel=1e-6, 
	             maxeval=999,
	             scaleTo=1.0,
	             scaleObjective=0,
	             normType=c("constant","rescale2", "mean", "rescale", "std", "len"),
	             scaleType=c("norm","nlmixr", "mult", "multAdd"),
	             scaleCmax=1e5,
	             scaleCmin=1e-5,
	             scaleC=1, #NULL,
	             scaleC0=1e5,
	             transformRUV=c("boxCox", "YeoJohnson")),
	rxControl=list(atol=1e-08,rtol=1e-06)
	)
{
# Control Argument Defaluts -----------------------------------------------
if (is.null(control$scaleTo)) {control$scaleTo = 1.0}
#if (is.null(control$transformRUV)) {control$transformRUV = 1.0}

# Error model Handling -------------------------------------------------------------
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

# Objective Function ------------------------------------------------------
	transformRUV <- 1
	
	obj = function(th)
	{
	  # unscale ------------------- #
	  unscaled.th <- numeric(length(th))
	  for (i in 1:length(th)) {th[i] <- unscalePar(th,i)}

	   # define parameters used for simulation, all parameters except the error terms
		.ixpar = npar
		theta = th[1:npar] 
		names(theta) = names(inits)[1:npar]
		theta = c(theta, fixPars)
    
		s = do.call(RxODE :: rxSolve, c(list(system, theta, evTable, rxControl)))
    #rxNorm(system) # add error piece, use the error parameters as parameters
    #nlmixr_err=boxCox(add.err^2+prop.err^(2*pow)*pred, lambda) for the predictions not the errors
    #model.test <- unlist(model)
    #print(model.test[3]=="prop")
    
		l = lapply(model, function(x) {
		  
		  err.combo = (x["err"]=="combo")+0  # returns 1 for combintion (add + prop) or 0 for prop or add alone
		  .ixpar <<- .ixpar+1 # create global parameter, that increaes the number of parameters (parameters being estimation, without err) by 1
		  sig = th[.ixpar:(.ixpar+err.combo)] #if "combo" is used, it grabs both terms, or else it grabs just add or prop
		  sig = if (x["err"]=="add") {   # changes sigma from a 1D vector into a 2D vector, with orientation sig[1] = add, sig[2] = prop
				c(sig, 0)
			} else if (x["err"]=="prop") {
				c(0, sig)
			} else {
				.ixpar <<- .ixpar+1
				sig
			}
			yp = s[rows,x["pred"]]
			sgy = thresh(sig[1]+yp*sig[2]) # variance of Y for FOCEi, threshold used to prevent 0?
			yo = data[, x["dv"]]
			

			
			# assign sigma for the error model that are used in var[Y]
			
      # Transform both sides (used for non-normal residuals)
			# may need to implement nlmixr::coxBox
			boxCox = function (x, lambda) {
			   if(lambda == 0) {
			     .h.x <- log(x)
			     }
			   else {
			     .h.x <- (x^lambda-1)/lambda
			   }
			   return(.h.x)
			}
			
			# yeoJohnson = need function or figure out how to use nlmixr::yeoJohnson
			yeoJohnson = function (x, lambda){
			 options(warn=-1)
			.h.x <- ifelse (x >= 0,
        ifelse (lambda != 0 & x >= 0, 
                ((x + 1)^lambda - 1)/lambda, 
                log(x + 1)),
			  ifelse (lambda != 2 & x < 0, 
			          -((-x+1)^(2-lambda)-1)/(2-lambda), 
			          -log(-x+1))
			)
			options(warn=0)
			 return(.h.x)
			}

			# need to define power sigma coefficient

      if (transformRUV == "boxCox") {
			   .prop.sig = 1 # need to reassign above
			   .add.sig = 0 # need to reassign above
			   .power = 1 # need to reassign above
			   
			   .h.x <- boxCox(yo, lambda) # obs
			   .h.y <- boxCox(yp, lambda) # pred
			   .h.y.var <- yp^(2*.power)*.prop.sig^2 + .add.sig^2  # variance of pred

			   # boxCox transformed -2 log-likelihood
			   .boxCox.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
			   
			   # back-transformed  -2 log-likelihood function, with penalty added
			   .n2ll = .boxCox.n2ll - 2*(lambda-1)*log(yo) -2*log(2*pi)
			   
			   # negative log-likelihood function for output
			   ll = .5*(.n2ll)
			   
			 } else if(transformRUV == "yeoJohnson"){
			   .prop.sig = 1 # need to reassign above
			   .add.sig = 0 # need to reassign above
			   .power.sig = 1 # need to reassign above
			   
			   .h.x <- yeoJohnson(yo, lambda) #obs
			   .h.y <- yeoJohnson(yp, lambda) #pred
			   .h.y.var <- yp^(2*.power)*.prop.sig^2 + .add.sig^2  # variance of pred
			     
			   # yeoJohnson transformed -2 log-likelihood
			   .yeoJohnson.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
			   
			   # back-transformed  -2 log-likelihood function, with penalty added
			   .n2ll <- ifelse(x >= 0, 
			                   .yeoJohnson.n2ll -2*((lambda-1)*log(x+1) -2*log(2*pi)),
			                   .yeoJohnson.n2ll -2*((1-lambda)*log(-x+1) -2*log(2*pi))
			                   )
			   
			   # negative log-likelihood function for output
			   ll = .5*(.n2ll)
			   
			 }
			 else{
			   ll = .5*((yo - yp)^2/sgy^2 + 2*log(sgy) + log(2*pi)) # negative log likelihood, as a vector
			 }

		sum(ll)
		})

		
		do.call("sum", l)  # same as return(as.numeric(l))
		
	}


# Parameter Normalization and Scaling -----------------------------------------------------------------------
	
	# normType assignment for scaling
 	normType <- control$normType
	if (normType == "constant") {
	  print("constant used")
	  C1 = 0
	  C2 = 1
	} else if (normType == "rescale2") {
	  print("rescale2 used")
	  C1 = (max(inits) + min(inits))/2
	  C2 = (max(inits) - min(inits))/2
	} else if (normType == "mean") {
	  print("mean used")
	  C1 = mean(inits)
	  C2 = max(inits) - min(inits)
	} else if (normType == "rescale") {
	  print("rescale used")
	  C1 = min(inits)
	  C2 = max(inits) - min(inits)
	} else if (normType == "std") {
	  print("std used")
	  C1 = mean(inits)
	  C2 = sd(inits)
	} else if (normType == "len") {
	  print("len used")
	  C1 = 0
	  C2 = sqrt(sum(inits*inits))
	}
	
	# scaleC assignment for scaling (adopted from foceFIT.R)
	scaleC <- double(length(inits));
	if (is.null(control$scaleC)){
	  scaleC <- rep(1, length(inits))
	} else {
	  scaleC <- as.double(control$scaleC);
	  if (length(inits) > length(scaleC)){
	    scaleC <- c(scaleC, rep(1, length(inits) - length(scaleC)));
	  } else if (length(inits) < length(scaleC)){
	    scaleC <- scaleC[seq(1, length(inits))];
	    warning("scaleC control option has more options than estimated parameters, please check.")
	  }
	}	
	
  # Function for scaling parameters based on scaleType
	  scaleType <- control$scaleType
	  scaleTo <- control$scaleTo
		scalePar <- function(x, i){
  	# simple scaling
  	if (scaleType == "norm") {
  	  print("norm used!")
  	  return((x[i]-C1)/C2)
  	}
  	# nlmixr
  	else if (scaleType == "nlmixr") {
  	  scaleTo = (inits[i] - C1)/C2
  	  print("nlmixr used!")
  	  return((x[i] - inits[i])/scaleC[i] + scaleTo)
  	}
  	# simple multiplicatice scaling
  	else if (scaleType == "mult") {
  	    if (scaleTo > 0) {
  	      print("mult 1 used!")
  	      return(x[i]/inits[i]*scaleTo)
  	      
  	    } else {
  	      print("mult 2 used!")
  	      return(x[i])
  	    }
  	}
  	# log non-log multiplicative scaling
  	else if (scaleType == "multAdd") {
  	  if(scaleTo > 0) {
  	    print("multAdd 1 used!")
  	    return((x[i]-inits[i]) + scaleTo)
  	  } else {
  	    print("multAdd 2 used!")
  	    return(x[i]/inits[i]*scaleTo)
  	  }
  	}
  	
# When should this be used? "norm" scaling is essentially no scaling when normType specified.
  	else {
  	  if(scaleTo > 0) {
  	    print("Other 1 used!")
  	    return((x[i] - inits[i]) + scaleTo)
  	  } else {
  	    print("Other 2 used!")
  	    return(x[i])
  	  }
  	}
	}
	
  # Function for unscaling parameters based on scaleType
	unscalePar <- function(x, i){
	  scaleType <- match.arg(control$scaleType, c("norm","nlmixr", "mult", "multAdd"))
	  # simple scaling
	  if (scaleType == "norm") {
	    return(x[i]*C2+C1);
	  }
	  # nlmixr
	  else if (scaleType == "nlmixr") {
	    scaleTo = (inits[i] - C1)/C2
	    return((x[i] - scaleTo)*scaleC[i] + inits[i])
	  }
	  # simple multiplicatice scaling
	  else if (scaleType == "mult") {
	    if (scaleTo > 0) {
	      return(x[i]*inits[i]/scaleTo)
	    } else {
	      return(x[i])
	    }
	  }
	  # log non-log multiplicative scaling
	  else if (scaleType == "multAdd") {
	    if(scaleTo > 0) {
	      return((x[i]- scaleTo) + inits[i])
	    } else {
	      return(x[i]*inits[i]/scaleTo)
	    }
	  }
	  
# When should this be used? "norm" scaling is essentially no scaling when normType specified.
	  else {
	    if(scaleTo > 0) {
	      return((x[i]-scaleTo)*1 + inits[i])
	    } else {
	      return(x[i])
	    }
	  }
	}
	
# Optimization Method -----------------------------------------------------------------------
	method <- match.arg(method)
	if (method =="bobyqa") {
	  control <- control[names(control) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]

	  # change control and see if it works better, default to these
	  #npt = npar*2 + 1
	  #rhobeg = 0.2
	  #rhoend = 1E-04
	  
	  # scaled the initial conditions
	  scaled.inits <- numeric(length(inits))
	  for (i in 1:length(inits)) {scaled.inits[i] <- scalePar(inits,i)}

	  # scaled lower boundary
	  if (is.null(lower)==FALSE){
  	  if (any(lower == -Inf)) {
  	    lower[which(lower==-Inf)]=0
  	    warning("Lower boundary of -Inf set to 0.")
  	  } 
  	  scaled.lower <-  numeric(length(lower))
  	  for (i in 1:length(lower)) {scaled.lower[i] <- scalePar(lower,i)}
	  } else {
	    scaled.lower <- NULL
	    scaled.upper <- NULL
	  }
	  
	  # scaled upper boundary
	  if (is.null(upper)==FALSE){
  	  scaled.upper <- numeric(length(upper))
  	  for (i in 1:length(upper)) {scaled.upper[i] <- scalePar(upper,i)}
	  } else {scaled.upper <- NULL}

	  fit = minqa::bobyqa(scaled.inits, obj, scaled.lower, scaled.upper, control=control)
	  fit$value <- fit$fval
    
	  print(inits)
	  print(fit$par)
	  
	  # unscale the output
	  par.temp <- numeric(length(fit$par))
	  for (i in 1:length(par.temp)) {par.temp[i] <- unscalePar(fit$par,i)}
	  fit$par <- par.temp
	  
	  print(fit$par)

	} else if (method=="lbfgsb3c") {
	  # remove any control elements nor required in lbfgsb3c()
	  lbfgsb3c.control <- list(trace=NULL,factr=NULL,pgtol=NULL,abstol=NULL,reltol=NULL,lmm=NULL,maxit=NULL,iprint=NULL)
	  if (any(names(lbfgsb3c.control) %in% names(control))) {
	    control = control[match(names(lbfgsb3c.control),names(control))[!is.na(match(names(lbfgsb3c.control),names(control)))]]
	  } else{control=NULL}
		fit = lbfgsb3c::lbfgsb3c(par = as.vector(inits), fn=obj, lower=lower, upper=upper, control=control)
	} else if (method=="Nelder-Mead") {
	  fit = mymin(as.vector(inits), obj, control=control)
	  fit$message=c("NON-CONVERGENCE", "NELDER_FTOL_REACHED")[1+fit$convergence]
	} else {
		if ("ftol_rel" %in% names(control)) {
			control$rel.tol = control$ftol_rel
			control$ftol_rel = NULL
		}
		if ("maxeval" %in% names(control)) {
			control$eval.max = control$maxeval
			control$maxeval = NULL
		}
		fit = nlminb(as.vector(inits), obj, lower, upper, control=control)
	}

# Hessian Calculation -----------------------------------------------------------------------
	fit$hessian = try(optimHess(fit$par, obj, control=control) , silent=TRUE)
	
	if(inherits(fit$hessian,"try-error")){
	  se = rep(NA, length(fit$par))
	  warning("standard error of the Hessian has failed")
	} else {
	  se = sqrt(diag(solve(fit$hessian)))
	}
	
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

####

uni_slice = function(x0, fr, rho=NULL, w=1, m=1000, lower=-1.0e20, upper=1.0e20)
{
	if (is.null(rho)) rho = environment(fr)
	.Call(slice_wrap, fr, rho, x0, w, as.integer(m), lower, upper, PACKAGE = 'nlmixr')$x1
}


genobj = function(system, model, evTable, inits, data, fixPars=NULL
	)
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
	
	## is this necessary ##
	have_zero = min(data$time) <= 0
	rows = if(have_zero) T else -1 # used in line 304 in obj()
	## ---------------- ##


	obj = function(th, do.ode.solving=T, negation=F)
	{
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
##' @param lst When \code{TRUE} return a list with
##'     list(object,warnings) instead of issuing the warnings.
##'     Otherwise, when \code{FALSE} issue the warnings and return the
##'     object.
##' @return The value of the expression or a list with the value of
##'     the expression and a list of warning messages
##' @author Matthew L. Fidler
##' @noRd
.collectWarnings <- function(expr,lst=FALSE){
    ws <- c();
    this.env <- environment()
    ret <- suppressWarnings(withCallingHandlers(expr,warning=function(w){assign("ws", unique(c(w$message, ws)), this.env)}))
    if (lst){
        return(list(ret, ws));
    } else {
        for (w in ws){
            warning(w)
        }
        return(ret);
    }
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

##' Generalized Cholesky Matrix Decomposition
##'
##'  Performs a (modified) Cholesky factorization of the form
##'
##'   t(P) \%*\% A \%*\% P  + E = t(R) \%*\% R
##'
##'  As detailed in Schnabel/Eskow (1990)
##'
##' @param matrix Matrix to be Factorized.
##' @param tol Tolerance; Algorithm suggests (.Machine$double.eps) ^ (1 / 3), default
##' @return Generalized Cholesky decomposed matrix.
##' @author Matthew L. Fidler (translation), Johannes Pfeifer, Robert
##'     B. Schnabel and Elizabeth Eskow
##'
##' @references
##'
##' matlab source: http://www.dynare.org/dynare-matlab-m2html/matlab/chol_SE.html; Slightly different return values
##'
##' Robert B. Schnabel and Elizabeth
##' Eskow. 1990. "A New Modified Cholesky Factorization," SIAM Journal
##' of Scientific Statistical Computing, 11, 6: 1136-58.
##'
##' Elizabeth Eskow and Robert B. Schnabel
##' 1991. "Algorithm 695 - Software for a New Modified Cholesky Factorization,"
##' ACM Transactions on Mathematical Software, Vol 17, No 3: 306-312
##'
##' @note
##'
##' This version does not pivot or return the E matrix
##'
##' @export
cholSE <- function(matrix, tol=(.Machine$double.eps) ^ (1 / 3)){
    .Call(`_nlmixr_cholSE_`, matrix, tol);
}

.setRoot <- function(){
    setwd("c:/");
}
