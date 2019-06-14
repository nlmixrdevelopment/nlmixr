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

dynmodel = function(system, model, evTable, inits, data, ...)
{

# dynmodel control -----------------------------------------------
dynmodelControl <- function(
                            fixPars=NULL,
                            lower = -Inf,
                            upper = Inf,
                            method=c("bobyqa", "Nelder-Mead", "lbfgsb3c", "PORT"),
                            ftol_rel=1e-6,
                            maxeval=999,
                            scaleTo=1.0,
                            scaleObjective=0,
                            normType=c("constant","rescale2", "mean", "rescale", "std", "len"),
                            scaleType=c("norm","nlmixr", "mult", "multAdd"),
                            scaleCmax=1e5,
                            scaleCmin=1e-5,
                            scaleC=1, #NULL,
                            scaleC0=1e5,
                            transformRUV=c("boxCox", "YeoJohnson"),
                            atol=1e-08,
                            rtol=1e-06
                          ) {
 if (is.null(lower)){
   lower = -Inf
 }
  }

print(dynmodel$control)

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

}


{
##' Control Options for dynmodel
##'
##' @param sigdig Optimization significant digits. This controls:
##'
##' \itemize{
##'
##'  \item The tolerance of the inner and outer optimization is \code{10^-sigdig}
##'
##'  \item The tolerance of the ODE solvers is
##'  \code{0.5*10^(-sigdig-2)}; For the sensitivity equations the
##'  default is \code{0.5*10^(-sigdig-1.5)} (only applicable for liblsoda)
##'
##'  \item The tolerance of the boundary check is \code{5 * 10 ^ (-sigdig + 1)}
##'
##'  \item The significant figures that some tables are rounded to.
##' }
##'
##' @param atolSens Sensitivity atol, can be different than atol with
##'     liblsoda.  This allows a less accurate solve for gradients (if desired)
##'
##' @param rtolSens Sensitivity rtol, can be different than rtol with
##'     liblsoda.  This allows a less accurate solve for gradients (if desired)
##'
##' @param epsilon Precision of estimate for n1qn1 optimization.
##'
##' @param maxstepsOde Maximum number of steps for ODE solver.
##'
##' @param print Integer representing when the outer step is
##'     printed. When this is 0 or do not print the iterations.  1 is
##'     print every function evaluation (default), 5 is print every 5
##'     evaluations.
##'
##' @param scaleTo Scale the initial parameter estimate to this value.
##'     By default this is 1.  When zero or below, no scaling is performed.
##'
##' @param scaleObjective Scale the initial objective function to this
##'     value.  By default this is 1.
##'
##' @param derivEps Forward difference tolerances, which is a
##'     vector of relative difference and absolute difference.  The
##'     central/forward difference step size h is calculated as:
##'
##'         \code{h = abs(x)*derivEps[1] + derivEps[2]}
##'
##' @param derivMethod indicates the method for calculating
##'     derivatives of the outer problem.  Currently supports
##'     "switch", "central" and "forward" difference methods.  Switch
##'     starts with forward differences.  This will switch to central
##'     differences when abs(delta(OFV)) <= derivSwitchTol and switch
##'     back to forward differences when abs(delta(OFV)) >
##'     derivSwitchTol.
##'
##' @param derivSwitchTol The tolerance to switch forward to central
##'     differences.
##'
##' @param covDerivMethod indicates the method for calculating the
##'     derivatives while calculating the covariance components
##'     (Hessian and S).
##'
##' @param covMethod Method for calculating covariance.  In this
##'     discussion, R is the Hessian matrix of the objective
##'     function. The S matrix is the sum of individual
##'     gradient cross-product (evaluated at the individual empirical
##'     Bayes estimates).
##'
##' \itemize{
##'
##'  \item "\code{r,s}" Uses the sandwich matrix to calculate the
##'  covariance, that is: \code{solve(R) \%*\% S \%*\% solve(R)}
##'
##'  \item "\code{r}" Uses the Hessian matrix to calculate the
##'  covariance as \code{2 \%*\% solve(R)}
##'
##'  \item "\code{s}" Uses the cross-product matrix to calculate the
##'  covariance as \code{4 \%*\% solve(S)}
##'
##'  \item "" Does not calculate the covariance step.
##' }
##'
##' @param covTryHarder If the R matrix is non-positive definite and
##'     cannot be corrected to be non-positive definite try estimating
##'     the Hessian on the unscaled parameter space.
##'
##' @param hessEps is a double value representing the epsilon for the Hessian calculation.
##'
##' @param centralDerivEps Central difference tolerances.  This is a
##'     numeric vector of relative difference and absolute difference.
##'     The central/forward difference step size h is calculated as:
##'
##'         \code{h = abs(x)*derivEps[1] + derivEps[2]}
##'
##' @param lbfgsLmm An integer giving the number of BFGS updates
##'     retained in the "L-BFGS-B" method, It defaults to 7.
##'
##' @param lbfgsPgtol is a double precision variable.
##'
##'     On entry pgtol >= 0 is specified by the user.  The iteration
##'     will stop when:
##'
##'        \code{max(\| proj g_i \| i = 1, ..., n) <= lbfgsPgtol}
##'
##'     where pg_i is the ith component of the projected gradient.
##'
##'     On exit pgtol is unchanged.  This defaults to zero, when the
##'     check is suppressed.
##'
##' @param lbfgsFactr Controls the convergence of the "L-BFGS-B"
##'     method.  Convergence occurs when the reduction in the
##'     objective is within this factor of the machine
##'     tolerance. Default is 1e10, which gives a tolerance of about
##'     \code{2e-6}, approximately 4 sigdigs.  You can check your
##'     exact tolerance by multiplying this value by
##'     \code{.Machine$double.eps}
##'
##' @param diagXform This is the transformation used on the diagonal
##'     of the \code{chol(solve(omega))}. This matrix and values are the
##'     parameters estimated in FOCEi. The possibilities are:
##'
##' \itemize{
##'  \item \code{sqrt} Estimates the sqrt of the diagonal elements of \code{chol(solve(omega))}.  This is the default method.
##'
##'  \item \code{log} Estimates the log of the diagonal elements of \code{chol(solve(omega))}
##'
##'  \item \code{identity} Estimates the diagonal elements without any transformations
##' }
##' @param sumProd Is a boolean indicating if the model should change
##'     multiplication to high precision multiplication and sums to
##'     high precision sums using the PreciseSums package.  By default
##'     this is \code{FALSE}.
##'
##'
##' @param optExpression Optimize the RxODE expression to speed up
##'     calculation. By default this is turned on.
##'
##' @param ci Confidence level for some tables.  By default this is
##'     0.95 or 95\% confidence.
##'
##' @param useColor Boolean indicating if focei can use ASCII color codes
##'
##' @param boundTol Tolerance for boundary issues.
##'
##' @param calcTables This boolean is to determine if the foceiFit
##'     will calculate tables. By default this is \code{TRUE}
##'
##' @param ... Ignored parameters
##'
##' @param maxInnerIterations Number of iterations for n1qn1
##'     optimization.
##'
##' @param maxOuterIterations Maximum number of L-BFGS-B optimization
##'     for outer problem.
##'
##' @param n1qn1nsim Number of function evaluations for n1qn1
##'     optimization.
##'
##' @param eigen A boolean indicating if eigenvectors are calculated
##'     to include a condition number calculation.
##'
##' @param addPosthoc Boolean indicating if posthoc parameters are
##'     added to the table output.
##'
##' @param printNcol Number of columns to printout before wrapping
##'     parameter estimates/gradient
##'
##' @param noAbort Boolean to indicate if you should abort the FOCEi
##'     evaluation if it runs into troubles.  (default TRUE)
##'
##' @param interaction Boolean indicate FOCEi should be used (TRUE)
##'     instead of FOCE (FALSE)
##'
##' @param cholSEOpt Boolean indicating if the generalized Cholesky
##'     should be used while optimizing.
##'
##' @param cholSECov Boolean indicating if the generalized Cholesky
##'     should be used while calculating the Covariance Matrix.
##'
##' @param fo is a boolean indicating if this is a FO approximation routine.
##'
##' @param cholSEtol tolerance for Generalized Cholesky
##'     Decomposition.  Defaults to suggested (.Machine$double.eps)^(1/3)
##'
##' @param cholAccept Tolerance to accept a Generalized Cholesky
##'     Decomposition for a R or S matrix.
##'
##' @param outerOpt optimization method for the outer problem
##'
##' @param innerOpt optimization method for the inner problem (not
##'     implemented yet.)
##'
##' @param stateTrim Trim state amounts/concentrations to this value.
##'
##' @param resetEtaP represents the p-value for reseting the
##'     individual ETA to 0 during optimization (instead of the saved
##'     value).  The two test statistics used in the z-test are either
##'     chol(omega^-1) \%*\% eta or eta/sd(allEtas).  A p-value of 0
##'     indicates the ETAs never reset.  A p-value of 1 indicates the
##'     ETAs always reset.
##'
##' @param resetThetaP represents the p-value for reseting the
##'     population mu-referenced THETA parameters based on ETA drift
##'     during optimization, and resetting the optimization.  A
##'     p-value of 0 indicates the THETAs never reset.  A p-value of 1
##'     indicates the THETAs always reset and is not allowed.  The
##'     theta reset is checked at the beginning and when nearing a
##'     local minima.  The percent change in objective function where
##'     a theta reset check is initiated is controlled in
##'     \code{resetThetaCheckPer}.
##'
##' @param resetThetaCheckPer represents objective function
##'     \% percentage below which resetThetaP is checked.
##'
##' @param resetThetaFinalP represents the p-value for reseting the
##'     population mu-referenced THETA parameters based on ETA drift
##'     during optimization, and resetting the optimization one final time.
##'
##' @param resetHessianAndEta is a boolean representing if the
##'     individual Hessian is reset when ETAs are reset using the
##'     option \code{resetEtaP}.
##'
##' @param diagOmegaBoundUpper This represents the upper bound of the
##'     diagonal omega matrix.  The upper bound is given by
##'     diag(omega)*diagOmegaBoundUpper.  If
##'     \code{diagOmegaBoundUpper} is 1, there is no upper bound on
##'     Omega.
##'
##' @param diagOmegaBoundLower This represents the lower bound of the
##'     diagonal omega matrix.  The lower bound is given by
##'     diag(omega)/diagOmegaBoundUpper.  If
##'     \code{diagOmegaBoundLower} is 1, there is no lower bound on
##'     Omega.
##'
##' @param rhobeg Beginning change in parameters for bobyqa algorithm
##'     (trust region).  By default this is 0.2 or 20% of the initial
##'     parameters when the parameters are scaled to 1. rhobeg and
##'     rhoend must be set to the initial and final values of a trust
##'     region radius, so both must be positive with 0 < rhoend <
##'     rhobeg. Typically rhobeg should be about one tenth of the
##'     greatest expected change to a variable.  Note also that
##'     smallest difference abs(upper-lower) should be greater than or
##'     equal to rhobeg*2. If this is not the case then rhobeg will be
##'     adjusted.
##'
##' @param rhoend The smallest value of the trust region radius that
##'     is allowed. If not defined, then 10^(-sigdig-1) will be used.
##'
##' @param npt The number of points used to approximate the objective
##'     function via a quadratic approximation for bobyqa. The value
##'     of npt must be in the interval [n+2,(n+1)(n+2)/2] where n is
##'     the number of parameters in par. Choices that exceed 2*n+1 are
##'     not recommended. If not defined, it will be set to 2*n + 1
##' @param eval.max Number of maximum evaluations of the objective function
##'
##' @param iter.max Maximum number of iterations allowed.
##'
##' @param rel.tol Relative tolerance before nlminb stops.
##'
##' @param x.tol X tolerance for nlmixr optimizers
##'
##' @param abstol Absolute tolerance for nlmixr optimizer
##'
##' @param reltol  tolerance for nlmixr
##'
##' @param gillK The total number of possible steps to determine the
##'     optimal forward/central difference step size per parameter (by
##'     the Gill 1983 method).  If 0, no optimal step size is
##'     determined.  Otherwise this is the optimal step size
##'     determined.
##'
##' @param gillRtol The relative tolerance used for Gill 1983
##'     determination of optimal step size.
##'
##' @param scaleType The scaling scheme for nlmixr.  The supported types are:
##'
##' \itemize{
##' \item \code{nlmixr}  In this approach the scaling is performed by the following equation:
##'
##'    v_{scaled} = (v_{current} - v_{init})/scaleC[i] + scaleTo
##'
##' The \code{scaleTo} parameter is specified by the \code{normType},
##' and the scales are specified by \code{scaleC}.
##'
##' \item \code{norm} This approach uses the simple scaling provided
##'     by the \code{normType} argument.
##'
##' \item \code{mult} This approach does not use the data
##' normalization provided by \code{normType}, but rather uses
##' multiplicitve scaling to a constant provided by the \code{scaleTo}
##' argument.
##'
##'   In this case:
##'
##'   v_{scaled} = v_{current}/v_{init}*scaleTo
##'
##' \item \code{multAdd} This approach changes the scaling based on
##' the parameter being specified.  If a parameter is defined in an
##' exponenital block (ie exp(theta)), then it is scaled on a
##' linearly, that is:
##'
##'   v_{scaled} = (v_{current}-v_{init}) + scaleTo
##'
##' Otherwise the parameter is scaled multiplicatively.
##'
##'    v_{scaled} = v_{current}/v_{init}*scaleTo
##'
##' }
##'
##' @param scaleC The scaling constant used with
##'     \code{scaleType=nlmixr}.  When not specified, it is based on
##'     the type of parameter that is estimated.  The idea is to keep
##'     the derivatives similar on a log scale to have similar
##'     gradient sizes.  Hence parameters like log(exp(theta)) would
##'     have a scaling factor of 1 and log(theta) would have a scaling
##'     factor of ini_value (to scale by 1/value; ie
##'     d/dt(log(ini_value)) = 1/ini_value or scaleC=ini_value)
##'
##'    \itemize{
##'
##'    \item For parameters in an exponential (ie exp(theta)) or
##'    parameters specifying powers, boxCox or yeoJohnson
##'    transformations , this is 1.
##'
##'    \item For additive, proportional, lognormal error structures,
##'    these are given by 0.5*abs(initial_estimate)
##'
##'    \item Factorials are scaled by abs(1/digamma(inital_estimate+1))
##'
##'    \item parameters in a log scale (ie log(theta)) are transformed
##'    by log(abs(initial_estimate))*abs(initial_estimate)
##'
##'    }
##'
##'    These parameter scaling coefficients are chose to try to keep
##'    similar slopes among parameters.  That is they all follow the
##'    slopes approximately on a log-scale.
##'
##'    While these are chosen in a logical manner, they may not always
##'    apply.  You can specify each parameters scaling factor by this
##'    parameter if you wish.
##'
##' @param scaleC0 Number to adjust the scaling factor by if the initial
##'     gradient is zero.
##'
##' @param scaleCmax Maximum value of the scaleC to prevent overflow.
##'
##' @param scaleCmin Minimum value of the scaleC to prevent underflow.
##'
##' @param normType This is the type of parameter
##'     normalization/scaling used to get the scaled initial valuse
##'     for nlmixr.  These are used with \code{scaleType} of.
##'
##'     With the exception of \code{rescale2}, these come
##'     from
##'     \href{https://en.wikipedia.org/wiki/Feature_scaling}{Feature
##'     Scaling}. The \code{rescale2} The rescaling is the same type
##'     described in the
##'     \href{http://apmonitor.com/me575/uploads/Main/optimization_book.pdf}{OptdesX}
##'     software manual.
##'
##'     In general, all all scaling formula can be described by:
##'
##'     v_{scaled} = (v_{unscaled}-C_{1})/C_{2}
##'
##'     Where
##'
##'
##'     The other data normalization approaches follow the following formula
##'
##'     v_{scaled} = (v_{unscaled}-C_{1})/C_{2};
##'
##' \itemize{
##'
##' \item \code{rescale2} This scales all parameters from (-1 to 1).
##'     The relative differences between the parameters are preserved
##'     with this approach and the constants are:
##'
##'     C_{1} = (max(all unscaled values)+min(all unscaled values))/2
##'
##'     C_{2} = (max(all unscaled values) - min(all unscaled values))/2
##'
##'
##' \item \code{rescale} or min-max normalization. This rescales all
##'     parmeters from (0 to 1).  As in the \code{rescale2} the
##'     relative differences are preserved.  In this approach:
##'
##'     C_{1} = min(all unscaled values)
##'
##'     C_{2} = max(all unscaled values) - min(all unscaled values)
##'
##'
##' \item \code{mean} or mean normalization.  This rescales to center
##'     the parameters around the mean but the parameters are from 0
##'     to 1.  In this approach:
##'
##'     C_{1} = mean(all unscaled values)
##'
##'     C_{2} = max(all unscaled values) - min(all unscaled values)
##'
##' \item \code{std} or standardization.  This standardizes by the mean
##'      and standard deviation.  In this approach:
##'
##'     C_{1} = mean(all unscaled values)
##'
##'     C_{2} = sd(all unscaled values)
##'
##' \item \code{len} or unit length scaling.  This scales the
##'    parameters to the unit length.  For this approach we use the Euclidean length, that
##'    is:
##'
##'     C_{1} = 0
##'
##'     C_{2} = sqrt(v_1^2 + v_2^2 + ... + v_n^2)
##'
##'
##' \item \code{constant} which does not perform data normalization. That is
##'
##'     C_{1} = 0
##'
##'     C_{2} = 1
##'
##' }
##'
##' @param gillStep When looking for the optimal forward difference
##'     step size, this is This is the step size to increase the
##'     initial estimate by.  So each iteration the new step size =
##'     (prior step size)*gillStep
##'
##' @param gillFtol The gillFtol is the gradient error tolerance that
##'     is accepable before issuing a warning/error about the gradient estimates.
##'
##' @param gillKcov The total number of possible steps to determine
##'     the optimal forward/central difference step size per parameter
##'     (by the Gill 1983 method) during the covariance step.  If 0,
##'     no optimal step size is determined.  Otherwise this is the
##'     optimal step size determined.
##'
##' @param gillStepCov When looking for the optimal forward difference
##'     step size, this is This is the step size to increase the
##'     initial estimate by.  So each iteration during the covariance
##'     step is equalt new step size = (prior step size)*gillStepCov
##'
##' @param gillFtolCov The gillFtol is the gradient error tolerance
##'     that is acceptable before issuing a warning/error about the
##'     gradient estimates during the covariance step.
##'
##' @param rmatNorm A parameter to normalize gradient step size by the
##'     parameter value during the calculation of the R matrix
##'
##' @param smatNorm A parameter to normalize gradient step size by the
##'     parameter value during the calculation of the S matrix
##'
##' @param covGillF Use the Gill calculated optimal Forward difference
##'     step size for the instead of the central difference step size
##'     during the central difference gradient calculation.
##'
##' @param optGillF Use the Gill calculated optimal Forward difference
##'     step size for the instead of the central difference step size
##'     during the central differences for optimization.
##'
##' @param covSmall The covSmall is the small number to compare
##'     covariance numbers before rejecting an estimate of the
##'     covariance as the final estimate (when comparing sandwich vs
##'     R/S matrix estimates of the covariance).  This number controls
##'     how small the variance is before the covariance matrix is
##'     rejected.
##'
##' @param adjLik In nlmixr, the objective function matches NONMEM's
##'     objective function, which removes a 2*pi constant from the
##'     likelihood calculation. If this is TRUE, the likelihood
##'     function is adjusted by this 2*pi factor.  When adjusted this
##'     number more closely matches the likelihood approximations of
##'     nlme, and SAS approximations.  Regardless of if this is turned
##'     on or off the objective function matches NONMEM's objective
##'     function.
##'
##' @param gradTrim The parameter to adjust the gradient to if the
##'     |gradient| is very large.
##'
##' @param gradCalcCentralSmall A small number that represents the value
##'     where |grad| < gradCalcCentralSmall where forward differences
##'     switch to central differences.
##'
##' @param gradCalcCentralLarge A large number that represents the value
##'     where |grad| > gradCalcCentralLarge where forward differences
##'     switch to central differences.
##'
##' @param etaNudge By default initial ETA estimates start at zero;
##'     Sometimes this doesn't optimize appropriately.  If this value
##'     is non-zero, when the n1qn1 optimization didn't perform
##'     appropriately, reset the Hessian, and nudge the ETA up by this
##'     value; If the ETA still doesn't move, nudge the ETA down by
##'     this value.  Finally if it doesn't move, reset it to zero and
##'     do not perform the optimization again.  This ETA nudge is only
##'     done on the first ETA optimization.
##'
##' @param maxOdeRecalc Maximum number of times to reduce the ODE
##'     tolerances and try to resolve the system if there was a bad
##'     ODE solve.
##'
##' @param odeRecalcFactor The factor to increase the rtol/atol with
##'     bad ODE solving.
##'
##' @param repeatGillMax If the tolerances were reduced when
##'     calculating the initial Gill differences, the Gill difference
##'     is repeated up to a maximum number of times defined by this
##'     parameter.
##'
##' @param stickyRecalcN The number of bad ODE solves before reducing
##'     the atol/rtol for the rest of the problem.
##'
##' @param nRetries If FOCEi doesn't fit with the current parameter
##'     estimates, randomly sample new parameter estimates and restart
##'     the problem.  This is similar to 'PsN' resampling.
##'
##' @inheritParams RxODE::rxSolve
##' @inheritParams minqa::bobyqa
##' @inheritParams foceiFit
##'
##' @details
##'
##' Note this uses the R's L-BFGS-B in \code{\link{optim}} for the
##' outer problem and the BFGS \code{\link[n1qn1]{n1qn1}} with that
##' allows restoring the prior individual Hessian (for faster
##' optimization speed).
##'
##' However the inner problem is not scaled.  Since most eta estimates
##' start near zero, scaling for these parameters do not make sense.
##'
##' This process of scaling can fix some ill conditioning for the
##' unscaled problem.  The covariance step is performed on the
##' unscaled problem, so the condition number of that matrix may not
##' be reflective of the scaled problem's condition-number.
##'
##' @author Matthew L. Fidler
##'
##' @seealso \code{\link{optim}}
##' @seealso \code{\link[n1qn1]{n1qn1}}
##' @seealso \code{\link[RxODE]{rxSolve}}
##' @export
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
