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

# thresh() and err.msg() --------------------------------------------------
thresh = function(x, cut=.Machine$double.xmin)
{
  x = abs(x)
	ifelse(x>cut, x, cut)
}

err.msg = function(x, pre="", post="")
{
  msg = paste0(x, collapse=", ")
  paste0(pre, msg, post)
}
# #########################################################################

# plot.dyn.ID() -----------------------------------------------------------
#' Plot of a non-population dynamic model fit
#'
#' Plot of a non-population dynamic model fit
#'
#' @param x a dynamodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
plot.dyn.ID = gof
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
# #########################################################################

# print.dyn.ID() ----------------------------------------------------------
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
# #########################################################################

# summary.dyn.ID() --------------------------------------------------------
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
# #########################################################################

# nmsimplex() and mymin() -------------------------------------------------
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

# Redundant function?
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

# #########################################################################

# nlmixrDynmodelConvert ---------------------------------------------------
#' Converting nlmixr Objects to dynmodel Objects
#' 
#' @param .nmf nlmixr object
#' @return list containing inputs for the dynmodel()
#' \itemize{
#' \item \code{$fixPars} - fixed parameters defined as \code{fixed()} in the nlmixr object
#' \item \code{$sigma} - error model parameters 
#' \item \code{$inits} - initial estimates for parameters in the model
#' \item \code{$lower} - lower boundaries for estimated parameters
#' \item \code{$upper} - upper boundaries for estimated parameters
#' \item \code{$system} - RxODE object that defines the structural model
#' \item \code{$model} - error model
#' }
#' 
#' @export


nlmixrDynmodelConvert <- function(.nmf){
  # Notes:
  # Description - function used to convert nlmixr() object to initial conditions and model used for dynmodel
  # input - nlmixr function (.nmf)
  # output - .return
  # TO DO:
  # add lower and upper outputs for inits
  
  # iniitalize list for output
  .return <- list()
  
  # convert nlmixr model to data.frame
  .nmf.original <- .nmf
  .nmf <- as.data.frame(.nmf$ini)
  
  # assign fixed terms
  .fix.index <- if (length(which(.nmf$fix==TRUE))==0) {NULL} else {which(.nmf$fix==TRUE)} # obtain row location for fixed terms
  .fixPars <- if (length(.nmf$est[.fix.index])==0) {NULL} else {exp(.nmf$est[.fix.index])}
  names(.fixPars) <- substring(.nmf$name[.fix.index],2)
  .return <- c(.return,fixPars=list(.fixPars))
  
  # assign theta terms(estimated terms excluding error terms)
  .theta.index <- 
    if (is.null(.fix.index)){
      which(!is.na(.nmf$ntheta) & is.na(.nmf$err),TRUE)
    } else {
      which(!is.na(.nmf$ntheta) & is.na(.nmf$err),TRUE)[-which(.nmf$fix == TRUE)]
    } # row location for theta values
  
  
  .theta <- exp(.nmf$est[.theta.index]) # initial estimate for theta values, back-transformed
  names(.theta) <- substring(.nmf$name[.theta.index],2)
  
  # assign sigma terms (estimated)
  .sigma.index <- 
    if (is.null(.fix.index)){
      which(!is.na(.nmf$ntheta) & !is.na(.nmf$err),TRUE)
    } else {
      which(!is.na(.nmf$ntheta) & !is.na(.nmf$err),TRUE)[-which(.nmf$fix == TRUE)]
    } # row location for theta values
  
  .sigma <- .nmf$est[.sigma.index] # initial estimate for theta values, back-transformed
  
  names(.sigma) <- .nmf$err[.sigma.index]
  .return <- c(.return,sigma=list(.sigma))
  
  # assign "inits", vector of theta and sigma terms # (will be used when the likelihood function is changed)
  .inits <- c(.theta)#,.sigma)
  .return$inits <- .inits
  
  # assign boundaries
  .lower <- exp(.nmf[,5][!is.na(.nmf["ntheta"]) & .nmf["fix"]==FALSE & is.na(.nmf["err"])]) # theta terms
  .lower <- c(.lower, .nmf[,5][!is.na(.nmf["ntheta"]) & .nmf["fix"]==FALSE & !is.na(.nmf["err"])]) # error terms
  
  .upper <- exp(.nmf[,7][!is.na(.nmf["ntheta"]) & .nmf["fix"]==FALSE & is.na(.nmf["err"])]) # theta terms
  .upper <- c(.upper, .nmf[,7][!is.na(.nmf["ntheta"]) & .nmf["fix"]==FALSE & !is.na(.nmf["err"])]) # error terms
  
  .return <- c(.return, lower = list(.lower), upper = list(.upper))
  
  # obtain system
  .system <- RxODE(.nmf.original$rxode.pred) # (use nlmixr_prde)
  .system$stateExtra <- NULL # remove extraState, the error model term shoudl not be inclcuded
  .system$lhs <- .system$lhs[-length(.system$lhs)] # remove the error model term
  .return <- c(.return,system=.system)
  
  # create error model
  .DV <- .nmf$condition[!is.na(.nmf$condition) & .nmf$condition != "ID"]
  .PRED <- "nlmixr_pred" # need to obtain from data? id dont know
  
  .formula <- list()
  for(i in 1:length(.sigma.index)) {
    if (i == 1) {
      .temp <- paste(.nmf$err[.sigma.index[i]],"(", .nmf$est[.sigma.index[i]], ")",sep="")
    } else {
      .temp <- paste("+ ",.nmf$err[.sigma.index[i]],"(", .nmf$est[.sigma.index[i]], ")",sep="")
    }
    .formula <- paste(.formula,.temp)
  }
  
  .model <- as.formula(paste(.DV,"~",.PRED,"+", .formula))
  .return <- c(.return,model=.model)
  
  # Output
  return(.return)
}
# #########################################################################

# dynmodelControl() -----------------------------------------------
#' Control Options for dynmodel
#' 
#' @inheritParams RxODE::rxSolve
#' @inheritParams foceiControl
#' @param ... other arguments apply to dynmodelControl
#' @export

# run devtools::document() to update documentation
# go through and make sure
# devtools::check_man() used to identify missing and problems

dynmodelControl <- function(...,
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
  scaleC=NULL,
  scaleC0=1e5,
  # RxODE
  atol=1e-08,
  rtol=1e-06,
  # bobyqaControl
  npt = NULL,
  rhobeg = 0.2,
  rhoend = 1E-04,
  iprint = 0,
  maxfun = NULL,
  # lbfgsb3c
  trace=0,
  factr=NULL,
  pgtol=NULL,
  abstol=NULL,
  reltol=NULL,
  lmm=NULL,
  maxit=NULL,
    #,iprint=NULL repreated above
  # nlminb (PORT)
  eval.max=NULL,
  iter.max=NULL,
    #trace=NULL,
  abs.tol=NULL,
  rel.tol=NULL,
  x.tol=NULL,
  xf.tol=NULL,
  step.min=NULL,
  step.max=NULL,
  sing.tol=NULL,
  scale.init=NULL,
  diff.g=NULL
  ) {
  if (missing(method)){method = "bobyqa"}
  if (missing(normType)){normType = "constant"}
  if (missing(scaleType)){scaleType = "norm"}
  
  .ret <- list(
    fixPars=fixPars,
    lower=lower,
    upper=upper,
    method=method,
    ftol_rel=ftol_rel,
    maxeval=maxeval,
    scaleTo=scaleTo,
    scaleObjective=scaleObjective,
    normType=normType,
    scaleType=scaleType,
    scaleCmax=scaleCmax,
    scaleCmin=scaleCmin,
    scaleC=scaleC, #NULL,
    scaleC0=scaleC0,
    # RxODE
    atol=atol,
    rtol=rtol,
    # bobyqa
    npt = npt,
    rhobeg = rhobeg,
    rhoend = rhoend,
    iprint = iprint, #also used in lbfgsb3c
    maxfun = maxfun,
    # lbfgsb3c
    trace=trace,
    factr=factr,
    pgtol=pgtol,
    abstol=abstol,
    reltol=reltol,
    lmm=lmm,
    maxit=maxit,
      #,iprint = iprint,
    # nlminb (PORT)
    eval.max=eval.max,
    iter.max=iter.max,
      #trace=NULL,
    abs.tol=abs.tol,
    rel.tol=rel.tol,
    x.tol=x.tol,
    xf.tol=xf.tol,
    step.min=step.min,
    step.max=step.max,
    sing.tol=sing.tol,
    scale.init=scale.init,
    diff.g=diff.g
  )
  
  class(.ret) <- "dynmodelControl"
  return(.ret)
  
}

# #########################################################################

# dynmodel()  ---------------------------------------------------------------
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

dynmodel = function(system, model, evTable, inits, data, control=list(), ...){
  # dynmodelControl Handling ------------------------------------------------
  if (!RxODE::rxIs(control, "dynmodelControl")){
    control <- do.call(dynmodelControl, control)
  }
  
  # reassign contorl names
  for (i in 1:length(control)){
    assign(names(control[i]),control[[i]])
  }
  
  # Error model  -------------------------------------------------------------
   if (class(model)=="formula") {
    model = list(model)
  }
  inits.err = NULL
  model = lapply(model, function(.model) {
    .model = unlist(lapply(attr(terms(.model),"variables"), as.list))
    .model = sapply(.model, deparse)
    

    # assign error terms    
    .sigma.add = if("add" %in% .model) {
      as.numeric(.model[which(.model=="add")+1])
    } else{NULL}
    
    .sigma.prop = if("prop" %in% .model) {
      as.numeric(.model[which(.model=="prop")+1])
    } else{NULL}
    
    .sigma.pow = if("pow" %in% .model) {
      as.numeric(.model[which(.model=="pow")+1])
    } else{NULL}
    
    .sigma.pow2 = if("pow2" %in% .model) {
      as.numeric(.model[which(.model=="pow2")+1])
    } else{NULL}
    
    .sigma.yeoJohnson = if("yeoJohnson" %in% .model) {
      as.numeric(.model[which(.model=="yeoJohnson")+1])
    } else{NULL}
    
    .sigma.boxCox = if("boxCox" %in% .model) {
      as.numeric(.model[which(.model=="boxCox")+1])
    } else{NULL}

    .sigma.norm = if("norm" %in% .model) {
      as.numeric(.model[which(.model=="norm")+1])
    } else{NULL}
    
    .sigma.dnorm = if("dnorm" %in% .model) {
      as.numeric(.model[which(.model=="dnorm")+1])
    } else{NULL}
    
    .sigma.logn = if("logn" %in% .model) {
      as.numeric(.model[which(.model=="logn")+1])
      .logn <<- T
    } else{NULL
      .logn <<- NULL}
    
    .sigma.dlnorm = if("dlnorm" %in% .model) {
      as.numeric(.model[which(.model=="dlnorm")+1])
      .dlnorm <<- T
    } else{NULL
      .dlnorm <<- NULL}
    
    .sigma.tbs = if("tbs" %in% .model) {
      as.numeric(.model[which(.model=="tbs")+1])
    } else{NULL}
    
    .sigma.tbsYj = if("tbsYj" %in% .model) {
      as.numeric(.model[which(.model=="tbsYj")+1])
    } else{NULL}
    
    # keep error model terms
    inits.err <- c(add=.sigma.add, prop=.sigma.prop, pow=.sigma.pow, pow2=.sigma.pow2, yeoJohnson=.sigma.yeoJohnson, boxCox=.sigma.boxCox,
                   norm=.sigma.norm, dnorm=.sigma.dnorm, tbs=.sigma.tbs, tbsYj=.sigma.tbs)
    
    inits.err <- inits.err[which(names(inits.err) %in% intersect(names(inits.err),.model))]
    inits.err <<- inits.err
    .model <- c("dv" = .model[2], "pred" = .model[3], inits.err)
    
    # error message for using power function
  })
  
  inits = c(inits, inits.err)
  if("pow2" %in% names(inits) & !("pow" %in% names(inits))){stop("Error Model: pow must be defined when using pow2")}


  # Check dynmodel() inputs, Define vars, modelVars, pars,  ------------
    # Check to make sure all there is consistency between error model, data. inits, and ODE model
    
    # Error "model" contains "data" variables?
    # get column names of data (Time and Observation)
    vars = names(data)
    # check to see if there is a discrepency between error model names and data
    nodef = setdiff(sapply(model, function(x) x["dv"]), vars)
    # print error message
    if (length(nodef)) {
      msg = err.msg(nodef, pre="var(s) not found in data: ")
      stop(msg)
    }
    
    # "system" variables contain error "model" variables?
    # obtain all variables from the system 
    modelVars = system$cmpMgr$get.modelVars()
    # reassign vars to combine state and lhs variables
    vars = c(modelVars$state, modelVars$lhs)
    # Check to see if the prediction term is in the error model
    nodef = setdiff(sapply(model, function(x) x["pred"]), vars)
    # print error message
    if (length(nodef)) {
      msg = err.msg(nodef, pre="modelVar(s) not found in model: ")
      stop(msg)
    }
    
    #  "system" variables contain estimated "init" variables and fixed "fixPars" variables?
    # obtain fixed and estimated parameters
    pars = modelVars$params
    # Check to see if there are values in pars, that are not in the initial conditions and fixed parameters
    nodef = setdiff(pars, c(names(inits), names(fixPars)))
    # print error message
    if (length(nodef)) {
      msg = err.msg(nodef, pre="par(s) not found: ")
      stop(msg)
    }
  
  # Additional assignment ---------------------------------------------------
  # number of estimated parameters, excluding the error terms
  npar = length(pars) - length(fixPars)
  
  # if the time of the observed "data" starts at zero, rows = T, else rows = -1 ??????
  #have_zero = min(data$time) <= 0
  #rows = if(have_zero) T else -1
  
  
  # Objective Function ------------------------------------------------------
  obj = function(th)
  {
    # unscale
    unscaled.th <- numeric(length(th))
    for (i in 1:length(th)) {th[i] <- unscalePar(th,i)}
    
    # define parameters used for simulation, all parameters except the error terms
    .ixpar = npar
    theta = th[1:npar] 
    names(theta) = names(inits)[1:npar]
    theta = c(theta, fixPars)
    
    # call rxODE for simulation
    rxControl <- c(atol=atol, rtol=rtol)
    s = do.call(RxODE :: rxSolve, c(list(object=system, params=theta, events=evTable), rxControl))
    
    # NOTES:
    #returnType="data.frame" - add to rxControl
    #rxNorm(system) # add error piece, use the error parameters as parameters
    #nlmixr_err=boxCox(add.err^2+prop.err^(2*pow)*pred, lambda) for the predictions not the errors
    
    # sum of log-likelihood function:
    l = lapply(model, function(x) {
      # name the inputs
      names(th) <- names(inits)
      
      # initialize sigmas for objective function
      add <- 0
      prop <- 0
      pow <- 1
      pow2 <- 1
      norm <- NULL
      dnorm <- NULL
      boxCox <- NULL
      tbs <- NULL
      yeoJohnson <- NULL
      tbsYJ <- NULL

      # assign names sigma names
      if(any(names(th) %in% names(model[[1]]))){
        for(i in 1:sum((names(th) %in% names(model[[1]])))) {
          assign(names(th[names(th) %in% names(model[[1]])])[i], as.numeric(th[names(th) %in% names(model[[1]])])[i])
        } 
      }
      
      if (!is.null(norm)) add <- norm
      if (!is.null(dnorm)) add <- dnorm
      if (!is.null(.logn)) {
        lambda <- 0
        th <- th[names(th)!="logn"]
      }
      if (!is.null(dlnorm)) {
        lambda <- 0
        th <- th[names(th)!="dlnorm"]
      }
      if (!is.null(boxCox)) lambda <- boxCox
      if (!is.null(tbs)) lambda <- tbs
      if (!is.null(yeoJohnson)) lambda <- yeoJohnson
      if (!is.null(tbsYJ)) lambda <- tbsYJ

     # predictted and observed values from RxODE
      #yp = s[rows,x["pred"]]
      yp = s[,x["pred"]]
      yo = data[, x["dv"]]
      
      # Transform both sides (used for non-normal residuals) ######### may need to implement nlmixr::coxBox and nlmixr::yeoJohnson
      boxCox = function (x, lambda) {
        if(lambda == 0) {
          .h.x <- log(x)
        }
        else {
          .h.x <- (x^lambda-1)/lambda
        }
        return(.h.x)
      }
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
      
      # log normal transformation 
      if (!is.null(.logn) | !is.null(.dlnorm)){
        .h.x <- log(yo) #boxCox(yo, lambda) # obs
        .h.y <- log(yp) #boxCox(yp, lambda) # pred
        
        if("pow" %in% names(model[[1]])) {
          .h.y.var <- yp^(2*pow2)*thresh(pow)^2 + thresh(add)^2  # variance of pred
        } else {
          .h.y.var <- yp^(2*pow2)*thresh(prop)^2 + thresh(add)^2  # variance of pred
        }
        
        # boxCox transformed -2 log-likelihood
        .logn.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
        
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll = .logn.n2ll - 2*(0-1)*log(yo) -2*log(2*pi) # lambda is zero here
        
        # negative log-likelihood function for output
        ll = .5*(.n2ll)
      }
      # boxCox Transform
      else if ("boxCox" %in% names(model[[1]]) | "tbs" %in% names(model[[1]])) {
        
        .h.x <- boxCox(yo, lambda) # obs
        .h.y <- boxCox(yp, lambda) # pred
        
        if("pow" %in% names(model[[1]])) {
          .h.y.var <- yp^(2*pow2)*thresh(pow)^2 + thresh(add)^2  # variance of pred
        } else {
          .h.y.var <- yp^(2*pow2)*thresh(prop)^2 + thresh(add)^2  # variance of pred
        }
        
        # boxCox transformed -2 log-likelihood
        .boxCox.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
        
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll = .boxCox.n2ll - 2*(lambda-1)*log(yo) -2*log(2*pi)
        
        # negative log-likelihood function for output
        ll = .5*(.n2ll)
      }
      # yeoJohnson Transform
      else if("yeoJohnson" %in% names(model[[1]]) | "tbsYJ" %in% names(model[[1]])) {
        
        .h.x <- yeo.johnson(yo, lambda) #obs
        .h.y <- yeo.johnson(yp, lambda) #pred

        # assign correct variance according to defined error model
        if (!any("prop" %in% names(model[[1]]))  & !any("pow" %in% names(model[[1]])) ) {
          .h.y.var <- 1
print('Here1')
        } else if (any("prop" %in% names(model[[1]]))  & !any("pow" %in% names(model[[1]])) ) {
          .h.y.var <- yp^(2*pow2)*thresh(prop)^2 + thresh(add)^2  # variance of pred

          
          cat("pow2",pow2)
print('Here2')

        } else {
          .h.y.var <- yp^(2*pow2)*thresh(pow)^2 + thresh(add)^2  # variance of pred
print('Here3')

        }

        # yeoJohnson transformed -2 log-likelihood
        .yeoJohnson.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
        
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll <- ifelse(yo >= 0, 
                        .yeoJohnson.n2ll -2*(lambda-1)*log(yo+1) -2*log(2*pi),
                        .yeoJohnson.n2ll -2*(1-lambda)*log(-yo+1) -2*log(2*pi)
                        
                        
        )
        if (yo>=0) print("y0>=0") else (print("yo<0"))

        # negative log-likelihood function for output
        ll = .5*(.n2ll)
      }
      # power model
      else if ("pow2" %in% names(model[[1]])) {
        sgy = thresh(add) + thresh(pow)*yp^(pow2)
        ll = .5*((yo - yp)^2/(sgy^2) + log(sgy^2) + log(2*pi))
      }
      # all other error models
      else {
        if (!any("add" %in% names(model[[1]])) & !any("prop" %in% names(model[[1]]))){
          sgy = 1  
        }else{
          sgy = thresh(add) + thresh(prop)*yp
        }
        ll = .5*((yo - yp)^2/(sgy^2) + log(sgy^2) + log(2*pi))
      }
      sum(ll)
    })
    do.call("sum", l)  # same as return(as.numeric(l)), l is a list for each value in the model?
  }

  
  # Scaling functions -----------------------------------------------------------------------
  # normType assignment for scaling (normalization type)
  if (normType == "constant") {
    C1 = 0
    C2 = 1
  } else if (normType == "rescale2") {
    C1 = (max(inits) + min(inits))/2
    C2 = (max(inits) - min(inits))/2
  } else if (normType == "mean") {
    C1 = mean(inits)
    C2 = max(inits) - min(inits)
  } else if (normType == "rescale") {
    C1 = min(inits)
    C2 = max(inits) - min(inits)
  } else if (normType == "std") {
    C1 = mean(inits)
    C2 = sd(inits)
  } else if (normType == "len") {
    C1 = 0
    C2 = sqrt(sum(inits*inits))
  }
  
  # scaleC assignment for scaling (adopted from foceFIT.R)
  if (is.null(scaleC) | length(scaleC) < length(inits)){
    scaleC <- rep(1, length(inits))
  } else {
    scaleC <- as.double(scaleC);
    if (length(inits) > length(scaleC)){
      scaleC <- c(scaleC, rep(1, length(inits) - length(scaleC)));
    } else if (length(inits) < length(scaleC)){
      scaleC <- scaleC[seq(1, length(inits))];
      warning("scaleC control option has more options than estimated parameters, please check.")
    }
  }	
  
  # Function for scaling parameters based on scaleType
  scalePar <- function(x, i){
    if (scaleType == "norm") { # simple scaling
      return((x[i]-C1)/C2)
    }	else if (scaleType == "nlmixr") { 	# nlmixr
      scaleTo = (inits[i] - C1)/C2
      return((x[i] - inits[i])/scaleC[i] + scaleTo)
    }	else if (scaleType == "mult") { # simple multiplicatice scaling
      if (scaleTo > 0) {
        return(x[i]/inits[i]*scaleTo)
      } else {
        return(x[i])
      }
    } else if (scaleType == "multAdd") { # log non-log multiplicative scaling
      if(scaleTo > 0) {
        return((x[i]-inits[i]) + scaleTo)
      } else {
        return(x[i]/inits[i]*scaleTo)
      }
    } else { # When should this be used? "norm" scaling is essentially no scaling when normType specified.
      if(scaleTo > 0) {
        return((x[i] - inits[i]) + scaleTo)
      } else {
        return(x[i])
      }
    }
  }
  
  # Function for unscaling parameters based on scaleType
  unscalePar <- function(x, i){
    if (scaleType == "norm") { # simple scaling
      return(x[i]*C2+C1);
    } else if (scaleType == "nlmixr") { 	  # nlmixr
      scaleTo = (inits[i] - C1)/C2
      return((x[i] - scaleTo)*scaleC[i] + inits[i])
    } else if (scaleType == "mult") { 	  # simple multiplicatice scaling
      if (scaleTo > 0) {
        return(x[i]*inits[i]/scaleTo)
      } else {
        return(x[i])
      }
    } else if (scaleType == "multAdd") { # log non-log multiplicative scaling
      if(scaleTo > 0) {
        return((x[i]- scaleTo) + inits[i])
      } else {
        return(x[i]*inits[i]/scaleTo)
      }
    } else { # When should this be used? "norm" scaling is essentially no scaling when normType specified.
      if(scaleTo > 0) {
        return((x[i]-scaleTo)*1 + inits[i])
      } else {
        return(x[i])
      }
    }
  }
  
  # Scale --------------
  if(normType != "constant" & scaleType != "norm"){
    # scaled the initial conditions
    inits.temp <- numeric(length(inits))
    for(i in 1:length(inits)){inits.temp[i] <- scalePar(inits,i)}
    .inits <- inits.temp
    # scaled lower boundary
    if (is.null(lower)==FALSE){
      if (any(lower == -Inf)) {
        lower[which(lower==-Inf)]=0
        warning("Lower boundary of -Inf set to 0.")
      } 
      lower.temp <- numeric(length(lower))
      for(i in 1:length(lower)){lower.temp[i] <- scalePar(lower,i)}
      .lower <- lower.temp
    } else {.lower <- NULL}
    # scaled upper boundary
    if (is.null(upper)==FALSE){
      upper.temp <- numeric(length(upper))
      for(i in 1:length(upper)){upper.temp[i] <- scalePar(upper,i)}
      .upper <- upper.temp
    } else {.upper <- NULL}
  } else {
    .inits <- inits
    .lower <- lower
    .upper <- upper
  }
  
  # Optimization -----------------------------------------------------------------------
  if (method =="bobyqa") {
    # remove any control elements not required in bobyqa()
    if (is.null(control$npt)) control$npt <- c(npt = length(.inits)*2 + 1)
    .control <- control[names(control) %in% c("npt", "rhobeg", "rhoend", "iprint")] #"maxfun"
    # Run bobyqa optimization
    fit = minqa::bobyqa(par=.inits, fn=obj, lower=.lower, upper=.upper, control=.control)
    fit$value <- fit$fval
  } else if (method=="lbfgsb3c") {
    # remove any control elements not required in lbfgsb3c()
    #lbfgsb3c.control <- list(trace=-1,factr=1e7,pgtol=0,abstol=0,reltol=0,lmm=5,maxit=99999,iprint=-1)
    #lbfgsb3c.control  <- c("trace","factr","pgtol","abstol","reltol","lmm","maxit","iprint")
    #if (any(names(lbfgsb3c.control) %in% names(control))) {
    #   .control = control[match(names(lbfgsb3c.control),names(control))[!is.na(match(names(lbfgsb3c.control),names(control)))]]
    # } else{.control=NULL}
    #.control <- control[names(control) %in% c("trace","factr","pgtol","abstol","reltol","lmm","maxit","iprint")] #"maxfun"
    
    # Run lbfbsb3c optimization
    fit = lbfgsb3c::lbfgsb3c(par = as.vector(.inits), fn=obj, lower=.lower, upper=.upper, gr=NULL)#, control=.control, gr=NULL)
    
  } else if (method=="Nelder-Mead") {
    # Run Nelder-Mead optimization
    if (.lower > 0 | .upper >0){warning("Optimization: Boundaries not used in Nelder-Mead")}
    fit = mymin(as.vector(.inits), obj, control=control)
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
    # Run nlminb (PORT) optimization
    nlminb.control <- list(eval.max=NULL,iter.max=NULL,trace=NULL,abs.tol=NULL,rel.tol=NULL,x.tol=NULL,
                           xf.tol=NULL,step.min=NULL, step.max=NULL,sing.tol=NULL,scale.init=NULL,diff.g=NULL)
    if (any(names(nlminb.control) %in% names(control))) {
      .control = control[match(names(nlminb.control),names(control))[!is.na(match(names(nlminb.control),names(control)))]]
    } else{.control=NULL}
    fit = nlminb(start = as.vector(.inits), objective = obj, gradient = NULL, hessian = NULL, 
                 scale = 1, control=.control, lower = .lower, upper = .upper)
  }
  
  # Hessian -----------------------------------------------------------------------
    fit$hessian = try(optimHess(fit$par, obj, control=control) , silent=TRUE)
    
    if(inherits(fit$hessian,"try-error")){
      se = rep(NA, length(fit$par))
      warning("standard error of the Hessian has failed")
    } else {
      se = sqrt(diag(solve(fit$hessian)))
    }
    
  # reassign the negative values to positive for add, prop/pow since they are standard deviations
    if (!is.na(match("add",names(inits)))) fit$par[match("add",names(inits))] = abs(fit$par[match("add",names(inits))])
    if (!is.na(match("prop",names(inits)))) fit$par[match("prop",names(inits))] = abs(fit$par[match("prop",names(inits))])
    if (!is.na(match("pow",names(inits)))) fit$par[match("pow",names(inits))] = abs(fit$par[match("pow",names(inits))])

 
  
  # dynmodel() Output -------------------------------------------------------
  # unscale optmized parameters here if scaling was used:
  if(normType != "constant" & scaleType != "norm"){
    par.temp <- numeric(length(fit$par))
    for (i in 1:length(par.temp)) {par.temp[i] <- unscalePar(fit$par,i)}
    fit$par <- par.temp
  }
  
  # create table for output
  res = cbind(fit$par, abs(se), abs(se/fit$par*100))
  dimnames(res) = list(names(inits), c("est", "se", "%cv"))
  
  # ??
  nobs = 0
  l = lapply(model, function(x) {
    yo =  data[,model[[1]]["dv"][[1]]]#data[, x["dv"]]
    nobs <<- nobs + length(yo)
  })
  if (!is.null(fit$objective)) fit$value = fit$objective
  
  # Output
  res = c(list(res=res, obj=obj, npar=length(fit$par), nobs=nobs, data=data), fit)
  class(res) = "dyn.ID"
  res

}

# #########################################################################





# ####################################################################### #
#
## MCMC Section
#
# ####################################################################### #
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

# dynmodel.mcmc() ---------------------------------------------------------
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
# #########################################################################


# ####################################################################### #
#
## Utilities for building nlmixr
#
# ####################################################################### #

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


# .collectWarnings --------------------------------------------------------
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
# #########################################################################

# nlmixrPrint() -----------------------------------------------------------
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
# #########################################################################

.dontRun <- function(...){
    ## This is for r checks, though they need to be loaded...
    vpc::vpc(...)
    dparser::dparse(...)
}

# cholSE() ----------------------------------------------------------------
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
# #########################################################################

.setRoot <- function(){
    setwd("c:/");
}
