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

#' @export
#' @rdname gof
plot.dyn.ID = gof
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

# as.focei.dynmodel() -----------------------------------------------------------
#' Output nlmixr format for dynmodel
#'
#' @param .dynmodelObject return object from
#' @export

# run devtools::document() to update documentation
# go through and make sure
# devtools::check_man() used to identify missing and problems
# Will outside of dynmodel, and called within dynmodel. If nlmixr input used, output nlmixr, if dynmodel input used, output dynmodel.

as.focei.dynmodel <- function(.dynmodelObject, .nlmixrObject, .data, .time, .theta, .fit, .message, .inits.err, .cov, .sgy, .dynmodelControl, .nobs2=0, .pt=proc.time(), .rxControl = RxODE::rxControl()){

  # setup ----
  .env <- new.env(parent=emptyenv()); # store information for attaching to fit
  .model <- RxODE::rxSymPySetupPred(.nlmixrObject$rxode.pred,
                                    function(){return(nlmixr_pred)},
                                    .nlmixrObject$theta.pars,
                                    .nlmixrObject$error,
                                    grad=FALSE,
                                    pred.minus.dv=TRUE, sum.prod=FALSE, #control$sumProd,
                                    theta.derivs=FALSE, optExpression=TRUE, #control$optExpression,
                                    run.internal=TRUE, only.numeric=TRUE)

  .nlmixrObject.df <- as.data.frame(.nlmixrObject$ini)
  # ####

  ## .parFixedDf ----
  .ci <- .dynmodelControl$ci
  .Estimates <- .dynmodelObject$res[,1]

  backTransformed <- function(.model){
    .temp.model <- .model
    # assign fixed terms
    .fix.index <- if (length(which(.nlmixrObject.df$fix==TRUE))==0) {NULL} else {which( .nlmixrObject.df$fix==TRUE)} # obtain row location for fixed terms
    .ref.fix <- substring(.nlmixrObject.df$name[.fix.index],2)

    .temp.log.fixPars.index <- intersect(.fix.index,.temp.model$log.etas)
    .temp.log.fixPars <-   .nlmixrObject.df $est[.temp.log.fixPars.index]
    names(.temp.log.fixPars) <- substring(  .nlmixrObject.df$name[.temp.log.fixPars.index],2)

    .temp.nonlog.fixPars.index <- setdiff(.fix.index,.temp.model$log.etas)
    .temp.nonlog.fixPars <-   .nlmixrObject.df$est[.temp.nonlog.fixPars.index]
    names( .temp.nonlog.fixPars ) <- substring(  .nlmixrObject.df$name[.temp.nonlog.fixPars.index],2)

    .fixPars <- c(.temp.log.fixPars,.temp.nonlog.fixPars)
    .fixPars <- .fixPars[order(factor(names(.fixPars), levels=.ref.fix))]

    # assign theta terms(estimated terms excluding error terms)
    .theta.index <<-
      if (is.null(.fix.index)){
        which(!is.na(  .nlmixrObject.df$ntheta) & is.na(  .nlmixrObject.df$err),TRUE)
      } else {
        which(!is.na(  .nlmixrObject.df$ntheta) & is.na(  .nlmixrObject.df$err),TRUE)[-which(  .nlmixrObject.df$fix == TRUE)]
      } # row location for theta values

    #.ref.theta <-   .nlmixrObject.df$name[.theta.index]
    #.temp.log.theta.index <- intersect(.theta.index,.temp.model$log.etas)

    return(.theta.index)
  }

  .Back.Transformed.Y.N <- (names(.Estimates) %in% .Estimates[backTransformed(.model)])
  .Back.Transformed <- ifelse(.Back.Transformed.Y.N, exp(.Estimates), .Estimates)
  .SE <- .dynmodelObject$res[,2]
  .RSE <- .dynmodelObject$res[,3]
  .z.score <- qnorm((1-.ci)/2, lower.tail = FALSE)
  .ci.lower <- ifelse(.Back.Transformed.Y.N, exp(.dynmodelObject$res[,1] - .z.score*.dynmodelObject$res[,2]), .dynmodelObject$res[,1] - .z.score*.dynmodelObject$res[,2])
  .ci.upper <- ifelse(.Back.Transformed.Y.N, exp(.dynmodelObject$res[,1] + .z.score*.dynmodelObject$res[,2]), .dynmodelObject$res[,1] + .z.score*.dynmodelObject$res[,2])

  .parFixedDf <- data.frame(Estimate = .Estimates, SE = .SE, RSE = .RSE,
                            Back.Transformed = .Back.Transformed, ci.lower = .ci.lower, ci.upper = .ci.upper)

  .env$parFixedDf <- .parFixedDf #as.list(.env)
  # ####

  ## .parFixed ----
  .digs <- .dynmodelControl$digs

  .Back.Transformed.label <-  paste0("Back-Transformed(",.ci*100,"%CI)")
  .Back.Transformed.df <- paste0(
    formatC(signif(.Back.Transformed, digits=.digs), digits=.digs, format="fg", flag="#"),
    rep(" (", length(.Back.Transformed)),
    formatC(signif(.ci.lower, digits=.digs), digits=.digs, format="fg", flag="#"),
    rep(", ", length(.Back.Transformed)),
    formatC(signif(.ci.upper, digits=.digs), digits=.digs, format="fg", flag="#"),
    rep(")", length(.Back.Transformed))
  )

  # assign labels to parFixed
  .parameters.df <- .nlmixrObject.df$label[.theta.index]
  if (any(is.na(.parameters.df))) {
    .parameters.df[is.na(.parameters.df)] <- ""
  } else {
    .parameters.df <- .nlmixrObject.df$label[.theta.index]
  }
  if (length(.parameters.df)<length(.Estimates)){
    .parameters.df <- c(.parameters.df, rep("", sum(!is.na(.nlmixrObject.df$err))))
  }



  .parFixed <- data.frame(
    .parameters.df,
    formatC(signif(.Estimates, digits=.digs), digits=.digs, format="fg", flag="#"),
    formatC(signif(.SE, digits=.digs), digits=.digs, format="fg", flag="#"),
    formatC(signif(.RSE, digits=.digs), digits=.digs, format="fg", flag="#"),
    .Back.Transformed.df)

  names(.parFixed) <- c("Parameter", "Est.", "SE", "%RSE", .Back.Transformed.label)

  .env$parFixed <- .parFixed
  # ####

  ## $covMethod ----
  .covMethod <- .dynmodelControl$covMethod;
  .env$covMethod <- .covMethod
  # ####

  ## $fixef ----
  .fixef <- .dynmodelObject$res[,1]
  .env$fixef <- .fixef
  # ####

  ## $nlmixrObject ----
  .temp.theta.index <- c(1:sum(as.data.frame(.nlmixrObject$ini)$fix == F & is.na(as.data.frame(.nlmixrObject$ini)$err)))
  .temp.replacements <- .dynmodelObject$res[.theta.index,1]
  ini <- as.data.frame(.nlmixrObject$ini)
  ini$est[.temp.theta.index] <- c(.temp.replacements)
  class(ini) <- c("nlmixrBounds", "data.frame")
  .nlmixrObject$ini <- ini

  .env$uif <- .nlmixrObject
  # ####

  ## $dynmodelObject ----
  .env$dynmodelObject <- .dynmodelObject
  # ####

  ## $cov ----
  labels <- names(.Estimates)
  dimnames(.cov) <- list(labels, labels)
  .env$cov <- .cov
  # ####

  # ####

  ## $nobs ----
  .env$nobs <- .dynmodelObject$nobs
  # ####

  ## $model ----
  .env$model <- .model
  # ####

  ## $objf ----
  .objf <-  .dynmodelObject$value-0.5*.dynmodelObject$nobs*log(2*pi)
  .env$objf <- .objf
  .env$objective <- .objf
  # ####


  ## $logLik ----
  ## returns value = -ll
  logLik <- -.dynmodelObject$value
  attr(logLik, "df") = .dynmodelObject$npar
  attr(logLik,"nobs") = .dynmodelObject$nobs
  class(logLik) = "logLik"
  .env$logLik <- logLik

  ## $objDf ----
  .aic <- AIC(logLik)
  .bic <- BIC(logLik)

  .objDf <- data.frame(OBJF = .objf, AIC = .aic, BIC = .bic, "Log-likelihood"=.dynmodelObject$value,
                       check.names=FALSE)
  .env$objDf <- .objDf
  # ####

  ## Required checks
  ## plot(obj)
  ## setOfv(obj) ## fail
  ## ranef(obj) ## fail

  ## omega ----
  .env$omega <- matrix(numeric(0), 0, 0)
  .env$omegaR <- matrix(numeric(0), 0, 0)
  # ####

  ## sigma ----
  .sigma <- .inits.err
  .env$sigma <- if (length(.sigma)>1) diag(.sigma) else .sigma
  # ####

  ## message ----
  .env$message <- .message
  # ####

  ## method ----
  .env$method <- "dynmodel"
  # ####

  ## extra ----
  .estimation.method <- .dynmodelControl$method
  .env$extra <- paste0(" (Estimation with " ,crayon::bold$yellow(.estimation.method),")")
  # ####

  ## fit ----
  .env$fit <- .fit
  # ####

  ## Additioanl output ----
  #.data <- .dynmodelObject$data # might need to change name do to input

  .temp <- nlmixrDynmodelConvert(.nlmixrObject)
  .temp.inits <- .nlmixrObject$dynmodel.fun(.temp$inits)
  .parameters <- c(.temp.inits, .temp$fixPars) #change parameters to final estimates. Should not be inits.
  #.parameters <- c(.temp$inits, .temp$fixPars) #change parameters to final estimates. Should not be inits.
  .system <- .temp$system

  .rxControl$returnType <- "data.frame"
  .rxControl$addDosing=TRUE

  .nlmixr.sim = do.call(RxODE :: rxSolve, c(list(object=.system, params=.parameters, events=.data), .rxControl))
  .ID <- if (is.null(.nlmixr.sim$id)) {rep(1, nrow(.data))} else {.nlmixr.sim$id}
  .TIME <- .nlmixr.sim$time
  .DV = RxODE::etTrans(.data,.system,addCmt=TRUE,dropUnits=TRUE,allTimeVar=TRUE)
  .DV <- .DV$DV
  .EVID <- .nlmixr.sim$evid
  .PRED <- .nlmixr.sim$nlmixr_pred
  .RES <- .DV - .PRED
  .sgy <- rep(.sgy[1],length(.DV))
  .WRES <- (1/sqrt(.sgy))*.RES

  .nlmixr.pred <- data.frame(ID = .ID, TIME = .TIME, DV = .DV, EVID = .EVID, PRED = .PRED, RES = .RES, WRES = .WRES)
  .nlmixr.pred <- cbind(.nlmixr.pred, .nlmixr.sim[ , -which(names(.nlmixr.sim) %in% c("time","evid","nlmixr_pred"))])

  class(.env) <- "nlmixrFitCoreSilent"
  .nlmixr.pred.temp <- c("nlmixrDynmodel", "nlmixrFitData", "nlmixrFitCore", "tbl_df", "tbl", "data.frame")
  attr(.nlmixr.pred.temp,".foceiEnv") <- .env
  class(.nlmixr.pred) <- .nlmixr.pred.temp

  ## $time ----
  .time$tableTime <- (proc.time() - .pt)["elapsed"]
  names(.time) <- c("setup", "scaling", "optimization", "Hessian", "run total", "table")
  .env$time <- .time
  # ####



  return(.nlmixr.pred)

  # ####

  ## optional
  ## $seed for mcmc
  ## $parHist
  ## $parHistStacked
}
# #########################################################################

# nlmixrDynmodelConvert() ---------------------------------------------------
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
#'

# run devtools::document() to update documentation
# go through and make sure
# devtools::check_man() used to identify missing and problems

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
  .temp.model <- RxODE::rxSymPySetupPred(.nmf.original$rxode.pred,
                                         function(){return(nlmixr_pred)},
                                         .nmf.original$theta.pars,
                                         .nmf.original$error,
                                         grad=FALSE,
                                         pred.minus.dv=TRUE, sum.prod=FALSE, #control$sumProd,
                                         theta.derivs=FALSE, optExpression=TRUE, #control$optExpression,
                                         run.internal=TRUE, only.numeric=TRUE)


  # assign fixed terms
  .fix.index <- if (length(which(.nmf$fix==TRUE))==0) {NULL} else {which(.nmf$fix==TRUE)} # obtain row location for fixed terms
  .ref.fix <- substring(.nmf$name[.fix.index],2)

  .temp.log.fixPars.index <- intersect(.fix.index,.temp.model$log.etas)
  .temp.log.fixPars <- .nmf$est[.temp.log.fixPars.index]
  names(.temp.log.fixPars) <- substring(.nmf$name[.temp.log.fixPars.index],2)

  .temp.nonlog.fixPars.index <- setdiff(.fix.index,.temp.model$log.etas)
  .temp.nonlog.fixPars <- .nmf$est[.temp.nonlog.fixPars.index]
  names( .temp.nonlog.fixPars ) <- substring(.nmf$name[.temp.nonlog.fixPars.index],2)

  .fixPars <- c(.temp.log.fixPars,.temp.nonlog.fixPars)
  .fixPars <- .fixPars[order(factor(names(.fixPars), levels=.ref.fix))]

  .return <- c(.return,fixPars=list(.fixPars))

  # assign theta terms(estimated terms excluding error terms)
  .theta.index <-
    if (is.null(.fix.index)){
      which(!is.na(.nmf$ntheta) & is.na(.nmf$err),TRUE)
    } else {
      which(!is.na(.nmf$ntheta) & is.na(.nmf$err),TRUE)[-which(.nmf$fix == TRUE)]
    } # row location for theta values

  .ref.theta <- .nmf$name[.theta.index]

  .temp.log.theta.index <- intersect(.theta.index,.temp.model$log.etas)
  .temp.log.theta <- .nmf$est[.temp.log.theta.index]
  names(.temp.log.theta) <- .nmf$name[.temp.log.theta.index]

  .temp.nonlog.theta.index <- setdiff(.theta.index,.temp.model$log.etas)
  .temp.nonlog.theta <- .nmf$est[.temp.nonlog.theta.index]
  names(.temp.nonlog.theta) <- .nmf$name[.temp.nonlog.theta.index]

  .theta <- c(.temp.nonlog.theta,.temp.log.theta)
  .theta <- .theta[order(factor(names(.theta), levels=.ref.theta))]

  # .theta.name.index <-  .theta.index + 2 + length(.theta.index) + length(.fix.index)

  # names(.theta) <- .temp.model$pred.only$lhs[.theta.name.index]

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
  .system <- RxODE(.nmf.original$rxode.pred) # (use nlmixr_pred)
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
#' @param CI Confidence interval range from 0-1, default is 0.95
#' @export

# run devtools::document() to update documentation
# go through and make sure
# devtools::check_man() used to identify missing and problems

# ADD RxODE::rxControl()

dynmodelControl <- function(...,
                            fixPars=NULL,
                            ci=0.95,
                            nlmixrOutput=FALSE,
                            digs=3,
                            lower = -Inf,
                            upper = Inf,
                            ## mma doesn't work
                            ## lbfgsbLG
                            ## slsqp
                            method=c("bobyqa", "Nelder-Mead", "lbfgsb3c", "L-BFGS-B", "PORT",
                                     "mma", "lbfgsbLG", "slsqp", "Rvmmin"),
                            ftol_rel=1e-6,
                            maxeval=999,
                            scaleTo=1.0,
                            scaleObjective=0,
                            normType=c("rescale2", "constant", "mean", "rescale", "std", "len"),
                            scaleType=c("nlmixr", "norm", "mult", "multAdd"),
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
                            print=1,
                            maxfun = NULL,
                            # lbfgsb3c
                            trace=0,
                            factr=NULL,
                            pgtol=NULL,
                            abstol=NULL,
                            reltol=NULL,
                            lmm=NULL,
                            maxit=100000L,
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
                            diff.g=NULL,
                            ## Sigdig
                            boundTol=NULL,
                            epsilon=NULL,
                            derivSwitchTol=NULL,
                            sigdig=4,
                            covMethod=c("nlmixrHess", "optimHess"),
                            # rxControl
                            gillK=10L,
                            gillStep=4,
                            gillFtol=0,
                            gillRtol=sqrt(.Machine$double.eps),
                            gillKcov=10L,
                            gillStepCov=2,
                            gillFtolCov=0,
                            rxControl = RxODE::rxControl()
                            ) {
  if (is.null(boundTol)){
    boundTol <- 5 * 10 ^ (-sigdig + 1)
  }
  if (is.null(epsilon)){
    epsilon <- 10 ^ (-sigdig - 1)
  }
  if (is.null(abstol)){
    abstol <- 10 ^ (-sigdig - 1)
  }
  if (is.null(reltol)){
    reltol <- 10 ^ (-sigdig - 1)
  }
  if (is.null(rhoend)){
    rhoend <- 10 ^ (-sigdig - 1);
  }
  if (is.null(factr)){
    factr <- 10 ^ (-sigdig - 1) / .Machine$double.eps;
  }
  if (is.null(atol)){
    atol <- 0.5 * 10 ^ (-sigdig - 2);
  }
  if (is.null(rtol)){
    rtol <- 0.5 * 10 ^ (-sigdig - 2);
  }
  if (is.null(rel.tol)){
    rel.tol <- 10 ^ (-sigdig - 1);
  }
  if (is.null(x.tol)){
    x.tol <- 10 ^ (-sigdig - 1);
  }
  if (is.null(derivSwitchTol)){
    derivSwitchTol <- 2 * 10 ^ (-sigdig-1);
  }

  if (missing(method)){method = "bobyqa"}
  if (missing(normType)){normType = "rescale2"}
  if (missing(scaleType)){scaleType = "nlmixr"}

  .ret <- list(
    fixPars=fixPars,
    ci=ci,
    nlmixrOutput=nlmixrOutput,
    digs=digs,
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
    print=print,
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
    diff.g=diff.g,
    covMethod=match.arg(covMethod),
    rxControl = rxControl,
    gillK=as.integer(gillK),
    gillKcov=as.integer(gillKcov),
    gillRtol=as.double(gillRtol),
    gillStep=as.double(gillStep),
    gillStepCov=as.double(gillStepCov)
  )
  .w <- which(sapply(.ret, is.null))
  .ret <- .ret[-.w];
  class(.ret) <- "dynmodelControl"
  return(.ret)

}

# #########################################################################

# dynmodel()  -------------------------------------------------------------
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
#' @author Wenping Wang, Mason McComb and Matt Fidler
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

dynmodel = function(system, model, inits, data, nlmixrObject=NULL, control=list(), ...){

  # Timing and environment --------------------------------------------------
  .pt <- proc.time()
  .time <- c()

  .dynmodel.env <- new.env(parent=emptyenv())

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
      .norm <<- T
    } else{NULL
      .norm <<- NULL}

    .sigma.dnorm = if("dnorm" %in% .model) {
      as.numeric(.model[which(.model=="dnorm")+1])
      .dnorm <<- T
    } else{NULL
      .dnorm <<- NULL}

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
                   norm=.sigma.norm, dnorm=.sigma.dnorm, tbs=.sigma.tbs, tbsYj=.sigma.tbsYj)

    inits.err <- inits.err[which(names(inits.err) %in% intersect(names(inits.err),.model))]
    inits.err <<- inits.err
    .model <- c("dv" = .model[2], "pred" = .model[3], inits.err)
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
  if (length(nodef) & is.null(nlmixrObject)) {
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
  if (is.null(nlmixrObject)) {
    pars = modelVars$params
    } else {
    .temp<-nlmixrDynmodelConvert(nlmixrObject)
    pars = names(.temp$inits)}

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

  # Objective Function ------------------------------------------------------

  yo = RxODE::etTrans(data,system,addCmt=TRUE,dropUnits=TRUE,allTimeVar=TRUE)
  yo = yo$DV[yo$EVID==0]
  nobs <- length(yo)

  .time$setupTime <- (proc.time() - .pt)["elapsed"]

  sgy<-c()
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

# function that translates the nmf$dynmodel.fun() parameters if function is null, dont apply translation.
    # RxODE(nmf$rxode.pred)
    # write function outside so it is faster

    if (!is.null(nlmixrObject)) {
      theta <- nlmixrObject$dynmodel.fun(theta)
    }
    # call rxODE for simulation

    s = do.call(RxODE :: rxSolve, c(list(object=system, params=theta, events=data), rxControl))

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
      tbsYj <- NULL

      # assign names sigma names
      if(any(names(th) %in% names(model[[1]]))){
        for(i in 1:sum((names(th) %in% names(model[[1]])))) {
          assign(names(th[names(th) %in% names(model[[1]])])[i], as.numeric(th[names(th) %in% names(model[[1]])])[i])
        }
      }

      if (!is.null(.norm)) add <- norm
      if (!is.null(.dnorm)) add <- dnorm
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
      if (!is.null(tbsYj)) lambda <- tbsYj

      # predictted and observed values from RxODE
      #yp = s[rows,x["pred"]]
      yp = s[,x["pred"]]
      #yo = data[, x["dv"]]

      # log normal transformation ----
      if (!is.null(.logn) | !is.null(.dlnorm)){
        .h.x <- boxCox(yo, lambda) #log(yo) # obs
        .h.y <- boxCox(yp, lambda) #log(yp)  # pred
        if("pow" %in% names(model[[1]])) {
          .h.y.var <- yp^(2*pow2)*thresh(pow)^2 + thresh(add)^2  # variance of pred
        } else {
          .h.y.var <- yp^(2*pow2)*thresh(prop)^2 + thresh(add)^2  # variance of pred
        }
        # boxCox transformed -2 log-likelihood
        .logn.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll = .logn.n2ll - 2*(0-1)*log(yo) -2*log(2*pi) # lambda is zero here
        sgy <<- h.y.var
        # negative log-likelihood function for output
        ll = .5*(.n2ll)
      }
      # boxCox Transform ----
      else if ("boxCox" %in% names(model[[1]]) | "tbs" %in% names(model[[1]])) {
        .h.x <- boxCox(yo, lambda) # obs
        .h.y <- boxCox(yp, lambda) # pred
        if("pow" %in% names(model[[1]])) {
          .h.y.var <- (yp^(2*pow2))*thresh(pow)^2 + thresh(add)^2  # variance of pred
        } else {
          .h.y.var <- yp^(2*pow2)*thresh(prop)^2 + thresh(add)^2  # variance of pred
        }
        # boxCox transformed -2 log-likelihood
        .boxCox.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll = .boxCox.n2ll - 2*(lambda-1)*log(yo) -2*log(2*pi)
        sgy <<- h.y.var
        # negative log-likelihood function for output
        ll = .5*(.n2ll)
      }
      # yeoJohnson Transform ----
      else if("yeoJohnson" %in% names(model[[1]]) | "tbsYj" %in% names(model[[1]])) {
        .h.x <- yeoJohnson(yo, lambda) #obs
        .h.y <- yeoJohnson(yp, lambda) #pred
        if("pow" %in% names(model[[1]])) {
          .h.y.var <- (yp^(2*pow2))*thresh(pow)^2 + thresh(add)^2  # variance of pred
        } else {
          .h.y.var <- yp^(2*pow2)*thresh(prop)^2 + thresh(add)^2  # variance of pred
        }
        # yeoJohnson transformed -2 log-likelihood
        .yeoJohnson.n2ll = log(.h.y.var) + ((.h.x - .h.y)^2)/.h.y.var
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll <- ifelse(yo >= 0,
                        .yeoJohnson.n2ll -2*(lambda-1)*log(yo+1) -2*log(2*pi),
                        .yeoJohnson.n2ll -2*(1-lambda)*log(-yo+1) -2*log(2*pi)
        )
        sgy <<- h.y.var
        # negative log-likelihood function for output
        ll = .5*(.n2ll)
      }
      # power model ----
      else if ("pow2" %in% names(model[[1]])) {
        sgy = thresh(add) + thresh(pow)*yp^(pow2)
        assign("sgy",sgy,envir = .dynmodel.env)
        sgy <<- sgy
        ll = .5*((yo - yp)^2/(sgy^2) + log(sgy^2) + log(2*pi))
      }
      # all other error models ----
      else {
        #  if (identical(c("dv","pred"),names(model[[1]]))){
        if (length(names(model[[1]]))==2){
          sgy = 1
          sgy <<- sgy
        }else{
          sgy = thresh(add) + thresh(prop)*yp
          assign("sgy",sgy,envir = .dynmodel.env)
        }
        ll = .5*((yo - yp)^2/(sgy^2) + log(sgy^2) + log(2*pi))
      }
      sgy <<- sgy
      sum(ll)
    })
    sgy <<- sgy
    do.call("sum", l)  # same as return(as.numeric(l)), l is a list for each value in the model?
  }

  # FIXME: Put options from control here gillK etc
  .funs <- nlmixrGradFun(obj, print=control$print,
                         gillRtol=control$gillRtol,
                         gillK=control$gillK,
                         gillStep=control$gillStep,
                         gillFtol=control$gillFtol,
                         thetaNames=names(inits))


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
  scaleC <- control$scaleC;
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
    .st <- proc.time()
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
    .time$scalingTime<- (proc.time() - .st)["elapsed"]
  } else {
    .st <- proc.time()
    .inits <- inits
    .lower <- lower
    .upper <- upper
    .time$scalingTime<-(proc.time() - .st)["elapsed"]
  }

  # Optimization -----------------------------------------------------------------------
  if (method == "bobyqa"){
    .optFun <- .bobyqa;
  } else if (any(method == c("nlminb", "PORT"))){
    .optFun <- .nlminb;
  } else if (method == "mma"){
    .optFun <- .nloptr;
  } else if (method == "slsqp"){
    .optFun <- .slsqp;
  } else if (method == "lbfgsbLG"){
    .optFun <- .lbfgsbLG;
  } else if (method == "Rvmmin"){
    .optFun <- .Rvmmin;
  } else if (method=="Nelder-Mead"){
    .optFun <- .mymin;
  } else if (method == "lbfgsb3c"){
    .optFun <- .lbfgsb3c
  } else if (method == "L-BFGS-B"){
    .optFun <- .lbfgsxbO
  } else {
    stop("Optimization method unknown.");
  }
  .ot <- proc.time()
  fit <- .optFun(as.vector(.inits), fn=.funs$eval, gr=.funs$gr, lower=.lower, upper=.upper, control=control);
  .message <- fit$message
  ## fit$value
  assign("fit",fit,envir = .dynmodel.env)
  .time$optimizationTime <- (proc.time() - .ot)["elapsed"]

  # Hessian -----------------------------------------------------------------------
  .ht <- proc.time()

  if (control$covMethod == "optimHess"){
    fit$hessian = try(optimHess(fit$par, obj, control=control) , silent=TRUE)
  } else {
    # FIXME: Put options from control here gillK etc
    fit$hessian = try(nlmixrHess(fit$par, obj,
                                 gillRtol=control$gillRtol,
                                 gillK=control$gillKcov,
                                 gillStep=control$gillStepCov,
                                 gillFtol=control$gillFtolCov) , silent=TRUE)
  }


  if(inherits(fit$hessian,"try-error")){
    se = rep(NA, length(fit$par))
    warning("standard error of the Hessian has failed")
  } else {
    cov.matrix = solve(fit$hessian)
    se = sqrt(diag(solve(fit$hessian)))
  }

  # reassign the negative values to positive for add, prop/pow since they are standard deviations
  if (!is.na(match("add",names(inits)))) fit$par[match("add",names(inits))] = abs(fit$par[match("add",names(inits))])
  if (!is.na(match("prop",names(inits)))) fit$par[match("prop",names(inits))] = abs(fit$par[match("prop",names(inits))])
  if (!is.na(match("pow",names(inits)))) fit$par[match("pow",names(inits))] = abs(fit$par[match("pow",names(inits))])
  if (!is.na(match("norm",names(inits)))) fit$par[match("norm",names(inits))] = abs(fit$par[match("norm",names(inits))])
  if (!is.na(match("dnorm",names(inits)))) fit$par[match("dnorm",names(inits))] = abs(fit$par[match("dnorm",names(inits))])

  .time$hessianTime <- (proc.time() - .ht)["elapsed"]


  # dynmodel Output -------------------------------------------------------
  # unscale optmized parameters here if scaling was used:
  if(normType != "constant" & scaleType != "norm"){
    par.temp <- numeric(length(fit$par))
    for (i in 1:length(par.temp)) {par.temp[i] <- unscalePar(fit$par,i)}
    fit$par <- par.temp
  }

  # create table for output
  res = cbind(fit$par, abs(se), abs(se/fit$par*100))
  dimnames(res) = list(names(inits), c("est", "se", "%cv"))


  # nobs = 0
  # l = lapply(model, function(x) {
  #   yo =  data[,model[[1]]["dv"][[1]]]#data[, x["dv"]]
  #   nobs <<- nobs + length(yo)
  # })
  if (!is.null(fit$objective)) fit$value = fit$objective

  # Output
  res = c(list(res=res, obj=obj, npar=length(fit$par), nobs=nobs, data=data), fit)
  class(res) = "dyn.ID"
  # Final Output ----------------------------------------------------------
  .time$totalTime <- (proc.time() - .pt)["elapsed"]

  .time <- as.data.frame(.time)
  names(.time) <- c("setup", "scaling", "optimization", "Hessian", "total")

  if (!is.null(nlmixrObject) & control$nlmixrOutput){
    nlmixr.ouptut <- as.focei.dynmodel(.dynmodelObject = res, .nlmixrObject = nlmixrObject, .data = data, .time = .time, .fit = fit, .message = .message, .inits.err = inits.err, .cov = cov.matrix, .sgy = sgy,
                                       .dynmodelControl = control, .nobs2=0, .pt=proc.time(), .rxControl = RxODE::rxControl())
    .hist <- .funs$hist();
    assign("parHistData", .hist, nlmixr.ouptut$env);
    return(nlmixr.ouptut)
  }
  else {
    return(res)
  }
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

genobj = function(system, model, evTable, inits, data, fixPars=NULL,
squared=T
){

  # Error model  -------------------------------------------------------------
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

  # Check dynmodel() inputs, Define vars, modelVars, pars,  ------------
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


  # Additional assignment ---------------------------------------------------
  ## is this necessary ##
  have_zero = min(data$time) <= 0
  rows = if(have_zero) T else -1 # used in line 304 in obj()
  ## ---------------- ##

  # Objective Function ------------------------------------------------------
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


# slice sampling  ---------------------------------------------------------
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

# pars.cp as a data frame, run in do all the ode solving at the end in parallel
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

  # Objective Function ------------------------------------------------------
  l = genobj(system, model, evTable, inits, data, fixPars, squared)
  rho = environment()
  pars = l$inits
  fr0 = l$obj

  if (is.null(seed)) seed=99
  set.seed(seed)

  # progress
  on.exit(RxODE::rxProgressAbort("Aborted MCMC Calculation"))
  RxODE::rxProgress(nsim)

  # slice sampling
  s = t(sapply(1:nsim, function(k,rho) {
    pars = do.slice(get("pars", rho), fr0)
    RxODE::rxTick()
    assign("pars", pars, rho)
  }, rho=rho))

  if (squared) s = s*s
  attr(s, "calls") <- calls
  attr(s, "obj") <- fr0
  attr(s, "class") <- "dyn.mcmc"
  RxODE::rxProgressStop()
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


# Utilities for nlmixr ----------------------------------------------------


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
# ########################################################################

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

##'@export
plot.nlmixrDynmodel <- function(x, y, ...){
    .lst  <- list();
    .tp  <- traceplot(x)
    if (!is.null(.tp)) .lst[[length(.lst)+1]] <- .tp;
    .dat <- as.data.frame(x);
    .p1 <- ggplot2::ggplot(.dat, ggplot2::aes_string("PRED", "DV")) +
           ggplot2::geom_abline(slope=1, intercept=0, col="red", size=1.2) +
            ggplot2::geom_point() + xlab("Predictions") +
            ggplot2::ggtitle("DV vs PRED")
    .lst[[length(.lst)+1]] <- .p1

    .p0 <- ggplot2::ggplot(.dat, ggplot2::aes_string(x="PRED", y="RES")) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(slope=0, intercept=0, col="red") +
        ggplot2::ggtitle("PRED vs RES")
    .lst[[length(.lst)+1]] <- .p0
    .p0 <- ggplot2::ggplot(.dat, ggplot2::aes_string(x="PRED", y="WRES")) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(slope=0, intercept=0, col="red") +
            ggplot2::ggtitle("PRED vs WRES")
    .lst[[length(.lst)+1]] <- .p0
    .p0 <- ggplot2::ggplot(.dat, ggplot2::aes_string(x="TIME", y="RES")) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(slope=0, intercept=0, col="red") +
        ggplot2::ggtitle("TIME vs RES")
    .lst[[length(.lst)+1]] <- .p0
    .p0 <- ggplot2::ggplot(.dat, ggplot2::aes_string(x="TIME", y="WRES")) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(slope=0, intercept=0, col="red") +
            ggplot2::ggtitle("TIME vs WRES")
    .lst[[length(.lst)+1]] <- .p0
    .ids <- unique(.dat$ID)
    .s <- seq(1, length(.ids), by=16)
    .j <- 0;
    for (i  in .s){
        .j <- .j + 1
        .tmp <- .ids[seq(i, i + 15)]
        .tmp <- .tmp[!is.na(.tmp)];
        .d1 <- .dat[.dat$ID %in% .tmp, ];

        .p3 <- ggplot2::ggplot(.d1, aes(x=TIME, y=DV)) +
            ggplot2::geom_point() +
            ggplot2::geom_line(aes(x=TIME, y=PRED), col="red", size=1.2) +
            ggplot2::facet_wrap(~ID) +
            ggplot2::ggtitle(sprintf("Individual Plots (%s of %s)", .j, length(.s)))
        .lst[[length(.lst)+1]] <- .p3
    }
    class(.lst)  <- "nlmixrPlotList"
    return(.lst)
}

##' Convert fit to classic dynmodel object
##'
##' @param x nlmixr object to convert to dynmodel object
##' @return dynmodel
##' @author Matthew Fidler
##' @export
as.dynmodel <- function(x){
    .ret <- try(x$dynmodelObject, silent=TRUE)
    if (inherits(.ret, "try-error") || is.null(.ret))
        stop("Cannot convert to dynmodel object");
    return(.ret)
}
