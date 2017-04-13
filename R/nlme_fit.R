## nlme_fit.R: population PK/PD modeling library
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

fmt_infusion_data <- function(dat) {
	x1 <- subset(dat, dat$EVID>10000 & dat$AMT>0)
	x2 <- subset(dat, dat$EVID>10000 & dat$AMT<0)
	x1$DUR <- x2$TIME - x1$TIME
	x1$AMT <- x1$AMT * x1$DUR
	x <- subset(dat, dat$EVID < 10000)
	x$DUR <- NA
	x <- rbind(x, x1)
	ord <- order(x$ID, x$TIME, -x$EVID)
	x[ord,]
}

cmt_fn_templ <- "
require(parallel)

user_fn <- function(<%=arg1%>, TIME, ID)
{
    unlist(parallel::mclapply(as.character(unique(ID)), function(subj)
    {
        sel.d <- ..ModList$ds$ID==as.integer(subj)
        dose<-..ModList$ds[sel.d, \"AMT\"]
        dstm<-..ModList$ds[sel.d, \"TIME\"]
        Tinf<-..ModList$ds[sel.d, \"DUR\"]

        sel <- ID==subj
        time.subj <- TIME[sel]
        s <- ..ModList$PKpars(<%=arg2%>)
        pkpars <- toupper(names(s))
        names(s) <- pkpars

        if (..ModList$oral) {
            if (is.element(\"TLAG\", pkpars)) {
                theta <- s[..ModList$refpars]
            } else {
                theta <- c(s[..ModList$refpars[1:(2*..ModList$ncmt+1)]], 0)    #no TLAG, set to 0
            }
        } else {
            theta <- c(s[..ModList$refpars[1:(2*..ModList$ncmt)]], 0, 0)
        }

        cp <- lin_cmt(time.subj, dstm, dose, Tinf, theta, ..ModList$oral, ..ModList$infusion, ..ModList$ncmt, ..ModList$parameterization)
        cp
    }, mc.cores=<%=mc.cores%>))
}
"


#' Fit nlme-based linear compartment mixed-effect model using closed form solution
#'
#' 'nlme_lin_cmpt' fits a linear one to three compartment model with
#' either first order absorption, or i.v. bolus, or i.v. infusion.  A
#' user specifies the number of compartments, route of drug
#' administrations, and the model parameterization. `nlmixr` supports
#' the clearance/volume parameterization and the micro constant
#' parameterization, with the former as the default.  Specification of
#' fixed effects, random effects and intial values follows the standard
#' nlme notations.
#'
#' @param dat data to be fitted
#' @param par_model list: model for fixed effects, randoms effects and initial values using nlme-type syntax.
#' @param ncmt numerical: number of compartments: 1-3
#' @param oral logical
#' @param infusion logical
#' @param tlag logical
#' @param parameterization numerical: type of parameterization, 1=clearance/volume, 2=micro-constants
#' @param par_trans function: calculation of PK parameters
#' @param mc.cores number of cores used in fitting (only for Linux)
#' @param ... additional nlme options
#' @return NULL
#' @author Wenping Wang
#' @examples
#' \dontrun{
#' library(nlmixr)
#'
#' dat <- read.table(system.file("examples/theo_md.txt", package = "nlmixr"), head=TRUE)
#' specs <- list(fixed=lKA+lCL+lV~1, random = pdDiag(lKA+lCL~1), start=c(lKA=0.5, lCL=-3.2, lV=-1))
#' fit <- nlme_lin_cmpt(dat, par_model=specs, ncmt=1, verbose=TRUE)
#' plot(augPred(fit,level=0:1))
#' summary(fit)
#'
#' }
#' @export
nlme_lin_cmpt <- function(dat, par_model,
	ncmt, oral=TRUE, infusion=FALSE, tlag=FALSE, parameterization=1,
	par_trans=get.parfn(oral, ncmt, parameterization, tlag),
	mc.cores=1, ...)
{
    if(oral*infusion) {
        msg <- "oral and infusion cannot be all TRUE"
        stop(msg)
    }

    #prep PKpars
    PKpars <- par_trans
	x <- deparse(body(PKpars))
	len <- length(x)
	x[len] <- "unlist(as.list(environment()))"
	x <- paste(c(x, "}"), collapse = "\n")
	body(PKpars) <- parse(text=x)

    #a new env with a ref in .GlobalEnv, holding model components
    #a hack due to non-std call by nlme
	assign("..ModList", new.env(), envir=.GlobalEnv)

    #master par list
    pm <- list(
        c("CL", "V", "KA", "TLAG"),
        c("CL", "V", "CLD", "VT", "KA", "TLAG"),
        c("CL", "V", "CLD", "VT", "CLD2", "VT2", "KA", "TLAG"),
        c("KE", "V", "KA", "TLAG"),
        c("KE", "V", "K12", "K21", "KA", "TLAG"),
        c("KE", "V", "K12", "K21", "K13", "K31", "KA", "TLAG")
    )
    dim(pm)<-c(3,2)


    #gen user_fn
    s <- formals(PKpars)
    arg1 <- paste(names(s), collapse=", ")
    arg2 <- paste(unlist(lapply(names(s), function(x) paste(x,"=",x,"[sel][1]", sep=""))), collapse=", ")
    arg3 <- sprintf("list(%s)", paste(names(s), "=.1", collapse=", "))
    brew::brew(text=cmt_fn_templ, output="fn.txt")
    source("fn.txt", local=..ModList)


    refpars <- pm[[ncmt, parameterization]]
    npars <- length(refpars)
    s <- do.call(PKpars, eval(parse(text=arg3)))
    pkpars <- names(s)

    if (!oral) {
        refpars <- refpars[1:(npars-2)]
        pkpars.aug <- pkpars
    } else {
        pkpars.aug <- pkpars
        if(!is.element("TLAG", pkpars)) {
            pkpars.aug <- c(pkpars, "TLAG")
        }
    }
    notdefed <- setdiff(refpars, pkpars.aug)
    if(length(notdefed)) {
        msg <- paste(c("undefined PK pars:", notdefed, collapse=" "))
        stop(msg)
    }

    #data prep
    ..ModList$oral<-oral
    ..ModList$infusion <- infusion
    ..ModList$ncmt<-ncmt
    ..ModList$parameterization<-parameterization
    if (infusion) {
    	dat <- fmt_infusion_data(dat)
    } else {
    	dat$DUR <- -1
    }
    ..ModList$ds <- dat[dat$EVID>0, c("ID", "TIME", "AMT", "DUR")]
    dat$DUR <- NULL
    dat <- dat[dat$EVID==0,]
    ..ModList$dat.g <- nlme::groupedData(DV~TIME|ID, dat)
    ..ModList$PKpars <- PKpars
	..ModList$refpars=refpars

    mod.specs <- list(model=as.formula(sprintf("DV ~ ..ModList$user_fn(%s, TIME, ID)", arg1)),
	    data = ..ModList$dat.g, fixed=par_model$fixed, random = par_model$random,
	    start=par_model$start, ...)
    do.call(getFromNamespace("nlme","nlme"), mod.specs)
}


ode_fn_templ = "
require(parallel)

user_fn <- function(<%=arg1%>, TIME, ID)
{
	z <- parallel::mclapply(as.character(unique(ID)), function(subj)
	{
		ev <- eventTable()
		ev$import.EventTable(subset(..ModList$dat.o, id==as.integer(subj)))
		#obs.rec <- ev$get.obs.rec()

		sel <- ID==subj
		plist <- ..ModList$PKpars(<%=arg2%>)
		inits <- plist$initCondition
		plist$initCondition <- NULL
		theta <- unlist(plist)

		if (any(theta>1e38)) {
		    warning('large parameter values. may rewrite par_trans.')
		    print(theta)
		}
		if (..ModList$debugODE) {
		    print(subj)
		    print(theta)
		}

		m <- ..ModList$m1$run(theta, ev, inits, transit_abs=<%=transit_abs%>, atol=<%=atol%>, rtol=<%=rtol%>)
		if (is.null(dim(m))) m = t(as.matrix(m))
		den <- if(is.null(..ModList$response.scaler)) 1 else theta[..ModList$response.scaler]
		m[, ..ModList$response]/den
	}, mc.cores=<%=mc.cores%>)
	unlist(z)
}
"

#' Fit nlme-based mixed-effect model using ODE implementation
#'
#' 'nlme_ode' fits a mixed-effect model described using ordinary differential
#' equation (ODEs). The ODE-definition follows RxODE syntax.
#' Specification of fixed effects, random effects and intial values follows
#' the standard nlme notations.
#'
#' @param dat.o data to be fitted
#' @param model a string containing the set of ordinary differential equations (ODE) and other expressions defining the changes in the dynamic  system. For details, see the sections \dQuote{Details} and  \dQuote{\code{RxODE Syntax}} below.
#' @param par_model list: model for fixed effects, randoms effects and initial values.
#' @param par_trans function: calculation of PK parameters
#' @param response names of the response variable
#' @param response.scaler optional response variable scaler. default is NULL
#' @param transit_abs a logical if transit absorption model is enabled
#' @param atol atol (absolute tolerance for ODE-solver)
#' @param rtol rtol (relative tolerance for ODE-solver)
#' @param debugODE a logical if debugging is enabled
#' @param mc.cores number of cores used in fitting (only for Linux)
#' @param ... additional nlme options
#' @return NULL
#' @details
#'    The ODE-based model specification may be coded inside a character
#'    string or in a text file, see Section \emph{RxODE Syntax} below for
#'    coding details.  An internal \code{RxODE} compilation manager object
#'    translates the ODE system into C, compiles it, and dynamically loads the
#'    object code into the current R session.  The call to \code{RxODE}
#'    produces an object of class \code{RxODE} which consists of a list-like
#'    structure (closure) with various member functions (see Section
#'    \emph{Value} below).
#'
#' @section RxODE Syntax:
#'
#'    An \code{RxODE} model specification consists of one or more
#'    statements terminated by semi-colons, \sQuote{\code{;}}, and
#'    optional comments (comments are delimited by \code{#} and an
#'    end-of-line marker).  \strong{NB:} Comments are not allowed
#'    inside statements.
#'
#'    A block of statements is a set of statements delimited by
#'    curly braces, \sQuote{\code{\{ ... \}}}.
#'    Statements can be either assignments or conditional \code{if}
#'    statements. Assignment statements can be either \dQuote{simple}
#'    assignments, where the left hand is an identifier (i.e., variable), or
#'    special \dQuote{time-derivative} assignments, where the left hand
#'    specifies the change of that variable with respect to time
#'    e.g., \code{d/dt(depot)}.
#'
#'    Expressions in assignment and \sQuote{\code{if}} statements can be
#'    numeric or logical (no character expressions are currently supported).
#'    Numeric expressions can include the following numeric operators
#'    (\sQuote{\code{+}}, \sQuote{\code{-}}, \sQuote{\code{*}},
#'    \sQuote{\code{/}}, \sQuote{\code{^}}),   and
#'    those mathematical functions defined in the C or the
#'    R math libraries (e.g., \code{fabs}, \code{exp}, \code{log}, \code{sin}).
#'    (Note that the modulo operator \sQuote{\code{\%}} is currently
#'    not supported.)
#'
#'    Identifiers in an \code{RxODE} model specification can refer to:
#'    \itemize{
#'       \item state variables in the dynamic system (e.g., compartments in a
#'       pharmacokinetic/pharmacodynamic model);
#'       \item implied input variable, \code{t} (time),
#'       \code{podo} (oral dose, for absorption models), and
#'       \code{tlast} (last time point);
#'       \item model parameters, (\code{ka} rate of absorption, \code{CL}
#'       clearance, etc.);
#'       \item others, as created by assignments as part of the model
#'       specification.
#'    }
#'
#'    Identifiers consist of case-sensitive alphanumeric characters,
#'    plus the underscore \sQuote{_} character.  \strong{NB:} the
#'    dot \sQuote{.} character is \strong{not} a valid character
#'    identifier.
#'
#'    The values of these variables at pre-specified time points are
#'    saved as part of the fitted/integrated/solved model (see
#'    \code{\link{eventTable}}, in particular its member function
#'    \code{add.sampling} that defines a set of time points at which
#'    to capture a snapshot of the system via the values of these variables).
#'
#'    The ODE specification mini-language is parsed with the help of
#'    the open source tool \emph{DParser}, Plevyak (2015).
#' @author Wenping Wang
#' @examples
#' \dontrun{
#' library(nlmixr)
#' ode <- "
#' d/dt(depot) =-KA*depot;
#' d/dt(centr) = KA*depot - KE*centr;
#' "
#' dat <- read.table(system.file("examples/theo_md.txt", package = "nlmixr"), head=TRUE)
#' mypar <- function(lKA, lKE, lCL)
#' {
#'     KA=exp(lKA)
#'     KE=exp(lKE)
#'     CL=exp(lCL)
#'     V = CL/KE
#' }
#'
#' specs <- list(fixed=lKA+lKE+lCL~1, random = pdDiag(lKA+lCL~1),
#' 	start=c(lKA=0.5, lKE=-2.5, lCL=-3.2))
#' fit <- nlme_ode(dat, model=ode, par_model=specs, par_trans=mypar,
#' 	response="centr", response.scaler="V")
#'
#' }
#' @export
nlme_ode <- function(dat.o, model, par_model, par_trans,
	response, response.scaler=NULL,
	transit_abs = FALSE,
	atol=1.0e-8, rtol=1.0e-8,
	debugODE=FALSE, mc.cores=1, ...)
{
    if (any(dat.o$EVID[dat.o$EVID>0]<101))
    	stop("incompatible EVID values")

    #a new env with a ref in .GlobalEnv, holding model components
    #a hack due to non-std call by nlme
	assign("..ModList", new.env(), envir=.GlobalEnv)

    #prep ode
    if (class(model)=="RxODE") ..ModList$m1 = model
    else if (class(model)=="character") {
		obj <- basename(tempfile())
		..ModList$m1 <- RxODE::RxODE(model = model, modName = obj)
	}
	else {
		stop('invalid model input')
	}

    #prep PKpars
    PKpars <- par_trans
	x <- deparse(body(PKpars))
	len <- length(x)
	x[len] <- "as.list(environment())"
	x <- paste(c(x, "}"), collapse = "\n")
	body(PKpars) <- parse(text=x)

    #gen user_fn
    s <- formals(PKpars)
    arg1 <- paste(names(s), collapse=", ")
    arg2 <- paste(unlist(lapply(names(s), function(x) paste(x,"=",x,"[sel][1]", sep=""))), collapse=", ")
    brew::brew(text=ode_fn_templ, output="fn.txt")
    source("fn.txt", local=..ModList)

    #data prep
	dat.g <- nlme::groupedData(DV~TIME|ID, subset(dat.o, dat.o$EVID==0))
	names(dat.o) <- tolower(names(dat.o))
	..ModList$response = response
	..ModList$response.scaler = response.scaler
	..ModList$dat.g = dat.g
	..ModList$dat.o = dat.o
	..ModList$PKpars = PKpars
	..ModList$debugODE = debugODE

    mod.specs <- list(model=as.formula(sprintf("DV ~ ..ModList$user_fn(%s, TIME, ID)", arg1)),
	    data = ..ModList$dat.g, fixed=par_model$fixed, random = par_model$random,
	    start=par_model$start, ...)

    do.call(nlme, mod.specs)
}

