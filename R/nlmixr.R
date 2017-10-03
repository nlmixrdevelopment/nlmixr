
.onAttach <- function(libname, pkgname){ ## nocov start
    ## Setup RxODE.prefer.tbl
    nlmixrSetupMemoize()
}

nlmixrSetupMemoize <- function(){
    reSlow <- rex::rex(".slow",end)
    f <- sys.function(-1)
    ns <- environment(f)
    .slow <- ls(pattern=reSlow,envir=ns);
    for (slow in .slow){
        fast <- sub(reSlow, "", slow);
        if (!memoise::is.memoised(get(fast, envir=ns)) && is.null(get(slow, envir=ns))){
            utils::assignInMyNamespace(slow, get(fast, envir=ns))
            utils::assignInMyNamespace(fast, memoise::memoise(get(slow, envir=ns)))
        }

    }
}

##' Clear memoise cache for nlmixr
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrForget <- function(){
    reSlow <- rex::rex(".slow",end)
    f <- sys.function(-1)
    ns <- environment(f)
    .slow <- ls(pattern=reSlow,envir=ns);
    for (slow in .slow){
        fast <- sub(reSlow, "", slow);
        memoise::forget(get(fast, envir=ns));
    }
}


##' @importFrom stats predict logLik na.fail pchisq
##' @importFrom brew brew
##' @importFrom lattice xyplot
##' @importFrom lattice trellis.par.get
##' @importFrom nlme nlme fixed.effects random.effects
##' @importFrom nlme groupedData
##' @importFrom nlme getData
##' @importFrom nlme pdDiag
##' @importFrom RxODE RxODE
##' @importFrom graphics abline lines matplot plot points title
##' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var AIC BIC asOneSidedFormula coef end fitted resid setNames start simulate
##' @importFrom utils assignInMyNamespace getFromNamespace head stack sessionInfo tail
##' @importFrom parallel mclapply
##' @importFrom lbfgs lbfgs
##' @importFrom methods is
##' @importFrom Rcpp evalCpp
##' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_line geom_abline xlab geom_smooth
##' @useDynLib nlmixr, .registration=TRUE


rex::register_shortcuts("nlmixr");
## GGplot use and other issues...
utils::globalVariables(c("DV", "ID", "IPRED", "IRES", "PRED", "TIME", "grp", "initCondition", "values", "nlmixr_pred", "iter", "val", "EVID"));

nlmixr.logo <- "         _             _             \n        | | %9s (_) %s\n  _ __  | | _ __ ___   _ __  __ _ __\n | '_ \\ | || '_ ` _ \\ | |\\ \\/ /| '__|\n | | | || || | | | | || | >  < | |\n |_| |_||_||_| |_| |_||_|/_/\\_\\|_|\n"

##' Messages the nlmixr logo...
##'
##' @param str String to print
##' @param version Version information (by default use package version)
##' @author Matthew L. Fidler
nlmixrLogo <- function(str="", version=sessionInfo()$otherPkgs$nlmixr$Version){
    message(sprintf(nlmixr.logo, str, version));
}
##' Dispaly nlmixr's version
##'
##' @author Matthew L. Fidler
##' @export
nlmixrVersion <- function(){
    nlmixrLogo()
}

##' nlmixr fits population PK and PKPD non-linear mixed effects models.
##'
##' nlmixr is an R package for fitting population pharmacokinetic (PK)
##' and pharmacokinetic-pharmacodynamic (PKPD) models.
##'
##' The nlmixr generalized function allows common access to the nlmixr
##' estimation routines.
##'
##' @param object Fitted object or function specifying the model.
##' @inheritParams nlmixr.fit
##' @param ... Other parameters
##' @return Either a nlmixr model or a nlmixr fit object
##' @author Matthew L. Fidler
##' @export
nlmixr <- function(object, data, est="nlme", control=list(), ...){
    UseMethod("nlmixr")
}

##' @rdname nlmixr
##' @export
nlmixr.function <- function(object, data, est="nlme", control=list(), ...){
    uif <- nlmixrUI(object);
    if (missing(data) && missing(est)){
        return(uif)
    } else {
        nlmixr.fit(uif, data, est, control=control, ...);
    }
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixrUI <- function(object, data, est="nlme", control=list(), ...){
    uif <- object
    if (missing(data) && missing(est)){
        return(uif)
    } else {
        nlmixr.fit(uif, data, est, control=control, ...);
    }
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixr.ui.nlme <- function(object, data, est="nlme", ...){
    env <- attr(object, ".focei.env")
    uif <- env$uif.new;
    if (missing(data) && missing(est)){
        return(uif)
    } else {
        if (missing(data)){
            data <- getData(object);
        }
        nlmixr.fit(uif, data, est, ...);
    }
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixr.ui.focei.fit <- nlmixr.nlmixr.ui.nlme

##' @rdname nlmixr
##' @export
nlmixr.nlmixr.ui.saem <- nlmixr.nlmixr.ui.nlme

##' Fit a nlmixr model
##'
##' @param uif Parsed nlmixr model (by \code{nlmixr(mod.fn)}).
##' @param data Dataset to estimate.  Needs to be RxODE compatible in
##'     EVIDs.
##' @param est Estimation method
##' @param control Estimation control options.  They could be
##'     \code{\link[nlme]{nlmeControl}}, \code{\link{saemControl}}
##' @param ... Parameters passed to estimation method.
##' @param sum.prod Take the RxODE model and use more precise
##'     products/sums.  Increases solving accuracy and solving time.
##' @param focei.translate Translate the model to FOCEi and then run
##'     the tables and objective function so that different estimation
##'     methodologies are comparable through OBJF.
##' @return nlmixr fit object
##' @author Matthew L. Fidler
##' @export
nlmixr.fit <- function(uif, data, est="nlme", control=list(), ...,
                       sum.prod=FALSE, focei.translate=TRUE){
    dat <- data;
    uif$env$infusion <- .Call(`_nlmixr_chkSolvedInf`, as.double(dat$EVID), as.integer(!is.null(uif$nmodel$lin.solved)));
    bad.focei <- "Problem calculating residuals, returning fit without residuals.";
    if (est == "saem"){
        pt <- proc.time()
        args <- as.list(match.call(expand.dots=TRUE))[-1]
        default <- saemControl();
        if (any(names(args) == "mcmc")){
            mcmc <- args$mcmc;
        } else if (any(names(control) == "mcmc")){
            mcmc <- control$mcmc;
        } else {
            mcmc <- default$mcmc;
        }
        if (any(names(args) == "ODEopt")){
            ODEopt <- args$ODEopt;
        } else if (any(names(control) == "ODEopt")){
            ODEopt <- control$ODEopt;
        } else {
            ODEopt <- default$ODEopt;
        }
        if (any(names(args) == "seed")){
            seed <- args$seed;
        } else if (any(names(control) == "seed")){
            seed <- control$seed;
        } else {
            seed <- default$seed;
        }
        if (any(names(args) == "print")){
            print <- args$print;
        } else if (any(names(control) == "print")){
            print <- control$print;
        } else {
            print <- default$print;
        }
        uif$env$mcmc <- mcmc;
        uif$env$ODEopt <- ODEopt;
        uif$env$sum.prod <- sum.prod
        model <- uif$saem.model
        cfg   = configsaem(model=model, data=dat, inits=uif$saem.init,
                           mcmc=mcmc, ODEopt=ODEopt, seed=seed);
        if (print > 1){
            cfg$print <- as.integer(print)
        }
        fit <- model$saem_mod(cfg);
        if (focei.translate){
            ret <- try(as.focei(fit, uif, pt, data=dat));
            if (inherits(ret, "try-error")){
                warning(bad.focei)
                return(fit)
            } else {
                return(ret)
            }
        } else  {
            return(fit);
        }
    } else if (est == "nlme"){
        pt <- proc.time()
        fun <- uif$nlme.fun;
        specs <- uif$nlme.specs;
        grp.fn <- uif$grp.fn
        dat$nlmixr.grp <- factor(apply(data, 1, function(x){
            cur <- x;
            names(cur) <- names(dat);
            with(as.list(cur), {
                return(grp.fn())
            })
        }));
        dat$nlmixr.num <- seq_along(dat$nlmixr.grp)
        print(any(dat$EVID[dat$EVID>0]<101))
        weight <- uif$nlme.var
        if (!is.null(uif$nmodel$lin.solved)){
            fit <- nlme_lin_cmpt(dat, par_model=specs,
                                 ncmt=uif$nmodel$lin.solved$ncmt,
                                 oral=uif$nmodel$lin.solved$oral,
                                 tlag=uif$nmodel$lin.solved$tlag,
                                 infusion=uif$env$infusion,
                                 parameterization=uif$nmodel$lin.solved$parameterization,
                                 par_trans=fun,
                                 weight=weight,
                                 verbose=TRUE,
                                 control=control,
                                 ...);
        } else {
            if (sum.prod){
                rxode <- RxODE::rxSumProdModel(uif$rxode.pred);
            } else {
                rxode <- uif$rxode.pred;
            }
            fit <- nlme_ode(dat,
                            model=rxode,
                            par_model=specs,
                            par_trans=fun,
                            response="nlmixr_pred",
                            weight=weight,
                            verbose=TRUE,
                            control=control,
                            ...);
        }
        ## Run FOCEi using same ETAs and THETA estimates to get
        ## comparable OBJFs and also extract table entries like
        ## CWRES.
        ## return(fit)
        if (focei.translate){
            ret <- as.focei(fit, uif, pt, data=dat)
            ## ret <- try(as.focei(fit, uif, pt, data=dat))
            ## if (inherits(ret, "try-error")){
            ##     warning(bad.focei)
            ##     return(fit);
            ## } else {
            ##     return(ret)
            ## }
            return(ret)
        } else  {
            return(fit);
        }
    } else if (est == "focei"){
        fit <- focei.fit(dat,
                         inits=uif$focei.inits,
                         PKpars=uif$theta.pars,
                         ## par_trans=fun,
                         model=uif$rxode.pred,
                         pred=function(){return(nlmixr_pred)},
                         err=uif$error,
                         lower=uif$focei.lower,
                         upper=uif$focei.upper,
                         theta.names=uif$focei.names,
                         eta.names=uif$eta.names,
                         ...)
        env <- attr(fit, ".focei.env")
        env$uif <- uif;
        uif.new <- uif;
        ns <- names(fit$theta);
        for (n in ns){
            uif.new$ini$est[uif.new$ini$name == n] <- fit$theta[n];
        }
        ome <- fit$omega;
        w <- which(!is.na(uif.new$ini$neta1))
        for (i in w){
            uif.new$ini$est[i] <- ome[uif.new$ini$neta1[i], uif.new$ini$neta2[i]];
        }
        env$uif.new <- uif.new;
        class(fit) <- c("nlmixr.ui.focei.fit", class(fit));
        return(fit);
    }
}
##' Control Options for SAEM
##'
##' @param seed Random Seed for SAEM step.  (Needs to be set for
##'     reproducibility.)  By default this is 99.
##' @param sa.iter Number of iterations in the Stochastic
##'     Approximation (SA) Step
##' @param em.iter Number of iterations in the
##'     Expectation-Maximization (EM) Step
##' @param nmc Number of Markov Chains?
##' @param nu Nu! A new parameter....?
##' @inheritParams RxODE::rxSolve
##' @param print The number it iterations that are completed before
##'     anything is printed to the console.  By default, this is 1.
##' @param ... Other arguments to control SAEM.
##' @return List of options to be used in \code{\link{nlmixr}} fit for
##'     SAEM.
##' @author Wenping Wang & Matthew L. Fidler
##' @export
saemControl <- function(seed=99,
                        sa.iter=200, em.iter=300,
                        nmc=3,
                        nu=c(2,2,2),
                        atol = 1e-08,
                        rtol = 1e-06,
                        stiff = TRUE,
                        transit_abs = FALSE,
                        print=1,
                        ...){
    list(mcmc=list(niter=c(sa.iter, em.iter=em.iter), nmc=nmc, nu=nu),
         ODEopt=list(atol=atol, rtol=rtol,
                     stiff=as.integer(stiff),
                     transit_abs = as.integer(transit_abs)),
         seed=seed,
         print=print, ...)

}
