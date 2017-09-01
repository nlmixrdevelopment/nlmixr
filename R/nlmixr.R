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
##' @param object Fitted object or
##' @param data Data for fit.
##' @param est Estimation routine.
##' @param ... Other parameters
##' @return Either a nlmixr model or a nlmixr fit object
##' @author Matthew L. Fidler
##' @export
nlmixr <- function(object, data, est="nlme", ...){
    UseMethod("nlmixr")
}

##' @rdname nlmixr
##' @export
nlmixr.function <- function(object, data, est="nlme", ...){
    uif <- nlmixrUI(object);
    if (missing(data) && missing(est)){
        return(uif)
    } else {
        nlmixr.fit(uif, data, est, ...);
    }
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixrUI <- function(object, data, est="nlme", ...){
    uif <- object
    if (missing(data) && missing(est)){
        return(uif)
    } else {
        nlmixr.fit(uif, data, est, ...);
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

##' Fit a nlmixr model
##'
##' @param uif Parsed nlmixr model (by \code{nlmixr(mod.fn)}).
##' @param data Dataset to estimate
##' @param est Estimation method
##' @param ... Parameters passed to estimation method.
##' @param focei.translate Translate the model to FOCEi and then run
##'     the tables and objective function so that different estimation
##'     methodologies are comparable through OBJF.
##' @return nlmixr fit object
##' @author Matthew L. Fidler
##' @export
nlmixr.fit <- function(uif, data, est="nlme", ..., focei.translate=TRUE){
    dat <- data;
    bad.focei <- "Problem calculating residuals, returning fit without residuals";
    if (est == "saem"){
        pt <- proc.time()
        model <- uif$saem.model
        cfg   = configsaem(model=model, data=dat, inits=uif$saem.init, ...);
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
        weight <- uif$nlme.var
        if (!is.null(uif$nmodel$lin.solved)){
            fit <- nlme_lin_cmpt(dat, par_model=specs,
                                 ncmt=uif$nmodel$lin.solved$ncmt,
                                 oral=uif$nmodel$lin.solved$oral,
                                 tlag=uif$nmodel$lin.solved$tlag,
                                 parameterization=uif$nmodel$lin.solved$parameterization,
                                 par_trans=fun,
                                 weight=weight,
                                 verbose=TRUE, ...);
        } else {
            fit <- nlme_ode(dat,
                            model=uif$rxode.pred,
                            par_model=specs,
                            par_trans=fun,
                            response="nlmixr_pred",
                            weight=weight,
                            verbose=TRUE,
                            ...);
        }
        ## Run FOCEi using same ETAs and THETA estimates to get
        ## comparable OBJFs and also extract table entries like
        ## CWRES.
        ## return(fit)
        if (focei.translate){
            ret <- try(as.focei(fit, uif, pt, data=dat))
            if (inherits(ret, "try-error")){
                warning(bad.focei)
                return(fit);
            } else {
                return(ret)
            }
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
