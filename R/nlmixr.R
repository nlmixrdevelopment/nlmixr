##' nlmixr fits population PK and PKPD non-linear mixed effects models.
##'
##' nlmixr is an R package for fitting population pharmacokinetic (PK)
##' and pharmacokinetic-pharmacodynamic (PKPD) models.
##' importFrom(Rcpp, evalCpp)
##' @importFrom brew brew
##' @importFrom lattice xyplot
##' @importFrom lattice trellis.par.get
##' @importFrom nlme nlme
##' @importFrom nlme groupedData
##' @importFrom nlme getData
##' @importFrom nlme pdDiag
##' @importFrom RxODE RxODE
##' @importFrom graphics abline lines matplot plot points title
##' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var AIC BIC asOneSidedFormula coef end fitted resid setNames start simulate
##' @importFrom utils assignInMyNamespace getFromNamespace head stack sessionInfo
##' @importFrom parallel mclapply
##' @importFrom lbfgs lbfgs
##' @importFrom methods is
##' @importFrom Rcpp evalCpp
##' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_line geom_abline
##' @useDynLib nlmixr, .registration=TRUE
"_PACKAGE"

rex::register_shortcuts("nlmixr");
## GGplot use and other issues...
utils::globalVariables(c("DV", "ID", "IPRED", "IRES", "PRED", "TIME", "grp", "initCondition", "values"));

nlmixr.logo <- "         _             _             \n        | | %9s (_) %s\n  _ __  | | _ __ ___   _ __  __ _ __\n | '_ \\ | || '_ ` _ \\ | |\\ \\/ /| '__|\n | | | || || | | | | || | >  < | |\n |_| |_||_||_| |_| |_||_|/_/\\_\\|_|\n"

##' Messages the nlmixr logo...
##'
##' @param str String to print
##' @author Matthew L. Fidler
nlmixrLogo <- function(str="", version=sessionInfo()$otherPkgs$nlmixr$Version){
    message(sprintf(nlmixr.logo, str, version));
}
##' nlmixr generalized function
##'
##' @param object Fitted object or
##' @param data Data for fit.
##' @param est Estimation routine.
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

##' Fit a nlmixr model
##'
##' @param uif Parsed nlmixr model (by \code{nlmixr(mod.fn)}).
##' @param data Dataset to estimate
##' @param est Estimation method
##' @param ... Parameters passed to estimation method.
##' @return nlmixr fit object
##' @author Matthew L. Fidler
##' @export
nlmixr.fit <- function(uif, data, est="nlme", ...){
    if (est == "nlme"){
        fun <- uif$nlme.fun;
        specs <- uif$nlme.specs;
        grp.fn <- uif$grp.fn
        dat <- data;
        dat$nlmixr.grp <- factor(apply(data, 1, function(x){
            cur <- x;
            names(cur) <- names(dat);
            with(as.list(cur), {
                return(grp.fn())
            })
        }));
        weight <- uif$nlme.var
        if (!is.null(uif$nmodel$lin.solved)){
            fit <- nlme_lin_cmpt(dat, par_model=specs,
                                 ncmt=uif$nmodel$lin.solved$ncmt,
                                 oral=uif$nmodel$lin.solved$oral,
                                 tlag=uif$nmodel$lin.solved$tlag,
                                 parameterization=uif$nmodel$lin.solved$parameterization,
                                 par_trans=fun,
                                 weight=weight,
                                 verbose=TRUE, ...)
            return(fit)
        } else {
            fit <- nlme_ode(dat,
                            model=uif$rxode.pred,
                            par_model=specs,
                            par_trans=fun,
                            response="nlmixr_pred",
                            weight=weight,
                            verbose=TRUE,
                            control = nlmeControl(pnlsTol = .01, msVerbose = TRUE),
                            ...)
            return(fit)
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
                         ...)
    }
}
