.onAttach <- function(libname, pkgname){ ## nocov start
    ## Setup RxODE.prefer.tbl
    nlmixrSetupMemoize()
    options(keep.source=TRUE)
}

nlmixrSetupMemoize <- function(){
    reSlow <- rex::rex(".slow", end)
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
##' @importFrom n1qn1 n1qn1
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
##' @importFrom dparser dparse
##' @importFrom vpc vpc
##' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_line geom_abline xlab geom_smooth
##' @importFrom RcppArmadillo armadillo_version
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

armaVersion <- function(){
    nlmixrLogo(str="RcppArmadiilo", RcppArmadillo::armadillo_version())
}

##' nlmixr fits population PK and PKPD non-linear mixed effects models.
##'
##' nlmixr is an R package for fitting population pharmacokinetic (PK)
##' and pharmacokinetic-pharmacodynamic (PKPD) models.
##'
##' The nlmixr generalized function allows common access to the nlmixr
##' estimation routines.
##'
##' @template uif
##'
##' @param object Fitted object or function specifying the model.
##' @inheritParams nlmixr_fit
##' @param ... Other parameters
##' @return Either a nlmixr model or a nlmixr fit object
##' @author Matthew L. Fidler, Rik Schoemaker
##' @export
nlmixr <- function(object, data, est="nlme", control=list(), calc.resid=TRUE, ...){
    UseMethod("nlmixr")
}

##' @rdname nlmixr
##' @export
nlmixr.function <- function(object, data, est="nlme", control=list(), calc.resid=TRUE, ...){
    uif <- nlmixrUI(object);
    class(uif) <- "list";
    uif$nmodel$model.name <- deparse(substitute(object))
    if (missing(data) && missing(est)){
        class(uif) <- "nlmixrUI"
        return(uif)
    } else {
        uif$nmodel$data.name <- deparse(substitute(data))
        class(uif) <- "nlmixrUI"
        nlmixr_fit(uif, data, est, control=control, calc.resid=calc.resid, ...);
    }
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixrUI <- function(object, data, est="nlme", control=list(), ...){
    uif <- object
    if (missing(data) && missing(est)){
        return(uif)
    } else {
        class(uif) <- "list";
        uif$nmodel$data.name <- deparse(substitute(data))
        class(uif) <- "nlmixrUI"
        nlmixr_fit(uif, data, est, control=control, ...);
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
        nlmixr_fit(uif, data, est, ...);
    }
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixr.ui.focei.fit <- nlmixr.nlmixr.ui.nlme

##' @rdname nlmixr
##' @export
nlmixr.nlmixr.ui.saem <- nlmixr.nlmixr.ui.nlme
##' Convert/Format the data appropriately for nlmixr
##'
##' @param data is the name of the data to convert.  Can be a csv file
##'     as well.
##' @return Appropriately formatted data
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrData <- function(data){
    UseMethod("nlmixrData");
}
##' @export
##' @rdname nlmixrData
nlmixrData.character <- function(data){
    if (!file.exists(data)){
        stop(sprintf("%s does not exist.", data))
    }
    if (regexpr(rex::rex(".csv", end), data) != -1){
        return(nlmixrData.default(utils::read.csv(data, na.strings=c(".", "NA", "na", ""))))
    } else {
        stop(sprintf("Do not know how to read in %s", data));
    }
}
##' @export
##' @rdname nlmixrData
nlmixrData.default <- function(data){
    dat <- as.data.frame(data);
    nm1 <- toupper(names(dat));
    for (n in c("ID", "EVID", "TIME", "DV", "AMT")){
        w <- which(nm1 == n)
        if (length(w) == 1L){
            names(dat)[w] <- n;
        } else if (length(w) == 0L){
            stop(sprintf("Need '%s' data item in dataset.", n))
        } else {
            stop(sprintf("Multiple '%s' columns in dataset.", n))
        }
    }
    if (is(dat$ID, "factor")){
        dat$ID <- paste(dat$ID);
    }
    idSort <- .Call(`_nlmixr_chkSortIDTime`, as.integer(dat$ID), as.double(dat$TIME), as.integer(dat$EVID));
    if (idSort == 3L){
        warning("NONMEM-style data converted to nlmixr/RxODE-style data.");
        return(nlmixrData.default(nmDataConvert(dat)));
    }
    backSort <- c()
    backSort2 <- c();
    if (is(dat$ID, "character")){
        lvl <- unique(dat$ID);
        lab <- paste(lvl)
        dat$ID <- factor(dat$ID, levels=lvl, labels=lab);
        backSort <- levels(dat$ID);
        backSort2 <- seq_along(backSort)
        dat$ID <- as.integer(dat$ID);
    } else {
        if (idSort == 2L){
            lvl <- unique(dat$ID);
            lab <- paste(lvl)
            dat$ID <- factor(dat$ID, levels=lvl, labels=lab);
            backSort <- levels(dat$ID);
            backSort2 <- seq_along(backSort)
            dat$ID <- as.integer(dat$ID);
        } else if (idSort == 0L){
            warning("Sorting by ID, TIME; Output fit may not be in the same order as input dataset.")
            dat <- dat[order(dat$ID, dat$TIME), ];
            lvl <- unique(dat$ID);
            lab <- paste(lvl)
            dat$ID <- factor(dat$ID, levels=lvl, labels=lab);
            backSort <- levels(dat$ID);
            backSort2 <- seq_along(backSort)
            dat$ID <- as.integer(dat$ID);
        }
    }
    attr(dat, "backSort") <- backSort;
    attr(dat, "backSort2") <- backSort2;
    return(dat);
}

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
##' @param calc.resid Translate the model to FOCEi and then run
##'     the tables and objective function so that different estimation
##'     methodologies are comparable through OBJF.
##' @return nlmixr fit object
##' @author Matthew L. Fidler
##' @export
nlmixr_fit <- function(uif, data, est="nlme", control=list(), ...,
                       sum.prod=FALSE, calc.resid=TRUE){
    start.time <- Sys.time();
    dat <- nlmixrData(data);
    up.covs <- toupper(uif$all.covs);
    up.names <- toupper(names(dat))
    for (i in seq_along(up.covs)){
        w <- which(up.covs[i] == up.names)
        if (length(w) == 1){
            names(dat)[w] = uif$all.covs[i];
        }
    }
    backSort <- attr(dat, "backSort");
    backSort2 <- attr(dat, "backSort2");
    attr(dat, "backSort") <- NULL;
    attr(dat, "backSort2") <- NULL;
    uif$env$infusion <- .Call(`_nlmixr_chkSolvedInf`, as.double(dat$EVID), as.integer(!is.null(uif$nmodel$lin.solved)));
    bad.focei <- "Problem calculating residuals, returning fit without residuals.";
    fix.dat <- function(x){
        if (length(backSort) > 0){
            cls <- class(x);
            class(x) <- "data.frame";
            x$ID <- factor(x$ID, backSort2, labels=backSort);
            class(x) <- cls;
            return(x);
        } else {
            return(x);
        }
    }
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
        cfg <- configsaem(model=model, data=dat, inits=uif$saem.init,
                           mcmc=mcmc, ODEopt=ODEopt, seed=seed);
        if (print > 1){
            cfg$print <- as.integer(print)
        }
        fit <- model$saem_mod(cfg);
        if (calc.resid){
            ret <- try(as.focei(fit, uif, pt, data=dat));
            if (inherits(ret, "try-error")){
                warning(bad.focei)
                return(fit)
            } else {
                ret <- fix.dat(ret);
                env <- attr(ret, ".focei.env")
                assign("start.time", start.time, env);
                assign("est", est, env);
                assign("stop.time", Sys.time(), env);
                ## Now remove dlls
                ## sf <- ret$uif$env$saem.fit
                ## rx <- attr(sf, "rx")
                ## if (is(rx, "RxODE")){
                ##     rxDelete(rx);
                ## }
                ## cpp <- attr(sf, "saem.cpp")
                ## if (file.exists(cpp)){
                ##     try(unlink(cpp))
                ## }
                ## dll <- attr(sf, "saem.dll")
                ## try({dyn.unload(dll)}, silent=TRUE);
                ## if (file.exists(dll))
                ##     unlink(dll);
                ## rxDelete(ret$model$pred.only)
                ## rxDelete(ret$model$ebe)
                ## rxDelete(ret$model$inner)
                return(ret)
            }
        } else {
            return(fit);
        }
    } else if (est == "nlme" || est == "nlme.mu" || est == "nlme.mu.cov" || est == "nlme.free"){
        pt <- proc.time()
        est.type <- est;
        if (est == "nlme.free"){
            fun <- uif$nlme.fun;
            specs <- uif$nlme.specs;
        } else if (est == "nlme.mu"){
            fun <- uif$nlme.fun.mu;
            specs <- uif$nlme.specs.mu;
        } else if (est == "nlme.mu.cov"){
            fun <- uif$nlme.fun.mu.cov
            specs <- uif$nlme.specs.mu.cov;
        } else {
            if (!is.null(uif$nlme.fun.mu.cov)){
                est.type <- "nlme.mu.cov"
                fun <- uif$nlme.fun.mu.cov
                specs <- uif$nlme.specs.mu.cov;
            } else if (!is.null(uif$nlme.mu.fun)){
                est.type <- "nlme.mu"
                fun <- uif$nlme.fun.mu
                specs <- uif$nlme.specs.mu;
            } else {
                est.type <- "nlme.free"
                fun <- uif$nlme.fun
                specs <- uif$nlme.fun.specs;
            }
        }
        grp.fn <- uif$grp.fn;
        dat$nlmixr.grp <- factor(apply(dat, 1, function(x){
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
        class(fit) <- c(est.type, class(fit));

        ## Run FOCEi using same ETAs and THETA estimates to get
        ## comparable OBJFs and also extract table entries like
        ## CWRES.
        ## return(fit)
        if (calc.resid){
            ret <- as.focei(fit, uif, pt, data=dat)
            ## ret <- try(as.focei(fit, uif, pt, data=dat))
            ## if (inherits(ret, "try-error")){
            ##     warning(bad.focei)
            ##     return(fit);
            ## } else {
            ##     return(ret)
            ## }
            ret <- fix.dat(ret);
            env <- attr(ret, ".focei.env")
            assign("start.time", start.time, env);
            assign("est", est, env);
            assign("stop.time", Sys.time(), env);
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
                         control=control,
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
        fit <- fix.dat(fit);
        env$uif.new <- uif.new;
        class(fit) <- c("nlmixr.ui.focei.fit", class(fit));
        if (uif$.clean.dll){
            focei.cleanup(fit)
        }
        assign("start.time", start.time, env);
        assign("est", est, env);
        assign("stop.time", Sys.time(), env);
        return(fit);
    } else if (est == "posthoc"){
        control.posthoc <- list(NOTRUN=TRUE,
                                ## inits.mat=mat,
                                cores=1,
                                find.best.eta=TRUE,
                                numeric=(!is.null(uif$nmodel$lin.solved)),
                                sum.prod=uif$env$sum.prod);
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
                         control=control.posthoc,
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
        fit <- fix.dat(fit);
        env$uif.new <- uif.new;
        class(fit) <- c("nlmixr.ui.focei.posthoc", class(fit));
        assign("start.time", start.time, env);
        assign("est", est, env);
        assign("stop.time", Sys.time(), env);
        return(fit);
    }
}
##' Control Options for SAEM
##'
##' @param seed Random Seed for SAEM step.  (Needs to be set for
##'     reproducibility.)  By default this is 99.
##'
##' @param nBurn Number of iterations in the Stochastic Approximation
##'     (SA), or burn-in step. This is equivalent to Monolix's \code{K_0} or \code{K_b}.
##'
##' @param nEm Number of iterations in the Expectation-Maximization
##'     (EM) Step. This is equivalent to Monolix's \code{K_1}.
##'
##' @param nmc Number of Markov Chains. By default this is 3.  When
##'     you increase the number of chains the numerical integration by
##'     MC method will be more accurate at the cost of more
##'     computation.  In Monolix this is equivalent to \code{L}
##'
##' @param nu This is a vector of 3 integers. They represent the
##'     numbers of transitions of the three different kernels used in
##'     the Hasting-Metropolis algorithm.  The default value is \code{c(2,2,2)},
##'     representing 40 for each transition initially (each value is
##'     multiplied by 20).
##'
##'     The first value represents the initial number of multi-variate
##'     Gibbs samples are taken from a normal distribution.
##'
##'     The second value represents the number of uni-variate, or multi-
##'     dimensional random walk Gibbs samples are taken.
##'
##'     The third value represents the number of bootstrap/reshuffling or
##'     uni-dimensional random samples are taken.
##'
##' @inheritParams RxODE::rxSolve
##' @param print The number it iterations that are completed before
##'     anything is printed to the console.  By default, this is 1.
##' @param ... Other arguments to control SAEM.
##' @return List of options to be used in \code{\link{nlmixr}} fit for
##'     SAEM.
##' @author Wenping Wang & Matthew L. Fidler
##' @export
saemControl <- function(seed=99,
                        nBurn=200, nEm=300,
                        nmc=3,
                        nu=c(2, 2, 2),
                        atol = 1e-06,
                        rtol = 1e-04,
                        stiff = TRUE,
                        transitAbs = FALSE,
                        print=1,
                        ...){
    .xtra <- list(...);
    .rm <- c();
    if (missing(transitAbs) && !is.null(.xtra$transit_abs)){
        transitAbs <- .xtra$transit_abs;
        .rm <- c(.rm, "transit_abs")
    }
    if (missing(nBurn) && !is.null(.xtra$n.burn)){
        nBurn <- .xtra$n.burn;
        .rm <- c(.rm, "n.burn")
    }
    if (missing(nEm) && !is.null(.xtra$n.em)){
        nEm <- .xtra$n.em;
        .rm <- c(.rm, "n.em")
    }
    .ret <- list(mcmc=list(niter=c(nBurn, nEm), nmc=nmc, nu=nu),
         ODEopt=list(atol=atol, rtol=rtol,
                     stiff=as.integer(stiff),
                     transitAbs = as.integer(transitAbs)),
         seed=seed,
         print=print, ...)
    if (length(.rm) > 0){
        .ret <- .ret[!(names(.ret) %in% .rm)]
    }
    .ret
}
