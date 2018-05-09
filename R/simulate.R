## Add RxODE THETA/ETA replacement mini DSL
repThetaEta <- function(x, theta=c(), eta=c()){
    ret <- eval(parse(text=sprintf("quote({%s})", x)));
    f <- function(x){
        if (is.atomic(x)) {
            return(x)
        } else if (is.name(x)) {
            return(x)
        } else if (is.pairlist(x)){
            return(x)
        } else if (is.call(x)) {
            if (identical(x[[1]], quote(`[`))){
                type <- tolower(as.character(x[[2]]))
                if (type == "theta"){
                    return(eval(parse(text=sprintf("quote(%s)", theta[as.numeric(x[[3]])]))))
                } else if (type == "eta"){
                    return(eval(parse(text=sprintf("quote(%s)", eta[as.numeric(x[[3]])]))))
                }
                stop("Only theta/eta translation supported.");
            } else {
                return(as.call(lapply(x, f)))
            }
        } else {
            stop("Don't know how to handle type ", typeof(x),
                 call. = FALSE)
        }
    }
    ret <- deparse(f(ret))[-1];
    ret <- paste(ret[-length(ret)], collapse="\n");
    return(RxODE::rxNorm(RxODE::rxGetModel(ret)))
}



##' Simulate a nlmixr solved system
##'
##' This takes the uncertainty in the model parameter estimates and to
##' simulate a number of "studies".  Each study simulates a
##' realization of the parameters from the uncertainty in the fixed
##' parameter estimates.  In addition the omega and sigma matrices are
##' simulated from the uncertainty in the Omega/Sigma matrices based
##' on the number of subjects and observations the model was based on.
##'
##' @inheritParams RxODE::rxSolve
##'
##' @export
nlmixrSim <- function(object, events=NULL, inits = NULL, scale = NULL,
                      covs = NULL, method = c("liblsoda", "lsoda", "dop853"),
                      transit_abs = NULL, atol = 1.0e-6, rtol = 1.0e-4,
                      maxsteps = 5000L, hmin = 0L, hmax = NULL, hini = 0L, maxordn = 12L, maxords = 5L, ...,
                      cores, covs_interpolation = c("linear", "locf", "nocb", "midpoint"),
                      add.cov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                      nCoresRV = 1L, sigmaIsChol = FALSE, nDisplayProgress=10000L,
                      amountUnits = NA_character_, timeUnits = "hours", stiff,
                      theta = NULL, eta = NULL, addDosing=FALSE, update.object=FALSE,do.solve=TRUE,
                      omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                      nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                      nStud = 1L, dfSub=0.0, dfObs=0.0, return.type=c("rxSolve", "matrix", "data.frame"),
                      seed=NULL, nsim=NULL){
    ## Get around lazy evaluation issues.
    force(events)
    method <- match.arg(method)
    covs_interpolation <- match.arg(covs_interpolation)
    return.type <- match.arg(return.type)
    mod <- gsub("rx_pred_~", "ipred=",
                gsub("rx_r_~.*;", "",
                     gsub(rex::rex("(0)~"), "(0)=",
                          gsub("=", "~", RxODE::rxNorm(object$model$pred.only), perl=TRUE))));
    nlmixr.data <- nlmixrData(getData(object))
    if (!RxODE::rxIs(events, "rx.event")){
        events <- nlmixr.data;
    }
    params <- fixed.effects(object);
    theta.n <- names(params);
    sim <- "\nsim=ipred"
    err <- object$uif$err;
    w <- which(!is.na(object$uif$err))
    mat <- diag(length(w))
    dimn <- character(length(w))
    for (i in seq_along(w)){
        ntheta <- object$uif$ini$ntheta[w[i]]
        cur <- theta.n[ntheta];
        dimn[i] <- cur;
        if (err[w[i]] == "add"){
            sim <- paste0(sim, "+", cur);
            mat[i, i] <- params[ntheta] ^ 2;
            params[ntheta] <- NA_real_;
        } else if (err[w[i]] == "prop"){
            sim <- paste0(sim, "+ipred*", cur);
            mat[i, i] <- params[ntheta] ^ 2;
            params[ntheta] <- NA_real_;
        }
    }
    params <- params[!is.na(params)]
    dimnames(mat) <- list(dimn, dimn);
    sigma <- mat;
    mod <- paste0(mod, "\n", sim);
    if (is.null(omega)){
        omega <- object$omega
        dm <- dim(object$omega)[1];
        eta.n <- dimnames(omega)[[1]]
        omega <- omega;
    }
    new.mod <- repThetaEta(mod, theta=theta.n, eta=eta.n);
    message("Compiling model...", appendLF=FALSE)
    newobj <- RxODE::RxODE(new.mod);
    on.exit({RxODE::rxUnload(newobj)});
    message("done");
    params <- params[-w];
    if (dfObs == 0.0){
        dfObs <- nobs(object);
    }
    if (dfSub == 0.0){
        dfSub <- length(unique(nlmixr.data$ID))
    }
    if (is.null(thetaMat) && nStud > 1L){
        thetaMat <- object$varFix;
    }
    ## if (!any(names(lst) == "return.type")){
    ##     lst$return.type <- "data.frame"
    ## }
    RxODE::rxSolve(object=newobj, params=params, events=events, inits=inits,
                   scale=scale, covs=covs, method=method, transit_abs=transit_abs,
                   atol=atol, rtol=rtol, maxsteps=maxsteps, hmin=hmin, hmax=hmax,
                   hini = hini, maxordn = maxordn, maxords = maxords, ...,
                   cores=cores, covs_interpolation = covs_interpolation, add.cov = add.cov,
                   matrix = matrix,
                   sigma = sigma, sigmaDf = sigmaDf, nCoresRV = nCoresRV, sigmaIsChol = sigmaIsChol,
                   nDisplayProgress = nDisplayProgress, amountUnits = amountUnits, timeUnits = timeUnits,
                   stiff=stiff, theta = theta, eta = eta, addDosing = addDosing, update.object = update.object,
                   do.solve = do.solve, omega = omega, omegaDf = omegaDf, omegaIsChol = omegaIsChol,
                   nSub = nSub, thetaMat = thetaMat, thetaDf = thetaDf, thetaIsChol = thetaIsChol,
                   nStud = nStud, dfSub = dfSub, dfObs = dfObs, return.type = return.type, seed = seed, nsim = nsim)
}
##' Predict a nlmixr solved system
##'
##' This takes uncertainty in the model and sets it to 0 and solve the
##' event information supplied.  When \code{ipred} is \code{TRUE}, set
##' residual predictions to 0.  When \code{ipred} is \code{FALSE}, set
##' residual predations.  When \code{ipred} is \code{NA}, do both
##' individual and population predictions.
##'
##' @inheritParams RxODE::rxSolve
##'
##' @export
nlmixrPred <- function(object, params=NULL, events=NULL, inits = NULL, scale = NULL,
                       covs = NULL, method = c("liblsoda", "lsoda", "dop853"),
                       transit_abs = NULL, atol = 1.0e-6, rtol = 1.0e-4,
                       maxsteps = 5000L, hmin = 0L, hmax = NULL, hini = 0L, maxordn = 12L, maxords = 5L, ...,
                       cores, covs_interpolation = c("linear", "locf", "nocb", "midpoint"),
                       add.cov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                       nCoresRV = 1L, sigmaIsChol = FALSE, nDisplayProgress=10000L,
                       amountUnits = NA_character_, timeUnits = "hours", stiff,
                       theta = NULL, eta = NULL, addDosing=FALSE, update.object=FALSE,do.solve=TRUE,
                       omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                       nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                       nStud = 1L, dfSub=0.0, dfObs=0.0, return.type=c("data.frame", "rxSolve", "matrix"),
                       seed=NULL, nsim=NULL,
                       ipred=FALSE){
    lst <- as.list(match.call()[-1]);
    if (RxODE::rxIs(lst$params, "rx.event")){
        if (!is.null(lst$events)){
            tmp <- lst$events;
            lst$events <- lst$params;
            lst$params <- tmp;
        } else {
            lst$events <- lst$params;
            lst$params <- NULL;
        }
    }
    if (!RxODE::rxIs(lst$events, "rx.event")){
        lst$events <- nlmixrData(getData(object));
    }
    message("Compiling model...", appendLF=FALSE)
    lst$object <- RxODE::RxODE(gsub(rex::rex("(0)~"), "(0)=",
                                    paste0(gsub("=", "~", RxODE::rxNorm(object$model$pred.only), perl=TRUE),
                                           "\npred=rx_pred_")));
    message("done")
    params <- fixed.effects(object);
    names(params) <- sprintf("THETA[%d]", seq_along(params))
    do.ipred <- FALSE
    do.pred <- FALSE
    if (is.na(ipred)) {
        do.ipred <- TRUE
        do.pred <- TRUE
    } else if (ipred){
        do.ipred <- TRUE
    } else {
        do.pred <- TRUE
    }
    if (do.ipred){
        re <- random.effects(object)[,-1];
        names(re) <- sprintf("ETA[%d]", seq_along(names(re)));
        ipred.par <- data.frame(re,t(params),
                                rx_err_=0, check.names=FALSE)
    }
    if (do.pred){
        neta <- dim(object$omega)[1]
        pred.par <- c(params, setNames(rep(0, neta + 1), c(sprintf("ETA[%d]", seq(1, neta)), "rx_err_")));
    }
    on.exit({RxODE::rxUnload(lst$object)});
    lst$return.type <- match.arg(return.type);
    if (!is.na(ipred)){
        if (do.pred){
            lst$params <- pred.par
        } else {
            lst$params <- ipred.par
        }
        ret <- do.call(getFromNamespace("rxSolve","RxODE"), lst, envir=parent.frame(2));
        if (do.ipred){
            names(ret) <- sub("pred", "ipred", names(ret))
        }
        return(ret)
    } else {
        lst$params <- pred.par
        ret.pred <- do.call(getFromNamespace("rxSolve","RxODE"), lst, envir=parent.frame(2));
        lst$params <- ipred.par
        ret.pred$ipred <- do.call(getFromNamespace("rxSolve","RxODE"), lst, envir=parent.frame(2))$pred
        return(ret.pred)
    }
}
##' @rdname nlmixrPred
##' @export
predict.focei.fit <- function(object, ...){
    nlmixrPred(object, ...)
}

##' Augmented Prediction for nlmixr fit
##'
##'
##' @param object Nlmixr fit object
##' @inheritParams nlme::augPred
##' @inheritParams RxODE::rxSolve
##' @return Stacked data.frame with observations, individual/population predictions.
##' @author Matthew L. Fidler
nlmixrAugPred <- function(object, ..., covs_interpolation = c("linear", "locf", "nocb", "midpoint"),
                          primary=NULL, minimum = NULL, maximum = NULL, length.out = 51L){
    if (!inherits(object, "focei.fit")){
        stop("Need a nlmixr fit object")
    }
    uif <- object$uif
    dat <- nlmixrData(getData(object))
    up.covs <- toupper(uif$all.covs);
    up.names <- toupper(names(dat))
    for (i in seq_along(up.covs)){
        w <- which(up.covs[i] == up.names)
        if (length(w) == 1){
            names(dat)[w] = uif$all.covs[i];
        }
    }
    r <- range(dat$TIME)
    if (is.null(minimum)){
        minimum <- r[1];
    }
    if (is.null(maximum)){
        maximum <- r[2];
    }
    new.time <- sort(unique(c(seq(minimum, maximum, length.out=length.out), dat$TIME)));
    ids <- unique(dat$ID)
    new.pts <- expand.grid(TIME=new.time, ID=ids);
    ## Add covariates in the augmented prediction
    covsi <- match.arg(covs_interpolation)
    all.covs <- uif$all.covs
    if (length(all.covs) > 0){
        fs <- c(locf=0, nocb=1, midpoint=0.5, linear=0)
        new.cov<- lapply(all.covs, function(cov){
            unlist(lapply(ids, function(id){
                dat.id <- dat[dat$ID == id, ];
                fun <- approxfun(dat.id$TIME, dat.id[[cov]], method=ifelse(covsi == "linear", "linear", "constant"),
                                 rule=2,
                                 f=fs[covsi]);
                return(fun(new.time))
            }))})
        names(new.cov) <- all.covs;
        new.pts <- cbind(new.pts, all.covs);
    }
    new.pts$EVID <- 0
    new.pts$AMT <- 0
    dat.old <- dat;
    dat <- rbind(dat[, names(new.pts)], new.pts);
    dat <- dat[order(dat$ID, dat$TIME), ];
    dat <- dat[!duplicated(paste(dat$ID, dat$TIME)), ];
    lst <- as.list(match.call()[-1])
    lst <- lst[!(names(lst) %in% c("primary", "minimum", "maximum", "length.out"))]
    lst$ipred <- NA
    lst$events <- dat
    lst$params <- NULL;
    dat.new <- do.call("nlmixrPred", lst, envir=parent.frame(2))
    dat.new <- data.frame(dat.new[, 1:2], stack(dat.new[,-(1:2)]))
    levels(dat.new$ind) <- gsub("pred", "Population", gsub("ipred", "Individual", levels(dat.new$ind)))
    names(dat.old) <- tolower(names(dat.old))
    dat.old <- dat.old[dat.old$evid == 0, ];
    dat.old <- data.frame(id=dat.old$id, time=dat.old$time, values=dat.old$dv, ind="Observed");
    return(rbind(dat.new, dat.old))
}

##' @rdname nlmixrAugPred
##' @export
augPred.focei.fit <- function(object, primary = NULL, minimum = min(primary), maximum = max(primary),
                              length.out = 51, ...){
    lst <- as.list(match.call()[-1])
    ret <- do.call("nlmixrAugPred", lst, envir=parent.frame(2))
    class(ret) <- c("nlmixrAugPred", "data.frame")
    return(ret)
}

##' @rdname nlmixrAugPred
##' @export
plot.nlmixrAugPred <- function(x, y, ...){
    ids <- unique(x$id)
    for (i  in seq(1, length(ids), by=16)){
        tmp <- ids[seq(i, i + 15)]
        tmp <- tmp[!is.na(tmp)];
        d1 <- x[x$id %in% tmp, ];
        dobs <- d1[d1$ind == "Observed", ];
        dpred <- d1[d1$ind != "Observed", ];
        p3 <- ggplot(d1,aes(time, values,col=ind)) + geom_line(data=dpred, size=1.2) +
            geom_point(data=dobs) + facet_wrap(~id)
        print(p3)
    }
}

##' @rdname nlmixrSim
##' @export
rxSolve.focei.fit <- function(object, events=NULL, inits = NULL, scale = NULL,
                              covs = NULL, method = c("liblsoda", "lsoda", "dop853"),
                              transit_abs = NULL, atol = 1.0e-6, rtol = 1.0e-4,
                              maxsteps = 5000L, hmin = 0L, hmax = NULL, hini = 0L, maxordn = 12L, maxords = 5L, ...,
                              cores, covs_interpolation = c("linear", "locf", "nocb", "midpoint"),
                              add.cov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                              nCoresRV = 1L, sigmaIsChol = FALSE, nDisplayProgress=10000L,
                              amountUnits = NA_character_, timeUnits = "hours", stiff,
                              theta = NULL, eta = NULL, addDosing=FALSE, update.object=FALSE,do.solve=TRUE,
                              omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                              nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                              nStud = 1L, dfSub=0.0, dfObs=0.0, return.type=c("rxSolve", "matrix", "data.frame"),
                              seed=NULL, nsim=NULL){
    do.call("nlmixrSim", as.list(match.call()[-1]), envir=parent.frame(2));
}

##' @rdname nlmixrSim
##' @export
simulate.focei.fit <- function(object, nsim=1, seed=NULL, ...){
    nlmixrSim(object, ..., nsim=nsim, seed=seed)
}

##' @rdname nlmixrSim
##' @export
solve.focei.fit <- function(a, b, ...){
    lst <- as.list(match.call()[-1])
    n <- names(lst)
    if (!missing(a)){
        n[n == "a"] <- "";
    }
    if (!missing(b)){
        n[n == "b"] <- "";
    }
    names(lst) <- n
    do.call("nlmixrSim", lst, envir=parent.frame(2))
}

##' @importFrom RxODE rxSolve
##' @export
RxODE::rxSolve
