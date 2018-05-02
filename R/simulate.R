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
nlmixrSim <- function(object, params=NULL, events=NULL, inits = NULL, scale = NULL,
                      covs = NULL, method = c("liblsoda", "lsoda", "dop853"),
                      transit_abs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
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
    lst <- as.list(match.call()[-1]);
    mod <- gsub("rx_r_~.*;", "",
                gsub(rex::rex("(0)~"), "(0)=",
                     gsub("=", "~", RxODE::rxNorm(object$model$pred.only), perl=TRUE)));
    if (RxODE::rxIs(params, "rx.event")){
        if (!is.null(events)){
            tmp <- events;
            events <- params;
            params <- tmp;
        } else {
            events <- params;
            params <- NULL;
        }
    }
    nlmixr.data <- nlmixrData(getData(object))
    if (!RxODE::rxIs(events, "rx.events")){
        events <- nlmixr.data;
    }
    lst$events <- events;
    if (is.null(params)){
        params <- fixed.effects(object);
        names(params) <- sprintf("THETA[%d]", seq_along(params))
    }
    sim <- "ipred=rx_pred_;\nsim=ipred"
    err <- object$uif$err;
    w <- which(!is.na(object$uif$err))
    mat <- diag(length(w))
    dimn <- character(length(w))
    for (i in seq_along(w)){
        cur <- sprintf("THETA[%d]", w[i]);
        dimn[i] <- cur;
        if (err[w[i]] == "add"){
            sim <- paste0(sim, "+", cur);
            mat[i, i] <- params[i] ^ 2;
        } else if (err[w[i]] == "prop"){
            sim <- paste0(sim, "+ipred*", cur);
            mat[i, i] <- params[i] ^ 2;
        }
    }
    dimnames(mat) <- list(dimn, dimn);
    lst$object <- RxODE::RxODE(paste0(mod, "\n", sim));
    params <- params[-w];
    lst$params <- params;
    if (dfObs == 0.0){
        lst$dfObs <- nobs(object);
    }
    if (dfSub == 0.0){
        lst$dfSub <- length(unique(nlmixr.data$ID))
    }
    omega <- object$omega
    dm <- dim(object$omega)[1];
    n <- sprintf("ETA[%d]", seq(1, dm))
    dimnames(omega) <- list(n, n)
    lst$omega <- omega;
    sigma <- mat;
    lst$sigma <- sigma;
    thetaMat <- object$varFix;
    d <- dimnames(thetaMat)[[1]];
    n <- names(fixed.effects(object))
    for (i in seq_along(n)){
        d <- gsub(n[i], sprintf("THETA[%d]", i), d, fixed=TRUE)
    }
    dimnames(thetaMat) <- list(d, d);
    lst$thetaMat <- thetaMat;
    on.exit({RxODE::rxUnload(lst$object)});
    if (!any(names(lst) == "return.type")){
        lst$return.type <- "data.frame"
    }
    do.call(getFromNamespace("rxSolve","RxODE"), lst)
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
                       transit_abs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
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
        ret <- do.call(getFromNamespace("rxSolve","RxODE"), lst);
        if (do.ipred){
            names(ret) <- sub("pred", "ipred", names(ret))
        }
        return(ret)
    } else {
        lst$params <- pred.par
        ret.pred <- do.call(getFromNamespace("rxSolve","RxODE"), lst);
        lst$params <- ipred.par
        ret.pred$ipred <- do.call(getFromNamespace("rxSolve","RxODE"), lst)$pred
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
    dat.new <- do.call("nlmixrPred", lst)
    dat.new <- data.frame(dat.new[, 1:2], stack(dat.new[,-(1:2)]))
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
    do.call("nlmixrAugPred", lst)
}

##' @rdname nlmixrSim
##' @export
rxSolve.focei.fit <- function(object, params=NULL, events=NULL, inits = NULL, scale = NULL,
                              covs = NULL, method = c("liblsoda", "lsoda", "dop853"),
                              transit_abs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
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
    do.call("nlmixrSim", as.list(match.call()[-1]));
}

##' @rdname nlmixrSim
##' @export
simulate.focei.fit <- function(object, nsim=1, seed=NULL, ...){
    nlmixrSim(object, ..., seed=seed, nsim=nsim);
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
    do.call("nlmixrSim", lst, envir=parent.frame(1))
}

##' @importFrom RxODE rxSolve
##' @export
RxODE::rxSolve
