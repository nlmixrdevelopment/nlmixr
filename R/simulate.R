##' Simulate a nlmixr solved system
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
    lst$object <- RxODE::RxODE(gsub(rex::rex("(0)~"), "(0)=",
                                    paste0(gsub("=", "~", RxODE::rxNorm(object$model$pred.only), perl=TRUE),
                                           "\nipred=rx_pred_\nsim=ipred+sqrt(rx_r_)*rx_err_")));
    if (RxODE::rxIs(params, "rx.event")){
        if (!is.null(events)){
            tmp <- events;
            events <- params;
            params <- events;
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
    sigma <- matrix(1, dimnames=list("rx_err_", "rx_err_"))
    lst$sigma <- sigma;
    thetaMat <- object$varFix;
    d <- dimnames(thetaMat)[[1]];
    n <- names(fixed.effects(object))
    for (i in seq_along(n)){
        d <- gsub(n[i], sprintf("THETA[%d]", i), d, fixed=TRUE)
    }
    dimnames(thetaMat) <- list(d, d);
    lst$thetaMat <- thetaMat;
    on.exit({RxODE::rxDelete(lst$object)});
    do.call(getFromNamespace("rxSolve","RxODE"), lst)
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
