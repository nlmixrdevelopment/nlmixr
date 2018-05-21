##' VPC based on ui model
##'
##' @param fit nlmixr fit object
##' @param data this is the data to use to augment the VPC fit.  By
##'     default is the fitted data, (can be retrieved by
##'     \code{\link[nlme]{getData}}), but it can be changed by specifying
##'     this argument.
##' @param n Number of VPC simulations.  By default 100
##' @inheritParams vpc::vpc
##' @inheritParams RxODE::rxSolve
##' @param ... Args sent to \code{\link[RxODE]{rxSolve}}
##' @return Simulated dataset (invisibly)
##' @author Matthew L. Fidler
##' @export
vpc_ui <- function(fit, data=NULL, n=100, bins = "jenks",
                   n_bins = "auto", bin_mid = "mean",
                   show = NULL, stratify = NULL, pred_corr = FALSE,
                   pred_corr_lower_bnd = 0, pi = c(0.05, 0.95), ci = c(0.05, 0.95),
                   uloq = NULL, lloq = NULL, log_y = FALSE, log_y_min = 0.001,
                   xlab = NULL, ylab = NULL, title = NULL, smooth = TRUE, vpc_theme = NULL,
                   facet = "wrap", labeller = NULL, vpcdb = FALSE, verbose = FALSE, ...,
                   nStud=NULL,
                   method=NULL){
    if (is.null(nStud)){
        nStud <- n;
    }
    if (is.numeric(data) || is.integer(data)){
        nStud <- n
        data <- NULL
    }
    tmp <- list(...)
    if (!is.null(tmp$nsim)){
        nStud <- tmp$nsim
    }
    con <- fit$fit$con;
    pt <- proc.time();
    message("Compiling VPC model...", appendLF=FALSE)
    mod <- gsub(rex::rex("(0)~"), "(0)=", paste0(gsub("=", "~", RxODE::rxNorm(fit$model$pred.only), perl=TRUE),
                                                 "\ndv=rx_pred_+sqrt(rx_r_)*rx_err_"))
    mod <- RxODE(mod);
    diff <- proc.time() - pt;
    message(sprintf("done (%.2f sec)", diff["elapsed"]));
    pt <- proc.time();
    theta <- fixed.effects(fit);
    names(theta) <- sprintf("THETA[%d]", seq_along(theta))
    omega <- fit$omega
    dm <- dim(fit$omega)[1];
    n <- sprintf("ETA[%d]", seq(1, dm))
    dimnames(omega) <- list(n, n)
    sigma <- matrix(1, dimnames=list("rx_err_", "rx_err_"))
    if (is.null(data)){
        dat <- nlmixrData(getData(fit));
    } else {
        dat <- data
    }
    if (is.null(method)){
        meth <- c("dop853", "lsoda", "liblsoda");
        meth <- meth[con$stiff + 1];
    } else {
        meth <- method;
    }
    sim <- RxODE::rxSolve(mod, params=theta, events=dat, omega=omega, nStud=nStud, sigma=sigma, add.cov=TRUE, return.type="data.frame",
                          atol=con$atol.ode, rtol=con$rtol.ode, maxsteps=con$maxsteps.ode,
                          hmin = con$hmin, hmax = con$hmax, hini = con$hini, transit_abs = con$transit_abs,
                          maxordn = con$maxordn, maxords = con$maxords, method=meth);
    diff <- proc.time() - pt;
    message(sprintf("done (%.2f sec)", diff["elapsed"]));
    onames <- names(dat)
    names(dat) <- tolower(onames)
    w <- which(duplicated(names(dat)));
    if (length(w) > 0){
        warning(sprintf("Dropping duplicate columns (case insensitive): %s", paste(onames, collapse=", ")))
        dat <- dat[, -w];
    }
    if (!is.null(stratify)){
        cols <- c(tolower(stratify), "dv")
        stratify <- tolower(stratify);
    }  else {
        cols <- c("dv");
    }
    dat <- dat[dat$evid == 0, ];
    ## Assume this is in the observed dataset. Add it to the current dataset
    if(!all(names(sim) %in% cols)){
        w <- cols[!(cols %in% names(sim))]
        if (length(w) >= 1){
            n <- names(sim)
            sim <- cbind(sim, dat[, w, drop = FALSE]);
            names(sim) <- c(n, w);
        }
    }
    RxODE::rxUnload(mod);
    ns <- loadNamespace("vpc");
    if (exists("vpc_vpc",ns)){
        vpcn <- "vpc_vpc"
    } else {
        vpcn <- "vpc"
    }
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    call <- call[names(call) %in% methods::formalArgs(getFromNamespace(vpcn,"vpc"))]
    call$obs_cols = list(id="id", dv="dv", idv="time")
    call$sim_cols = list(id="id", dv="dv", idv="time")
    call$stratify = stratify
    p = do.call(getFromNamespace(vpcn,"vpc"), c(list(sim=sim, obs=dat), call), envir = parent.frame(1))
    print(p);
    return(invisible(sim));
}


##' @rdname vpc_ui
##' @export
vpc.nlmixr.ui.focei <- function(sim, ...){
    vpc_ui(fit=sim, ...);
}

##' @rdname vpc_ui
##' @export
vpc.nlmixr.ui.saem <- function(sim, ...){
    vpc_ui(fit=sim, ...);
}

##' @rdname vpc_ui
##' @export
vpc.nlmixr.ui.nlme <- function(sim, ...){
    vpc_ui(fit=sim, ...);
}

##' @rdname vpc_ui
##' @S3method vpc ui
##' @export vpc.ui
vpc.ui <- function(sim, ...){
    vpc_ui(fit=sim, ...);
}
