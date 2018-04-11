##' VPC based on ui model
##'
##' @param fit nlmixr fit object
##' @param n Number of VPC simulations.  By default 100
##' @inheritParams vpc::vpc
##' @inheritParams RxODE::rxSolve
##' @param ... Args sent to rxSolve
##' @author Matthew L. Fidler
##' @export
vpc.ui <- function(fit, n=100, bins = "jenks",
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
    con <- fit$fit$con;
    pt <- proc.time();
    message("Compiling VPC model...", appendLF=FALSE)
    mod <- gsub(rex::rex("(0)~"), "(0)=", paste0(gsub("=", "~", rxNorm(fit$model$pred.only), perl=TRUE),
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
    dat <- nlmixrData(getData(fit));
    if (is.null(method)){
        meth <- c("dop853", "lsoda", "liblsoda");
        meth <- meth[con$stiff + 1];
    } else {
        meth <- method;
    }
    sim <- rxSolve(mod, params=theta, events=dat, omega=omega, nStud=nStud, sigma=sigma, add.cov=TRUE, return.type="data.frame",
                   atol=con$atol.ode, rtol=con$rtol.ode, maxsteps=con$maxsteps.ode,
                   hmin = con$hmin, hmax = con$hmax, hini = con$hini, transit_abs = con$transit_abs,
                   maxordn = con$maxordn, maxords = con$maxords, method=meth);
    diff <- proc.time() - pt;
    message(sprintf("done (%.2f sec)", diff["elapsed"]));
    names(dat) <- tolower(names(dat))
    w <- which(duplicated(names(dat)));
    if (length(w) > 0){
        warning("Dropping duplicate columns (case insensitive)")
        dat <- dat[, -w];
    }
    if (!is.null(stratify)){
        cols <- c(tolower(stratify), "dv")
    }  else {
        cols <- c("dv");
    }
    dat <- dat[dat$evid == 0, ];
    ## Assume this is in the observed dataset. Add it to the current dataset
    if(!all(names(sim) %in% cols)){
        w <- cols[!(cols %in% names(sim))]
        n <- names(sim)
        sim <- cbind(sim, dat[, w]);
        names(sim) <- c(n, w);
        diff <- proc.time() - pt;
        pt <- proc.time();
    }
    rxDelete(mod);
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    call <- call[names(call) %in% formalArgs(getFromNamespace("vpc","vpc"))]
    call$obs_cols = list(id="id", dv="dv", idv="time")
    call$sim_cols = list(id="id", dv="dv", idv="time")
    call$stratify = stratify
    do.call(getFromNamespace("vpc","vpc"), c(list(sim=sim, obs=dat), call), envir = parent.frame(1))
}
