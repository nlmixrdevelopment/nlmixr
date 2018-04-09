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
                   nStud=NULL, method="lsoda"){
    if (is.null(nStud)){
        nStud <- n;
    }
    message("Compiling VPC model...", appendLF=FALSE)
    mod <- gsub(rex::rex("(0)~"), "(0)=", paste0(gsub("=", "~", rxNorm(fit$model$pred.only), perl=TRUE),
           "\ndv=rx_pred_+sqrt(rx_r_)*rx_err_"))
    mod <- RxODE(mod);
    message("done!");
    theta <- fixed.effects(fit);
    names(theta) <- sprintf("THETA[%d]", seq_along(theta))
    omega <- fit$omega
    dm <- dim(fit$omega)[1];
    n <- sprintf("ETA[%d]", seq(1, dm))
    dimnames(omega) <- list(n, n)
    sigma <- matrix(1, dimnames=list("rx_err_", "rx_err_"))
    dat <- getData(fit);
    sim <- rxSolve(mod, params=theta, events=dat, omega=omega, nStud=nStud, sigma=sigma, add.cov=TRUE, return.type="data.frame",
                   method=method);
    names(dat) <- tolower(names(dat))
    ## rxDelete(mod);
    vpc::vpc(sim = sim, obs = dat,
             bins = bins, n_bins = n_bins, bin_mid = bin_mid, obs_cols = list(id="id", dv="dv", idv="time"),
             sim_cols = list(id="id", dv="dv", idv="time"),
             software = "auto", show = show, stratify = stratify, pred_corr = pred_corr,
             pred_corr_lower_bnd = pred_corr_lower_bnd, pi = pi, ci = ci,
             uloq = uloq, lloq = lloq, log_y = log_y, log_y_min = log_y_min,
             xlab = xlab, ylab = ylab, title = title, smooth = smooth, vpc_theme = vpc_theme,
             facet = facet, labeller = labeller, vpcdb = vpcdb, verbose = verbose)
}
