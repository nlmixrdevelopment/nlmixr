nner.cppcalc.resid.fit <- function(fit, data, con){
    etas <- random.effects(fit);
    thetas <- fixed.effects(fit);
    pars <- .Call(`_nlmixr_nlmixrParameters`, thetas, etas);
    meth <- c("dop853", "lsoda", "liblsoda");
    meth <- meth[con$stiff + 1];
    preds <- list(ipred=RxODE::rxSolve(fit$env$model$inner, pars$ipred,data,return.type="data.frame",
                                atol=con$atol.ode, rtol=con$rtol.ode, maxsteps=con$maxsteps.ode,
                                hmin = con$hmin, hmax = con$hmax, hini = con$hini, transit_abs = con$transit_abs,
                                maxordn = con$maxordn, maxords = con$maxords, method=meth),
                  pred=RxODE::rxSolve(fit$env$model$inner, pars$pred,data,return.type="data.frame",
                               atol=con$atol.ode, rtol=con$rtol.ode, maxsteps=con$maxsteps.ode,
                               hmin = con$hmin, hmax = con$hmax, hini = con$hini, transit_abs = con$transit_abs,
                               maxordn = con$maxordn, maxords = con$maxords, method=meth));
    lst <- .Call(`_nlmixr_nlmixrResid`, preds,fit$omega, fit$DV, etas, pars$eta.lst);
    ## Add Empirical Bayes Estimates
    df <- RxODE::rxSolve(fit$env$model$ebe, pars$ipred,data,return.type="data.frame",
                  hmin = con$hmin, hmax = con$hmax, hini = con$hini, transit_abs = con$transit_abs,
                  maxordn = con$maxordn, maxords = con$maxords, method=(c("lsoda", "dop853", "liblsoda"))[con$stiff])[, -(1:2)];
    df <- df[, !(names(df) %in% (c("nlmixr_pred", names(etas), names(thetas))))]
    if (any(names(df) %in% names(lst[[1]]))){
        warning("Calculated residuals like IPRED are masked by nlmixr calculated values");
        df <- df[, !(names(df) %in% names(lst[[1]]))];
    }
    if (any(names(df) %in% names(lst[[3]]))){
        df <- df[, !(names(df) %in% names(lst[[3]]))];
    }
    lst[[4]] <- df;
    return(lst);
}
