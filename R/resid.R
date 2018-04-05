calc.resid.fit <- function(fit, data){
    etas <- random.effects(fit);
    thetas <- fixed.effects(fit);
    pars <- .Call(`_nlmixr_nlmixrParameters`, thetas, etas);
    print(pars$ipred)
    preds <- list(ipred=rxSolve(fit$env$model$inner, pars$ipred,data,return.type="data.frame"),
                  pred=rxSolve(fit$env$model$inner, pars$pred,data,return.type="data.frame"));
    lst <- .Call(`_nlmixr_nlmixrResid`, preds,fit$omega, fit$DV, etas, pars$eta.lst);
    ## Add Empirical Bayes Estimates
    df <- rxSolve(fit$env$model$ebe, pars$ipred,data,return.type="data.frame")[, -(1:2)];
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
