calc.resid.fit <- function(fit, data){
    etas <- random.effects(fit);
    thetas <- fixed.effects(fit);
    pars <- .Call(`_nlmixr_nlmixrParameters`, fixed.effects(fit), etas);
    preds <- list(ipred=rxSolve(fit$env$model$inner, pars$ipred,data,return.type="data.frame"),
                  pred=rxSolve(fit$env$model$inner, pars$pred,data,return.type="data.frame"));
    lst <- .Call(`_nlmixr_nlmixrResid`, preds,fit$omega, fit$DV, etas, pars$eta.lst);
    ## Add Empirical Bayes Estimates
    df <- rxSolve(fit$env$model$ebe, pars$ipred,data,return.type="data.frame")[, -(1:2)];
    df <- df[, !(names(df) %in% (c("nlmixr_pred", names(etas), names(thetas))))]
    lst[[3]] <- df;
    return(lst);
}
