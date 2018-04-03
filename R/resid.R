calc.resid.fit <- function(fit, data){
    etas <- random.effects(fit);
    pars <- .Call(`_nlmixr_nlmixrParameters`, fixed.effects(fit), etas);
    preds <- list(ipred=rxSolve(fit$env$model$inner, pars$ipred,data,return.type="data.frame"),
                  pred=rxSolve(fit$env$model$inner, pars$pred,data,return.type="data.frame"));
    return(.Call(`_nlmixr_nlmixrResid`, preds,fit$omega, fit$DV, etas, pars$eta.lst));
}
