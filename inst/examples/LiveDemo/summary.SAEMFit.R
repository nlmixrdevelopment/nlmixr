summary.saemFit<-function (x, ...) 
{
    fit = x
    th = fit$Plambda
    nth = length(th)
    H = solve(fit$Ha[1:nth, 1:nth])
    se = sqrt(diag(H))
    m = cbind(exp(th), th, se)
    lhsVars = scan("LHS_VARS.txt", what = "", quiet = T)
    if (length(lhsVars) == nth) 
        dimnames(m)[[1]] = lhsVars
    dimnames(m)[[2]] = c("th", "log(th)", "se(log_th)")
    cat("THETA:\n")
    print(m)
    cat("\nOMEGA:\n")
    print(fit$Gamma2_phi1)
    cat("\nIIV:\n")
    print(sqrt(diag(fit$Gamma2_phi1)))
    if (any(fit$sig2 == 0)) {
        cat("\nSIGMA:\n")
        print(max(fit$sig2^2))
        cat("\nSIGMA (SD):\n")
        print(max(fit$sig2))
    }
    else {
        cat("\nARES & BRES:\n")
        print(fit$sig2)
    }
    invisible(list(theta = th, se = se, H = H, omega = fit$Gamma2_phi1, 
        eta = fit$mpost_phi))
   THETA<<-m
}
