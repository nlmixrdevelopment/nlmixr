print.lme <- function (x, ...) 
{
    dd <- x$dims
    if (inherits(x, "nlme")) {
        cat("Nonlinear mixed-effects model fit by ")
        cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
        cat("  Model:", deparse(x$call$model), "\n")
    }
    else {
        cat("Linear mixed-effects model fit by ")
        cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
    }
    #cat("  Data:", deparse(x$call$data), "\n")
    if (!is.null(x$call$subset)) {
        cat("  Subset:", deparse(asOneSidedFormula(x$call$subset)[[2L]]), 
            "\n")
    }
    cat("  Log-", ifelse(x$method == "REML", "restricted-", ""), 
        "likelihood: ", format(x$logLik), "\n", sep = "")
    fixF <- x$call$fixed
    if (inherits(fixF, "formula") || is.call(fixF) || is.name(fixF)) {
        cat("  Fixed:", deparse(x$call$fixed), "\n")
    }
    else {
        cat("  Fixed:", deparse(lapply(fixF, function(el) as.name(deparse(el)))), 
            "\n")
    }
    print(fixef(x))
    cat("\n")
    print(summary(x$modelStruct), sigma = x$sigma)
    cat("Number of Observations:", dd[["N"]])
    cat("\nNumber of Groups: ")
    Ngrps <- dd$ngrps[1:dd$Q]
    if ((lNgrps <- length(Ngrps)) == 1) {
        cat(Ngrps, "\n")
    }
    else {
        sNgrps <- 1:lNgrps
        aux <- rep(names(Ngrps), sNgrps)
        aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps, 
            lNgrps))[!lower.tri(diag(lNgrps))])
        names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
        cat("\n")
        print(rev(Ngrps))
    }
    invisible(x)
}


print.summary.lme <- function (x, verbose = FALSE, ...) 
{
    dd <- x$dims
    verbose <- verbose || attr(x, "verbose")
    if (inherits(x, "nlme")) {
        cat("Nonlinear mixed-effects model fit by ")
        cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
        cat("  Model:", deparse(x$call$model), "\n")
    }
    else {
        cat("Linear mixed-effects model fit by ")
        cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
    }
    #cat(" Data:", deparse(x$call$data), "\n")
    if (!is.null(x$call$subset)) {
        cat("  Subset:", deparse(asOneSidedFormula(x$call$subset)[[2L]]), 
            "\n")
    }
    print(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = c(x$logLik), 
        row.names = " "))
    if (verbose) {
        cat("Convergence at iteration:", x$numIter, "\n")
    }
    cat("\n")
    print(summary(x$modelStruct), sigma = x$sigma, reEstimates = x$coef$random, 
        verbose = verbose)
    cat("Fixed effects: ")
    fixF <- x$call$fixed
    if (inherits(fixF, "formula") || is.call(fixF)) {
        cat(deparse(x$call$fixed), "\n")
    }
    else {
        cat(deparse(lapply(fixF, function(el) as.name(deparse(el)))), 
            "\n")
    }
    xtTab <- as.data.frame(x$tTable)
    wchPval <- match("p-value", names(xtTab))
    for (i in names(xtTab)[-wchPval]) {
        xtTab[, i] <- format(zapsmall(xtTab[, i]))
    }
    xtTab[, wchPval] <- format(round(xtTab[, wchPval], 4))
    if (any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 
        0))) {
        levels(xtTab[, wchPval])[wchLv] <- "<.0001"
    }
    row.names(xtTab) <- dimnames(x$tTable)[[1L]]
    print(xtTab)
    if (nrow(x$tTable) > 1) {
        corr <- x$corFixed
        class(corr) <- "correlation"
        print(corr, title = " Correlation:", ...)
    }
    cat("\nStandardized Within-Group Residuals:\n")
    print(x$residuals)
    cat("\nNumber of Observations:", x$dims[["N"]])
    cat("\nNumber of Groups: ")
    Ngrps <- dd$ngrps[1:dd$Q]
    if ((lNgrps <- length(Ngrps)) == 1) {
        cat(Ngrps, "\n")
    }
    else {
        sNgrps <- 1:lNgrps
        aux <- rep(names(Ngrps), sNgrps)
        aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps, 
            lNgrps))[!lower.tri(diag(lNgrps))])
        names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
        cat("\n")
        print(rev(Ngrps))
    }
    invisible(x)
}
