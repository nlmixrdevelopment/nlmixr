##' @export
anova.nlmixr.ui.focei.fit <- function(object, ..., test = TRUE, type = c("sequential", "marginal"),
                                      adjustSigma = TRUE, Terms, L, verbose = FALSE){
    ancall <- sys.call()
    ancall$verbose <- ancall$test <- ancall$type <- NULL
    object <- list(object, ...)
    termsClass <- vapply(object, data.class, "")
    valid.cl <- c("nlmixr.ui.saem", "nlmixr.ui.nlme", "nlmixr.ui.focei.fit");
    if (!all(match(termsClass, valid.cl, 0))) {
        valid.cl <- paste0("\"", valid.cl, "\"")
        stop(gettextf("objects must inherit from classes %s, or %s",
                      paste(head(valid.cl, -1), collapse = ", "), tail(valid.cl,
                                                                       1)), domain = NA)
    }
    rt <- length(object)
    ## FIXME different responses...?
    ## resp <- vapply(object, function(el) deparse(getResponseFormula(el)[[2L]]),
    ##                "")
    ## subs <- as.logical(match(resp, resp[1L], FALSE))
    ## if (!all(subs))
    ##     warning("some fitted objects deleted because response differs from the first model")
    ## if (sum(subs) == 1)
    ##     stop("first model has a different response from the rest")
    ## object <- object[subs]
    ##
    ## termsModel <- lapply(object, function(el) formula(el)[-2])
    ## Fixme estMeth
    ## estMeth <- vapply(object, function(el) if (is.null(val <- el[["method"]]))
    ##                                            NA_character_
    ##                                        else val, "")
    ## FIXME: different estimation methods?
    ## if (length(uEst <- unique(estMeth[!is.na(estMeth)])) >
    ##     1) {
    ##     stop("all fitted objects must have the same estimation method")
    ## }
    ## estMeth[is.na(estMeth)] <- uEst
    ## REML <- uEst == "REML"
    ## if (REML) {
    ##     aux <- vapply(termsModel, function(el) {
    ##         tt <- terms(el)
    ##         val <- paste(sort(attr(tt, "term.labels")), collapse = "&")
    ##         if (attr(tt, "intercept") == 1)
    ##             paste(val, "(Intercept)", sep = "&")
    ##         else val
    ##     }, ".")
    ##     if (length(unique(aux)) > 1) {
    ##         warning("fitted objects with different fixed effects. REML comparisons are not meaningful.")
    ##     }
    ## }
    ## termsCall <- lapply(object, function(el) {
    ##     if (is.null(val <- el$call) && is.null(val <- attr(el,
    ##                                                        "call")))
    ##         stop("objects must have a \"call\" component or attribute")
    ##     val
    ## })
    ## termsCall <- vapply(termsCall, function(el) paste(deparse(el),
    ##                                                   collapse = ""), "")
    aux <- lapply(object, logLik)
    if (length(unique(sapply(object, nobs))) > 1) {
        stop("all fitted objects must use the same number of observations")
    }
    dfModel <- vapply(aux, attr, 1, "df")
    logLik <- vapply(aux, c, 1.1)
    aod <- data.frame(Model = 1:rt, df = dfModel,
                      AIC = vapply(aux, AIC, 1),
                      BIC = vapply(object, BIC, 1),
                      logLik = logLik,
                      check.names = FALSE)
    if (test) {
        ddf <- diff(dfModel)
        if (sum(abs(ddf)) > 0) {
            effects <- rep("", rt)
            for (i in 2:rt) {
                if (ddf[i - 1] != 0) {
                    effects[i] <- paste(i - 1, i, sep = " vs ")
                }
            }
            pval <- rep(NA, rt - 1)
            ldf <- as.logical(ddf)
            lratio <- 2 * abs(diff(logLik))
            lratio[!ldf] <- NA
            pval[ldf] <- pchisq(lratio[ldf], abs(ddf[ldf]),
                                lower.tail = FALSE)
            aod <- data.frame(aod, Test = effects,
                              L.Ratio = c(NA, lratio),
                              `p-value` = c(NA, pval),
                              check.names = FALSE)
        }
    }
    attr(aod, "rt") <- rt
    attr(aod, "verbose") <- verbose
    class(aod) <- c("anova.nlmixr", "data.frame")
    return(aod)
}

##' @export
anova.nlmixr.ui.nlme <- anova.nlmixr.ui.focei.fit

##' @export
anova.nlmixr.ui.saem <- anova.nlmixr.ui.focei.fit

print.anova.nlmixr <- function (x, verbose = attr(x, "verbose"), ...)
{
    ox <- x
    if ((rt <- attr(x, "rt")) == 1) {
        if (!is.null(lab <- attr(x, "label"))) {
            cat(lab)
            if (!is.null(L <- attr(x, "L"))) {
                print(zapsmall(L), ...)
            }
        }
        pval <- format(round(x[, "p-value"], 4))
        pval[as.double(pval) == 0] <- "<.0001"
        x[, "F-value"] <- format(zapsmall(x[, "F-value"]))
        x[, "p-value"] <- pval
        print(as.data.frame(x), ...)
    }
    else {
        ## if (verbose) {
        ##     cat("Call:\n")
        ##     objNams <- row.names(x)
        ##     for (i in 1:rt) {
        ##         cat(" ", objNams[i], ":\n", sep = "")
        ##         cat(" ", as.character(x[i, "call"]), "\n")
        ##     }
        ##     cat("\n")
        ## }
        x <- as.data.frame(x[, -1])
        for (i in names(x)) {
            xx <- x[[i]]
            if (i == "p-value") {
                xx <- round(xx, 4)
                xna <- is.na(xx)
                xx[!xna] <- format(xx[!xna])
                xx[as.double(xx) == 0] <- "<.0001"
                xx[xna] <- ""
            }
            else {
                if (match(i, c("AIC", "BIC", "logLik", "L.Ratio"),
                          0)) {
                    xna <- is.na(xx)
                    xx <- zapsmall(xx)
                    xx[xna] <- 0
                    xx <- format(xx)
                    xx[xna] <- ""
                }
            }
            x[[i]] <- xx
        }
        print(as.data.frame(x), ...)
    }
    invisible(ox)
}

focei.cleanup <- function(obj){
    if (is(obj, "focei.fit")){
        model <- obj$model
        for (m in c("ebe", "pred.only", "inner", "outer")){
            if (is(model[[m]], "RxODE")) {RxODE::rxUnload(model[[m]])}
        }
        if (file.exists(model$cache.file)){
            unlink(model$cache.file)
        }
    }
}

##  LocalWords:  focei nlme linCmt solveC SolvedC UI saem VarCorr AIC
##  LocalWords:  nlmixr REML FOCEi dplyr logLik
