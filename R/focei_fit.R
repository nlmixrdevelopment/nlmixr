regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                         group(any_of("0":"9"), ".", some_of("0":"9"))),
                      maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9"))));
regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"));
regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
regNum <- rex::rex(maybe("-"), or(regDecimalint, regFloat1, regFloat2))

##' Construct RxODE linCmt function
##'
##' @param fun function to convert to solveC syntax
##' @return SolvedC RxODE object
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
constructLinCmt <- function(fun){
    pars <- nlmixrfindLhs(body(fun));
    lines <- deparse(body(fun))[-1];
    lines <- lines[-length(lines)];
    ret <- RxODE::rxLinCmtTrans(sprintf("%s\nCentral=linCmt(%s);\n", paste(lines, collapse="\n"), paste(pars, collapse=", ")));
    ret <- strsplit(ret, "\n")[[1]];
    ret <- paste(ret[-seq_along(lines)], collapse="\n");
    return(ret)
}

is.focei <- function(x){
    env <- attr(x, ".focei.env");
    fit <- env$fit;
    if (length(names(x)) == length(fit$data.names)){
        return(all(names(x) == fit$data.names) &&
               length(x[, 1]) == fit$data.len);
    } else {
        return(FALSE)
    }
}

##' Return composite nlme/focei to nlme
##'
##' @param x nlme/focei object from common UI.
##' @return nlme object or NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
as.nlme <- function(x){
    if (is.focei(x)){
        env <- attr(x, ".focei.env");
        fit <- env$fit;
        nlme <- NULL;
        uif <- NULL
        if (is(x, "nlmixr.ui.nlme")){
            nlme <- fit$nlme;
            return(nlme);
        }
    }
    return(NULL);
}

##' Return composite saem/focei to saem
##'
##' @param x saem/focei object from common UI.
##' @return saem object or NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
as.saem <- function(x){
    if (is.focei(x)){
        env <- attr(x, ".focei.env");
        fit <- env$fit;
        saem <- NULL;
        uif <- NULL
        if (is(x, "nlmixr.ui.saem")){
            saem <- fit$saem;
            return(saem);
        }
    }
    return(NULL);
}


##' @importFrom nlme VarCorr
##' @export
VarCorr.nlmixr.ui.nlme <- function(x, sigma = 1, ...){
    VarCorr(as.nlme(x), sigma=sigma, ...)
}

use.utf <- function() {
    opt <- getOption("cli.unicode", NULL)
    if (! is.null(opt)) {
        isTRUE(opt)
    } else {
        l10n_info()$`UTF-8` && !is.latex()
    }
}

is.latex <- function() {
    if (!("knitr" %in% loadedNamespaces())) return(FALSE)
    get("is_latex_output", asNamespace("knitr"))()
}

#' Print a focei fit
#'
#' Print a first-order conditional non-linear mixed effect model with
#' interaction fit
#'
#' @param x a focei model fit
#' @param ... additional arguments
#' @return NULL
#' @export
print.focei.fit <- function(x, ...) {
    if (is.focei(x)){
        parent <- parent.frame(2);
        bound <- do.call("c", lapply(ls(parent), function(cur){
            if (identical(parent[[cur]], x)){
                return(cur)
            }
            return(NULL);
        }));
        if (length(bound) > 1) bound <- bound[1];
        if (length(bound) == 0){
            bound  <- ""
        }
        env <- attr(x, ".focei.env");
        fit <- env$fit;
        args <- as.list(match.call(expand.dots = TRUE));
        if (any(names(args) == "n")){
            n <- args$n;
        } else {
            n <- 3L;
        }
        if (any(names(args) == "width")){
            width <- args$width;
        } else {
            width <- NULL;
        }
        saem <- NULL;
        nlme <- NULL;
        uif <- NULL
        if (is(x, "nlmixr.ui.nlme")){
            nlme <- fit$nlme;
            uif <- env$uif;
            if (is(nlme, "nlme.free")){
                text <- "Free-form"
            } else {
                if (use.utf()){
                    mu <- "\u03BC";
                } else {
                    mu <- "mu"
                }
                if (is(nlme, "nlme.mu")){
                    text <- sprintf("%s-ref", mu)
                } else if (is(nlme, "nlme.mu.cov")){
                    text <- sprintf("%s-ref & covs", mu)
                }
            }
            message(cli::rule(paste0(crayon::bold$blue("nlmix"), crayon::bold$red("r"), " ", crayon::bold$yellow("nlme"), " fit by ",
                                     crayon::bold$yellow(ifelse(nlme$method == "REML", "REML", "maximum likelihood"))," (",
                                     crayon::italic(paste0(ifelse(is.null(uif$nmodel$lin.solved), "ODE", "Solved"),
                                                           "; ", text)), ")")))
        } else if (is(x, "nlmixr.ui.saem")){
            saem <- fit$saem;
            uif <- env$uif;
            message(cli::rule(paste0(crayon::bold$blue("nlmix"), crayon::bold$red("r"), " ", crayon::bold$yellow("SAEM"), " fit (",
                             crayon::italic(ifelse(is.null(uif$nmodel$lin.solved), "ODE", "Solved")), "); ",
                             crayon::blurred$italic("OBJF calculated from FOCEi approximation"))))
        } else {
            ## sprintf("nlmixr FOCEI fit (%s)\n", ifelse(fit$focei.control$grad, "with global gradient", "without global gradient")
            if (use.utf()){
                eta.name <- "\u03B7";
            } else {
                eta.name <- "eta"
            }
            type <- crayon::italic(ifelse(is.null(uif$nmodel$lin.solved), "ODE", "Solved"));
            if (is(x, "nlmixr.ui.focei.posthoc")){
                message(cli::rule(paste0(crayon::bold$blue("nlmix"), crayon::bold$red("r"), " ", crayon::bold$yellow("FOCEi"), " posthoc ",
                                         eta.name, " estimation (",type, ")")))
            } else {
                message(cli::rule(paste0(crayon::bold$blue("nlmix"), crayon::bold$red("r"), " ", crayon::bold$yellow("FOCEi"), " fit ",
                                         ifelse(fit$focei.control$grad, "with global gradient", "without global gradient"),
                                         " (",type , ")")));
            }
        }
        if (any(names(fit) == "condition.number")){
            df.objf <- data.frame(OBJF=fit$objective, AIC=AIC(x), BIC=BIC(x), "Log-likelihood"=as.numeric(logLik(x)),
                                  "Condition Number"=fit$condition.number,
                                  row.names="", check.names=FALSE)
        } else {
            df.objf <- data.frame(OBJF=fit$objective, AIC=AIC(x), BIC=BIC(x),"Log-likelihood"=as.numeric(logLik(x)),
                                  row.names="", check.names=FALSE)
        }
        if (!is.null(nlme)){
            message(paste0(crayon::bold$yellow("FOCEi"), "-based goodness of fit metrics:"))
        }
        nlmixrPrint(df.objf)
        if (!is.null(nlme)){
            message(paste0("\n",crayon::bold$yellow("nlme"), "-based goodness of fit metrics:"))
            df.objf <- data.frame(AIC=AIC(as.nlme(x)), BIC=BIC(as.nlme(x)),"Log-likelihood"=as.numeric(logLik(as.nlme(x))),
                                  row.names="", check.names=FALSE)
            nlmixrPrint(df.objf)
        }
        message(paste0("\n", cli::rule(paste0(crayon::bold("Time"), " (sec; ", crayon::yellow(bound), crayon::bold$blue("$time"), "):"))));
        print(fit$time);
        message(paste0("\n", cli::rule(paste0(crayon::bold("Parameters"), " (", crayon::yellow(bound), crayon::bold$blue("$par.fixed"), "):"))));
        print(x$par.fixed)
        tmp <- fit$omega
        diag(tmp) <- 0;
        if (all(tmp == 0)){
            message("\n  No correlations in between subject variability (BSV) matrix")
        } else {
            message("\n  Correlations in between subject variability (BSV) matrix:")
            rs <- fit$omega.R
            lt <- lower.tri(rs);
            dn1 <- dimnames(fit$omega.R)[[2]]
            nms <- apply(which(lt,arr.ind=TRUE),1,function(x){sprintf("R(%s)",paste(dn1[x],collapse=", "))});
            lt <- structure(rs[lt], .Names=nms)
            lt <- lt[lt != 0]
            digs <- 3;
            lts <- sapply(lt, function(x){
                x <- abs(lt);
                ret <- "<"
                if (x > 0.7){
                    ret <- ">" ## Strong
                } else if (x > 0.3){
                    ret <- "=" ## Moderate
                }
                return(ret)
            })
            nms <- names(lt);
            lt <- sprintf("%s%s", formatC(signif(lt, digits=digs),digits=digs,format="fg", flag="#"), lts)
            names(lt) <- nms;
            lt <- gsub(rex::rex("\""), "", paste0("    ", R.utils::captureOutput(print(lt))));
            if (crayon::has_color()){
                lt <- gsub(rex::rex(capture(regNum), ">"), "\033[1m\033[31m\\1 \033[39m\033[22m", lt, perl=TRUE)
                lt <- gsub(rex::rex(capture(regNum), "="), "\033[1m\033[32m\\1 \033[39m\033[22m", lt, perl=TRUE)
                lt <- gsub(rex::rex(capture(regNum), "<"), "\\1 ", lt, perl=TRUE)
            } else {
                lt <- gsub(rex::rex(capture(regNum), or(">", "=", "<")), "\\1 ", lt, perl=TRUE)
            }
            message(paste(lt, collapse="\n"), "\n")
        }
        message(paste0("  Full BSV covariance (", crayon::yellow(bound), crayon::bold$blue("$omega"), ") or correlation (", crayon::yellow(bound), crayon::bold$blue("$omega.R"), "; diagonals=SDs)"));
        message(paste0("  Distribution stats (mean/skewness/kurtosis/p-value) available in ",
                       crayon::yellow(bound), crayon::bold$blue("$shrink")));

        message(paste0("\n", cli::rule(paste0(crayon::bold("Fit Data"), " (object", ifelse(bound == "", "", " "),
                                              crayon::yellow(bound),
                                              " is a modified ", crayon::blue("data.frame"), "):"))))

        is.dplyr <- requireNamespace("dplyr", quietly = TRUE);
        if (!is.dplyr){
            is.data.table <- requireNamespace("data.table", quietly = TRUE);
            if (is.data.table){
                print(data.table::data.table(x), topn=n);
            } else {
                print(head(as.matrix(x), n = n));
            }
        } else {
            print(dplyr::as.tbl(x), n = n, width = width);
        }
    } else {
        print(as.data.frame(x));
    }
}

##' Extract log likelihood for focei fit
##'
##' @param object focei fit object
##' @param ... Additional arguments (currently ignored)
##' @return logLik object
##' @author Matthew L. Fidler
##' @export
logLik.focei.fit <- function(object, ...){
    ## Gives AIC as side-effect
    env <- attr(object, ".focei.env");
    fit <- env$fit;
    num <- -fit$objective/2;
    attr(num,"df") <- length(fit$par);
    class(num) <- "logLik";
    return(num);
}
##' Extract the number of observations for the focei fit
##'
##' @param object focei fit object
##' @param ... other parameters (ignored)
##' @return Number of observations
##' @author Matthew L. Fidler
##' @importFrom stats nobs
##' @export
nobs.focei.fit <- function(object,...){
    ## nobs and logLik are needed for BIC
    data <- object
    return(length(data[, 1]))
}
##' Extract variance/covariance matrix for FOCEI fit
##'
##' @param object Focei fit object
##' @param ... ignored parameters
##' @param type covariance type == blank for current. "r.s" for the
##'     sandwich matrix, "s" for the S matrix, "r" for the R matrix.
##'     When requesting another matrix, the model is updated assuming
##'     that covariance matrix is correct.
##' @return variance/covariance matrix
##' @author Matthew L. Fidler
##' @export
vcov.focei.fit <- function(object, ..., type=c("", "r.s", "s", "r")){
    type <- match.arg(type);
    if (type == ""){
        env <- attr(object, ".focei.env");
        fit <- env$fit;
        return(fit$cov);
    } else if (type == "r.s"){
        env <- attr(object, ".focei.env");
        fit <- env$fit;
        fit$cov <- fit$cov.r.s;
        fit$se <- sqrt(diag(fit$cov))
    } else if (type == "s"){
        env <- attr(object, ".focei.env");
        fit <- env$fit;
        fit$cov <- fit$cov.s;
        fit$se <- sqrt(diag(fit$cov))
    } else if (type == "r"){
        env <- attr(object, ".focei.env");
        fit <- env$fit;
        fit$cov <- fit$cov.r;
        fit$se <- sqrt(diag(fit$cov));
    }
    nms <- names(fit$theta);
    w <- seq_along(nms)
    fit$par.data.frame <- data.frame(est=fit$theta, se=fit$se[w], "%cv"=abs(fit$se[w] / fit$theta * 100),
                                     check.names=FALSE, row.names=nms);
    if (env$con$eigen){
        eig <- try(eigen(fit$cov,TRUE,TRUE)$values, silent=TRUE);
        if (!inherits(eig, "try-error")){
            fit$eigen <- eig
            tmp <- sapply(fit$eigen, abs)
            fit$condition.number <- max(tmp) / min(tmp);
        }
    }
    assign("fit", fit, env);
    return(fit$cov);
}

##' Extract the fitted values from the model
##'
##' @param object Fit object
##' @param ... other parameters
##' @param population when true (default false), calculate/extract the
##'     population predictions; When false, calculate/extract the individual
##'     predictions.
##'
##'     If this is a matrix with the number columns equal to the
##'     number of ETAs in model, and the number of rows equals to the
##'     number of subjects in the dataset, then these etas will be
##'     used to fit the data.
##' @param type The type of fitted object to be extracted.  When the
##'     value is "fitted", this gives the individual or population
##'     fitted values. When "posthoc", this extracts the posthoc
##'     deviations from the typical values or ETAs.
##' @return Individual/population predictions
##' @author Matthew L. Fidler
##' @export
fitted.focei.fit <- function(object, ..., population=FALSE,
                             type=c("fitted", "posthoc")){
    type <- match.arg(type);
    if (type == "fitted"){
        if (population){
            return(object$CPRED);
        } else {
            return(object$PRED);
        }
    }
    env <- attr(object, ".focei.env");
    fit <- env$fit;
    if (is(population, "logical")){
        tmp <- fit$etas.df
        if (!is.null(tmp)){
            return(tmp)
        }
    }
    ofv.FOCEi <- env$ofv.FOCEi;
    dat <- object
    env <- environment(ofv.FOCEi);
    if (class(population) == "matrix"){
        if (nrow(population) == env$nSUB && ncol(population) == env$nETA){
            old.mat <- env$inits.mat
            assign("inits.mat", population, env)
            on.exit({assign("inits.mat", old.mat, env)});
        } else {
            stop(sprintf("The matrix of etas specified by population needs to be %s rows by %s columns", env$nSUB, env$nETA));
        }
    }
    tmp <- tempfile()
    on.exit({sink();unlink(tmp)}, add=TRUE)
    sink(tmp);
    x <- ofv.FOCEi(fit$par)
    if (type == "posthoc"){
        d1 <- data.frame(do.call("rbind", lapply(attr(x,"subj"), function(s) {
                                              matrix(as.vector(attr(s, type)), nrow=1)
                                          })));
        names(d1) <- paste0("ETA", seq_along(d1[1, ]))
        eta.names <- fit$eta.names;
        if (!is.null(eta.names)){
            if (length(names(d1)) == length(eta.names)){
                names(d1) <- eta.names;
            }
        }
        d1 <- data.frame(ID=unique(object$ID), d1)
        return(d1)
    }
}

## nlme/lme objects
## getCovariate
## getGroups
## getResponse
## getVarCov
## logDet
## pairs
## pdConstruct?
## qqnorm
## VarCorr
## augPred
## comparePred

##' @importFrom nlme getData
##' @export
getData.focei.fit <- function(object){
    env <- attr(object, ".focei.env");
    return(env$orig.data)
}

##' @importFrom nlme ranef
##' @export
ranef.focei.fit <- function(object, ...){
    fitted.focei.fit(object, ..., type="posthoc")
}

##' @importFrom nlme fixef
##' @export
fixef.focei.fit <- function(object, ...){
    args <- as.list(match.call(expand.dots = TRUE))[-1];
    if (any(names(args) == "full")){
        if (args$full){
            uif <- object$uif;
            if (any(names(args) == "ci")){
                ci <- args$ci
            } else {
                ci <- 0.95;
            }
            qn <- stats::qnorm(1-(1-ci)/2);
            if (!is.null(uif)){
                saem <- object$saem;
                nlme <- object$nlme;
                lab <- paste(uif$ini$label[!is.na(uif$ini$ntheta)]);
                lab[lab == "NA"] <- "";
                lab <- gsub(" *$", "", gsub("^ *", "", lab));
                if (!is.null(nlme)){
                    ttab <- data.frame(summary(nlme)$tTable);
                    row.names(ttab) <- .fixNlmeNames(row.names(ttab), uif);
                    names(ttab) <- c("Estimate", "SE", "DF", "t-value", "p-value")
                    ttab <- ttab[, 1:2]
                    tmp <- object$par.data.frame;
                    tmp$SE <- NA
                    tmp <- tmp[!(row.names(tmp) %in% row.names(ttab)), ,drop = FALSE]
                    df <- data.frame(Parameter=lab, rbind(ttab,tmp))
                } else if (!is.null(saem)){
                    nth <- length(uif$saem.theta.name)
                    se <- structure(sqrt(diag(RxODE::rxInv(saem$Ha[1:nth,1:nth]))), .Names=uif$saem.theta.name)
                    df <- object$par.data.frame
                    df$SE <- NA
                    for (v in names(se)){
                        df[v, "SE"] <- se[v]
                    }
                    df <- data.frame(Parameter=lab, df);
                }
                if (!is.null(saem) || !is.null(nlme)){
                    df <- data.frame(df, "CV"=abs(df$SE / df$Estimate * 100));
                } else {
                    df <- object$par.data.frame;
                    df <- data.frame(Parameter=lab, df)
                }

                df <- data.frame(df, Untransformed=df$Estimate,
                                 Lower.ci=df$Estimate - qn * df$SE,
                                 Upper.ci=df$Estimate + qn * df$SE,
                                 check.names=FALSE);
                ## Now get the ETA information
                mu.ref <- unlist(uif$mu.ref);
                if (length(mu.ref) > 0){
                    n.mu.ref <- names(mu.ref)
                    ome <- object$omega
                    mu.ref <- structure(as.vector(n.mu.ref), .Names=as.vector(mu.ref));
                    log.eta <- uif$log.eta;
                    digs <- 3;
                    cvp <- sapply(row.names(df), function(x){
                        y <- mu.ref[x];
                        if (is.na(y)) return(" ");
                        v <- ome[y, y];
                        if (any(y == log.eta)){
                            sprintf("%s%%", formatC(signif(sqrt(exp(v) - 1) * 100, digits=digs),digits=digs,format="fg", flag="#"));
                        } else {
                            sprintf("%s", formatC(signif(sqrt(v),digits=digs), digits=digs,format="fg", flag="#"))
                        }
                    })
                    shrink <- object$shrink;
                    errs <- as.data.frame(uif$ini);
                    errs <- paste(errs[which(!is.na(errs$err)), "name"]);
                    sh <- sapply(row.names(df), function(x){
                        y <- mu.ref[x];
                        if (is.na(y)) {
                            if (any(x == errs)){
                                v <- shrink[7, "IWRES"];
                                if (is.na(v)){
                                    return(" ")
                                }
                            } else {
                                return(" ")
                            }
                        } else {
                            v <- shrink[7, y];
                        }
                        t <- ">"
                        if (v < 0){
                        } else  if (v < 20){
                            t <- "<"
                        } else if (v < 30){
                            t <- "="
                        }
                        sprintf("%s%%%s", formatC(signif(v, digits=digs),digits=digs,format="fg", flag="#"), t);
                    })
                    df <- data.frame(df, BSV.cv.sd=cvp, "Shrink(SD)%"=sh, check.names=FALSE);
                }
                if (!is(object, "nlmixr.ui.focei.posthoc")){
                    df <- data.frame(df,
                                     Lower.ci=df$Estimate - qn * df$SE,
                                     Upper.ci=df$Estimate + qn * df$SE,
                                     check.names=FALSE);

                }
                log.theta <- which(row.names(df) %in% uif$log.theta);
                if (length(log.theta) > 0){
                    df$Untransformed[log.theta] <- exp(df$Untransformed[log.theta])
                    if (!is(object, "nlmixr.ui.focei.posthoc")){
                        df$Lower.ci[log.theta] <- exp(df$Lower.ci[log.theta])
                        df$Upper.ci[log.theta] <- exp(df$Upper.ci[log.theta])
                    }
                }
                w <- which(uif$ini$err == "prop")
                if (length(w) > 0){
                    name <- uif$ini$name[w];
                    w <- which(row.names(df) %in% name);
                    if (length(w) > 0){
                        df$Untransformed[w] <- df$Untransformed[w] * 100
                        if (!is(object, "nlmixr.ui.focei.posthoc")){
                            df$Lower.ci[w] <- df$Lower.ci[w] * 100
                            df$Upper.ci[w] <- df$Upper.ci[w] * 100
                        }
                    }
                }
                attr(df, "ci") <- ci;
                class(df) <- c("nlmixr.par.fixed", "data.frame");
                return(df);
            }
            return(object$par.data.frame);
        }
    }
    return(object$theta)
}

##'@export
print.nlmixr.par.fixed <- function(x, ...){
    df <- x;
    class(df) <- "data.frame";
    digs <- 3;
    digs.cv <- 1;
    ## Rik's formating functions (modified)
    F2<-function(x,digits){
        ifelse(is.na(x), " ",
               formatC(signif(x,digits=digits), digits=digits,format="fg", flag="#"))}
    FFP<-function(x,digits){ifelse(is.na(x), " ",
                                   paste0(formatC(x,digits=digits,format="f",flag="#"), "%"))}
    df$Estimate <- F2(df$Estimate, digs);
    if (is(x, "nlmixr.ui.focei.posthoc")){
        df$SE <- F2(df$SE, digs);
        df$CV <- paste0(FFP(df$CV, digs.cv))
    }
    ci <- attr(x, "ci")
    df$Untransformed <- sprintf("%s%s %s%s%s%s%s", F2(df$Untransformed, digs),
                                ifelse(x$Estimate * 100 == x$Untransformed, "%", ""),
                                ifelse(is.na(df$Lower.ci) | is.na(df$Lower.ci), "", "("),
                                F2(df$Lower.ci, digs),
                                ifelse(is.na(df$Lower.ci) | is.na(df$Lower.ci), "", ", "),
                                F2(df$Upper.ci, digs),
                                ifelse(is.na(df$Lower.ci) | is.na(df$Lower.ci), "", ")"));
    df$CV <- gsub("NA", "", sprintf("%s", F2(df$CV, digs)));
    names(df) <- gsub("CV", "%RSE", names(df))
    names(df) <- gsub("Untransformed", sprintf("Back-transformed(%s%%CI)", ci * 100), names(df))
    if (any(names(df) == "BSV.cv.sd")){
        tmp <- df$BSV.cv.sd
        tmp <- tmp[tmp != " "];
        if (all(regexpr(rex::rex("%", end), tmp) != -1)){
            names(df) <- gsub("BSV.cv.sd", "BSV(CV%)", names(df));
        } else if (all(regexpr(rex::rex("%", end), tmp) == -1)){
            names(df) <- gsub("BSV.cv.sd", "BSV(SD)", names(df));
        } else {
            names(df) <- gsub("BSV.cv.sd", "BSV(CV% or SD)", names(df));
        }
    }
    df$SE <- gsub("NA", "", sprintf("%s", F2(df$SE, digs)));
    df <- df[, regexpr("[.]ci", names(df)) == -1]
    if (all(df$Parameter == "")){
        df <- df[, -1];
    } else {
        df$Parameter <- paste(df$Parameter);
        w <- which(df$Parameter == "");
        if (length(w) > 0){
            df$Parameter[w] <- row.names(df)[w];
        }
    }
    pdf <- R.utils::captureOutput(print(df));
    if (crayon::has_color()){
        pdf <- gsub(rex::rex(capture(regNum), "%>"), "\033[1m\033[31m\\1%\033[39m\033[22m ", pdf, perl=TRUE)
        pdf <- gsub(rex::rex(capture(regNum), "%="), "\033[1m\033[32m\\1%\033[39m\033[22m ", pdf, perl=TRUE)
        pdf <- gsub(rex::rex(capture(regNum), "%<"), "\\1% ", pdf, perl=TRUE)
        pdf <- gsub(rex::rex(capture(or(c(row.names(df), names(df))))), "\033[1m\\1\033[22m", pdf, perl=TRUE);
    } else {
        pdf <- gsub(rex::rex(capture(regNum), "%", or(">", "=", "<")), "\\1% ", pdf, perl=TRUE)
    }
    pdf <- paste(pdf, collapse="\n");
    message(pdf)
    return(invisible(df))
}

##'@export
print.nlmixr.shrink <- function(x, ...){
    tmp <- x;
    class(tmp) <- NULL;
    tmp <- (sprintf("%s%%", round(tmp, 3)));
    names(tmp) <- names(x);
    tmp <- data.frame(t(tmp))
    rownames(tmp) <- "";
    print(tmp)
}

##' Extract residuals from the FOCEI fit
##'
##' @param object focei.fit object
##' @param ... Additional arguments
##' @param type Residuals type fitted.
##' @return residuals
##' @author Matthew L. Fidler
##' @export
residuals.focei.fit <- function(object, ..., type=c("ires", "res", "iwres", "wres", "cwres", "cpred", "cres")){
    return(object[, toupper(match.arg(type))]);
}
##' Convert focei fit to a data.frame
##'
##' @param x focei fit object
##' @param row.names row names for the data frame.
##' @param optional If TRUE setting row names and column names is
##'     optional.
##' @param ... Other parameters
##' @return A data frame with no information about the focei fit
##' @author Matthew L. Fidler
##' @export
as.data.frame.focei.fit <- function(x, row.names = NULL, optional = FALSE, ...){
    tmp <- x;
    attr(tmp, ".focei.env") <- NULL;
    cls <- class(tmp);
    cls <- cls[cls != "focei.fit"];
    class(tmp) <- cls;
    as.data.frame(tmp, row.names=row.names, optional=optional, ...);
}

## FIXME: simulate function

## FIXME: show?
## FIXME: qqnorm?
## FIXME: family?


parseOM <- function(OMGA){
    re = "\\bETA\\[(\\d+)\\]\\b"
    .offset = as.integer(0)
    lapply(1:length(OMGA), function(k) {
        s = OMGA[[k]]
        f = eval(parse(text=(sprintf("y~%s", deparse(s[[2]])))))
        r = unlist(lapply(attr(terms(f),"variables"), deparse))[-(1:2)]
        nr = length(r)

        ix = grep(re, r)
        if(nr-length(ix)) stop("invalid OMGA specs")

        ix = as.integer(sub(re, "\\1", r))
        if (any(ix - (.offset+1:nr))) stop("invalid OMGA specs")
        .offset <<- .offset + nr
        eval(s[[3]])
    })
}

genOM <- function(s)
{
    getNR = function(a) round(sqrt(2 * length(a) + 0.25) - 0.1)
    nr = sum(sapply(s, getNR))
    .mat <- matrix(0, nr, nr)
    .offset = as.integer(0)
    j = lapply(1:length(s), function(k) {
        a = s[[k]]
        p <- getNR(a)
        starts = row(.mat) > .offset  & col(.mat) > .offset
        .mat[col(.mat) >= row(.mat) & col(.mat) <= .offset+p & starts] <<- a
        .offset <<- .offset+p
    })
    a = .mat[col(.mat) >= row(.mat)]
    .mat <- t(.mat)
    .mat[col(.mat) >= row(.mat)] <- a
    .mat
}
##' Parameter history extraction
##'
##' @param x object to extract parameter history from
##' @param stacked Should the data frame be stacked (default FALSE)
##' @param ... Other parameters (currently ignored)
##' @return A data frame with the parameter history.  NULL if no
##'     information is available
##' @author Matthew L. Fidler & Wenping Wang
##' @export
par.hist <- function(x, stacked=FALSE, ...){
    UseMethod("par.hist")
}
##' @rdname par.hist
##' @export
par.hist.nlmixr.ui.saem <- function(x, stacked=FALSE, ...){
    uif <- x$uif;
    theta.names <- c(uif$saem.theta.name, uif$saem.omega.name, uif$saem.res.name);
    m <- x$saem$par_hist
    if (stacked){
        df <- data.frame(
            val=as.vector(m),
            par=rep(theta.names, each=nrow(m)),
            iter=rep(1:nrow(m), ncol(m)))
        return(df)
    } else {
        dimnames(m) <- list(NULL, theta.names);
        df <- data.frame(iter=rep(1:nrow(m)), as.data.frame(m));
        return(df)
    }
}
##' @rdname par.hist
##' @export
par.hist.nlmixr.ui.focei.fit <- function(x, stacked=FALSE, ...){
    m <- x$fit.df
    if (stacked){
        m <- data.frame(stack(m[,-1]), iter=m$iter);
        names(m) <- c("val", "par", "iter")
    }
    return(m)
}

par.hist.default <- function(x, stacked=FALSE, ...){
    return(NULL)
}
##' Produce trace-plot for fit if applicable
##'
##' @param x fit object
##' @param ... other parameters
##' @return Fit traceplot or nothing.
##' @author Rik Schoemaker, Wenping Wang & Matthew L. Fidler
##' @export
traceplot <- function(x, ...){
    UseMethod("traceplot");
}

##' @rdname traceplot
##' @export
traceplot.focei.fit <- function(x, ...){
    m <- par.hist(x, stacked=TRUE);
    if (!is.null(m)){
        p0 <- ggplot(m, aes(iter, val)) +
            geom_line() +
            facet_wrap(~par, scales = "free_y")
        if (!is.null(x$mcmc)){
            p0 <- p0 + ggplot2::geom_vline(xintercept=x$mcmc$niter[1], col="blue", size=1.2);
        }
        print(p0)
    }
}

#' Plot a focei.fit plot
#'
#' Plot some standard goodness of fit plots for the focei fitted object
#'
#' @param x a focei fit object
#' @param ... additional arguments
#' @return NULL
#' @author Wenping Wang & Matthew Fidler
#' @export
plot.focei.fit <- function(x, ...) {
    traceplot(x);
    dat <- as.data.frame(x);
    d1 <- data.frame(DV=dat$DV, stack(dat[, c("PRED", "IPRED")]))

    p1 <- ggplot(d1,aes(values,DV)) +
        facet_wrap(~ind) +
        geom_abline(slope=1, intercept=0, col="red", size=1.2) +
        geom_smooth(col="blue", lty=2, formula=DV ~ values + 0, size=1.2) +
        geom_point() +
        xlab("Predictions");

    print(p1)

    p2 <- ggplot(dat, aes(x=IPRED, y=IRES)) +
        geom_point() +
        geom_abline(slope=0, intercept=0, col="red")
    print(p2)

    p2 <- ggplot(dat, aes(x=TIME, y=IRES)) +
        geom_point() +
        geom_abline(slope=0, intercept=0, col="red")
    print(p2)


    ids <- unique(dat$ID)
    for (i  in seq(1, length(ids), by=16)){
        tmp <- ids[seq(i, i + 15)]
        tmp <- tmp[!is.na(tmp)];
        d1 <- dat[dat$ID %in% tmp, ];

        p3 <- ggplot(d1, aes(x=TIME, y=DV)) +
            geom_point() +
            geom_line(aes(x=TIME, y=IPRED), col="red", size=1.2) +
            geom_line(aes(x=TIME, y=PRED), col="blue", size=1.2) +
            facet_wrap(~ID)
        print(p3)
    }
}

                                        #
##TODOs
##- NM opt
##- AGQ

##' FOCEI Fit for nlmixr
##'
##' @param data Data to fit
##' @param inits Initialization list
##' @param PKpars Pk Parameters
##' @param diag.xform The form of the diagonal of the Omega matrix to
##'     be estimated.  This can be "sqrt" "log" or "identity"
##' @param model The RxODE model to use
##' @param pred The Prediction function
##' @param err The Error function
##' @param lower Lower bounds
##' @param upper Upper Bounds
##' @param control Control list
##' @param theta.names Names of the thetas to be used in the final object
##' @param eta.names Eta names to be used in the final object
##' @param ... Ignored parameters
##' @return A focei fit object
##' @author Matthew L. Fidler and Wenping Wang
##' @export
focei.fit <- function(data,
                      inits,
                      PKpars,
                      diag.xform=c("sqrt", "log", "identity"),
                      model=NULL,
                      pred=NULL,
                      err=NULL,
                      lower= -Inf,
                      upper= Inf,
                      control=list(),
                      theta.names=NULL,
                      eta.names=NULL,
                      ...){
    UseMethod("focei.fit");
}
##' @export
##' @rdname focei.fit
focei.fit.nlmixr.ui.nlme <- function(data, inits, ...){
    name.data <- TRUE;
    if (missing(inits)){
        inits <- getData(data);
        name.data <- FALSE
    } else if (!is(inits, "data.frame")){
        stop("The second argument needs to be missing, data, or a prior fit.");
    }
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    names(call)[1] <- "object"
    if (name.data){
        names(call)[2] <- "data";
    } else {
        call$data <- inits;
    }
    call$est <- "focei";
    return(do.call(getFromNamespace("nlmixr","nlmixr"), call, envir = parent.frame(1)));
}
##' @export
##' @rdname focei.fit
focei.fit.nlmixr.ui.focei.fit <- focei.fit.nlmixr.ui.nlme
##' @export
##' @rdname focei.fit
focei.fit.nlmixrUI <- function(data, inits, ...){
    name.data <- TRUE;
    if (!is(inits, "data.frame")){
        stop("The second argument needs to a data frame to fit from the UI function.");
    }
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    names(call)[1] <- "object"
    if (name.data){
        names(call)[2] <- "data";
    } else {
        call$data <- inits;
    }
    call$est <- "focei";
    return(do.call(getFromNamespace("nlmixr","nlmixr"), call, envir = parent.frame(1)));
}
##' @export
##' @rdname focei.fit
focei.fit.function <- focei.fit.nlmixrUI

##' @export
##' @rdname focei.fit
focei.fit.data.frame <- function(...){
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    return(.collectWarnings(do.call(focei.fit.data.frame0, call, envir=parent.frame(1))))
}

focei.fit.data.frame0 <- function(data,
                                  inits,
                                  PKpars,
                                  diag.xform=c("sqrt", "log", "identity"),
                                  model=NULL,
                                  pred=NULL,
                                  err=NULL,
                                  lower= -Inf,
                                  upper= Inf,
                                  control=list(),
                                  theta.names=NULL,
                                  eta.names=NULL){

    orig.data <- data;
    do.table <- FALSE;
    sink.file <- tempfile();
    orig.sink.number <- sink.number(type="output");
    do.sink <- TRUE;
    fit.df <- NULL;
    sink.close <- function(n=orig.sink.number){
        if (do.sink){
            while(sink.number(type="output") > n){
                sink(type="output");
            }
        }
        ## message("Closed sink");
        ## unlink(sink.file);
    }
    sink.start <- function(do.it=TRUE){
        ## sink.close();
        if (do.sink){
            if (do.it){
                ## message("Starting Sink to", sink.file);
                sink(sink.file, type="output");
            }
        }
    }
    sink.get <- function(){
        lines <- NULL;
        if (do.sink){
            if (file.exists(sink.file)){
                if (file.size(sink.file) != 0){
                    lines <- readLines(sink.file);
                    return(lines);
                }
            }
        }
        return(lines);
    }
    on.exit({sink.close();
        lines <- sink.get();
        unlink(sink.file);
        if (!is.null(lines)){
            message("After Error:");
            message(paste(paste("##", lines), collapse="\n"));
        }
        running <- FALSE}, add=TRUE)
    ##data = dat; PKpars=mypars; diag.xform="sqrt"; model=list(); control=list()
    ##model options

    ##options
    con <- list(
        DEBUG.ODE = F,
        DEBUG = F, RESET.INITS.MAT = T,
        TRACE.INNER=FALSE,
        TOL.INNER=1e-4,
        trace = 0,
        atol.ode=1e-6,
        rtol.ode=1e-6,
        hmin = 0L,
        hmax = NULL,
        hini = 0L,
        maxordn = 12L,
        maxords = 5L,
        stiff=1L,
        atol.outer=1e-6,
        rtol.outer=1e-6,
        maxsteps.ode = 99999,
        reltol.outer = 1e-4,
        absltol.outer = 1e-4,
        cores=1,
        transitAbs=FALSE,
        NONMEM=TRUE,
        NOTRUN=FALSE,
        PRINT.PARS=FALSE,
        ## NONMEM does PRED-DV
        ## Almquist does DV-PRED
        pred.minus.dv=TRUE,
        switch.solver=FALSE,
        cov.method="r,s",
        factr=1e10,
        grad=FALSE,
        accept.eta.size=1.5,
        sigdig=0,
        precision=0, ## 0 = No ridge penalty.
        ridge.decay=0, ## 0 = no decay; Inf = No ridge
        reset.precision=NULL,
        eigen=TRUE,
        rhobeg=.2,
        rhoend=1e-2,
        npt=NULL,
        ## FIXME bounds need to be changed.. They aren't working
        est.chol.omegaInv=FALSE, ##RxODE::rxSymPyVersion() >= 1,
        add.posthoc=TRUE,
        extra.output=TRUE, ## Display extra output on each iteration
        inits.mat=NULL,
        find.best.eta=TRUE,
        inner.opt="n1qn1",
        inner.opt="lbfgs",
        save.curve=TRUE,
        numDeriv.method1="simple",
        numDeriv.method2="simple",
        numDeriv.swap=2.3,
        sum.prod=TRUE,
        theta.grad=FALSE,
        scale.to=1,
        optim="L-BFGS-B"
    )

    curi <- 0;

    this.env <- environment()

    pt <- proc.time()
    nmsC <- names(con)
    if (length(control) > 0){
        control <- control[sapply(names(control), function(x){!is.null(control[[x]])})]
    }
    con[(namc <- names(control))] <- control
    if (!is.null(control$transit_abs)){
        con$transitAbs <- control$transit_abs;
    }
    if (length(noNms <- namc[!namc %in% c(nmsC, "transit_abs")]))
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

    running <- TRUE
    if (con$NOTRUN){
        running <- FALSE
    }

    numDeriv.method <- con$numDeriv.method1

    cl <- NULL;
    if (.Platform$OS.type == "windows" && con$cores > 1){
        cl <- parallel::makePSOCKcluster(con$cores)
        on.exit({parallel::stopCluster(cl)}, add=TRUE);
    }

    if (con$sigdig == 0 & con$ridge.decay == 0 && !con$extra.output){
        do.sink <- FALSE
    }
    optim <- con$optim;
    optim.method <- match.arg(optim, c( ## "n1qn1",
                                         "bobyqa",
                                         "L-BFGS-B",
                                         "BFGS",
                                         ## "lbfgs",
                                         "lbfgsb3",
                                         ## "newuoa",
                                         "nlminb"##,
                                         ## "uobyqa"
                                     ))
    grad.methods <- c("BFGS", "L-BFGS-B", "lbfgs", "lbfgsb3", "nlminb", "mma", "slsqp", "lbfgs-nlopt", "tnewton_precond_restart",
                      "tnewton_precond", "tnewton", "var1", "var2", "n1qn1")
    print.grad <- any(optim.method == grad.methods);
    if (con$grad && !print.grad){
        if (!con$NOTRUN){
            message("Warning; You selected a gradient method, but the optimization procedure doesn't require the gradient.\nIgnoring gradient")
        }
        con$grad <- FALSE;
        con$theta.grad <- FALSE
    }
    if (con$NOTRUN){
        con$theta.grad <- FALSE;
        print.grad <- FALSE;
    }
    if(is(model, "RxODE") || is(model, "character")) {
        ODEmodel <- TRUE
        if (class(pred) != "function"){
            stop("pred must be a function specifying the prediction variables in this model.")
        }
    }
    else {
        ODEmodel <- TRUE
        model <- constructLinCmt(PKpars);
        pred <- eval(parse(text="function(){return(Central);}"))
    }

    square <- function(x) x*x
    diag.xform <- match.arg(diag.xform)
    diag.xform.inv = c("sqrt"="square", "log"="exp", "identity"="identity")[diag.xform]

    ##process inits
    if (is.null(err)){
        err <-eval(parse(text=paste0("function(){err",paste(inits$ERROR[[1]],collapse=""),"}")));
    }
    ## print(th0.om)
    model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err, grad=con$grad,
                                     pred.minus.dv=con$pred.minus.dv, sum.prod=con$sum.prod,
                                     theta.derivs=con$theta.grad, run.internal=TRUE);
    ## rxCat(model$inner);
    if (!con$NOTRUN){
        message(sprintf("Original Compartments=%s", length(RxODE::rxState(model$obj))))
        message(sprintf("\t Inner Compartments=%s", length(RxODE::rxState(model$inner))))
        if (con$grad){
            message(sprintf("\t Outer Compartments=%s", length(RxODE::rxState(model$outer))))
        }
    }
    cov.names <- par.names <- RxODE::rxParams(model$pred.only);
    cov.names <- cov.names[regexpr(rex::rex(start, or("THETA", "ETA"), "[", numbers, "]", end), cov.names) == -1];
    colnames(data) <- sapply(names(data), function(x){
        if (any(x == cov.names)){
            return(x)
        } else {
            return(toupper(x))
        }
    })
    pcov <- c()
    lhs <- c(names(RxODE::rxInits(model$pred.only)), RxODE::rxLhs(model$pred.only))
    if (length(lhs) > 0){
        cov.names <- cov.names[regexpr(rex::rex(start, or(lhs), end), cov.names) == -1];
    }
    if (length(cov.names) > 0){
        if (!all(cov.names %in% names(data))){
            message("Model:")
            RxODE::rxCat(model$pred.only)
            message("Needed Covariates:")
            nlmixrPrint(cov.names)
            stop("Not all the covariates are in the dataset.")
        }
        message("Needed Covariates:")
        nlmixrPrint(cov.names)
    }

    ## RxODE(rxNorm(model$inner), modName="test");
    if (is.null(model$extra.pars)){
        nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)))
    } else {
        nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)),
                 sprintf("ERR[%s]", seq_along(model$extra.pars)))
    }
    if (!is.null(theta.names) && (length(inits$THTA) + length(model$extra.pars)) == length(theta.names)){
        nms <- theta.names;
    }
    if (length(lower) == 1){
        lower <- rep(lower, length(inits$THTA));
    } else if (length(lower) != length(inits$THTA)){
        print(inits$THTA)
        print(lower)
        stop("Lower must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
    }
    if (length(upper) == 1){
        upper <- rep(upper, length(inits$THTA));
    } else if (length(lower) != length(inits$THTA)){
        stop("Upper must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
    }

    extra.pars <- c();
    if (!is.null(model$extra.pars)){
        model$extra.pars <- eval(call(diag.xform, model$extra.pars))
        if (length(model$extra.pars) > 0){
            inits$THTA <- c(inits$THTA, model$extra.pars);
            lower.err <- rep(con$atol.ode * 10, length(model$extra.pars));
            upper.err <- rep(Inf, length(model$extra.pars));
            lower <-c(lower, lower.err);
            upper <- c(upper, upper.err);
        }
    }

    ##FIXME
    ##data
    if (is.null(data$ID)) stop('"ID" not found in data')
    if (is.null(data$DV)) stop('"DV" not found in data')
    if (is.null(data$EVID)) data$EVID = 0
    if (is.null(data$AMT)) data$AMT = 0
    ## Make sure they are all double amounts.
    for (v in c("TIME", "AMT", "DV", cov.names))
        data[[v]] <- as.double(data[[v]]);
    data.sav = data
    ds <- data[data$EVID > 0, c("ID", "TIME", "AMT", cov.names)]
    data <- data[data$EVID == 0, c("ID", "TIME", "DV", cov.names)]
    ## keep the covariate names the same as in the model
    w <- which(!(names(data.sav) %in% cov.names))
    names(data.sav)[w] <- tolower(names(data.sav[w]))         #needed in ev

    lh = parseOM(inits$OMGA)
    nlh = sapply(lh, length)
    osplt = rep(1:length(lh), nlh)
    lini = list(inits$THTA, unlist(lh));
    nlini = sapply(lini, length)
    nsplt = rep(1:length(lini), nlini)

    om0 = genOM(lh)
    rxSym <- RxODE::rxSymInvCreate(mat=om0, diag.xform=diag.xform, chol=con$est.chol.omegaInv);
    th0.om <- rxSym$th;

    w.th <- seq_along(inits$THTA);
    w.om <- seq_along(rxSym$th) + length(inits$THTA);

    ## Add lower bounds for omega matrix.
    w <- RxODE::rxSymDiag(rxSym);
    th0.om <- unlist(th0.om);
    lower.om <- rep(-Inf, length(th0.om));
    upper.om <- rep(Inf, length(th0.om));
    lower.om[w] <- con$atol.ode * 10;
    upper.om[w] <- th0.om[w] * 1e4;
    lower <- c(lower, lower.om);
    upper <- c(upper, upper.om);
    inits.vec = c(inits$THTA, th0.om)

    if (any(inits.vec == 0)){
        warning("Some of the initial conditions were 0, changing to 0.0001");
        inits.vec[inits.vec == 0] <- 0.0001;
    }
    if (!con$NOTRUN){
        message("Boundaries:");
        nlmixrPrint(data.frame(lower,inits.vec,upper));
    }
    names(inits.vec) = NULL
    if (con$scale.to == 0){
        con$scale.to <- NULL;
    }
    if (is.null(con$scale.to)){
        par.lower <- lower;
        par.upper <- upper
    } else {
        par.lower <- lower / inits.vec * con$scale.to;
        par.upper <- upper / inits.vec * con$scale.to;
    }

    if (!con$NOTRUN && !is.null(con$scale.to)){
        message("Scaled Boundaries:");
    }
    w.neg <- which(par.lower > par.upper);
    tmp <- par.lower[w.neg];
    par.lower[w.neg] <- par.upper[w.neg];
    par.upper[w.neg] <- tmp;
    if (!con$NOTRUN){
        if (!is.null(con$scale.to)){
            nlmixrPrint(data.frame(par.lower,scaled=rep(con$scale.to, length(inits.vec)),par.upper))
        }
        if (do.sink){
            message("\nKey:")
            if (!is.null(con$scale.to)){
                message(" S: Scaled Parameter values");
            }
            message(" U: Unscaled Parameter values");
            message(" X: eXp(Unscaled Parameter value)");
            if(print.grad){
                message(" G: Gradient");
            }
            message(" D: Significant Figures");
            message(" Optimization output displayed with comments, i.e. ##\n");
        } else {
            message("");
        }
    }
    nTHTA = nlini[1]
    nETA  = nrow(om0)
    ID.all = unique(data[,"ID"])        #FIXME
    ID.ord = order(ID.all)
    names(ID.ord) = ID.all
    nSUB  = length(ID.all)

    find.best.eta <- con$find.best.eta;
    first <- TRUE

    ofv.FOCEi.ind <- memoise::memoise(function(pars) {
                                      cur.diff <<- NULL;
                                      if(con$PRINT.PARS) print(pars)
                                      ## initial parameters are scale.to which translates to scale.to * inits.vec / scale.to = inits.vecs
                                      if (!is.null(con$scale.to)){
                                          pars <- pars * inits.vec / con$scale.to
                                      }
                                      if (!is.null(last.pars) && con$accept.eta.size != 0 && find.best.eta && con$grad){
                                          inits.mat.bak <- inits.mat;
                                          inits.mat <- rxUpdateEtas(last.dEta.dTheta, matrix(pars - last.pars, ncol=1), inits.mat,con$accept.eta.size);
                                      } else {
                                          inits.mat.bak <- NULL;
                                      }
                                      ## last.pars
                                      THETA <- pars[w.th];
                                      rxSymEnv <-  RxODE::rxSymInv(rxSym, pars[w.om]);
                                      llik.one <- function(subj, con){
                                          ev <- RxODE::eventTable()
                                          dati <- data.sav[data.sav$id==subj, ]
                                          if (length(cov.names) > 0){
                                              cur.cov <- dati[, cov.names, drop = FALSE];
                                          } else {
                                              cur.cov <- NULL;
                                          }
                                          ev$import.EventTable(dati)
                                          .wh = ID.ord[as.character(subj)]
                                          if(con$DEBUG>10 && .wh==1) print(inits.mat[.wh,, drop = FALSE])               #FIXME
                                          if (con$DEBUG.ODE) print("i'm here :)")
                                          c.hess <- NULL;
                                          if (con$save.curve && !first && con$inner.opt == "n1qn1") c.hess <- inits.c.hess[.wh, ];
                                          args <- list(model, ev, theta=THETA, eta=inits.mat[.wh,, drop = FALSE], c.hess=c.hess,
                                                       dv=dati$dv[ev$get.obs.rec()], inv.env=rxSymEnv,
                                                       nonmem=con$NONMEM, invisible=1-con$TRACE.INNER, epsilon=con$TOL.INNER,
                                                       id=subj, inits.vec=inits.vec, cov=cur.cov, estimate=find.best.eta,
                                                       atol=con$atol.ode, rtol=con$rtol.ode, maxsteps=con$maxsteps.ode,
                                                       atol.outer=con$atol.outer, rtol.outer=con$rtol.outer,
                                                       hmin = con$hmin, hmax = con$hmax, hini = con$hini, transitAbs = con$transitAbs,
                                                       maxordn = con$maxordn, maxords = con$maxords, stiff=con$stiff,
                                                       pred.minus.dv=con$pred.minus.dv, switch.solver=con$switch.solver,
                                                       inner.opt=con$inner.opt, add.grad=print.grad, numDeriv.method=numDeriv.method,
                                                       table=do.table);
                                          if (!is.null(con$scale.to)){
                                              args$scale.to <- con$scale.to
                                          }
                                          ##method=numDeriv.method
                                          if (!is.null(inits.mat.bak)){
                                              args$eta.bak=inits.mat.bak[.wh, ];
                                          }
                                          ret <- do.call(getFromNamespace("rxFoceiInner","nlmixr"), args)
                                          attr(ret, "wh") <- .wh; ## FIXME move to C?
                                          return(ret)
                                      }
                                      ##------------------------------
                                      ## parallelize
                                      if (con$cores > 1){
                                          if (.Platform$OS.type != "windows") {
                                              llik.subj <- parallel::mclapply(X = ID.all, FUN = llik.one, mc.cores = con$cores, con=con)
                                          } else {
                                              llik.subj <- parallel::parLapply(cl, X = ID.all, fun = llik.one, con=con)
                                          }
                                      } else {
                                          llik.subj <- lapply(ID.all, llik.one, con=con);
                                      }
                                      ## Use wenping's solution for inits.mat
                                      m = t(sapply(llik.subj, function(x) {
                                          c(attr(x, "wh"), attr(x, "posthoc"))
                                      }))
                                      if(con$RESET.INITS.MAT) inits.mat[m[,1],] <<- m[,-1]
                                      ## Save last curvature
                                      if (con$save.curve && con$inner.opt == "n1qn1" && !is.null(attr(llik.subj[[1]], "c.hess"))){
                                          m <- t(sapply(llik.subj, function(x){
                                              c(attr(x, "wh"), attr(x, "c.hess"))
                                          }))
                                          inits.c.hess[m[,1],] <<- m[,-1]
                                      }
                                      if (con$grad && con$accept.eta.size != 0){
                                          last.pars <<- pars;
                                          last.dEta.dTheta <<- lapply(llik.subj, function(x){attr(x, "dEta.dTheta")});
                                      }
                                      llik.subj
                                  })

    ofv.FOCEi.vec = function(pars) {
        llik.subj = ofv.FOCEi.ind(pars)
        unlist(llik.subj)
    }

    last.ofv <- NULL;
    last.step.ofv <- NULL;
    cur.diff <- NULL;

    last.pars <- NULL;
    last.step.pars <- NULL;
    last.dEta.dTheta <- NULL;

    last.penalty <- NULL;
    sigdig.fit <- NULL;
    last.sigdig <- NULL;

    if (con$sigdig != 0 && !any(optim.method == c("BFGS", "L-BFGS-B", "lbfgsb3"))){
        warning("Significant figures doesn't make sense with this optimization routine.  Resetting to no sigdig exit")
        con$sigdig <- 0;
    }

    regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                             group(any_of("0":"9"), ".", some_of("0":"9"))),
                          maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9"))));
    regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"));
    regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))));
    regNum <- rex::rex(maybe("-"), or(regDecimalint, regFloat1, regFloat2));

    digs <- 5;
    optim.obj <- function(lines, prefix="o"){
        if (class(lines) == "numeric"){
            ofv <- lines;
            ## if (any(optim.method == c("lbfgs"))){
            ##     return(paste0(prefix, ".", RxODE::rxCout(ofv)));
            ## } else
            if (any(optim.method == c("L-BFGS-B", "BFGS"))){
                return(sprintf("%s.%.6f", prefix, ofv))
            } else if (any(optim.method == c("nlminb", "newuoa", "bobyqa", "uobyqa", "n1qn1"))){
                return(paste0(prefix, ".", gsub("^ *", "", sprintf("%#14.8g",ofv))));
            } else if (optim.method == "lbfgsb3"){
                ## This uses dblepr("",)
                ##  call dblepr(" f=",-1, f, 1)
                ## According to r documentation this is the same as printing a number.
                ## I believe this is the same as format(#)
                return(sprintf("%s.%s", prefix, format(last.ofv)));
            } else {
                ## The R interface uses Rprintf("%f")
                return(sprintf("%s.%s", prefix, sprintf("%f", last.ofv)));
            }
        } else {
            if (any(optim.method == "lbfgs")){
                reg <- rex::rex("fx = ",capture(regNum));
            } else if (any(optim.method == c("L-BFGS-B", "BFGS"))){
                reg <- rex::rex("iter", any_spaces, any_numbers, any_spaces, "value", any_spaces, capture(regNum));
            } else if (any(optim.method == c("newuoa", "nlminb", "bobyqa", "uobyqa", "n1qn1"))){
                reg <- rex::rex(any_spaces, any_numbers, ":", any_spaces, capture(regNum), any_spaces, ":",anything);
            } else if  (optim.method == "lbfgsb3") {
                reg <- rex::rex(anything, "f", any_spaces, "=", any_spaces,
                                capture(regNum), any_spaces, end);
            } else {
                reg <- rex::rex(anything, "f(x)", any_spaces, "=", any_spaces, capture(regNum), any_spaces, end);
            }
            w <- which(regexpr(reg, lines, perl=TRUE) != -1)
            if (length(w) > 0){
                w <- w[length(w)];
                return(gsub(reg, "o.\\1", lines[w], perl=TRUE));
            } else {
                return("");
            }
        }
    }
    curi <- 0;
    ofv.cache <- new.env(parent=emptyenv())
    ofv.cache$last1 <- NULL;
    ofv.cache$last <- NULL;
    ofv.cache$first <- NULL;
    ofv.cache$echo.ridge <- TRUE;

    subset.lines <- function(lines){
        ## Warning: H2 seems singular; Using pseudo-inverse
        w <- which(regexpr(rex::rex("Warning: H2 seems singular; Using pseudo-inverse"), lines) != -1);
        if (length(w) > 0){
            bad.inv.H2 <- TRUE;
            lines <- lines[-w];
        } else {
            bad.inv.H2 <- FALSE;
        }
        ## error: inv(): matrix seems singular
        w <- which(regexpr(rex::rex("error: inv(): matrix seems singular"), lines) != -1);
        if (length(w) > 0){
            bad.inv <- TRUE;
            lines <- lines[-w];
        } else {
            bad.inv <- FALSE;
        }
        w <- which(regexpr(rex::rex("error: chol(): decomposition failed"), lines) != -1);
        ## error: chol(): decomposition failed
        if (length(w) > 0){
            bad.chol <- TRUE;
            lines <- lines[-w];
        } else {
            bad.chol <- FALSE;
        }
        cor.nearPD <- FALSE;
        ## Warning: The Hessian is non-positive definite, correcting with nearPD
        w <- which(regexpr(rex::rex("Warning: The Hessian is non-positive definite, correcting with nearPD"), lines) != -1);
        if (length(w) > 0){
            cor.nearPD <- TRUE;
            lines <- lines[-w];
        }
        w <- which(regexpr(rex::rex("Warning: The Hessian chol() decomposition failed, but we can use the (slower) determinant instead"), lines) != -1);
        cor.det <- FALSE
        if (length(w) > 0){
            cor.det <- TRUE;
            lines <- lines[-w];
        }

        w <- which(regexpr(rex::rex("Warning: Hessian (H) seems singular; Using pseudo-inverse"), lines) != -1);
        if (length(w) > 0){
            lines <- lines[-w];
            message("## Hessian(s)  inverted with pinv()");
        }
        lines <- lines[regexpr(rex::rex(start, any_spaces, end),
                               lines) == -1];

        lines <- lines[regexpr(rex::rex(start, any_spaces, end),
                               lines) == -1];
        if (cor.det || cor.nearPD){
            message("## Hessian(s) corrected by ", appendLF=FALSE);
            if (cor.det && cor.nearPD){
                message("both det and nearPD");
            } else if (cor.nearPD){
                message("nearPD");
            } else {
                message("det");
            }
        }
        if (bad.inv.H2){
            message("## Hessian(s) from second order derivatives inverted with pinv()");
        }
        message(paste(paste0("## ", lines), collapse="\n"));
        return(lines);
    }

    print.step <- function(lst=NULL, last.ofv.txt=NULL){
        reset <- FALSE
        last <- last.ofv.txt
        if (is.null(lst)){
            last.info <- get(last, envir=ofv.cache, inherits=FALSE);
        } else{
            last.info <- lst
        }
        message(sprintf("Step %s: %s", curi, substr(last, 3, nchar(last))), appendLF=FALSE);
        if (is.null(ofv.cache$last)){
            message("")
        } else if (ofv.cache$last$ofv > last.info$ofv){
            message(sprintf(" (%s reduction)", format(ofv.cache$last$ofv - last.info$ofv)))
            reset <- TRUE;
        } else {
            message(sprintf(" (%s increase)", format(last.info$ofv - ofv.cache$last$ofv)))
        }
        if (is.null(con$scale.to)){
            p <- rbind(data.frame(t(c(last.info$pars))),
                       data.frame(exp(t(c(last.info$pars)))));
            rn <- c("U:", "X:")
        } else {
            p <- rbind(data.frame(t(c(last.info$pars))),
                       data.frame(t(c(last.info$pars * inits.vec / con$scale.to))),
                       data.frame(exp(t(c(last.info$pars * inits.vec / con$scale.to)))))
            rn <- c("S:", "U:", "X:")
        }
        ## message(paste0("\n S: ", paste(sapply(last.info$pars, function(x){sprintf("%#8g", x)}), collapse=" ")))
        ## message(paste0(" U: ", paste(sapply(last.info$pars*inits.vec, function(x){sprintf("%#8g", x)}), collapse=" ")))
        ## message(paste0(" X: ", paste(sapply(last.info$pars*inits.vec, function(x){sprintf("%#8g", exp(x))}), collapse=" ")))
        grad.txt <- optim.obj(last.info$llik, "l");
        if (exists(grad.txt, envir=ofv.cache, inherits=FALSE)){
            last.grad <- get(grad.txt, envir=ofv.cache, inherits=FALSE);
            p <- rbind(p,
                       data.frame(t(c(last.grad))))
            rn <- c(rn, "G:")
        }
        sdig <- rep(-Inf, length(last.info$pars))
        if (is.null(ofv.cache$last)){
            reset <- TRUE;
            rownames(p) <- rn;
            if ((length(theta.names) + length(eta.names)) == length(names(p))){
                names(p) <- c(theta.names, paste0("ome(", eta.names, ")"));
            } else {
                names(p)[seq_along(nms)] <- nms
            }
            print(p)
        } else {
            last.pars <- ofv.cache$last$pars
            cur.pars <- last.info$pars
            if (all(last.pars == cur.pars) && !is.null(ofv.cache$last1)){
                message("## Parameters the same, use last iteration for sigdig calculation...");
                last.pars <- ofv.cache$last1$pars;
                ofv.cache$last <- ofv.cache$last1;
            }
            ##
            ## NONMEM's sigdigs as per http://cognigencorp.com/nonmem/current/2012-October/3654.html
            ##
            if (!is.null(con$scale.to)){
                sdig <- -log10(abs(cur.pars - last.pars) * 0.1 / con$scale.to);
                rn <- c(rn, "D:");
                p <- rbind(p, data.frame(t(sdig)));
            } else {
                sdig <- -log10(abs(cur.pars - last.pars) / abs(inits.vec));
                rn <- c(rn, "D:");
                p <- rbind(p, data.frame(t(sdig)));
            }
            rownames(p) <- rn;
            if ((length(theta.names) + length(eta.names)) == length(names(p))){
                names(p) <- c(theta.names, paste0("ome(", eta.names, ")"));
            } else {
                names(p)[seq_along(nms)] <- nms
            }
            print(p)
        }
        return(list(p, sdig, reset, last.info))
    }
    ofv.FOCEi <- function(pars) {
        llik.subj <- ofv.FOCEi.ind(pars)
        first <<- FALSE
        llik <- -2*PreciseSums::psSum(unlist(llik.subj));
        corrected <- do.call("sum", (lapply(llik.subj, function(x){attr(x, "corrected")})))
        ofv <- llik;
        reset <- FALSE;
        sigdig.exit <- FALSE
        cur.pars <- pars;
        if (running){
            ## Print and calculate sigdigs.  If sigdigs exit, then exit.
            ## Can possibly have earlier exit criterion.
            sink.close();
            lines <- sink.get();
            if (!is.null(lines)){
                lines <- subset.lines(lines);
                last <- optim.obj(lines);
                if (last != ""){
                    if (exists(last, envir=ofv.cache, inherits=FALSE)){
                        p <- print.step(last.ofv.txt=last)
                        sdig <- p[[2]];
                        reset <- p[[3]];
                        last.info <- p[[4]];
                        p <- p[[1]]
                        if (all(sdig > con$numDeriv.swap)){
                            numDeriv.method <- con$numDeriv.method2
                        }
                        if (con$sigdig > 0 && all(sdig > con$sigdig) && curi > 3){
                            sigdig.exit <- TRUE
                        }
                        if (is.null(ofv.cache$first)){
                            ofv.cache$first <- llik;
                        }
                        cur.diff <<- (ofv.cache$first - llik);
                        w <- which(row.names(p) == "U:");
                        curi <<- curi + 1;
                        tmp <- cbind(data.frame(iter=curi, objf=as.numeric(substr(last, 3, nchar(last))), p[w,, drop = TRUE]));
                        row.names(tmp) <- NULL;
                        fit.df <<- rbind(fit.df, tmp, check.names=FALSE);
                    } else {
                        message("Warning: Prior objective function not found...");
                        reset <- TRUE;
                    }
                    message("");
                }
            }
        }
        if (is.null(cur.diff)){
            extra <- 1
        } else if (con$ridge.decay != 0){
            extra <- exp(-con$ridge.decay * cur.diff);
        } else {
            extra <- 1;
        }
        ## ridge penalty llik - 0.5*(sum(pars)^2*con$precision):
        last.penalty <<- 0.5 * sum( pars * pars * con$precision) * extra;
        if (reset && con$ridge.decay != 0 && ofv.cache$echo.ridge){
            message(sprintf("Current ridge penalty: %s", last.penalty));
            if (sprintf("%s", last.penalty) == "0"){
                assign("echo.ridge", FALSE, envir=ofv.cache, inherits=FALSE)
            }
        }
        lst <- list(pars=pars, llik=llik);
        llik <- llik + last.penalty;
        attr(llik, "subj") = llik.subj
        last.ofv <<- as.numeric(llik);
        lst$ofv <- last.ofv;
        last.pars <<- pars;
        assign(optim.obj(last.ofv), lst, envir=ofv.cache, inherits=FALSE);
        if (reset){
            if (!is.null(ofv.cache$last)){
                ofv.cache$last1 <- ofv.cache$last
            }
            last.info <- get(last, envir=ofv.cache, inherits=FALSE);
            ofv.cache$last <- last.info;
        }
        if (running && any(optim == c("L-BFGS-B", "BFGS", "lbfgs")) &&
            ((is.null(con$scale.to) && sum(pars - inits.vec) == 0) ||
             (!is.null(con$scale.to) && all(pars == con$scale.to)))){
            last.ofv.txt <- optim.obj(last.ofv)
            p <- print.step(lst, last.ofv.txt)[[1]];
            if (is.null(ofv.cache$first)){
                ofv.cache$first <- llik;
                ofv.cache$last <- lst;
            }
            cur.diff <<- (ofv.cache$first - llik);
            w <- which(row.names(p) == "U:");
            curi <<- curi + 1;
            tmp <- cbind(data.frame(iter=curi, objf=as.numeric(substr(last.ofv.txt, 3, nchar(last.ofv.txt))), p[w,, drop = TRUE]));
            row.names(tmp) <- NULL;
            fit.df <<- rbind(fit.df, tmp, check.names=FALSE);
        }
        sink.start(running);
        if (sigdig.exit){
            message("Exit based on sigdig convergence");
            if (llik < last.info$ofv){
                cur.pars <- pars;
                ofv <- as.vector(llik)
            } else {
                ofv <- last.info$ofv;
            }
            fit <- list()
            fit$par = cur.pars;
            fit$objective = ofv;
            fit$convergence = 2
            fit$message = "sigdig convergence"
            sigdig.fit <<- fit;
            stop("sigidig exit");
        }
        if (running){
            return(as.vector(llik))
        } else {
            ## sink.close(0);
            return(llik)
        }
    }

    gr.FOCEi <- function(pars, envir=this.env){
        ret <- .Call(`_nlmixr_grFOCEi`, pars, envir);
        return(ret)
    }

    s.FOCEi <- function(pars, envir=this.env){
        .Call(`_nlmixr_sFOCEi`, pars, envir)
    }

    np <- length(inits.vec)
    meth <- c();

    if (is.null(con$scale.to)){
        par0 <- inits.vec
    } else {
        par0 <- rep(con$scale.to, length(inits.vec));
    }

    opt <- function(){
        fit <- NULL;
        if (optim.method == "lbfgs"){
            fit <- try({lbfgs::lbfgs(call_eval=ofv.FOCEi, call_grad=gr.FOCEi, vars=par0)});
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!(length(lines) == 1L && lines[1L] == ""))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            } else if (!is.null(fit) && any(names(fit) == "value")) {
                fit$objective <- fit$value;
                fit$value <- NULL;
            }
        } else if (any(optim.method == c("L-BFGS-B", "BFGS"))){
            fit <- try({stats::optim(par=par0, fn=ofv.FOCEi, gr=gr.FOCEi,
                                     method=optim.method,
                                     control=list(trace=1, REPORT=1),
                                     lower=ifelse(optim.method == "L-BFGS-B", par.lower, -Inf),
                                     upper=ifelse(optim.method == "L-BFGS-B", par.upper, Inf))});
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!(length(lines) == 1L && lines[1L] == ""))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            } else if (!is.null(fit) && any(names(fit) == "value")) {
                fit$objective <- fit$value;
                fit$value <- NULL;
            }
        } else if (optim.method=="lbfgsb3"){
            prm <- par0;
            if (requireNamespace("lbfgsb3", quietly = TRUE)){
                fit <- try({lbfgsb3::lbfgsb3(prm=prm,
                                             fn=ofv.FOCEi, gr=gr.FOCEi,
                                             control=list(trace=1## , pgtol = con$reltol.outer
                                                          ),
                                             lower=par.lower,
                                             upper=par.upper
                                             )})
            } else {
                stop("This requires lbfgsb3 to be installed;  Please install it.")
            }
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!(length(lines) == 1L && lines[1L] == ""))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            } else if (!is.null(fit) && any(names(fit) == "prm")) {
                ## pgtol should be able to be specified, Currently only trace and iprint are used.
                fit$par = fit$prm
                fit$prm = NULL
                fit$objective <- fit$f;
                fit$f <- NULL;
            }
        } else if (optim.method == "bobyqa"){
            if (!is.null(con$npt)) npt=con$npt
            else npt = 2*np+1
            fit <- minqa::bobyqa(par0,
                                 ofv.FOCEi, ## gradient=gr.FOCEi,
                                 control=list(rhobeg=con$rhobeg, rhoend=con$rhoend, npt=npt, iprint=3),
                                 lower=par.lower,
                                 upper=par.upper
                                 );
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!is.null(lines))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            }
        } else if (optim.method == "newuoa"){
            if (!is.null(con$npt)) npt=con$npt
            else npt = 2*np+1
            fit <- minqa::newuoa(par0,
                                 ofv.FOCEi, ## gradient=gr.FOCEi,
                                 control=list(rhobeg=con$rhobeg, rhoend=con$rhoend, npt=npt, iprint=3)
                                 ## lower=par.lower,
                                 ## upper=par.upper
                                 );
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!is.null(lines))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            }
        } else if (optim.method == "uobyqa"){
            fit <- minqa::uobyqa(par0,
                                 ofv.FOCEi, ## gradient=gr.FOCEi,
                                 control=list(rhobeg=con$rhobeg, rhoend=con$rhoend, npt=npt, iprint=3)
                                 ## lower=par.lower,
                                 ## upper=par.upper
                                 );
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!is.null(lines))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            }
        } else if (optim.method=="nlminb") {
            fit <- try({nlminb(par0,
                               ofv.FOCEi, gradient=gr.FOCEi,
                               control=list(trace=1, rel.tol=con$reltol.outer),
                               lower=par.lower,
                               upper=par.upper)});
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!is.null(lines))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            }
        }
        return(fit);
    }
    if (is.null(con$inits.mat)){
        inits.mat <- matrix(0, nSUB, nETA)
    } else if (sum(dim(con$inits.mat) - c(nSUB, nETA)) == 0){
        inits.mat <- con$inits.mat;
    } else {
        stop(sprintf("Mismatch of dimensions; ETA matrix supplied: %sx%s, expected %sx%s", dim(con$inits.mat)[1], dim(con$inits.mat)[2],
                     nSUB, nETA))
    }
    if (con$inner.opt == "n1qn1"){
        inits.c.hess <- matrix(0, nSUB, nETA * (nETA + 13) / 2)
    }
    if (.Platform$OS.type == "windows" && con$cores > 1){
        ## Copy over all of the objects within scope to
        ## all clusters.
        ## Taken from https://github.com/nathanvan/mcmc-in-irt/blob/master/post-10-mclapply-hack.R

        ## while( identical( this.env, globalenv() ) == FALSE ) {
        parallel::clusterExport(cl,
                                ls(all.names=TRUE, envir=this.env),
                                envir=this.env)
        ##     this.env <- parent.env(environment())
        ## }
        parallel::clusterExport(cl,
                                ls(all.names=TRUE, envir=globalenv()),
                                envir=globalenv())
    }
    setup.time <- proc.time() - pt;
    pt <- proc.time()
    if (con$NOTRUN) {
        sink(sink.file, type="output");
        pt <- proc.time()
        fit = list()
        fit$par = rep(con$scale.to, length(inits.vec))
        fit$value = as.vector(ofv.FOCEi(rep(con$scale.to, length(inits.vec))))
        fit$convergence = -99
        fit$message = "notrun"
        optim.time <- proc.time() - pt;
        fit$cov.time <- proc.time() - proc.time();
        fit$objective <- fit$value;
        sink(type="output");
        lines <- readLines(sink.file);
        unlink(sink.file);
        w <- which(regexpr(rex::rex("Warning: The Hessian is non-positive definite, correcting with nearPD"), lines) != -1);
        if (length(w) > 0){
            warning("Hessian correcting during objective function calculation.")
        }
    }
    else {
        if (con$grad){
            sink.start()
            fit <- opt();
            if (!is.null(fit)){
                last.pars <- NULL;
                last.dEta.dTheta <- NULL;
                if ((!is.null(con$reset.precision) && con$reset.precision != con$precision) ||
                    (con$precision > 0 && con$ridge.decay > 0 && !is.infinite(con$ridge.decay) && sprintf("%s", last.penalty) != "0")){
                    ofv.cache <- new.env(parent=emptyenv())
                    ofv.cache$last <- NULL;
                    ofv.cache$first <- NULL;
                    ofv.cache$echo.ridge <- TRUE;

                    last.ofv <- NULL;
                    last.step.ofv <- NULL;
                    cur.diff <- NULL;

                    last.pars <- NULL;
                    last.step.pars <- NULL;
                    last.dEta.dTheta <- NULL;

                    last.penalty <- NULL;
                    sigdig.fit <- NULL;

                    if (is.null(con$reset.precision)){
                        con$precision <- 0
                    } else {
                        con$precision <- con$reset.precision;
                    }
                    if (is.null(scale.to)){
                        inits.vec = fit$par;
                    } else {
                        inits.vec = fit$par*inits.vec / con$scale.to;
                    }
                    sink.start();
                    fit <- opt();
                }
            }
        } else {
            fit = opt();
            if (!is.null(fit)){
                last.pars <- NULL;
                last.dEta.dTheta <- NULL;

                if ((!is.null(con$reset.precision) && con$reset.precision != con$precision) ||
                    (con$precision > 0 && con$ridge.decay > 0 && !is.infinite(con$ridge.decay) && sprintf("%s", last.penalty) != "0")){
                    last.ofv <- NULL;
                    last.step.ofv <- NULL;
                    cur.diff <- NULL;

                    last.pars <- NULL;
                    last.step.pars <- NULL;
                    last.dEta.dTheta <- NULL;

                    last.penalty <- NULL;
                    sigdig.fit <- NULL;
                    if (is.null(con$reset.precision)){
                        con$precision <- 0
                    } else {
                        con$precision <- con$reset.precision;
                    }
                    if (is.null(scale.to)){
                        inits.vec = fit$par;
                    } else {
                        inits.vec = fit$par*inits.vec / con$scale.to;
                    }
                    sink.start();
                    fit = opt();
                }
            }
        }
        if (inherits(fit, "try-error"))
            fit <- NULL
        if (!is.null(fit)){
            running <- FALSE;
            fit$precision = con$precision;
            con$precision <- 0;
            fit$grad = con$grad;
            fit$ridge.decay <- con$ridge.decay;
            fit$sigdig <- con$sigdig;
            con$sigdig <- 0 ;
            fit$ridge.decay <- con$ridge.decay
            con$ridge.decay <- 0;
            fit$reset.precision <- con$reset.precision;
            con$reset.precision <- 0;
            optim.time <- proc.time() - pt;
            pt <- proc.time();
            fit$setup.time <- setup.time
            fit$optim.time <- optim.time
            ## if (con$fix.eta.for.grad){
            ##     find.best.eta <- FALSE; ## Keep etas.
            ## }
            ## Cov
            ## - invert second derivative of the llik
            ## - Lewis method
            ## - Lineaize model
            ## Easy to complete
            ## Can I approximate the model
            ## Stochastic estimation
            ## Bootstrap
            ## LLik are asymptotic.  Not likely asymptotic.
            ## Psim assumes the fixed effects and random effects
            ## - Assumed the random effects is zero.
            ## - Not sure what is accepted as an input
            if (!print.grad && !con$grad){
                print.grad <- TRUE; ## Calculate covariance a bit qicker
                memoise::forget(ofv.FOCEi.ind);
            }
            message("Calculate covariance...")
            sink.start();
            scale.to <- con$scale.to;
            if (is.null(scale.to)){
                fit$par.unscaled = fit$par;
            } else {
                fit$par.unscaled = fit$par*inits.vec / con$scale.to;
            }
            ## Perform on unscaled parameters
            con$scale.to <- NULL;
            ## Use First Order condition for covariance
            ## inits.mat <- matrix(0, nSUB, nETA)
            if (con$cov.method != "s"){
                R1 <- optimHess(fit$par.unscaled, ofv.FOCEi, gr.FOCEi);
                Rinv <- RxODE::rxInv(R1)
                if (det(R1) <= 0){
                    warning("Non positive-definite Hessian matrix when calculating the covariance; Falling back to S-matrix covariance")
                    fit$hessian.bad <- R1;
                    ## R1 <- as.matrix(Matrix::nearPD(R1)$mat);
                    ## Rinv <- RxODE::rxInv(R1);
                    con$cov.method <- "s"
                    fit$warning <- "Non positive-definite Hessian matrix when calculating the covariance; Falling back to S-matrix covariance"
                    Rinv <- NULL
                }
            } else {
                Rinv <- NULL
            }
            ## if (con$grad){
            s = s.FOCEi(fit$par.unscaled);
            con$scale.to <- scale.to;
            rm(scale.to)
            if (con$cov.method != "s"){
                fit$cov.r.s <- Rinv %*% s %*% Rinv
                fit$cov.r <- 2 * Rinv;
                fit$cov.s <- 4 * RxODE::rxInv(s);
                if (con$cov.method == "r"){
                    fit$cov <- fit$cov.r;
                } else {
                    fit$cov <- fit$cov.r.s;
                }
            } else {
                fit$cov.s <- 4 * RxODE::rxInv(s);
                fit$cov <- fit$cov.s;
            }
            if (any(diag(fit$cov) <= 0)){
                ## Can be semi-definite (i.e. some eigenvalues can be zero, but assume not.)
                w <- "Diagonal elements of covariance were negative; Finding nearest positive definate matrix to estimate covariance.";
                if (any(names(fit) == "warning")){
                    fit$warning <- c(fit$warning,w)
                } else {
                    fit$warning <- w;
                }
                warning(w);
                fit$cov.bad <- fit$cov
                fit$cov <- as.matrix(Matrix::nearPD(fit$cov)$mat);
            }
            lines <- sink.get();
            if (!is.null(lines)){
                cor.nearPD <- FALSE;
                ## Warning: The Hessian is non-positive definite, correcting with nearPD
                w <- which(regexpr(rex::rex("Warning: The Hessian is non-positive definite, correcting with nearPD"), lines) != -1);
                if (length(w) > 0){
                    cor.nearPD <- TRUE;
                }
                w <- which(regexpr(rex::rex("Warning: The Hessian chol() decomposition failed, but we can use the (slower) determinant instead"), lines) != -1);
                cor.det <- FALSE
                if (length(w) > 0)
                    cor.det <- TRUE
                if (cor.det || cor.nearPD){
                    message("## Hessian(s) corrected by ", appendLF=FALSE);
                    if (cor.det && cor.nearPD){
                        message("both det and nearPD");
                    } else if (cor.nearPD){
                        message("nearPD");
                    } else {
                        message("det");
                    }
                }
            }
            message("done\n");
            if (con$eigen){
                fit$eigen <- eigen(fit$cov,TRUE,TRUE)$values;
                tmp <- sapply(fit$eigen, abs)
                fit$condition.number <- max(tmp) / min(tmp);
            }
            fit$se <- sqrt(diag(fit$cov))
            fit$last.penalty <- last.penalty;
            fit$cov.time <- proc.time() - pt;
            fit$sigdig <- last.sigdig;
        }
    }
    if (!is.null(fit)){
        ## class(con) <- "focei.fit.con"
        fit$focei.control <- con;
        if (is.null(fit$par.unscaled)){
            if (is.null(con$scale.to)){
                fit$par.unscaled = fit$par;
            } else {
                fit$par.unscaled = fit$par*inits.vec / con$scale.to;
            }
        }
        tmp <- fit$par.unscaled[seq_along(nms)];
        names(tmp) <- nms;
        fit$theta <- tmp;
        for (n in c("fval", "value")){
            if (any(names(fit) == n)){
                fit$objective <- fit[[n]];
            }
        }
        lD <- fit$par.unscaled[-seq_along(nms)];
        rxSymEnv <-  RxODE::rxSymInv(rxSym, lD);
        fit$omega <- rxSymEnv$omega;
        d <- sqrt(diag(fit$omega))
        d1 <- d;
        d <- RxODE::rxInv(diag(length(d)) * d);
        omega.r <- d %*% fit$omega %*% d;
        diag(omega.r) <- d1;

        if (length(eta.names) == dim(fit$omega)[1]){
            dimnames(fit$omega) <- list(eta.names, eta.names);
            fit$eta.names <- eta.names;
            dimnames(omega.r) <- dimnames(fit$omega)
        }
        fit$omega.R <- omega.r
        fit$fit.df <- fit.df;
        w <- seq_along(nms)
        if (con$NOTRUN){
            fit$par.data.frame <- data.frame(Estimate=fit$theta, row.names=nms);
        } else {
            fit$par.data.frame <- data.frame(Estimate=fit$theta, SE=fit$se[w], CV=abs(fit$se[w] / fit$theta * 100),
                                             check.names=FALSE, row.names=nms);
        }
        env <- environment(ofv.FOCEi);
        attr(data, ".focei.env") <- env;
        class(data) <- c("focei.fit", "data.frame")
    } else {
        stop("Fitting data unsuccessful.");
    }
    memoise::forget(ofv.FOCEi.ind);
    ## Allow IPRED/PRED to be calculated by taking out the caching
    ## (which is per-parameter, not per parameter and eta)
    ## ofv.FOCEi.ind <- ofv.FOCEi.ind.slow;
    running <- FALSE;
    sink.close()
    con$cores <- 1; ## don't run parallel for extracting information
    find.best.eta <- FALSE;
    print.grad <- FALSE; ## Turn off slow gradient calculation
    do.table <- TRUE
    message("Calculating Table Variables...")
    tmp <- c()
    pt <- proc.time();
    res <- calc.resid.fit(data, orig.data, con);
    fit$shrink <- res[[2]];
    fit$con <- con;
    data <- cbind(as.data.frame(data), res[[1]], res[[3]], res[[4]]);
    table.time <- proc.time() - pt;
    fit$table.time <- table.time;
    fit$data.names <- names(data);
    fit$data.len <- length(data[, 1]);
    fit$time <- data.frame(setup=setup.time["elapsed"],
                           optimize=optim.time["elapsed"],
                           covariance=fit$cov.time["elapsed"],
                           table=fit$table.time["elapsed"],
                           row.names="");
    fit$model <- model;
    message("done")
    attr(data, ".focei.env") <- env
    class(data) <- c("focei.fit", "data.frame")
    data
}


##' @export
`$.focei.fit` <-  function(obj, arg, exact = TRUE){
    m <- as.data.frame(obj);
    ret <- m[[arg, exact = exact]];
    if (is.null(ret)){
        env <- attr(obj, ".focei.env");
        if (arg == "uif"){
            return(env$uif)
        } else if (arg == "par.hist.stacked"){
            return(par.hist(obj, stacked=TRUE))
        } else if (arg == "par.hist"){
            return(par.hist(obj))
        } else if (arg == "par.fixed"){
            return(fixed.effects(obj, full=TRUE));
        } else if (arg == "eta"){
            return(random.effects(obj));
        } else if (arg == "seed"){
            if (is(obj, "nlmixr.ui.saem")){
                return(attr(as.saem(obj),"saem.cfg")$seed)
            } else {
                return(NA);
            }
        } else if (arg == "model.name"){
            return(env$uif$model.name);
        } else if (arg == "data.name"){
            return(env$uif$data.name);
        } else if (arg == "simInfo"){
            return(.simInfo(obj));
        } else {
            fit <- env$fit;
            ret <- fit[[arg, exact = exact]]
            if (is.null(ret)){
                if (exists(arg, env, inherits=FALSE)){
                    return(get(arg, env, inherits=FALSE));
                } ## else if (exists(arg, env$uif$env, inherits=FALSE)){
                ##     return(get(arg, env$uif$env, inherits=FALSE));
                ## }
                else {
                    return(NULL)
                }
            } else {
                return(ret)
            }
        }
    } else {
        return(ret)
    }
}

##' @importFrom utils str
##' @export
str.focei.fit <- function(object, ...){
    uif <- object$uif;
    if (is.null(uif)){
        cat('FOCEi combined dataset and properties\n');
    } else {
        cat('nlmixr UI combined dataset and properties\n')
    }
    m <- as.data.frame(object);
    str(m)
    env <- attr(object, ".focei.env");
    fit <- env$fit;
    str(fit)
    str(as.list(env$fit))
    str(as.list(object$uif$env))
    cat(" $ par.hist         : Parameter history (if available)\n")
    cat(" $ par.hist.stacked : Parameter history in stacked form for easy plotting (if available)\n")
    cat(" $ par.fixed        : Fixed Effect Parameter Table\n")
    cat(" $ eta              : Individual Parameter Estimates\n")
    cat(" $ seed             : Seed (if applicable)\n");
    cat(" $ model.name       : Model name (from R function)\n");
    cat(" $ data.name        : Name of R data input\n");
    cat(" $ simInfo          : RxODE list for simulation\n");
}

##' @export
head.focei.fit <- function(x, ...){
    ## hack to be consistent with method... and satisfy CMD CHECK
    args <- as.list(match.call(expand.dots=TRUE))
    if (any(names(args) == "n")){
        n <- args$n
    } else {
        n <- 6
    }
    message('FOCEI combined dataset and list');
    env <- attr(x, ".focei.env");
    message('Head List:')
    fit <- env$fit;
    print(head(fit, n=floor(n / 2), ...))
    message('Head Data (Also accessible by $fit):')
    m <- as.data.frame(x);
    return(head(m, n=floor(n / 2), ...))
}

##' @export
tail.focei.fit <- function(x, ...){
    ## hack to be consistent with method... and satisfy CMD CHECK
    args <- as.list(match.call(expand.dots=TRUE))
    if (any(names(args) == "n")){
        n <- args$n
    } else {
        n <- 6
    }
    message('FOCEI combined dataset and list');
    env <- attr(x, ".focei.env");
    message('Tail List:')
    fit <- env$fit;
    print(tail(fit, n=floor(n / 2), ...))
    message('Tail Data (Also accessible by $fit):')
    m <- as.data.frame(x);
    return(tail(m, n=floor(n / 2), ...))
}


##' Convert fit to FOCEi style fit
##'
##' @param object Fit object to convert to FOCEi-style fit.
##' @param uif Unified Interface Function
##' @param pt Proc time object
##' @param ... Other Parameters
##' @param data The data to pass to the FOCEi translation.
##' @return A FOCEi fit style object.
##' @author Matthew L. Fidler
as.focei <- function(object, uif, pt=proc.time(), ..., data){
    UseMethod("as.focei");
}

##' Get the FOCEi theta or eta specification for model.
##'
##' @param object Fit object
##' @param uif User interface function or object
##' @param ... Other parameters
##' @return List for the OMGA list in FOCEi
##' @author Matthew L. Fidler
focei.eta <- function(object, uif, ...){
    UseMethod("focei.eta");
}

##' Get the FOCEi theta specification for the model
##'
##' @inheritParams focei.eta
##' @return Parameter estimates for Theta
focei.theta <- function(object, uif, ...){
    UseMethod("focei.theta");
}

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
