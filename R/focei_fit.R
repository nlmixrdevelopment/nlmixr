##' Construct RxODE linCmt function
##'
##' @param fun function to convert to solveC syntax
##' @return SolvedC RxODE object
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
constructLinCmt <- function(fun){
    pars <- nlmixrfindLhs(body(fun));
    ret <- RxODE::rxLinCmtTrans(sprintf("Central=linCmt(%s);\n", paste(pars, collapse=", ")));
    return(ret)
}

is.focei <- function(x){
    env <- attr(x, ".focei.env");
    fit <- env$fit;
    if (length(names(x)) == length(fit$data.names)){
        return(all(names(x) == fit$data.names) &&
               length(x[, 1]) == fit$data.len)
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
        env <- attr(x, ".focei.env");
        fit <- env$fit;
        args <- as.list(match.call(expand.dots = TRUE));
        if (any(names(args) == "n")){
            n <- args$n;
        } else {
            n <- 6L;
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
            message(sprintf("nlmixr nlme fit by %s (%s)\n", ifelse(nlme$method == "REML", "REML", "maximum likelihood"),
                            ifelse(is.null(uif$nmodel$lin.solved), "ODE", "Solved")));
        } else if (is(x, "nlmixr.ui.saem")){
            saem <- fit$saem;
            uif <- env$uif;
            message(sprintf("nlmixr SAEM fit (%s)\n", ifelse(is.null(uif$nmodel$lin.solved), "ODE", "Solved")))
        } else {
            message(sprintf("nlmixr FOCEI fit (%s)\n", ifelse(fit$control$grad, "with global gradient", "without global gradient")));
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
            message("FOCEi-based goodness of fit metrics:")
        }
        RxODE::rxPrint(df.objf)
        if (!is.null(nlme)){
            message("\nnlme-based goodness of fit metrics:")
            df.objf <- data.frame(AIC=AIC(as.nlme(x)), BIC=BIC(as.nlme(x)),"Log-likelihood"=as.numeric(logLik(as.nlme(x))),
                                  row.names="", check.names=FALSE)
            RxODE::rxPrint(df.objf)
        }
        message("\nTime (sec):");
        print(fit$time);
        message("\nParameters:")
        if (!is.null(nlme)){
            ttab <- summary(nlme)$tTable;
            print(ttab)
            tmp <- fit$par.data.frame
            message("\nResidual Errors")
            print(tmp[!(row.names(tmp) %in% row.names(ttab)), ,drop = FALSE]);
        } else if (!is.null(saem)){

            df <- fit$par.data.frame;
            th = saem$Plambda
            nth = length(th)
            H = solve(saem$Ha[1:nth,1:nth])
            se = sqrt(diag(H))
            m =  cbind(exp(th), th, se)
            dimnames(m)[[2]] = c("est", "log(est)", "se(log_est)")
            trans <- uif$saem.theta.trans;
            nm <- rep("", nth)
            dfr <- rownames(df)
            for (i in seq_along(trans)){
                ## i = focei trans[i] = saem
                nm[trans[i]] <- dfr[i];
            }
            dimnames(m)[[1]] <- nm
            message("From SAEM:")
            print(m)

            ## Augment with SEs from SAEM.
            ## message("SAEM!");
            message("FOCEi:")
            print(df);
        } else {
            print(fit$par.data.frame);
        }
        message("\nOmega:");
        print(fit$omega);
        is.dplyr <- requireNamespace("dplyr", quietly = TRUE);
        if (!is.dplyr){
            message("\nFit Data (head):")
            print(head(as.matrix(x), n = n));
        } else {
            message("\nFit Data:")
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
    fit$par.data.frame <- data.frame(est=fit$theta, se=fit$se[w], "%cv"=fit$se[w] / fit$theta * 100,
                                     check.names=FALSE, row.names=nms);
    if (env$con$eigen){
        fit$eigen <- eigen(fit$cov,TRUE,TRUE)$values;
        tmp <- sapply(fit$eigen, abs)
        print(tmp);
        fit$condition.number <- max(tmp) / min(tmp);
    }
    assign("fit", fit, env);
    return(fit$cov);
}

##' Extract the fitted values from the model
##'
##' @param object Fit object
##' @param ... other parameters
##' @param population when true (default false), calculate the
##'     population predictions; When false, calculate the individual
##'     predictions.
##'
##'     If this is a matrix with the number columns equal to the
##'     number of ETAs in model, and the number of rows equals to the
##'     number of subjects in the dataset, then these etas will be
##'     used to fit the data.
##' @param type The type of fitted object to be extracted.  When the
##'     value is "fitted", this gives the individual or population
##'     fitted values. The "Vi" option gives the variance estimate for
##'     the individual.  The "Vfo" gives the Variance under the fo
##'     assumption when population=FALSE, and the FOCE assumption when
##'     population=TRUE. "dErr_dEta" gives the df/deta*eta. When
##'     "posthoc", this extracts the posthoc deviations from the
##'     typical values or ETAs.
##' @return Individual/population predictions
##' @author Matthew L. Fidler
##' @export
fitted.focei.fit <- function(object, ..., population=FALSE,
                             type=c("fitted", "Vi", "Vfo", "dErr_dEta", "dR_dEta", "posthoc")){
    env <- attr(object, ".focei.env");
    fit <- env$fit;
    type <- match.arg(type);
    ofv.FOCEi <- env$ofv.FOCEi;
    dat <- object
    env <- environment(ofv.FOCEi);
    if (class(population) == "logical"){
        if (population && type != "posthoc"){
            old.mat <- env$inits.mat
            assign("inits.mat", matrix(0, env$nSUB, env$nETA), env)
            on.exit({assign("inits.mat", old.mat, env)});
        }
    } else if (class(population) == "matrix"){
        if (nrow(population) == env$nSUB && ncol(population) == env$nETA){
            old.mat <- env$inits.mat
            assign("inits.mat", population, env)
            on.exit({assign("inits.mat", old.mat, env)});
        } else {
            stop(sprintf("The matrix of etas specified by population needs to be %s rows by %s columns", env$nSUB, env$nETA));
        }
    }
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
    } else {
        return(unlist(lapply(attr(x,"subj"), function(s) attr(s, type))))
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
    return(object$theta)
}

##' Extract residuals from the FOCEI fit
##'
##' @param object focei.fit object
##' @param ... Additional arguments
##' @param type Residuals type fitted.
##' @param etas ETAs matrix to use for the calculation.
##' @return residuals
##' @author Matthew L. Fidler
##' @export
residuals.focei.fit <- function(object, ..., type=c("ires", "res", "iwres", "wres", "cwres", "cpred", "cres"),
                                etas=NULL){
    env <- attr(object, ".focei.env");
    if (!is.null(etas)){
        if (class(etas) == "matrix" && nrow(etas) == env$nSUB && ncol(etas) == env$nETA){
            old.mat <- env$inits.mat
            assign("inits.mat", etas, env)
            on.exit({assign("inits.mat", old.mat, env)});
        } else {
            stop(sprintf("The matrix of etas specified by etas needs to be %s rows by %s columns", env$nSUB, env$nETA))
        }
    }
    dat <- object;
    DV <- dat$DV;
    type <- match.arg(type);
    ## up.type <- toupper(type)
    ## if (any(up.type == names(object))){
    ##     return(dat[, ]);
    ## }
    if (any(type == c("res", "wres", "cwres", "cres", "cpred", "cpredi"))){
        ## population=TRUE => eta=0
        if (type == "wres"){
            ## Efo= f(|eta=0)
            ## cov = Vfo|eta=0+Vi|eta=0
            ## These WRES are calculated the same as Hooker 2007, but don't seem to match NONMEM.
            pred <- fitted(object, population=TRUE)
            ## In theory the decorrelation is done by eigenvector decomposition. However, it doens't seem to affect the WRES values.
            ## W <- fitted(object, population=TRUE, type="Vfo") + fitted(object, population=TRUE, type="Vi")
            ## tmp <- aggregate(W,by=list(object$ID),FUN=function(x){e <- eigen(diag(x)); return(diag(e$vectors %*% diag(sqrt(e$values)) %*% rxInv(e$vectors)))});
            ## tmp <- data.frame(ID=tmp[,1],stack(data.frame(tmp[,-1])));
            ## tmp <- tmp[order(tmp$ID,tmp$ind),];
            ## tmp <- sqrt(tmp$values)
            ## Even though it does not match NONMEM, it should be adequate metric for WRES...
            W <- fitted(object, population=TRUE, type="Vfo") + fitted(object, population=TRUE, type="Vi")
            return((DV - pred) / sqrt(W));
        } else if (type == "cwres"){
            ## According to Hooker 2007, the Vi (or the
            ## dh/deta*Sigma*dh/deta) should be calculated under eta=0
            ## conditions.  However, I don't think this is correct;
            ## Neither does NONMEM; They use Vi under the post-hoc
            ## predicted conditions
            pred <- fitted(object, population=FALSE) - fitted(object, population=FALSE, type="dErr_dEta");
            W <- sqrt(fitted(object, population=FALSE, type="Vfo") + fitted(object, population=FALSE, type="Vi"))
            return((DV - pred) / W);
        } else if (type == "cpred"){
            pred <- fitted(object, population=FALSE) - fitted(object, population=FALSE, type="dErr_dEta");
            return(pred);
        } else if (type == "cres"){
            pred <- fitted(object, population=FALSE) - fitted(object, population=FALSE, type="dErr_dEta");
            return(DV - pred);
        } else {
            pred <- fitted(object, population=TRUE);
            return(DV - pred);
        }
    } else if (any(type == c("ires", "iwres"))){
        ipred <- fitted(object, population=FALSE);
        if (type == "iwres"){
            W <- sqrt(fitted(object, population=FALSE, type="Vi"));
            return((DV - ipred) / W);
        } else {
            return(DV - ipred);
        }
    }
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

##' Simulate response based on FOCEi model's datathe fitted object's dataset
##'
##' @param object Focei object
##' @param nsim Number of simulated parameters
##' @param seed Seed to start with (if specified)
##' @param ... Other parameters
##' @return New dataset based on original dataset
##' @export
##' @author Matthew L. Fidler
simulate.focei.fit <- function(object, nsim=1, seed=NULL, ...){
    if(!is.null(seed)){
        set.seed(seed);
    }
    nsub <- length(unique(object$ID));
    nobs <- length(object$ID);
    do.call("rbind",
            lapply(seq(1, nsim),
                   function(x){
                eta.mat <- mvtnorm::rmvnorm(nsub, sigma = object$omega, ...)
                ipred <- fitted(object, population=eta.mat);
                dv <- ipred + sapply(fitted(object, population=eta.mat, type="Vi"),
                                     function(r){ return(rnorm(1, sd=sqrt(r)))});
                return(data.frame(SIMID=x, ID=object$ID, TIME=object$TIME, DV=dv, IPRED=ipred))
            }))
}
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
    m = x$saem$par_hist
    df = data.frame(
        val=as.vector(m),
        par=rep(1:ncol(m), each=nrow(m)),
        iter=rep(1:nrow(m), ncol(m))
    )
    if (!stacked){
        df <-  data.frame(iter=df$iter, unstack(df,val~par));
    }
    return(df)
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
    m = par.hist(x, stacked=TRUE);
    if (!is.null(m)){
        p0 = ggplot(m, aes(iter, val)) +
            geom_line() +
            facet_wrap(~par, scales = "free_y")

        print(p0)
    }

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
##'     be estimated.  This can be "sqrt" "log" or "identity".
##' @param optim The optimization method to be used
##' @param model The RxODE model to use
##' @param pred The Prediction function
##' @param err The Error function
##' @param lower Lower bounds
##' @param upper Upper Bounds
##' @param control Control list
##' @param calculate.vars is a list of variables that will be
##'     calculated after the FOCEI estimation is complete.
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
                      optim=c(
                          ## "n1qn1",
                          "lbfgs",
                          "L-BFGS-B",
                          "bobyqa",
                          "lbfgsb3",
                          ## "newuoa",
                          "nlminb"##,
                          ## "uobyqa"
                      ),
                      model=NULL,
                      pred=NULL,
                      err=NULL,
                      lower= -Inf,
                      upper= Inf,
                      control=list(),
                      calculate.vars=c("pred", "ipred", "ires", "res", "iwres", "wres", "cwres", "cpred", "cres"),
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
        stop("The second argument needs to be missing or data to piple a fit to FOCEi");
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
    return(collectWarnings(do.call(focei.fit.data.frame0, call, envir=parent.frame(1))))
}

focei.fit.data.frame0 <- function(data,
                                  inits,
                                 PKpars,
                                 diag.xform=c("sqrt", "log", "identity"),
                                 optim=c(
                                     ## "n1qn1",
                                     "L-BFGS-B", ## From optim
                                     "bobyqa",
                                     "lbfgsb3",
                                     ## "newuoa",
                                     "nlminb"##,
                                     ## "uobyqa"
                                 ),
                                 model=NULL,
                                 pred=NULL,
                                 err=NULL,
                                 lower= -Inf,
                                 upper= Inf,
                                 control=list(),
                                 calculate.vars=c("pred", "ipred", "ires", "res", "iwres", "wres", "cwres", "cpred", "cres"),
                                 theta.names=NULL,
                                 eta.names=NULL){
    orig.data <- data;
    colnames(data) <- toupper(names(data));
    do.table <- FALSE;
    sink.file <- tempfile();
    orig.sink.number <- sink.number()
    do.sink <- TRUE;
    fit.df <- NULL;
    sink.close <- function(n=orig.sink.number){
        if (do.sink){
            while(sink.number() > n){
                sink();
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
                sink(sink.file);
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
    on.exit({sink.close(0);
        lines <- sink.get();
        unlink(sink.file);
        if (!is.null(lines)){
            message("After Error:");
            message(paste(paste("##", lines), collapse="\n"));
        }
        running <- FALSE})
    running <- TRUE
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
        atol.outer=1e-6,
        rtol.outer=1e-6,
        maxsteps.ode = 99999,
        reltol.outer = 1e-4,
        absltol.outer = 1e-4,
        cores=2,
        transit_abs=FALSE,
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
        sigdig=5,
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
        save.curve=TRUE,
        numDeriv.method1="simple",
        numDeriv.method2="simple",
        numDeriv.swap=2.3,
        sum.prod=TRUE,
        theta.grad=TRUE
    )

    curi <- 0;

    pt <- proc.time()
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

    numDeriv.method <- con$numDeriv.method1

    cl <- NULL;
    if (.Platform$OS.type == "windows" && con$cores > 1){
        cl <- parallel::makePSOCKcluster(con$cores)
        on.exit({parallel::stopCluster(cl)}, add=TRUE);
    }

    if (con$sigdig == 0 & con$ridge.decay == 0 && !con$extra.output){
        do.sink <- FALSE
    }

    optim.method <- match.arg(optim);
    grad.methods <- c("L-BFGS-B", "lbfgs", "lbfgsb3", "nlminb", "mma", "slsqp", "lbfgs-nlopt", "tnewton_precond_restart",
                      "tnewton_precond", "tnewton", "var1", "var2", "n1qn1")
    print.grad <- any(optim.method == grad.methods);
    if (con$grad && !print.grad){
        if (!con$NOTRUN){
            message("Warning; You selected a gradient method, but the optimization procedure doesn't require the gradient.\nIgnoring gradient")
        }
        con$grad <- FALSE;
    }
    if (con$NOTRUN){
        con$theta.grad <- FALSE;
        print.grad <- FALSE;
    }
    if(is(model, "RxODE") || is(model, "character")) {
        ODEmodel = TRUE
        if (class(pred) != "function"){
            stop("pred must be a function specifying the prediction variables in this model.")
        }
    }
    else {
        ODEmodel = TRUE
        model <- constructLinCmt(PKpars);
        pred <- eval(parse(text="function(){return(Central);}"))
    }

    square = function(x) x*x
    diag.xform <- match.arg(diag.xform)
    diag.xform.inv = c("sqrt"="square", "log"="exp", "identity"="identity")[diag.xform]

    ##process inits
    if (is.null(err)){
        err <-eval(parse(text=paste0("function(){err",paste(inits$ERROR[[1]],collapse=""),"}")));
    }

    ## print(th0.om)
    model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err, grad=con$grad,
                                     pred.minus.dv=con$pred.minus.dv, sum.prod=con$sum.prod,
                                     theta.derivs=con$theta.grad);
    message(sprintf("Original Compartments=%s", length(RxODE::rxState(model$obj))))
    message(sprintf("\t Inner Compartments=%s", length(RxODE::rxState(model$inner))))
    if (con$grad){
        message(sprintf("\t Outer Compartments=%s", length(RxODE::rxState(model$outer))))
    }
    cov.names <- RxODE::rxParams(model$inner);
    cov.names <- cov.names[regexpr(rex::rex(start, or("THETA", "ETA"), "[", numbers, "]", end), cov.names) == -1];
    lhs <- c(names(RxODE::rxInits(model$inner)), RxODE::rxLhs(model$inner))
    if (length(lhs) > 0){
        cov.names <- cov.names[regexpr(rex::rex(start, or(lhs), end), cov.names) == -1];
    }
    if (length(cov.names) > 0){
        if (!all(cov.names %in% names(data))){
            message("Needed Covariates:")
            RxODE::rxPrint(cov.names)
            stop("Not all the covariates are in the dataset.")
        }
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
        RxODE::rxPrint(data.frame(lower,inits.vec,upper));
    }
    names(inits.vec) = NULL
    par.lower <- lower / inits.vec;
    par.upper <- upper / inits.vec;
    if (!con$NOTRUN){
        message("Scaled Boundaries:");
    }
    w.neg <- which(par.lower > par.upper);
    tmp <- par.lower[w.neg];
    par.lower[w.neg] <- par.upper[w.neg];
    par.upper[w.neg] <- tmp;
    if (!con$NOTRUN){
        RxODE::rxPrint(data.frame(par.lower,rep(1, length(inits.vec)),par.upper))
        if (do.sink){
            message("\nKey:")
            message(" S: Scaled Parameter values");
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

    ofv.FOCEi.ind.slow <- function(pars) {
        cur.diff <<- NULL;
        if(con$PRINT.PARS) print(pars)
        pars = pars*inits.vec
        if (!is.null(last.pars) && con$accept.eta.size != 0 && find.best.eta && con$grad){
            inits.mat.bak <- inits.mat;
            inits.mat <- RxODE::rxUpdateEtas(last.dEta.dTheta, matrix(pars - last.pars, ncol=1), inits.mat,con$accept.eta.size);
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
                         pred.minus.dv=con$pred.minus.dv, switch.solver=con$switch.solver,
                         inner.opt=con$inner.opt, add.grad=print.grad, numDeriv.method=numDeriv.method,
                         table=do.table);
            ##method=numDeriv.method
            if (!is.null(inits.mat.bak)){
                args$eta.bak=inits.mat.bak[.wh, ];
            }
            ret <- do.call(getFromNamespace("rxFoceiInner","RxODE"), args)
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
    }
    ofv.FOCEi.ind <- memoise::memoise(ofv.FOCEi.ind.slow)

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
            if (any(optim.method == c("lbfgs"))){
                return(paste0(prefix, ".", RxODE::rxCout(ofv)));
            } else if (any(optim.method == c("L-BFGS-B"))){
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
            } else if (any(optim.method == c("L-BFGS-B"))){
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
    curi <<- 0;
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
    ofv.FOCEi <- function(pars) {
        llik.subj <- ofv.FOCEi.ind(pars)
        first <<- FALSE
        llik <- -2*do.call("sum", llik.subj);
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
                        last.info <- get(last, envir=ofv.cache, inherits=FALSE);
                        message(sprintf("Step %s: %s", curi, substr(last, 3, nchar(last))), appendLF=FALSE);
                        fit.df <<- rbind(fit.df, data.frame(t(c(iter=curi, objf=as.numeric(substr(last, 3, nchar(last))), last.info$pars*inits.vec))))
                        curi <<- curi + 1;
                        if (is.null(ofv.cache$last)){
                            message("")
                        } else if (ofv.cache$last$ofv > last.info$ofv){
                            message(sprintf(" (%s reduction)", format(ofv.cache$last$ofv - last.info$ofv)))
                            reset <- TRUE;
                        } else {
                            message(sprintf(" (%s increase)", format(last.info$ofv - ofv.cache$last$ofv)))
                        }
                        p <- rbind(data.frame(t(c(last.info$pars))),
                                   data.frame(t(c(last.info$pars * inits.vec))),
                                   data.frame(exp(t(c(last.info$pars * inits.vec)))))
                        rn <- c("S:", "U:", "X:")
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
                            sdig <- -log10(abs(cur.pars - last.pars)) - 1;
                            rn <- c(rn, "D:");
                            p <- rbind(p,
                                       data.frame(t(sdig)));
                            rownames(p) <- rn;
                            if ((length(theta.names) + length(eta.names)) == length(names(p))){
                                names(p) <- c(theta.names, paste0("ome(", eta.names, ")"));
                            } else {
                                names(p)[seq_along(nms)] <- nms
                            }
                            print(p)
                            if (all(sdig > con$numDeriv.swap)){
                                numDeriv.method <- con$numDeriv.method2
                            }
                            if (con$sigdig > 0 && all(sdig > con$sigdig)){
                                sigdig.exit <- TRUE
                            }
                            if (is.null(ofv.cache$first)){
                                ofv.cache$first <- llik;
                            }
                            cur.diff <<- (ofv.cache$first - llik);
                        }
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

    gr.FOCEi <- function(pars){
        llik.subj <- ofv.FOCEi.ind(pars)
        llik <- -2*do.call("sum", llik.subj);
        ## if (con$grad){
            if (con$ridge.decay != 0 && is.null(cur.diff)){
                cur.diff <<- 0
                extra <- 1
            } else if (con$ridge.decay != 0){
                extra <- exp(-con$ridge.decay * cur.diff);
            } else {
                extra <- 1
            }
            gr <- -2 * Reduce("+", lapply(llik.subj, function(x){
                                       return(attr(x,"grad"))
                                   })) + pars * con$precision  * extra;
            assign(optim.obj(llik, "l"), gr, envir=ofv.cache, inherits=FALSE);
        return(gr);
        ## } else {
        ##     if (con$fix.eta.for.grad){
        ##         find.best.eta <<- FALSE; ## Keep etas.
        ##         on.exit({find.best.eta <<- FALSE;}, add=TRUE)
        ##     }
        ##     gr <- numDeriv::grad(ofv.FOCEi, pars)
        ##     assign(optim.obj(llik, "l"), gr, envir=ofv.cache, inherits=FALSE);
        ##     return(gr);
        ## }
    }

    s.FOCEi <- function(pars){
        llik.subj <- ofv.FOCEi.ind(pars);
        S <- Reduce('+', lapply(llik.subj, function(x){
                             ## FIXME: when grad=FALSE calculate individual numerical deriavative
                             ## Gradients on normal scale to form cross-products to be summed.
                             return(crossprod(matrix(as.vector(attr(x,"grad"))*inits.vec, nrow=1)));
                         }))
        return(S);
    }
    np <- length(inits.vec)
    meth <- c();

    opt <- function(){
        fit <- NULL;
        if (optim.method == "lbfgs"){
            par <- rep(1, length(inits.vec));
            fit <- try({lbfgs::lbfgs(call_eval=ofv.FOCEi, call_grad=gr.FOCEi, vars=par)});
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
        } else if (optim.method == "L-BFGS-B"){
            par <- rep(1, length(inits.vec));
            fit <- try({stats::optim(par=par, fn=ofv.FOCEi, gr=gr.FOCEi,
                                     method=optim.method,
                                     control=list(trace=1, REPORT=1),
                                     lower=par.lower,
                                     upper=par.upper,
                                     hessian=regexpr("r", con$cov.method) != -1)});
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
            prm <- rep(1, length(inits.vec));
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
            fit <- minqa::bobyqa(rep(1, length(inits.vec)),
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
            fit <- minqa::newuoa(rep(1, length(inits.vec)),
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
            fit <- minqa::uobyqa(rep(1, length(inits.vec)),
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
        } else if (optim.method=="n1qn1"){
            fit <- try({n1qn1(ofv.FOCEi, gr.FOCEi,
                              rep(1, length(inits.vec)),
                              invisible=0)});
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
            fit <- try({nlminb(rep(1, length(inits.vec)),
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
        stop("Mismatch of dimensions")
    }
    if (con$inner.opt == "n1qn1"){
        inits.c.hess <- matrix(0, nSUB, nETA * (nETA + 13) / 2)
    }
    if (.Platform$OS.type == "windows" && con$cores > 1){
        ## Copy over all of the objects within scope to
        ## all clusters.
        ## Taken from https://github.com/nathanvan/mcmc-in-irt/blob/master/post-10-mclapply-hack.R

        this.env <- environment()
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
        pt <- proc.time()
        fit = list()
        fit$par = rep(1, length(inits.vec))
        fit$value = as.vector(ofv.FOCEi(rep(1, length(inits.vec))))
        fit$convergence = -99
        fit$message = "notrun"
        optim.time <- proc.time() - pt;
        fit$cov.time <- proc.time() - proc.time();
        fit$objective <- fit$value;
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
            m = diag(inits.vec);
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
            ## Use First Order condition for covariance
            ## inits.mat <- matrix(0, nSUB, nETA)
            if (con$cov.method != "s"){
                if (any(names(fit) == "hessian")){
                    R1 <- fit$hessian
                } else {
                    R1 <- optimHess(fit$par, ofv.FOCEi, gr.FOCEi);
                }
                Rinv <- RxODE::rxInv(R1)
                if (det(R1) <= 0){
                    warning("Non positive-definite Hessian matrix when calculating the covariance; Falling back to S-matrix covariance")
                    fit$hessian.bad <- R1;
                    ## R1 <- as.matrix(Matrix::nearPD(R1)$mat);
                    ## Rinv <- RxODE::rxInv(R1);
                    con$cov.method <- "s"
                    fit$warning <- "Non positive-definite Hessian matrix when calculating the covariance; Falling back to S-matrix covariance"
                    Rinv <- NULL
                } else {
                    Rinv = m %*% Rinv %*% m
                }
            } else {
                Rinv <- NULL
            }
            ## if (con$grad){
            s = s.FOCEi(fit$par)
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
        fit$control <- con;
        fit$par.unscaled = fit$par*inits.vec
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
        if (length(eta.names) == dim(fit$omega)[1]){
            dimnames(fit$omega) <- list(eta.names, eta.names);
            fit$eta.names <- eta.names;
        }
        if (!is.null(fit.df)){
            names(fit.df) <- c("iter", "objf", nms, eta.names)
            fit$fit.df <- fit.df;
        }
        w <- seq_along(nms)
        if (con$NOTRUN){
            fit$par.data.frame <- data.frame(est=fit$theta, "exp(est)"=exp(fit$theta), row.names=nms);
        } else {
            fit$par.data.frame <- data.frame(est=fit$theta, "exp(est)"=exp(fit$theta), "se(est)"=fit$se[w], "%cv(est)"=fit$se[w] / fit$theta * 100,
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
    ofv.FOCEi.ind <- ofv.FOCEi.ind.slow;
    running <- FALSE;
    sink.close()
    con$cores <- 1; ## don't run parallel for extracting infomration
    find.best.eta <- FALSE;
    print.grad <- FALSE; ## Turn off slow gradient calculation
    old.model <- model;
    model$theta <- NULL;
    model$outer <- NULL;
    do.table <- TRUE
    message("Calculating Table Variables...")
    pt <- proc.time();
    if (any("ipred" == calculate.vars)){
        ## message("\tIPRED", appendLF=FALSE)
        data$IPRED <- fitted(data, population=FALSE);
        calculate.vars <- calculate.vars[calculate.vars != "ipred"];
        ## message("done.")
    }
    if (any("pred" == calculate.vars)){
        ## message("\tPRED", appendLF=FALSE)
        data$PRED <- fitted(data, population=TRUE)
        calculate.vars <- calculate.vars[calculate.vars != "pred"];
        ## message("done.")
    }
    for (v in calculate.vars){
        if (v != "pred"){
            ## message(sprintf("\t%s", v), appendLF=FALSE)
            data[, toupper(v)] <- resid(data, type=v);
            ## message("done.")
        }
    }
    if (con$add.posthoc){
        etas <- fitted(data, type="posthoc")
        data <- merge(data, etas);
        ## Drops the class/environment; Put back in.
        attr(data, ".focei.env") <- env;
        class(data) <- c("focei.fit", "data.frame")
    }
    ## data <- merge(data, m);
    table.time <- proc.time() - pt;
    fit$table.time <- table.time;
    fit$data.names <- names(data);
    fit$data.len <- length(data[, 1]);
    fit$time <- data.frame(setup=setup.time["elapsed"],
                           optimize=optim.time["elapsed"],
                           covariance=fit$cov.time["elapsed"],
                           table=fit$table.time["elapsed"],
                           row.names="");
    model <- old.model;
    fit$model <- model;
    message("done")
    data
}


##' @export
`$.focei.fit` <-  function(obj, arg, exact = TRUE){
    m <- as.data.frame(obj);
    ret <- m[[arg, exact = exact]];
    if (is.null(ret)){
        env <- attr(obj, ".focei.env");
        fit <- env$fit;
        return(fit[[arg, exact = exact]])
    } else {
        return(ret)
    }
}

##' @importFrom utils str
##' @export
str.focei.fit <- function(object, ...){
    message('FOCEI combined dataset and list');
    m <- as.data.frame(object);
    str(m)
    env <- attr(object, ".focei.env");
    fit <- env$fit;
    str(fit)
}

##' @export
head.focei.fit <- function(object, n=6, ...){
    message('FOCEI combined dataset and list');
    env <- attr(object, ".focei.env");
    message('Head List:')
    fit <- env$fit;
    print(head(fit, n=floor(n / 2), ...))
    message('Head Data (Also accessible by $fit):')
    m <- as.data.frame(object);
    return(head(m, n=floor(n / 2), ...))
}

##' @export
tail.focei.fit <- function(object, n=6, ...){
    message('FOCEI combined dataset and list');
    env <- attr(object, ".focei.env");
    message('Tail List:')
    fit <- env$fit;
    print(tail(fit, n=floor(n / 2), ...))
    message('Tail Data (Also accessible by $fit):')
    m <- as.data.frame(object);
    return(tail(m, n=floor(n / 2), ...))
}


##' Convert fit to FOCEi style fit
##'
##' @param object Fit object to convert to FOCEi-style fit.
##' @param uif Unified Interface Function
##' @param pt Proc time object
##' @param ...  Other Parameters
##' @return A FOCEi fit style object.
##' @author Matthew L. Fidler
as.focei <- function(object, uif, pt=proc.time(), ...){
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
            aod <- data.frame(aod, Test = effects, L.Ratio = c(NA,
                                                               lratio), `p-value` = c(NA, pval), check.names = FALSE)
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
