##' Find the assignments in R expression
##'
##' @param x R expression
##' @return list of assigned parameters
##' @author Hadley Wickham and Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrfindLhs <- function(x) {
    ## Modified from http://adv-r.had.co.nz/Expressions.html find_assign4
    if (is.atomic(x) || is.name(x)) {
        character()
    } else if (is.call(x)) {
        if ((identical(x[[1]], quote(`<-`)) ||
             identical(x[[1]], quote(`=`))) &&
            is.name(x[[2]])) {
            lhs <- as.character(x[[2]])
        } else {
            lhs <- character()
        }
        unique(c(lhs, unlist(lapply(x, nlmixrfindLhs))))
    } else if (is.pairlist(x)) {
        unique(unlist(lapply(x, nlmixrfindLhs)))
    } else {
        stop("Don't know how to handle type ", typeof(x),
             call. = FALSE)
    }
}
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
        cat(sprintf("nlmixr FOCEI fit (%s)\n\n", ifelse(fit$control$grad, "with global gradient", "without global gradient")));
        if (any(names(fit) == "condition.number")){
            print(data.frame(OBJF=fit$objective, AIC=AIC(x), BIC=BIC(x), "Condition Number"=fit$condition.number,
                             row.names="", check.names=FALSE))
        } else {
            print(data.frame(OBJF=fit$objective, AIC=AIC(x), BIC=BIC(x),
                             row.names="", check.names=FALSE))
        }
        cat("\nTime (sec):\n");
        print(fit$time);
        cat("\nParameters:\n")
        print(fit$par.data.frame);
        cat("\nOmega:\n");
        print(fit$omega);
        is.dplyr <- requireNamespace("dplyr", quietly = TRUE);
        if (!is.dplyr){
            cat("\nFit Data (head):\n")
            print(head(as.matrix(x), n = n));
        } else {
            cat("\nFit Data:\n")
            print(dplyr::as.tbl(x), n = n, width = width);
        }
    } else {
        print(as.data.frame(x));
    }
}

##' Extract log likelihood for focei fit
##'
##' @param fit focei fit object
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
##' @param ...  other parameters
##' @param population when true (default false), calculate the
##'     population predictions; When false, calculate the individual
##'     predictions.
##'
##'     If this is a matrix with the number columns equal to the
##'     number of ETAs in model, and the number of rows equals to the
##'     number of subjects in the dataset, then these etas will be
##'     used to fit the data.
##' @return Individual predictions
##' @author Matthew L. Fidler
##' @export
fitted.focei.fit <- function(object, ..., population=FALSE,
                             type=c("fitted", "Vi", "Vfo", "posthoc")){
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
    ## print(attr(x,"subj"))
    if (type == "posthoc"){
        d1 <- data.frame(do.call("rbind", lapply(attr(x,"subj"), function(s) attr(s, type))));
        names(d1) <- paste0("ETA", seq_along(d1[1, ]))
        d1 <- data.frame(ID=unique(object$ID), d1)
        return(d1)
    } else {
        return(unlist(lapply(attr(x,"subj"), function(s) attr(s, type))))
    }
}

##' Extract residuals from the FOCEI fit
##'
##' @param object focei.fit object
##' @param ... Additional arguments
##' @param type Type of resduals fitted.
##' @return residuals
##' @author Matthew L. Fidler
##' @export
residuals.focei.fit <- function(object, ..., type=c("ires", "res", "iwres", "wres", "cwres")){
    dat <- object;
    DV <- dat$DV;
    type <- match.arg(type);
    up.type <- toupper(type)
    if (any(up.type == names(object))){
        return(dat[, object]);
    }
    if (any(type == c("res", "wres", "cwres"))){
        pred <- fitted(object, population=TRUE);
        if (type == "wres"){
            ## population=TRUE => eta=0
            W <- sqrt(fitted(object, population=TRUE, type="Vfo") + fitted(object, population=TRUE, type="Vi"))
            return((DV - pred) / W);
        } else if (type == "cwres"){
            W <- sqrt(fitted(object, population=FALSE, type="Vfo") + fitted(object, population=TRUE, type="Vi"))
            return((DV - pred) / W);
        } else {
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
    dat <- as.data.frame(x);
    d1 <- data.frame(DV=dat$DV, stack(dat[, c("PRED", "IPRED")]))

    p1 <- ggplot2::ggplot(d1,ggplot2::aes(values,DV)) +
        ggplot2::facet_wrap(~ind) +
        ggplot2::geom_abline(slope=1, intercept=0, col="red", size=1.2) +
        ggplot2::geom_smooth(col="blue", lty=2, formula=DV ~ values + 0, size=1.2) +
        ggplot2::geom_point() +
        ggplot2::xlab("Predictions");

    print(p1)

    p2 <- ggplot2::ggplot(dat, ggplot2::aes(x=IPRED, y=IRES)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope=0, intercept=0, col="red")
    print(p2)

    ids <- unique(dat$ID)
    for (i  in seq(1, length(ids), by=16)){
        tmp <- ids[seq(i, i + 15)]
        tmp <- tmp[!is.na(tmp)];
        d1 <- dat[dat$ID %in% tmp, ];

        p3 <- ggplot2::ggplot(d1, ggplot2::aes(x=TIME, y=DV)) +
            ggplot2::geom_point() +
            ggplot2::geom_line(ggplot2::aes(x=TIME, y=IPRED), col="red", size=1.2) +
            ggplot2::geom_line(ggplot2::aes(x=TIME, y=PRED), col="blue", size=1.2) +
            ggplot2::facet_wrap(~ID)
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
##' @param inits Initilization list
##' @param PKpars Pk Parameters
##' @param diag.xform
##' @param optim
##' @param model
##' @param pred
##' @param err
##' @param lower
##' @param upper
##' @param control
##' @param calculate.vars
##' @return
##' @author Matthew L. Fidler
##' @export
focei.fit <- function(data,
                      inits,
                      PKpars,
                      diag.xform=c("sqrt", "log", "identity"),
                      optim=c(
                          "bobyqa",
                          "lbfgsb3",
                          "newuoa",
                          "nlminb",
                          "uobyqa",
                          ## These are nlopt functions...
                          "praxis",
                          "mma", ## Errors (Grad method)
                          "slsqp", ## Errors
                          "lbfgs-nlopt", ## Errors
                          "tnewton_precond_restart",
                          "tnewton_precond",
                          "tnewton",
                          "var2",
                          "var1",
                          "cobyla-nlopt",
                          "bobyqa-nlopt",
                          "newuoa-nlopt",
                          "nelder-mead-nlopt",
                          "sbplx"),
                      model=NULL,
                      pred=NULL,
                      err=NULL,
                      lower= -Inf,
                      upper= Inf,
                      control=list(),
                      calculate.vars=c("pred", "ipred", "ires", "res", "iwres", "wres", "cwres")){
    data <- data;
    colnames(data) <- toupper(names(data));
    sink.file <- tempfile();
    orig.sink.number <- sink.number()
    do.sink <- TRUE;
    sink.close <- function(n=orig.sink.number){
        if (do.sink){
            while(sink.number() > n){
                sink();
            }
        }
        ## cat("Closed sink\n");
        ## unlink(sink.file);
    }
    sink.start <- function(do.it=TRUE){
        ## sink.close();
        if (do.sink){
            if (do.it){
                ## cat("Starting Sink to", sink.file, "\n");
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
            cat("After Error:\n");
            cat(paste(paste("##", lines), collapse="\n"), "\n");
        }
        options(RxODE.warn.on.assign=oldAssign);
        running <- FALSE})
    oldAssign <- getOption("RxODE.warn.on.assign");
    options(RxODE.warn.on.assign=FALSE);
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
        atol.ode=1e-13,
        rtol.ode=1e-13,
        atol.outer=1e-8,
        rtol.outer=1e-8,
        maxsteps.ode = 99999,
        reltol.outer = 1e-4,
        cores=parallel::detectCores(),
        transit_abs=FALSE,
        NONMEM=TRUE,
        NOTRUN=FALSE,
        PRINT.PARS=FALSE,
        ## NONMEM does PRED-DV
        ## Almquist does DV-PRED
        pred.minus.dv=TRUE,
        switch.solver=FALSE,
        cov.method="grad",
        factr=1e10,
        grad=TRUE,
        accept.eta.size=1.5,
        sigdig=0,
        precision=0, ## 0 = No ridge penalty.
        ridge.decay=0, ## 0 = no decay; Inf = No ridge
        reset.precision=NULL,
        eigen=TRUE,
        fix.eta.for.grad=FALSE,
        rhobeg=.2,
        rhoend=1e-2,
        npt=NULL,
        est.chol.omegaInv=TRUE,
        add.posthoc=TRUE,
        extra.output=TRUE ## Display extra output on each iteration
    )

    curi <- 0;

    cl <- NULL;
    if (.Platform$OS.type == "windows" && con$cores > 1){
        cl <- parallel::makePSOCKcluster(con$cores)
        on.exit({parallel::stopCluster(cl)}, add=TRUE);
    }

    pt <- proc.time()
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

    if (con$sigdig == 0 & con$ridge.decay == 0 && !con$extra.output){
        do.sink <- FALSE
    }

    optim.method <- match.arg(optim);
    grad.methods <- c("lbfgs", "lbfgsb3", "nlminb", "mma", "slsqp", "lbfgs-nlopt", "tnewton_precond_restart",
                      "tnewton_precond", "tnewton", "var1", "var2")
    print.grad <- any(optim.method == grad.methods);
    if (con$grad && !print.grad){
        cat("Warning; You selected a gradient method, but the optimization procedure doesn't require the gradient.\nIgnoring gradient\n")
        con$grad <- FALSE;
    }
    if(class(model)=="RxODE") {
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
    model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err, grad=con$grad, pred.minus.dv=con$pred.minus.dv, run.internal=TRUE);
    cov.names <- RxODE::rxParams(model$inner);
    cov.names <- cov.names[regexpr(rex::rex(start, or("THETA", "ETA"), "[", numbers, "]", end), cov.names) == -1];
    if (length(cov.names) > 0){
        if (!all(cov.names %in% names(data))){
            stop("Not all the covariates are in the dataset.")
        }
    }

    ## RxODE(rxNorm(model$inner), modName="test");
    nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)),
             sprintf("ERR[%s]", seq_along(model$extra.pars)))
    if (length(lower) == 1){
        lower <- rep(lower, length(inits$THTA));
    } else if (length(lower) != length(inits$THTA)){
        stop("Lower must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
    }
    if (length(upper) == 1){
        upper <- rep(upper, length(inits$THTA));
    } else if (length(lower) != length(inits$THTA)){
        stop("Upper must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
    }

    extra.pars <- eval(call(diag.xform, model$extra.pars))
    if (length(model$extra.pars) > 0){
        inits$THTA <- c(inits$THTA, extra.pars);
        lower.err <- rep(con$atol.ode * 10, length(model$extra.pars));
        upper.err <- rep(Inf, length(model$extra.pars));
        lower <-c(lower, lower.err);
        upper <- c(upper, upper.err);
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
    lower <- c(lower, lower.om);
    upper <- c(upper, upper.om);
    inits.vec = c(inits$THTA, th0.om)
    if (any(inits.vec == 0)){
        warning("Some of the initial conditions were 0, chainging to 0.0001");
        inits.vec[inits.vec == 0] <- 0.0001;
    }
    cat("Boundaries:\n");
    print(data.frame(lower,inits.vec,upper))
    names(inits.vec) = NULL
    par.lower <- lower / inits.vec;
    par.upper <- upper / inits.vec;
    cat("Scaled Boundaries:\n");
    w.neg <- which(par.lower > par.upper);
    tmp <- par.lower[w.neg];
    par.lower[w.neg] <- par.upper[w.neg];
    par.upper[w.neg] <- tmp;
    print(data.frame(par.lower,rep(1, length(inits.vec)),par.upper))
    if (do.sink){
        cat("\nKey:\n")
        cat(" S: Scaled Parameter values\n");
        cat(" U: Unscaled Parameter values\n");
        if(print.grad){
            cat(" G: Gradient\n");
        }
        cat(" D: Significant Figures\n");
        cat(" Optimization output displayed with comments, i.e. ##\n\n");
    } else {
        cat("\n");
    }
    nTHTA = nlini[1]
    nETA  = nrow(om0)
    ID.all = unique(data[,"ID"])        #FIXME
    ID.ord = order(ID.all)
    names(ID.ord) = ID.all
    nSUB  = length(ID.all)

    find.best.eta <- TRUE;

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
            if(con$DEBUG>10 && .wh==1) print(inits.mat[.wh,])               #FIXME
            if (con$DEBUG.ODE) print("i'm here :)")
            args <- list(model, ev, theta=THETA, eta=inits.mat[.wh,],
                         dv=dati$dv[ev$get.obs.rec()], inv.env=rxSymEnv,
                         nonmem=con$NONMEM, invisible=1-con$TRACE.INNER, epsilon=con$TOL.INNER,
                         id=subj, inits.vec=inits.vec, cov=cur.cov, estimate=find.best.eta,
                         atol=con$atol.ode, rtol=con$rtol.ode, maxsteps=con$maxsteps.ode,
                         atol.outer=con$atol.outer, rtol.outer=con$rtol.outer,
                         pred.minus.dv=con$pred.minus.dv, switch.solver=con$switch.solver);
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
                gc();
                llik.subj <- parallel::mclapply(X = ID.all, FUN = llik.one, cores = con$cores, con=con)
            } else {
                llik.subj <- parallel::parLapply(cl, X = ID.all, fun = llik.one, con=con)
                ## nocov end
            }
        } else {
            llik.subj <- lapply(ID.all, llik.one, con=con);
        }
        ## Use wenping's solution for inits.mat
        m = t(sapply(llik.subj, function(x) {
            c(attr(x, "wh"), attr(x, "posthoc"))
        }))
        if(con$RESET.INITS.MAT) inits.mat[m[,1],] <<- m[,-1]


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
            if (any(optim.method == c("nlminb", "newuoa", "bobyqa", "uobyqa"))){
                return(paste0(prefix, ".", gsub("^ *", "", sprintf("%#14.8g",ofv))));
    } else if (optim.method == "lbfgsb3"){
        ## This uses dblepr("",)
        ##  call dblepr(" f=",-1, f, 1)
        ## According to r documentation this is the same as printing a number.
        ## I believe this is the same as format(#)
        return(sprintf("%s.%s", prefix, format(last.ofv)));
    } else {
        ## The R interface uses Rprintf("%f") for nloptr
        return(sprintf("%s.%s", prefix, sprintf("%f", last.ofv)));
    }
        } else {
            if (any(optim.method == c("newuoa", "nlminb", "bobyqa", "uobyqa"))){
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
            cat("## Hessian(s)  inverted with pinv()\n");
        }
        lines <- lines[regexpr(rex::rex(start, any_spaces, end),
                               lines) == -1];

        lines <- lines[regexpr(rex::rex(start, any_spaces, end),
                               lines) == -1];
        if (cor.det || cor.nearPD){
            cat("## Hessian(s) corrected by ");
            if (cor.det && cor.nearPD){
                cat("both det and nearPD\n");
            } else if (cor.nearPD){
                cat("nearPD\n");
            } else {
                cat("det\n");
            }
        }
        if (bad.inv.H2){
            cat("## Hessian(s) from second order derivatives inverted with pinv()\n");
        }
        cat(paste(paste(paste0("## ", lines), collapse="\n"), "\n"));
        return(lines);
    }

    ofv.FOCEi <- function(pars) {
        llik.subj <- ofv.FOCEi.ind(pars)
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
                        cat(sprintf("Step %s: %s", curi, substr(last, 3, nchar(last))));
                        curi <<- curi + 1;
                        if (is.null(ofv.cache$last)){
                        } else if (ofv.cache$last$ofv > last.info$ofv){
                            cat(sprintf(" (%s reduction)", format(ofv.cache$last$ofv - last.info$ofv)))
                            reset <- TRUE;
                        } else {
                            cat(sprintf(" (%s increase)", format(last.info$ofv - ofv.cache$last$ofv)))
                        }
                        cat(paste0("\n S: ", paste(sapply(last.info$pars, function(x){sprintf("%#8g", x)}), collapse=" "), "\n"))
                        cat(paste0(" U: ", paste(sapply(last.info$pars*inits.vec, function(x){sprintf("%#8g", x)}), collapse=" "), "\n"))
                            grad.txt <- optim.obj(last.info$llik, "l");
                        if (exists(grad.txt, envir=ofv.cache, inherits=FALSE)){
                            last.grad <- get(grad.txt, envir=ofv.cache, inherits=FALSE);
                            cat(paste0(" G: ", paste(sapply(last.grad, function(x){sprintf("%#8g", x)}), collapse=" "), "\n"))
                            }
                        if (is.null(ofv.cache$last)){
                            reset <- TRUE;
                        } else {
                            last.pars <- ofv.cache$last$pars
                            cur.pars <- last.info$pars
                            if (all(last.pars == cur.pars) && !is.null(ofv.cache$last1)){
                                cat("## Parameters the same, use last iteration for sigdig caluclation...\n");
                                last.pars <- ofv.cache$last1$pars;
                                ofv.cache$last <- ofv.cache$last1;
                            }
                            ##
                            ## NONMEM's sigdigs as per http://cognigencorp.com/nonmem/current/2012-October/3654.html
                            ##
                            sdig <- -log10(abs(cur.pars - last.pars)) - 1;
                            cat(paste0(" D: ", paste(sapply(sdig, function(x){sprintf("%f", x)}), collapse=" "), "\n"))
                            if (con$sigdig > 0 && all(sdig > con$sigdig)){
                                sigdig.exit <- TRUE
                            }
                            if (is.null(ofv.cache$first)){
                                ofv.cache$first <- llik;
                            }
                            cur.diff <<- (ofv.cache$first - llik);
                        }
                    } else {
                        cat("Warning: Prior objective function not found...\n");
                        reset <- TRUE;
                    }
                    cat("\n");
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
            cat(sprintf("Current ridge penalty: %s\n", last.penalty));
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
        if (con$grad){
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
        } else {
            if (con$fix.eta.for.grad){
                find.best.eta <<- FALSE; ## Keep etas.
                on.exit({find.best.eta <<- FALSE;}, add=TRUE)
            }
            gr <- numDeriv::grad(ofv.FOCEi, pars)
            assign(optim.obj(llik, "l"), gr, envir=ofv.cache, inherits=FALSE);
            return(gr);
        }
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
    nlopt.gradient <- c("mma", "slsqp", "nlopt_lbfgs",
                        "tnewton_precond_restart", "tnewton_precond", "nlopt_ld_tnewton_restart",
                        "tnewton", "var2", "var1")
    meth["mma"] <- "NLOPT_LD_MMA";
    meth["slsqp"] <- "NLOPT_LD_SLSQP";
    meth["lbfgs-nlopt"] <-  "NLOPT_LD_LBFGS";
    meth["tnewton_precond_restart"] <- "NLOPT_LD_TNEWTON_PRECOND_RESTART"
    meth["tnewton_precond"] <- "NLOPT_LD_TNEWTON_PRECOND"
    meth["tnewton_restart"] <- "NLOPT_LD_TNEWTON_RESTART"
    meth["tnewton"] <- "NLOPT_LD_TNEWTON"
    meth["var2"] <- "NLOPT_LD_VAR2"
    meth["var1"]  <- "NLOPT_LD_VAR1"
    nlopt.gradient.free <- c("cobyla-nlopt", "bobyqa-nlopt", "newuoa-nlopt", "praxis", "nelder-mead-nlopt", "sbplx");
    meth["cobyla-nlopt"] <- "NLOPT_LN_COBYLA";
    meth["newuoa-nlopt"] <- "NLOPT_LN_NEWUOA_BOUND";
    meth["bobyqa-nlopt"] <- "NLOPT_LN_BOBYQA";
    meth["praxis"] <- "NLOPT_LN_PRAXIS";
    meth["nelder-mead-nlopt"] <- "NLOPT_LN_NELDERMEAD";
    meth["sbplx"] <- "NLOPT_LN_SBPLX"
    opt <- function(){
        fit <- NULL;
        if (any(optim.method == c(nlopt.gradient.free, nlopt.gradient))){
            if (any(optim.method == nlopt.gradient)){
                fit <- try({nloptr::nloptr(rep(1, length(inits.vec)),
                                           ofv.FOCEi,
                                           gr.FOCEi,
                                           lb=par.lower,
                                           ub=par.upper,
                                           opts=list(algorithm=meth[optim.method],
                                                     "print_level"=2,
                                                     ftol_rel=con$reltol.outer))})
            } else {
                if (any(optim.method == c("bobyqa-nlopt"))){
                    local <- list("step1"=con$rhobeg);
                } else {
                    local <- c();
                }
                fit <- nloptr::nloptr(rep(1, length(inits.vec)),
                               ofv.FOCEi,
                               lb=par.lower,
                               ub=par.upper,
                               opts=list(algorithm=meth[optim.method],
                                         "print_level"=2,
                                         ftol_rel=con$reltol.outer,
                                         "local_opts"=local))
                ## fit <- try({});
            }
            if (inherits(fit, "try-error") && !is.null(sigdig.fit)){
                if (attr(fit, "condition")$message == "sigidig exit"){
                    fit <- sigdig.fit
                } else {
                    lines <- sink.get();
                    if (!is.null(lines))
                        message(paste0(lines, collapse="\n"), "\n");
                    fit <- NULL;
                }
            } else {
                fit$par <- fit$solution
            }
        } else if (optim.method=="lbfgsb3"){
            prm <- rep(1, length(inits.vec));
            fit <- try({lbfgsb3:::lbfgsb3(prm=prm,
                                          fn=ofv.FOCEi, gr=gr.FOCEi,
                                          control=list(trace=1## , pgtol = con$reltol.outer
                                                       ),
                                          lower=par.lower,
                                          upper=par.upper
                                          )})
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
        } else if (optim.method=="nlminb") {
            fit <- try({nlminb2(rep(1, length(inits.vec)),
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
    inits.mat <- matrix(0, nSUB, nETA)
    if (.Platform$OS.type == "windows"){
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
            if (con$fix.eta.for.grad){
                find.best.eta <- FALSE; ## Keep etas.
            }
            if (con$cov.method=="hessian" && any(names(fit) == "Hessian.inv")) {
                message("Calulate covariance...")
                sink.start();
                R2 <- fit$Hessian.inv;
                fit$Hessian.inv = NULL
            } else if (any(con$cov.method==c("hessian", "grad")) && con$grad){
                message("Calulate covariance...")
                sink.start();
                ## Use First Order condition for covariance
                old.mat <- inits.mat;
                inits.mat <- matrix(0, nSUB, nETA)
                R1 <- optimHess(fit$par, ofv.FOCEi, gr.FOCEi);
                inits.mat <- old.mat;
                Rinv <- RxODE::rxInv(R1)
            } else {
                message("Calulate covariance...")
                sink.start();
                R1 <- optimHess(fit$par, ofv.FOCEi);
                Rinv <- RxODE::rxInv(R1)
            }
            if (con$cov.method!="hessian"){
                if (det(R1) <= 0){
                    warning("Non positive-definite Hessian matrix when calculating the covariance; Correcting with nearPD")
                    fit$hessian.bad <- R1;
                    R1 <- as.matrix(Matrix::nearPD(R1)$mat);
                    Rinv <- RxODE::rxInv(R1);
                    fit$warning <- "Non positive-definite Hessian matrix when calculating the covariance; Correcting with nearPD"
                }
            }
            Rinv = m %*% Rinv %*% m
            if (con$grad){
                s = s.FOCEi(fit$par)
                fit$cov.r.s <- Rinv %*% s.FOCEi(fit$par) %*% Rinv
                fit$cov.r <- 2 * Rinv;
                fit$cov.s <- 4 * RxODE::rxInv(s);
                fit$cov <- fit$cov.r;
            } else {
                fit$cov <- 2 * Rinv
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
                    cat("## Hessian(s) corrected by ");
                    if (cor.det && cor.nearPD){
                        cat("both det and nearPD\n");
                    } else if (cor.nearPD){
                        cat("nearPD\n");
                    } else {
                        cat("det\n");
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
        for (n in c("fval")){
            if (any(names(fit) == n)){
                fit$objective <- fit[[n]];
            }
        }
        lD <- fit$par.unscaled[-seq_along(nms)];
        rxSymEnv <-  RxODE::rxSymInv(rxSym, lD);
        fit$omega <- rxSymEnv$omega;
        w <- seq_along(nms)
        if (con$NOTRUN){
            fit$par.data.frame <- data.frame(est=fit$theta, row.names=nms);
        } else {
            fit$par.data.frame <- data.frame(est=fit$theta, se=fit$se[w], "%cv"=fit$se[w] / fit$theta * 100,
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
    cat("Calculating Table Variables...\n")
    pt <- proc.time();
    if (any("ipred" == calculate.vars)){
        data$IPRED <- fitted(data, population=FALSE);
        calculate.vars <- calculate.vars[calculate.vars != "ipred"];
    }
    if (any("pred" == calculate.vars)){
        data$PRED <- fitted(data, population=TRUE);
        calculate.vars <- calculate.vars[calculate.vars != "pred"];
    }
    for (v in calculate.vars){
        if (any(v == c("ires", "res", "iwres", "wres", "cwres"))){
            data[, toupper(v)] <- resid(data, type=v);
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
    cat("done\n")
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
    cat('FOCEI combined dataset and list\n');
    m <- as.data.frame(object);
    str(m)
    env <- attr(object, ".focei.env");
    fit <- env$fit;
    str(fit)
}
