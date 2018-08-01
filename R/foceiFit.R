.regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                          group(any_of("0":"9"), ".", some_of("0":"9"))),
                       maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9"))));
.regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"));
.regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
.regNum <- rex::rex(maybe("-"), or(.regDecimalint, .regFloat1, .regFloat2))

##' Control Options for FOCEi
##'
##' @param sigdig Optimization significant digits. This controls:
##'
##'  \itemize{
##' \item Defaults for optimization and ODE solving
##' \itemize{
##' \item The tolerance of the inner and outer optimization is 10^-sigdig
##' \item The tolerance of the ODE solvers is 10^(-sigdig-1)
##' }
##' \item The significant figures that some tables are rounded to.
##' }
##'
##' @param epsilon Precision of estimate for n1qn1 optimization.

##' @param maxstepsOde Maximum number of steps for ODE solver.
##'
##'
##' @param printInner Integer representing when the inner step is
##'     printed. By default this is 0 or do not print.  1 is print
##'     every function evaluation, 5 is print every 5 evaluations.
##'
##' @param printOuter Integer representing when the outer step is
##'     printed. When this is 0 or do not print the iterations.  1 is
##'     print every function evaluation (default), 5 is print every 5
##'     evaluations.
##'
##' @param scaleTo Scale the initial parameter estimate to this value.
##'     By default this is 1.  When zero or below, no scaling is performed.
##'
##' @param scaleObjective Scale the initial objective function to this
##'     value.  By default this is 1.  When \code{scaleObjective} is
##'     greater than zero, this scaling is performed by:
##'
##'      scaledObj = (currentObj / |initialObj|) * scaleObjective
##'
##'     Therefore, if the initial objective function is negative, the
##'     initial scaled objective function would be negative as well.
##'     When \code{scaleObjective} is less than zero, no scaling is
##'     performed.
##'
##' @param derivEps Central/Forward difference tolerances, which is a
##'     vector of relative difference and absolute difference.  The
##'     central/forward difference step size h is calculated as:
##'
##'         h = abs(x)*derivEps[1]+derivEps[2]
##'
##' @param derivMethod indicates the method for calculating
##'     derivatives of the outer problem.  Currently supports
##'     "central" and "forward" difference methods.
##'
##' @param covDerivMethod indicates the method for calculating the
##'     derivatives while calculating the covariance components
##'     (Hessian and S).
##'
##' @param covMethod Method for calculating covariance.  In this
##'     discussion, R is the Hessian matrix of the objective
##'     function. The S matrix is the sum of each individual's
##'     gradient cross-product (evaluated at the individual empirical
##'     Bayes estimates).
##' \itemize{
##' \item "\code{r,s}" Uses the sandwich matrix to calculate the covariance, that is: \code{R^-1 * S * R^-1}
##'
##' \item "\code{r}" Uses the Hessian matrix to calculate the
##'      covariance as \code{2*R^-1}
##' \item "\code{s}" Uses the crossproduct matrix to calculate the covariance as \code{4*S^-1}
##' \item "" Does not calculate the covariance step.
##' }
##'
##' @param lbfgsLmm An integer giving the number of BFGS updates
##'     retained in the "L-BFGS-B" method, It defaults to 40.
##'
##' @param lbfgsPgtol is a double precision variable.
##'
##'     On entry pgtol >= 0 is specified by the user.  The iteration
##'     will stop when:
##'
##'        \code{max{|proj g_i | i = 1, ..., n} <= lbfgsPgtol}
##'
##'     where pg_i is the ith component of the projected gradient.
##'
##'     On exit pgtol is unchanged.  This defaults to zero, when the
##'     check is suppressed.
##'
##' @param lbfgsFactr Controls the convergence of the "L-BFGS-B"
##'     method.  Convergence occurs when the reduction in the
##'     objective is within this factor of the machine
##'     tolerance. Default is 1e10, which gives a tolerance of about
##'     \code{2e-6}, approximately 4 sigdigs.  You can check your
##'     exact tolerance by multiplying this value by
##'     \code{.Machine$double.eps}
##'
##' @param diagXform This is the transformation used on the diagonal
##'     of the \code{chol(inv(omega))}. This matrix and values are the
##'     parameters estimated in FOCEi. The possibilities are:
##'
##' \itemize{
##' \item sqrt Estimates the sqrt of the diagonal elements of \code{chol(inv(omega))}.  This is the default method.
##' \item log Estimates the log of the diagonal elements of \code{chol(inv(omega))}
##' \item identity Estimates the diagonal elements without any transformations
##' }
##'
##' c("sqrt", "log", "identity"),
##'
##' @param sumProd Is a boolean indicating if the model should change
##'     multiplication to high precision multiplication and sums to
##'     high precision sums using the PreciseSums package.  By default
##'     this is \code{FALSE}.
##'
##'
##' @param optExpression Optimize the RxODE expression to speed up
##'     calculation. By default this is turned on.
##'
##' @param ci Confidence level for some tables.  By default this is
##'     0.95 or 95% confidence.
##'
##' @param ...
##' @param maxInnerIterations Number of iterations for n1qn1
##'     optimization.
##'
##' @param maxOuterIterations Maximum number of L-BFGS-B optimization
##'     for outer problem.
##' @param n1qn1nsim Number of function evaluations for n1qn1
##'     optimization.
##'
##' @param eigen A boolean indicating if eigenvectors are calculated
##'     to include a condition number calculation.
##'
##' @param addPosthoc Boolean indicating if posthoc parameters are
##'     added to the table output.
##'
##' @inheritParams RxODE::rxSolve
##' @inheritParams RxODE::rxSymPySetupPred
##'
##' @details
##'
##' Note this uses the R's L-BFGS-B in \code{\link{optim}} for the
##' outer problem and the BFGS \code{\link[n1qn1]{n1qn1}} with that
##' allows restoring the prior individual Hessian (for faster
##' optimization speed).
##'
##' By default FOCEi scales the outer problem parameters to 1.0 for
##' the initial parameter estimates and scales the objective function
##' to 1.0, as suggested by the NAG library
##' (https://www.nag.com/numeric/fl/nagdoc_fl25/html/e04/e04intro.html)
##' and scipy
##' (https://www.scipy-lectures.org/advanced/mathematical_optimization/).
##'
##' However the inner problem is not scaled.  Since most eta estimates
##' start near zero, scaling for these parameters do not make sense.
##'
##' This process of scaling can fix some ill conditioning for the
##' unscaled problem.  The covariance step is performed on the
##' unscaled problem, so the condition number of that matrix may not
##' be reflective of the scaled problem's condition-number.
##'
##' @author Matthew L. Fidler
##'
##' @seealso \code{\link{optim}}
##' @seealso \code{\link[n1qn1]{n1qn1}}
##' @export
foceiControl <- function(sigdig=4,
                         epsilon=NULL, #1e-4,
                         maxInnerIterations=10000,
                         maxOuterIterations=50000,
                         n1qn1nsim=NULL,
                         method = c("liblsoda", "lsoda", "dop853"),
                         transitAbs = NULL, atol = NULL, rtol = NULL,
                         maxstepsOde = 5000L, hmin = 0L, hmax = NULL, hini = 0, maxordn = 12L, maxords = 5L, cores,
                         covsInterpolation = c("linear", "locf", "nocb", "midpoint"),
                         printInner=0L,
                         printOuter=1L,
                         scaleTo=1.0,
                         scaleObjective=1.0,
                         derivEps=c(1.0e-5, 1.0e-5),
                         derivMethod=c("forward", "central"),
                         covDerivMethod=c("central", "forward"),
                         covMethod=c("r,s", "r", "s", ""),
                         lbfgsLmm=40L,
                         lbfgsPgtol=0,
                         lbfgsFactr=NULL, #1e-4 / .Machine$double.eps, ## .Machine$double.eps*x=1e-5
                         eigen=TRUE,
                         addPosthoc=TRUE,
                         diagXform=c("sqrt", "log", "identity"),
                         sumProd=FALSE,
                         optExpression=TRUE,
                         ci=0.95,
                         ## outerOpt=c("lbfgsb", "qnbd"),
                         ..., stiff){
    if (is.null(epsilon)){
        epsilon <- 10 ^ (-sigdig)
    }
    if (is.null(lbfgsFactr)){
        lbfgsFactr <- 10 ^ (-sigdig - 1) / .Machine$double.eps;
    }
    if (is.null(atol)){
        atol <- 0.5 * 10 ^ (-sigdig - 1);
    }
    if (is.null(rtol)){
        rtol <- 0.5 * 10 ^ (-sigdig - 1);
    }
    .xtra <- list(...);
    if (is.null(transitAbs) && !is.null(.xtra$transit_abs)){  # nolint
        transitAbs <- .xtra$transit_abs;  # nolint
    }
    if (missing(covsInterpolation) && !is.null(.xtra$covs_interpolation)){  # nolint
        covsInterpolation <- .xtra$covs_interpolation; # nolint
    }
    if (missing(maxInnerIterations) && !is.null(.xtra$max_iterations)){  # nolint
        maxInnerIterations <- .xtra$max_iterations; # nolint
    }
    if (!missing(stiff) && missing(method)){
        if (RxODE::rxIs(stiff, "logical")){
            if (stiff){
                method <- "lsoda"
                warning("stiff=TRUE has been replaced with method = \"lsoda\".")
            } else {
                method <- "dop853"
                warning("stiff=FALSE has been replaced with method = \"dop853\".")
            }
        }
    }  else {
        method <- match.arg(method);
    }
    .methodIdx <- c("lsoda"=1L, "dop853"=0L, "liblsoda"=2L);
    method <- as.integer(.methodIdx[method]);
    derivMethod <- match.arg(derivMethod);
    .methodIdx <- c("forward"=0L, "central"=1L);
    derivMethod <- as.integer(.methodIdx[derivMethod])
    covDerivMethod <- .methodIdx[match.arg(covDerivMethod)];
    if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1];
    covsInterpolation <- tolower(match.arg(covsInterpolation,
                                           c("linear", "locf", "LOCF", "constant", "nocb", "NOCB", "midpoint")))
    if (covsInterpolation == "constant") covsInterpolation <- "locf";
    covsInterpolation  <- as.integer(which(covsInterpolation == c("linear", "locf", "nocb", "midpoint")) - 1);
    if (missing(cores)){
        cores <- RxODE::rxCores();
    }
    if (missing(n1qn1nsim)){
        n1qn1nsim <- 10 * maxInnerIterations + 1;
    }
    if (length(covMethod) == 1){
        if (covMethod == ""){
            covMethod <- 0L
        }
    }
    if (RxODE::rxIs(covMethod, "character")){
        covMethod <- match.arg(covMethod);
        .covMethodIdx <- c("r,s" = 1L, "r"=2L, "s"=3L);
        covMethod <- .covMethodIdx[match.arg(covMethod)];
    }
    ## outerOpt <- match.arg(outerOpt)
    ## .outerIdx <- c("lbfgsb"=0L, "qnbd"=1L)
    ## outerOpt <- as.integer(.outerIdx[outerOpt]);
    .ret <- list(maxOuterIterations=as.integer(maxOuterIterations),
                 maxInnerIterations=as.integer(maxInnerIterations),
                 method=method,
                 transitAbs=transitAbs,
                 atol=atol,
                 rtol=rtol,
                 maxstepsOde=maxstepsOde,
                 hmin=hmin,
                 hmax=hmax,
                 hini=hini,
                 maxordn=maxordn,
                 maxords=maxords,
                 cores=cores,
                 covsInterpolation=covsInterpolation,
                 n1qn1nsim=as.integer(n1qn1nsim),
                 printInner=as.integer(printInner),
                 printOuter=as.integer(printOuter),
                 lbfgsLmm=as.integer(lbfgsLmm),
                 lbfgsPgtol=as.double(lbfgsPgtol),
                 lbfgsFactr=as.double(lbfgsFactr),
                 scaleTo=scaleTo,
                 epsilon=epsilon,
                 derivEps=derivEps,
                 derivMethod=derivMethod,
                 covDerivMethod=covDerivMethod,
                 covMethod=covMethod,
                 eigen=as.integer(eigen),
                 addPosthoc=as.integer(addPosthoc),
                 diagXform=match.arg(diagXform),
                 sumProd=sumProd,
                 optExpression=optExpression,
                 outerOpt=0L,
                 ci=as.double(ci),
                 sigdig=as.double(sigdig),
                 scaleObjective=as.double(scaleObjective))
    class(.ret) <- "foceiControl"
    return(.ret);
}
.parseOM <- function(OMGA){
    .re = "\\bETA\\[(\\d+)\\]\\b"
    .offset = as.integer(0)
    lapply(1:length(OMGA), function(.k) {
        .s = OMGA[[.k]]
        .f = eval(parse(text=(sprintf("y~%s", deparse(.s[[2]])))))
        .r = unlist(lapply(attr(terms(.f),"variables"), deparse))[-(1:2)]
        .nr = length(.r)

        .ix = grep(.re, .r)
        if(.nr-length(.ix)) stop("invalid OMGA specs")

        .ix = as.integer(sub(.re, "\\1", .r))
        if (any(.ix - (.offset+1:.nr))) stop("invalid OMGA specs")
        .offset <<- .offset + .nr
        eval(.s[[3]])
    })
}

.genOM <- function(s)
{
    .getNR = function(.a) round(sqrt(2 * length(.a) + 0.25) - 0.1)
    .nr = sum(sapply(s, .getNR))
    .mat <- matrix(0, .nr, .nr)
    .offset = as.integer(0)
    j = lapply(1:length(s), function(k) {
        .a = s[[k]]
        .p <- .getNR(.a)
        .starts = row(.mat) > .offset  & col(.mat) > .offset
        .mat[col(.mat) >= row(.mat) & col(.mat) <= .offset+.p & .starts] <<- .a
        .offset <<- .offset+.p
    })
    .a = .mat[col(.mat) >= row(.mat)]
    .mat <- t(.mat)
    .mat[col(.mat) >= row(.mat)] <- .a
    .mat
}
##' FOCEi fit
##'
##' @param data Data to fit
##' @param inits Initialization list
##' @param PKpars Pk Parameters
##' @param model The RxODE model to use
##' @param pred The Prediction function
##' @param err The Error function
##' @param lower Lower bounds
##' @param upper Upper Bounds
##' @param fixed Boolean vector indicating what parameters should be
##'     fixed.
##' @param skipCov Boolean vector indicating what parameters should be
##'     fixed when calculating covariances
##' @param control FOCEi options Control list
##' @param thetaNames Names of the thetas to be used in the final
##'     object
##' @param etaNames Eta names to be used in the final object
##' @param etaMat Eta matrix for initial estimates or final estimates
##'     of the ETAs.
##' @param ... Ignored parameters
##' @return A focei fit object
##' @author Matthew L. Fidler and Wenping Wang
##' @return FOCEi fit object
##' @author Matthew L. Fidler & Wenping Wang
##' @export
##' @examples
##'
##'
##' mypar2 <- function ()
##' {
##'     ka <- exp(THETA[1] + ETA[1])
##'     cl <- exp(THETA[2] + ETA[2])
##'     v  <- exp(THETA[3] + ETA[3])
##' }
##'
##' mod <- RxODE({
##'     d/dt(depot) <- -ka * depot
##'     d/dt(center) <- ka * depot - cl / v * center
##'     cp <- center / v
##' })
##'
##' pred <- function() cp
##'
##' err <- function(){
##'     err <- add(0.1)
##' }
##'
##' inits <- list(THTA=c(0.5, -3.2, -1),
##'               OMGA=list(ETA[1] ~ 1, ETA[2] ~ 2, ETA[3] ~ 1));
##'
##' fit <- foceiFit(theo_sd, inits, mypar2,mod,pred,err)
##'
foceiFit <- function(data,
                     inits,
                     PKpars,
                     model=NULL,
                     pred=NULL,
                     err=NULL,
                     lower= -Inf,
                     upper= Inf,
                     fixed = NULL,
                     skipCov=NULL,
                     control=foceiControl(),
                     thetaNames=NULL,
                     etaNames=NULL,
                     etaMat=NULL,
                     ...,
                     env=NULL){
    .pt <- proc.time();
    loadNamespace("n1qn1");
    if (!RxODE::rxIs(control, "foceiControl")){
        control <- do.call(foceiControl, control);
    }
    if (is.null(env)) .ret <- new.env(parent=emptyenv())
    else .ret <- env;
    .ret$origData <- data;
    .ret$etaNames <- etaNames;
    .ret$lower <- lower;
    .ret$upper <- upper;
    .ret$thetaFixed <- fixed;
    .ret$skipCov <- skipCov;
    .ret$control <- control;
    if(is(model, "RxODE") || is(model, "character")) {
        .ret$ODEmodel <- TRUE
        if (class(pred) != "function"){
            stop("pred must be a function specifying the prediction variables in this model.")
        }
    } else {
        ## Fixme
        .ret$ODEmodel <- TRUE
        model <- constructLinCmt(PKpars);
        pred <- eval(parse(text="function(){return(Central);}"))
    }
    .square <- function(x) x*x
    .ret$diagXformInv = c("sqrt"=".square", "log"="exp", "identity"="identity")[control$diagXform]
    if (is.null(err)){
        err <-eval(parse(text=paste0("function(){err",paste(inits$ERROR[[1]],collapse=""),"}")));
    }
    .ret$model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err, grad=(control$derivMethod == 2L),
                                          pred.minus.dv=TRUE, sum.prod=control$sumProd,
                                          theta.derivs=FALSE, optExpression=control$optExpression);

    .covNames <- .parNames <- RxODE::rxParams(.ret$model$pred.only);
    .covNames <- .covNames[regexpr(rex::rex(start, or("THETA", "ETA"), "[", numbers, "]", end), .covNames) == -1];
    colnames(data) <- sapply(names(data), function(x){
        if (any(x == .covNames)){
            return(x)
        } else {
            return(toupper(x))
        }
    })
    .lhs <- c(names(RxODE::rxInits(.ret$model$pred.only)), RxODE::rxLhs(.ret$model$pred.only))
    if (length(.lhs) > 0){
        .covNames <- .covNames[regexpr(rex::rex(start, or(.lhs), end), .covNames) == -1];
    }
    if (length(.covNames) > 0){
        if (!all(.covNames %in% names(data))){
            message("Model:")
            RxODE::rxCat(model$pred.only)
            message("Needed Covariates:")
            nlmixrPrint(.covNames)
            stop("Not all the covariates are in the dataset.")
        }
        message("Needed Covariates:")
        print(.covNames);
    }
    if (is.null(.ret$model$extra.pars)){
        .nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)))
    } else {
        if (is.null(skipCov)){
            if (is.null(fixed)){
                .tmp <- rep(FALSE, length(inits$THTA))
            } else {
                .tmp <- c(fixed, rep(FALSE, length(inits$THTA) - length(fixed)))
            }
            .ret$skipCov <- c(.tmp,
                              rep(TRUE, length(.ret$model$extra.pars)))
        }
        .nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)),
                  sprintf("ERR[%s]", seq_along(.ret$model$extra.pars)))
    }
    if (!is.null(thetaNames) && (length(inits$THTA) + length(.ret$model$extra.pars)) == length(thetaNames)){
        .nms <- thetaNames;
    }
    .ret$thetaNames <- .nms;
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

    .extraPars <- c();
    if (!is.null(.ret$model$extra.pars)){
        .ret$model$extra.pars <- eval(call(control$diagXform, .ret$model$extra.pars))
        if (length(.ret$model$extra.pars) > 0){
            inits$THTA <- c(inits$THTA, .ret$model$extra.pars);
            .lowerErr <- rep(control$atol * 10, length(.ret$model$extra.pars));
            .upperErr <- rep(Inf, length(.ret$model$extra.pars));
            lower <-c(lower, .lowerErr);
            upper <- c(upper, .upperErr);
        }
    }

    if (is.null(data$ID)) stop('"ID" not found in data')
    if (is.null(data$DV)) stop('"DV" not found in data')
    if (is.null(data$EVID)) data$EVID = 0
    if (is.null(data$AMT)) data$AMT = 0
    ## Make sure they are all double amounts.
    for (.v in c("TIME", "AMT", "DV", .covNames))
        data[[.v]] <- as.double(data[[.v]]);
    .ret$dataSav = data;
    .ds <- data[data$EVID > 0, c("ID", "TIME", "AMT", .covNames)]
    data <- data[data$EVID == 0, c("ID", "TIME", "DV", .covNames)]
    ## keep the covariate names the same as in the model
    .w <- which(!(names(.ret$dataSav) %in% .covNames))
    names(.ret$dataSav)[.w] <- tolower(names(.ret$dataSav[.w]))         #needed in ev

    .lh = .parseOM(inits$OMGA)
    .nlh = sapply(.lh, length)
    .osplt = rep(1:length(.lh), .nlh)
    .lini = list(inits$THTA, unlist(.lh));
    .nlini = sapply(.lini, length)
    .nsplt = rep(1:length(.lini), .nlini)

    .om0 = .genOM(.lh)
    if (length(etaNames) == dim(.om0)[1]){
        .ret$etaNames <- .ret$etaNames
    } else {
        .ret$etaNames <- sprintf("ETA[%d]", seq(1, dim(.om0)[1]))
    }
    .ret$rxInv <- RxODE::rxSymInvCholCreate(mat=.om0, diag.xform=control$diagXform);

    .ret$thetaIni <- inits$THTA

    if (any(.ret$thetaIni == 0)){
        warning("Some of the initial conditions were 0, changing to 0.0001");
        .ret$thetaIni[.ret$thetaIni == 0] <- 0.0001;
    }
    names(.ret$thetaIni) <- sprintf("THETA[%d]", seq_along(.ret$thetaIni))
    .ret$etaMat <- etaMat
    .ret$setupTime <- (proc.time() - .pt)["elapsed"];
    .ret <- RxODE::foceiFitCpp_(.ret);
    .pt <- proc.time();
    .etas <- .ret$ranef
    .thetas <- .ret$fixef
    .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas);
    .preds <- list(ipred=RxODE::rxSolve(.ret$model$inner, .pars$ipred, .ret$dataSav,returnType="data.frame",
                                        atol=.ret$control$atol, rtol=.ret$control$rtol, maxsteps=.ret$control$maxstepsOde,
                                        hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini,
                                        transitAbs = .ret$control$TransitAbs,
                                        maxordn = .ret$control$maxordn, maxords = .ret$control$maxords,
                                        method=.ret$control$method),
                   pred=RxODE::rxSolve(.ret$model$inner, .pars$pred, .ret$dataSav, returnType="data.frame",
                                       atol=.ret$control$atol, rtol=.ret$control$rtol, maxsteps=.ret$control$maxstepsOde,
                                       hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini,
                                       transitAbs = .ret$control$transitAbs, maxordn = .ret$control$maxordn,
                                       maxords = .ret$control$maxords, method=.ret$control$method));
    .lst <- .Call(`_nlmixr_nlmixrResid`, .preds, .ret$omega, data$DV, .etas, .pars$eta.lst);
    .df <- RxODE::rxSolve(.ret$model$ebe, .pars$ipred,.ret$dataSav,returnType="data.frame",
                          hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini, transitAbs = .ret$control$transitAbs,
                          maxordn = .ret$control$maxordn, maxords = .ret$control$maxords,
                          Method=.ret$control$method)[, -(1:2)];
    .df <- .df[, !(names(.df) %in% c("nlmixr_pred", names(.etas), names(.thetas)))]
    if (any(names(.df) %in% names(.lst[[1]]))){
        warning("Calculated residuals like IPRED are masked by nlmixr calculated values");
        .df <- .df[, !(names(.df) %in% names(.lst[[1]]))];
    }
    if (any(names(.df) %in% names(.lst[[3]]))){
        .df <- .df[, !(names(.df) %in% names(.lst[[3]]))];
    }
    .lst[[4]] <- .df;
    .ret$shrink <- .lst[[2]];
    .df <- cbind(as.data.frame(data), .lst[[1]], .lst[[3]], .lst[[4]]);
    .ret$tableTime <- (proc.time() - .pt)["elapsed"];
    .ret$time <- data.frame(.ret$time, table=.ret$tableTime);
    .isDplyr <- requireNamespace("dplyr", quietly = TRUE);
    if (!.isDplyr){
        .isDataTable <- requireNamespace("data.table", quietly = TRUE);
        if (.isDataTable){
            .df <- data.table::data.table(.df)
        }
    } else {
        .df <- dplyr::as.tbl(.df);
    }
    .cls <- class(.df);
    .cls <- c("foceiFitData", "foceiFitCore", .cls)
    attr(.cls, ".foceiEnv") <- .ret;
    class(.df) <- .cls;
    return(.df)
}

##' @export
`$.foceiFitData` <-  function(obj, arg, exact = FALSE){
    .ret <- obj[[arg]]
    if (is.null(.ret)){
        .cls <- class(obj);
        .env <- attr(.cls, ".foceiEnv");
        if (arg == "omega.R") arg <- "omegaR"
        if (exists(arg, envir=.env)){
            .ret <- get(arg, envir=.env);
        }
    }
    return(.ret)
}

##'@export
logLik.foceiFitCore <- function(object, ...){
    object$logLik
}

##'@export
nobs.foceiFitCore <- function(object, ...){
    object$nobs
}

##'@export
vcov.foceiFitCore <- function(object, ...){
    object$cov
}

##'@export
getData.foceiFitCore <- function(object){
    object$origData
}

##'@export
ranef.foceiFitCore <- function(object, ...){
    object$ranef;
}

##'@export
fixef.foceiFitCore <- function(object, ...){
    object$fixef;
}

##'@export
print.foceiFitCore <- function(x, ...){
    .parent <- parent.frame(2);
    .bound <- do.call("c", lapply(ls(.parent), function(.cur){
                               if (identical(.parent[[.cur]], x)){
                                   return(.cur)
                               }
                               return(NULL);
                           }))
    message(cli::rule(paste0(crayon::bold$blue("nlmix"), crayon::bold$red("r"), " ", crayon::bold$yellow(x$method), " fit ",
                             x$extra)))
    print(x$objDf)
    message(paste0("\n", cli::rule(paste0(crayon::bold("Time"), " (sec; ", crayon::yellow(.bound), crayon::bold$blue("$time"), "):"))));
    print(x$time)
    message(paste0("\n", cli::rule(paste0(crayon::bold("Population Parameters"), " (", crayon::yellow(.bound), crayon::bold$blue("$fixedDf"), " & ", crayon::yellow(.bound), crayon::bold$blue("$fixedDfSig"), "):"))));
    print(x$fixedDfSig)
    ################################################################################
    ## Correlations
    .tmp <- x$omega
    diag(.tmp) <- 0;
    if (all(.tmp == 0)){
        message("\n  No correlations in between subject variability (BSV) matrix")
    } else {
        message("\n  Correlations in between subject variability (BSV) matrix:")
        .rs <- fit$omegaR
        .lt <- lower.tri(.rs);
        .dn1 <- dimnames(fit$omegaR)[[2]]
        .nms <- apply(which(.lt,arr.ind=TRUE),1,function(x){sprintf("R(%s)",paste(.dn1[x],collapse=", "))});
        .lt <- structure(.rs[.lt], .Names=.nms)
        .lt <- .lt[.lt != 0]
        .digs <- 3;
        lts <- sapply(.lt, function(x){
            x <- abs(.lt);
            .ret <- "<"
            if (x > 0.7){
                .ret <- ">" ## Strong
            } else if (x > 0.3){
                .ret <- "=" ## Moderate
            }
            return(.ret)
        })
        .nms <- names(.lt);
        .lt <- sprintf("%s%s", formatC(signif(.lt, digits=.digs),digits=.digs,format="fg", flag="#"), lts)
        names(.lt) <- .nms;
        .lt <- gsub(rex::rex("\""), "", paste0("    ", R.utils::captureOutput(print(.lt))));
        if (crayon::has_color()){
            .lt <- gsub(rex::rex(capture(.regNum), ">"), "\033[1m\033[31m\\1 \033[39m\033[22m", .lt, perl=TRUE)
            .lt <- gsub(rex::rex(capture(.regNum), "="), "\033[1m\033[32m\\1 \033[39m\033[22m", .lt, perl=TRUE)
            .lt <- gsub(rex::rex(capture(.regNum), "<"), "\\1 ", .lt, perl=TRUE)
        } else {
            .lt <- gsub(rex::rex(capture(.regNum), or(">", "=", "<")), "\\1 ", .lt, perl=TRUE)
        }
        message(paste(.lt, collapse="\n"), "\n")
    }
    message(paste0("  Full BSV covariance (", crayon::yellow(.bound), crayon::bold$blue("$omega"), ") or correlation (", crayon::yellow(.bound), crayon::bold$blue("$omegaR"), "; diagonals=SDs)"));
    message(paste0("  Distribution stats (mean/skewness/kurtosis/p-value) available in ",
                   crayon::yellow(.bound), crayon::bold$blue("$shrink")));
    if (RxODE::rxIs(x, "foceiFitData")){
        message(paste0("\n", cli::rule(paste0(crayon::bold("Fit Data"), " (object", ifelse(.bound == "", "", " "),
                                              crayon::yellow(.bound),
                                              " is a modified ", crayon::blue("data.frame"), "):"))))
        if (RxODE::rxIs(x, "tbl") || RxODE::rxIs(x, "data.table")){
            NextMethod(x)
        } else{
            print(head(x));
        }
    }
}
