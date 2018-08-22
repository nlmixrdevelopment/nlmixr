.regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                          group(any_of("0":"9"), ".", some_of("0":"9"))),
                       maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9"))));
.regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"));
.regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
.regNum <- rex::rex(maybe("-"), or(.regDecimalint, .regFloat1, .regFloat2))

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

##' Control Options for FOCEi
##'
##' @param sigdig Optimization significant digits. This controls:
##'
##'  Defaults for optimization and ODE solving
##'
##'  The tolerance of the inner and outer optimization is \code{10^-sigdig}
##'
##'  The tolerance of the ODE solvers is \code{10^(-sigdig-1)}
##'
##'  The significant figures that some tables are rounded to.
##'
##' @param epsilon Precision of estimate for n1qn1 optimization.
##'
##' @param maxstepsOde Maximum number of steps for ODE solver.
##'
##' @param printInner Integer representing when the inner step is
##'     printed. By default this is 0 or do not print.  1 is print
##'     every function evaluation, 5 is print every 5 evaluations.
##'
##' @param print Integer representing when the outer step is
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
##'      \code{scaledObj = currentObj / |initialObj| * scaleObjective}
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
##'         \code{h = abs(x)*derivEps[1] + derivEps[2]}
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
##'
##'  "\code{r,s}" Uses the sandwich matrix to calculate the covariance, that is: \eqn{R^-1 \times S \times R^-1}
##'
##'  "\code{r}" Uses the Hessian matrix to calculate the covariance as \eqn{2\times R^-1}
##'
##'  "\code{s}" Uses the crossproduct matrix to calculate the covariance as \eqn{4\times S^-1}
##'
##'  "" Does not calculate the covariance step.
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
##'  \code{sqrt} Estimates the sqrt of the diagonal elements of \code{chol(inv(omega))}.  This is the default method.
##'
##'  \code{log} Estimates the log of the diagonal elements of \code{chol(inv(omega))}
##'
##'  \code{identity} Estimates the diagonal elements without any transformations
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
##' @param ... Ignored parameters
##'
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
##' @details
##'
##' Note this uses the R's L-BFGS-B in \code{\link{optim}} for the
##' outer problem and the BFGS \code{\link[n1qn1]{n1qn1}} with that
##' allows restoring the prior individual Hessian (for faster
##' optimization speed).
##'
##' By default FOCEi scales the outer problem parameters to 1.0 for
##' the initial parameter estimates and scales the objective function
##' to 1.0, as suggested by the
##' \href{https://www.nag.com/numeric/fl/nagdoc_fl25/html/e04/e04intro.html}{NAG library}
##' and \href{scipy}{https://www.scipy-lectures.org/advanced/mathematical_optimization/}.
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
                         maxInnerIterations=1000,
                         maxOuterIterations=5000,
                         n1qn1nsim=NULL,
                         method = c("liblsoda", "lsoda", "dop853"),
                         transitAbs = NULL, atol = NULL, rtol = NULL,
                         maxstepsOde = 5000L, hmin = 0L, hmax = NULL, hini = 0, maxordn = 12L, maxords = 5L, cores,
                         covsInterpolation = c("linear", "locf", "nocb", "midpoint"),
                         printInner=0L,
                         print=1L,
                         printNcol=8L,
                         scaleTo=1.0,
                         scaleObjective=1.0,
                         derivEps=c(1.0e-5, 1.0e-5),
                         derivMethod=c("forward", "central"),
                         covDerivMethod=c("forward", "central"),
                         covMethod=c("r,s", "r", "s", ""),
                         lbfgsLmm=50L,
                         lbfgsPgtol=0,
                         lbfgsFactr=NULL, #1e-4 / .Machine$double.eps, ## .Machine$double.eps*x=1e-5
                         eigen=TRUE,
                         addPosthoc=TRUE,
                         diagXform=c("sqrt", "log", "identity"),
                         sumProd=FALSE,
                         optExpression=TRUE,
                         ci=0.95,
                         ## outerOpt=c("lbfgsb", "qnbd"),
                         useColor=crayon::has_color(),
                         boundTol=NULL,
                         calcTables=TRUE,
                         ..., stiff){
    if (is.null(boundTol)){
        boundTol <- 5 * 10 ^ (-sigdig + 1)
    }
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
                 print=as.integer(print),
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
                 scaleObjective=as.double(scaleObjective),
                 useColor=as.integer(useColor),
                 boundTol=as.double(boundTol),
                 calcTables=calcTables,
                 printNcol=as.integer(printNcol));
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
##' FOCEi fit
##'
##' @param data Data to fit; Needs to be RxODE compatible and have \code{DV},
##'     \code{AMT}, \code{EVID} in the dataset.
##' @param inits Initialization list
##' @param PKpars Pk Parameters function
##' @param model The RxODE model to use
##' @param pred The Prediction function
##' @param err The Error function
##' @param lower Lower bounds
##' @param upper Upper Bounds
##' @param fixed Boolean vector indicating what parameters should be
##'     fixed.
##' @param skipCov Boolean vector indicating what parameters should be
##'     fixed when calculating covariances
##' @param control FOCEi options Control list.  See
##'     \code{\link{foceiControl}}
##' @param thetaNames Names of the thetas to be used in the final
##'     object.
##' @param etaNames Eta names to be used in the final object.
##' @param etaMat Eta matrix for initial estimates or final estimates
##'     of the ETAs.
##' @param ... Ignored parameters
##' @return A focei fit object
##' @author Matthew L. Fidler and Wenping Wang
##' @return FOCEi fit object
##' @export
##' @examples
##' \dontrun{
##' ## Comparison to Wang2007 objective functions
##'
##' mypar2 = function ()
##' {
##'     k = theta[1] * exp(eta[1]);
##' }
##'
##' mod <- RxODE({
##'     ipre = 10 * exp(-k * t)
##' })
##' pred <- function() ipre
##'
##' errProp <- function(){
##'   return(prop(0.1))
##' }
##'
##' inits <- list(THTA=c(0.5),
##'               OMGA=list(ETA[1] ~ 0.04));
##'
##' w7 <- Wang2007
##'
##' w7$DV <- w7$Y
##' w7$EVID <- 0
##' w7$AMT <- 0
##'
##' ## Wang2007 prop error 39.458 for NONMEM, nlmixr matches.
##' fitP <- foceiFit(w7, inits, mypar2,mod,pred,errProp,
##'      control=foceiControl(maxOuterIterations=0,covMethod=""))
##'
##' ## Note if you have the etas you can evaluate the likelihood
##' ## of an arbitrary model.  It doesn't have to be solved by
##' ## FOCEi
##'
##' etaMat <- matrix(fitP$eta[,-1])
##'
##' fitP2 <- foceiFit(w7, inits, mypar2,mod,pred,errProp, etaMat=etaMat,
##'       control=foceiControl(maxOuterIterations=0,maxInnerIterations=0,
##'       covMethod=""))
##'
##'
##' errAdd <- function(){
##'   return(add(0.1))
##' }
##'
##' ## Wang2007 add error of -2.059 for NONMEM FOCE=NONMEM FOCEi;
##' ## nlmixr matches.
##' fitA <- foceiFit(w7, inits, mypar2,mod,pred,errAdd,
##'      control=foceiControl(maxOuterIterations=0,covMethod=""))
##'
##' ## Extending Wang2007 to add+prop with same dataset
##' errAddProp <- function(){
##'   return(add(0.1) + prop(0.1));
##' }
##'
##' fitAP <- foceiFit(w7, inits, mypar2,mod,pred,errAddProp,
##'      control=foceiControl(maxOuterIterations=0,covMethod=""))
##'
##' ## Checking lognormal
##'
##' errLogn <- function(){
##'    return(lnorm(0.1));
##' }
##'
##' ## First run the fit with the nlmixr lnorm error
##'
##' fitLN <- foceiFit(w7, inits, mypar2,mod,pred,errLogn,
##'      control=foceiControl(maxOuterIterations=0,covMethod=""))
##'
##'
##' ## Next run on the log-transformed space
##' w72 <- w7; w72$DV <- log(w72$DV)
##'
##' predL <- function() log(ipre)
##'
##' fitLN2 <- foceiFit(w72, inits, mypar2,mod,predL,errAdd,
##'      control=foceiControl(maxOuterIterations=0,covMethod=""))
##'
##' ## Correct the fitLN2's objective function to be on the normal scale
##' print(fitLN2$objective + 2*sum(w72$DV))
##'
##' ## Note the objective function of the lognormal error is on the normal scale.
##' print(fitLN$objective)
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
##'     return(add(0.1))
##' }
##'
##' inits <- list(THTA=c(0.5, -3.2, -1),
##'               OMGA=list(ETA[1] ~ 1, ETA[2] ~ 2, ETA[3] ~ 1));
##'
##' ## Remove 0 concentrations (should be lloq)
##'
##' d <- theo_sd[theo_sd$EVID==0 & theo_sd$DV>0 | theo_sd$EVID>0,];
##'
##' fit1 <- foceiFit(d, inits, mypar2,mod,pred,err)
##'
##' ## you can also fit lognormal data with the objective function on the same scale
##'
##' errl <- function(){
##'     return(lnorm(0.1))
##' }
##'
##' fit2 <- foceiFit(d, inits, mypar2,mod,pred,errl)
##'
##' ## You can also use the standard nlmixr functions to run FOCEi
##'
##' library(data.table);
##' datr <- Infusion_1CPT;
##' datr$EVID<-ifelse(datr$EVID==1,10101,datr$EVID)
##' datr<-data.table(datr)
##' datr<-datr[EVID!=2]
##' datro<-copy(datr)
##' datIV<-datr[AMT>0][,TIME:=TIME+AMT/RATE][,AMT:=-1*AMT]
##' datr<-rbind(datr,datIV)
##'
##' one.compartment.IV.model <- function(){
##'   ini({ # Where initial conditions/variables are specified
##'     # '<-' or '=' defines population parameters
##'     # Simple numeric expressions are supported
##'     lCl <- 1.6      #log Cl (L/hr)
##'     lVc <- log(90)  #log V (L)
##'     # Bounds may be specified by c(lower, est, upper), like NONMEM:
##'     # Residuals errors are assumed to be population parameters
##'     prop.err <- c(0, 0.2, 1)
##'     # Between subject variability estimates are specified by '~'
##'     # Semicolons are optional
##'     eta.Vc ~ 0.1   #IIV V
##'     eta.Cl ~ 0.1; #IIV Cl
##'   })
##'   model({ # Where the model is specified
##'     # The model uses the ini-defined variable names
##'     Vc <- exp(lVc + eta.Vc)
##'     Cl <- exp(lCl + eta.Cl)
##'     # RxODE-style differential equations are supported
##'     d / dt(centr) = -(Cl / Vc) * centr;
##'     ## Concentration is calculated
##'     cp = centr / Vc;
##'     # And is assumed to follow proportional error estimated by prop.err
##'     cp ~ prop(prop.err)
##'    })}
##'
##' fitIVp <- nlmixr(one.compartment.IV.model, datr, "focei");
##'
##' ## You can also use the Cox-Box Transform of both sides with
##' ## proportional error (Donse 2016)
##'
##' one.compartment.IV.model <- function(){
##' ini({ # Where initial conditions/variables are specified
##'     ## '<-' or '=' defines population parameters
##'     ## Simple numeric expressions are supported
##'     lCl <- 1.6      #log Cl (L/hr)
##'     lVc <- log(90)  #log V (L)
##'     ## Bounds may be specified by c(lower, est, upper), like NONMEM:
##'     ## Residuals errors are assumed to be population parameters
##'     prop.err <- c(0, 0.2, 1)
##'     zeta <- c(0.1, 1, 10)
##'     add.err <- c(0, 0.001)
##'     lambda <- c(-2, 1, 2)
##'     ## Between subject variability estimates are specified by '~'
##'     ## Semicolons are optional
##'     eta.Vc ~ 0.1   #IIV V
##'     eta.Cl ~ 0.1; #IIV Cl
##' })
##' model({ ## Where the model is specified
##'     ## The model uses the ini-defined variable names
##'     Vc <- exp(lVc + eta.Vc)
##'     Cl <- exp(lCl + eta.Cl)
##'     ## RxODE-style differential equations are supported
##'     d / dt(centr) = -(Cl / Vc) * centr;
##'     ## Concentration is calculated
##'     cp = centr / Vc;
##'     ## And is assumed to follow proportional error estimated by prop.err
##'     cp ~ pow(prop.err, zeta) + add(add.err) + tbs(lambda)
##' })}
##'
##' fitIVtbs <- nlmixr(one.compartment.IV.model, datr, "focei")
##'
##' ## If you want to use a variance normalizing distribution with
##' ## negative/positive data you can use the Yeo-Johnson transformation
##' ## as well.  This is implemented by the tbsYj(lambda) function.
##' one.compartment.IV.model <- function(){
##' ini({ # Where initial conditions/variables are specified
##'     ## '<-' or '=' defines population parameters
##'     ## Simple numeric expressions are supported
##'     lCl <- 1.6      #log Cl (L/hr)
##'     lVc <- log(90)  #log V (L)
##'     ## Bounds may be specified by c(lower, est, upper), like NONMEM:
##'     ## Residuals errors are assumed to be population parameters
##'     prop.err <- c(0, 0.2, 1)
##'     delta <- c(0.1, 1, 10)
##'     add.err <- c(0, 0.001)
##'     lambda <- c(-2, 1, 2)
##'     ## Between subject variability estimates are specified by '~'
##'     ## Semicolons are optional
##'     eta.Vc ~ 0.1   #IIV V
##'     eta.Cl ~ 0.1; #IIV Cl
##' })
##' model({ ## Where the model is specified
##'     ## The model uses the ini-defined variable names
##'     Vc <- exp(lVc + eta.Vc)
##'     Cl <- exp(lCl + eta.Cl)
##'     ## RxODE-style differential equations are supported
##'     d / dt(centr) = -(Cl / Vc) * centr;
##'     ## Concentration is calculated
##'     cp = centr / Vc;
##'     ## And is assumed to follow proportional error estimated by prop.err
##'     cp ~ pow(prop.err, delta) + add(add.err) + tbsYj(lambda)
##' })}
##'
##' fitIVyj <- nlmixr(one.compartment.IV.model, datr, "focei")
##'
##' }
foceiFit <- function(data, ...){
    UseMethod("foceiFit")
}
##'@rdname foceiFit
##'@export
foceiFit.data.frame <- function(...){
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    return(.collectWarnings(do.call(foceiFit.data.frame0, call, envir=parent.frame(1))))
}

.updateParFixed <- function(.ret){
    if (exists("uif", envir=.ret) & exists("shrink", envir=.ret)){
        .uif <- .ret$uif;
        .lab <- paste(.uif$ini$label[!is.na(.uif$ini$ntheta)]);
        .lab[.lab == "NA"] <- "";
        .lab <- gsub(" *$", "", gsub("^ *", "", .lab));
        .muRef <- unlist(.uif$mu.ref);
        if (length(.muRef) > 0){
            .nMuRef <- names(.muRef)
            .ome <- .ret$omega
            .muRef <- structure(as.vector(.nMuRef), .Names=as.vector(.muRef));
            .logEta <- .uif$log.eta;
            .digs <- .ret$control$sigdig;
            .cvOnly <- TRUE;
            .sdOnly <- TRUE;
            .cvp <- sapply(row.names(.ret$popDfSig), function(x){
                .y <- .muRef[x];
                if (is.na(.y)) return(" ");
                .v <- .ome[.y, .y];
                if (any(.y == .logEta)){
                    .sdOnly <<- FALSE;
                    sprintf("%s%%", formatC(signif(sqrt(exp(.v) - 1) * 100, digits=.digs),
                                            digits=.digs, format="fg", flag="#"));
                } else {
                    .cvOnly <<- FALSE;
                    sprintf("%s", formatC(signif(sqrt(.v),digits=.digs),
                                          digits=.digs, format="fg", flag="#"));
                }
            })

            .shrink <- .ret$shrink;
            .errs <- as.data.frame(.uif$ini);
            .errs <- paste(.errs[which(!is.na(.errs$err)), "name"]);
            .sh <- sapply(row.names(.ret$popDfSig), function(x){
                .y <- .muRef[x];
                if (is.na(.y)) {
                    if (any(x == .errs)){
                        .v <- .shrink[7, "IWRES"];
                        if (length(.v) != 0) return(" ")
                        if (is.null(.v)) return(" ")
                        if (is.na(.v)){
                            return(" ")
                        }
                    } else {
                        return(" ")
                    }
                } else {
                    .v <- .shrink[7, .y];
                }
                if (length(.v) != 1) return(" ");
                .t <- ">"
                if (.v < 0){
                } else  if (.v < 20){
                    .t <- "<"
                } else if (.v < 30){
                    .t <- "="
                }
                sprintf("%s%%%s", formatC(signif(.v, digits=.digs),
                                          digits=.digs, format="fg", flag="#"), .t);
            })
            .ret$parFixed <- data.frame(.ret$popDfSig, "BSD"=.cvp, "Shrink(SD)%"=.sh, check.names=FALSE);
            names(.ret$parFixed)[ifelse(exists("cov", envir=.ret), 5, 3)] <- ifelse(.sdOnly, "BSV(SD)", ifelse(.cvOnly, "BSV(CV%)", "BSV(CV% or SD)"))
            if (!all(.lab == "")){
                .ret$parFixed <- data.frame(Parameter=.lab, .ret$parFixed, check.names=FALSE)
            }
        } else {
            .ret$parFixed <- .ret$popDfSig
        }
    } else {
        .ret$parFixed <- .ret$popDfSig
    }
    class(.ret$parFixed) <- c("nlmixrParFixed", "data.frame");
}

foceiFit.data.frame0 <- function(data,
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
    on.exit({RxODE::rxSolveFree()});
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
    .covNames <- .parNames <- c();
    if (!exists("noLik", envir=.ret)){
        .ret$model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err, grad=(control$derivMethod == 2L),
                                              pred.minus.dv=TRUE, sum.prod=control$sumProd,
                                              theta.derivs=FALSE, optExpression=control$optExpression, run.internal=TRUE);

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
                RxODE::rxCat(.ret$model$pred.only)
                message("Needed Covariates:")
                nlmixrPrint(.covNames)
                stop("Not all the covariates are in the dataset.")
            }
            message("Needed Covariates:")
            print(.covNames);
        }
        .extraPars <- .ret$model$extra.pars
    } else {
        if (.ret$noLik){
            .ret$model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err, grad=(control$derivMethod == 2L),
                                                  pred.minus.dv=TRUE, sum.prod=control$sumProd,
                                                  theta.derivs=FALSE, optExpression=control$optExpression, run.internal=TRUE,
                                                  only.numeric=TRUE);
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
                    RxODE::rxCat(.ret$model$pred.only)
                    message("Needed Covariates:")
                    nlmixrPrint(.covNames)
                    stop("Not all the covariates are in the dataset.")
                }
                message("Needed Covariates:")
                print(.covNames);
            }
            .extraPars <- .ret$model$extra.pars
        } else {
            .extraPars <- NULL;
        }
    }
    .ret$skipCov <- skipCov;
    if (is.null(skipCov)){
        if (is.null(fixed)){
            .tmp <- rep(FALSE, length(inits$THTA))
        } else {
            .tmp <- c(fixed, rep(FALSE, length(inits$THTA) - length(fixed)))
        }
        if (exists("uif", envir=.ret)){
            .uifErr <- .ret$uif$ini$err[!is.na(.ret$uif$ini$ntheta)];
            .uifErr <- sapply(.uifErr, function(x){
                if (is.na(x)) return(FALSE);
                return(!any(x == c("pow2", "tbs", "tbsYj")))
            })
            .tmp <- (.tmp | .uifErr);
        }
        .ret$skipCov <- c(.tmp,
                          rep(TRUE, length(.extraPars)))
    }
    if (is.null(.extraPars)){
        .nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)))
    } else {
        .nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)),
                  sprintf("ERR[%s]", seq_along(.extraPars)))
    }
    if (!is.null(thetaNames) && (length(inits$THTA) + length(.extraPars)) == length(thetaNames)){
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

    if (!is.null(.extraPars)){
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

    if (any(.ret$thetaIni == 0 && control$scaleTo > 0)){
        warning("Some of the initial conditions were 0, changing to 0.0001");
        .ret$thetaIni[.ret$thetaIni == 0] <- 0.0001;
    }
    names(.ret$thetaIni) <- sprintf("THETA[%d]", seq_along(.ret$thetaIni))
    .ret$etaMat <- etaMat
    .ret$setupTime <- (proc.time() - .pt)["elapsed"];
    if (exists("uif", envir=.ret)){
        .uif <- .ret$uif
        .ret$logThetas <- as.integer(which(setNames(sapply(.uif$focei.names,function(x)any(x==.uif$log.theta)),NULL)))
    }
    if (exists("noLik", envir=.ret)){
        if (!.ret$noLik){
            .ret$.params <- c(sprintf("THETA[%d]", seq_along(.ret$thetaIni)),
                              sprintf("ETA[%d]", seq(1, dim(.om0)[1])));
            .ret$.thetan <- length(.ret$thetaIni);
            .ret$nobs <- sum(data$EVID == 0)
        }
    }
    .ret <- RxODE::foceiFitCpp_(.ret);
    if (!control$calcTables){
        return(.ret);
    }
    if (exists("noLik", envir=.ret)){
        if (.ret$noLik){
            message("Calculating residuals/tables")
            .pt <- proc.time();
            .etas <- .ret$ranef
            .thetas <- .ret$fixef
            .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas);
            .preds <- list(ipred=RxODE::rxSolve(.ret$model$pred.only, .pars$ipred, .ret$dataSav,returnType="data.frame.TBS",
                                                atol=.ret$control$atol, rtol=.ret$control$rtol, maxsteps=.ret$control$maxstepsOde,
                                                hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini,
                                                transitAbs = .ret$control$TransitAbs,
                                                maxordn = .ret$control$maxordn, maxords = .ret$control$maxords,
                                                method=.ret$control$method),
                           pred=RxODE::rxSolve(.ret$model$pred.only, .pars$pred, .ret$dataSav, returnType="data.frame",
                                               atol=.ret$control$atol, rtol=.ret$control$rtol, maxsteps=.ret$control$maxstepsOde,
                                               hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini,
                                               transitAbs = .ret$control$transitAbs, maxordn = .ret$control$maxordn,
                                               maxords = .ret$control$maxords, method=.ret$control$method),
                           cwres=FALSE);
        } else {
            .pt <- proc.time();
            .etas <- .ret$ranef
            .thetas <- .ret$fixef
            .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas);
            .ret$shrink <- .Call(`_nlmixr_nlmixrShrink`, .ret$omega, .etas, .pars$eta.lst[-(dim(.ret$omega)[1] + 1)]);
            .updateParFixed(.ret);
            return(.ret);
        }
    } else{
        message("Calculating residuals/tables")
        .pt <- proc.time();
        .etas <- .ret$ranef
        .thetas <- .ret$fixef
        .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas);
        .preds <- list(ipred=RxODE::rxSolve(.ret$model$inner, .pars$ipred, .ret$dataSav,returnType="data.frame.TBS",
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
    }
    .lst <- .Call(`_nlmixr_nlmixrResid`, .preds, .ret$omega, data$DV, .preds$ipred$rxLambda, .preds$ipred$rxYj, .etas, .pars$eta.lst);
    if (is.null(.preds$cwres)){
        .df <- RxODE::rxSolve(.ret$model$pred.only, .pars$ipred,.ret$dataSav,returnType="data.frame",
                              hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini, transitAbs = .ret$control$transitAbs,
                              maxordn = .ret$control$maxordn, maxords = .ret$control$maxords,
                              Method=.ret$control$method)[, -(1:4)];
    } else {
        .df <- .preds$ipred[, -c(1:4, length(names(.preds$ipred)) - 0:1), drop = FALSE];
    }
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
    .updateParFixed(.ret);
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
    .cls <- c("nlmixrFitData", "nlmixrFitCore", .cls)
    attr(.cls, ".foceiEnv") <- .ret;
    class(.df) <- .cls;
    message("done.")
    return(.df)
}

##' @rdname foceiFit
##' @export
focei.fit <- foceiFit

##' @export
`$.nlmixrFitCore` <- function(obj, arg, exact = FALSE){
    .env <- obj;
    if (arg == "par.hist") arg <- "parHist"
    if (arg == "par.hist.stacked") arg <- "parHistStacked"
    if (arg == "omega.R") arg <- "omegaR"
    if (arg == "par.fixed") arg <- "parFixed"
    if (arg == "eta") arg <- "ranef"
    if (arg == "theta") arg <- "fixef"
    if (arg == "varFix") arg <- "cov"
    if (arg == "thetaMat") arg <- "cov"
    if (exists(arg, envir=.env)){
        return(get(arg, envir=.env));
    }
    if (arg == "env"){
        return(.env)
    }
    if (exists("uif", .env)){
        .uif <- .env$uif;
        if (arg == "modelName") arg <- "model.name"
        if (arg == "dataName") arg <- "data.name"
        .ret <- `$.nlmixrUI`(.uif, arg);
        if (!is.null(.ret)) return(.ret)
        .env2  <- `$.nlmixrUI`(.uif, "env");
        if (exists(arg, envir=.env2)){
            return(get(arg, envir=.env2))
        }
    }
}

##' @export
`$.nlmixrFitData` <-  function(obj, arg, exact = FALSE){
    .ret <- obj[[arg]]
    if (is.null(.ret)){
        .cls <- class(obj);
        .env <- attr(.cls, ".foceiEnv");
        .ret <- `$.nlmixrFitCore`(.env, arg, exact);
        if (is.null(.ret)){
            if (arg == "simInfo"){
                return(.simInfo(obj))
            }
        }
    }
    return(.ret)
}

##' Return composite nlme/focei to nlme
##'
##' @param x nlme/focei object from common UI.
##' @return nlme object or NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
as.nlme <- function(object, ...){
    .ret <- object$nlme
    if (is.null(.ret)) stop("Cannot convert to nlme.");
    return(.ret)
}
##' Return composite saem/focei to saem
##'
##' @param x saem/focei object from common UI.
##' @return saem object or NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
as.saem <- function(x){
    .ret <- x$saem
    if (is.null(.ret)) stop("Cannot convert to saem.");
    return(.ret)
}

##' @importFrom nlme VarCorr
##' @export
VarCorr.nlmixrFitCore <- function(x, sigma = 1, ...){
    VarCorr(as.nlme(x), sigma=sigma, ...)
}

##' @export
str.nlmixrFitData <- function(object, ...){
    NextMethod(object);
    .env <- object$env
    ## cat(" $ par.hist         : Parameter history (if available)\n")
    ## cat(" $ par.hist.stacked : Parameter history in stacked form for easy plotting (if available)\n")
    cat(" $ omega            : Omega matrix\n")
    cat(" $ omegaR           : Omega Correlation matrix\n")
    cat(" $ shrink           : Shrinkage table, includes skewness, kurtosis, and eta p-values\n")
    cat(" $ parFixed         : Fixed Effect Parameter Table\n")
    cat(" $ theta            : Fixed Parameter Estimates\n")
    cat(" $ eta              : Individual Parameter Estimates\n")
    cat(" $ seed             : Seed (if applicable)\n");
    if (exists("uif", envir=object$env)){
        cat(" $ meta             : Model meta information environment\n");
        cat(" $ modelName        : Model name (from R function)\n");
        cat(" $ dataName         : Name of R data input\n");
        cat(" $ simInfo          : RxODE list for simulation\n");
    }

}


##' Extract residuals from the FOCEI fit
##'
##' @param object focei.fit object
##' @param ... Additional arguments
##' @param type Residuals type fitted.
##' @return residuals
##' @author Matthew L. Fidler
##' @export
residuals.nlmixrFitData <- function(object, ..., type=c("ires", "res", "iwres", "wres", "cwres", "cpred", "cres")){
    return(object[, toupper(match.arg(type))]);
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
plot.nlmixrFitData <- function(x, ...) {
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


    p2 <- ggplot(dat, aes(x=IPRED, y=IWRES)) +
        geom_point() +
        geom_abline(slope=0, intercept=0, col="red")
    print(p2)

    p2 <- ggplot(dat, aes(x=TIME, y=IWRES)) +
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

##'@export
logLik.nlmixrFitCore <- function(object, ...){
    object$logLik
}

##'@export
nobs.nlmixrFitCore <- function(object, ...){
    object$nobs
}

##'@export
vcov.nlmixrFitCore <- function(object, ...){
    object$cov
}

##'@export
getData.nlmixrFitCore <- function(object){
    object$origData
}

##'@export
ranef.nlmixrFitCore <- function(object, ...){
    object$ranef;
}

##'@export
fixef.nlmixrFitCore <- function(object, ...){
    object$fixef;
}

##'@export
print.nlmixrFitCore <- function(x, ...){
    .parent <- parent.frame(2);
    .bound <- do.call("c", lapply(ls(.parent), function(.cur){
                               if (identical(.parent[[.cur]], x)){
                                   return(.cur)
                               }
                               return(NULL);
                           }))
    if (length(.bound) == 0){
        .bound <- ""
    } else if (length(.bound) > 2){
        .bound <- .bound[order(sapply(.bound, nchar))];
        .bound <- .bound[1];
    }
    .posthoc <- (x$control$maxOuterIterations == 0L & x$control$maxInnerIterations > 0L)
    .posthoc <- ifelse(.posthoc, paste0(crayon::bold(" posthoc"), " estimation"), " fit");
    message(cli::rule(paste0(crayon::bold$blue("nlmix"), crayon::bold$red("r"), " ", crayon::bold$yellow(x$method),.posthoc,
                             x$extra)))
    print(x$objDf)
    message(paste0("\n", cli::rule(paste0(crayon::bold("Time"), " (sec; ", crayon::yellow(.bound), crayon::bold$blue("$time"), "):"))));
    print(x$time)
    message(paste0("\n", cli::rule(paste0(crayon::bold("Population Parameters"), " (", crayon::yellow(.bound), crayon::bold$blue("$parFixed"), "):"))));
    .pf <- R.utils::captureOutput(print(x$parFixed))
    if (crayon::has_color()){
        .pf <- gsub(rex::rex(capture(.regNum), "%>"), "\033[1;31m\\1%\033[0m ", .pf, perl=TRUE)
        .pf <- gsub(rex::rex(capture(.regNum), "%="), "\033[1;32m\\1%\033[0m ", .pf, perl=TRUE)
        .pf <- gsub(rex::rex(capture(.regNum), "="), "\033[1;32m\\1\033[0m ", .pf, perl=TRUE)
        .pf <- gsub(rex::rex(capture(.regNum), "%<"), "\\1% ", .pf, perl=TRUE)
        .tmp <- c(row.names(x$parFixed), names(x$parFixed))
        .tmp <- .tmp[order(-sapply(.tmp, nchar))]
        .pf <- gsub(rex::rex(boundary,capture(or(.tmp)), boundary), "\033[1m\\1\033[0m", .pf, perl=TRUE);
        .pf <- gsub(rex::rex(capture(or(.tmp))), "\033[1m\\1\033[0m", .pf, perl=TRUE);
        .pf <- gsub(rex::rex("FIXED"), "\033[1;32mFIXED\033[0m", .pf, perl=TRUE)
    } else {
        .pf <- gsub(rex::rex(capture(.regNum), "%", or(">", "=", "<")), "\\1% ", .pf, perl=TRUE)
        .pf <- gsub(rex::rex(capture(.regNum), "="), "\\1 ", .pf, perl=TRUE)
    }
    message(paste(.pf, collapse="\n"), "\n")
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
        .lts <- sapply(.lt, function(x){
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
        .lt <- sprintf("%s%s", formatC(signif(.lt, digits=.digs),digits=.digs,format="fg", flag="#"), .lts)
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
    if (RxODE::rxIs(x, "nlmixrFitData")){
        .dfName <- "data.frame";
        if (RxODE::rxIs(x, "tbl"))  .dfName <- "tibble"
        if (RxODE::rxIs(x, "data.table"))  .dfName <- "data.table"
        message(paste0("\n", cli::rule(paste0(crayon::bold("Fit Data"), " (object", ifelse(.bound == "", "", " "),
                                              crayon::yellow(.bound),
                                              " is a modified ", crayon::blue(.dfName), "):"))))
        if (RxODE::rxIs(x, "tbl") || RxODE::rxIs(x, "data.table")){
            .oldOps <- options();
            on.exit(options(.oldOps))
            options(tibble.print_max = 3, tibble.print_min = 3)
            NextMethod(x)
        } else{
            print(head(x));
        }
    }
}


## FIXME fitted?


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
traceplot.nlmixrFitCore <- function(x, ...){
    .m <- x$parHistStacked;
    if (!is.null(.m)){
        .p0 <- ggplot(.m, aes(iter, val)) +
            geom_line() +
            facet_wrap(~par, scales = "free_y")
        if (!is.null(x$mcmc)){
            .p0 <- .p0 + ggplot2::geom_vline(xintercept=x$mcmc$niter[1], col="blue", size=1.2);
        }
        print(.p0)
    }
}

##' @export
getVarCov.nlmixrFitCore <- function (obj, ...){
    .env <- obj;
    if (RxODE::rxIs(obj, "nlmixrFitData")){
        .env <- obj$env;
    } else {
        class(.env) <- NULL
    }
    if (exists("cov", envir=.env)) return(.env$cov);
    .pt <- proc.time();
    .args <- list(...);
    .control <- .env$control;
    .control$maxInnerIterations <- 0L;
    .control$maxOuterIterations <- 0L;
    .control$boundTol <- 0
    .control$calcTables <- FALSE;
    .dat <- getData(obj);
    .uif <- obj$uif;
    .mat <- as.matrix(nlme::random.effects(obj)[,-1]);
    .skipCov <- obj$skipCov;
    .inits <- list(THTA=as.vector(nlme::fixed.effects(obj)),
                   OMGA=focei.eta.nlmixrFitCore(obj))
    .fit2 <- foceiFit.data.frame(data=.dat,
                                 inits=.inits,
                                 PKpars=.uif$theta.pars,
                                 ## par_trans=fun,
                                 model=.uif$rxode.pred,
                                 pred=function(){return(nlmixr_pred)},
                                 err=.uif$error,
                                 lower=.uif$focei.lower,
                                 upper=.uif$focei.upper,
                                 thetaNames=.uif$focei.names,
                                 etaNames=.uif$eta.names,
                                 etaMat=.mat,
                                 skipCov=.skipCov,
                                 control=.control)
    .env$cov <- .fit2$cov;
    .env$popDf <- .fit2$popDf;
    .env$popDfSig <- .fit2$popDfSig;
    .updateParFixed(.env);
    .parent <- parent.frame(2);
    .bound <- do.call("c", lapply(ls(.parent), function(.cur){
                               if (identical(.parent[[.cur]], obj)){
                                   return(.cur)
                               }
                               return(NULL);
                           }))
    message(paste0("Updated original fit object ", ifelse(is.null(.bound), "", crayon::yellow(.bound))))
    .env$time$covariance <- (proc.time() - .pt)["elapsed"];
    return(.env$cov);
}

focei.eta.nlmixrFitCore <- function(object, ...){
    .uif <- object$uif;
    ## Reorder based on translation
    .df <- as.data.frame(.uif$ini);
    .eta <- .df[!is.na(.df$neta1), ];
    .len <- length(.eta$name)
    .curOme <- object$omega;
    .curLhs <- character()
    .curRhs <- numeric()
    .ome <- character()
    for (.i in seq_along(.eta$name)){
        .lastBlock <- FALSE;
        if (.i == .len){
            .lastBlock <- TRUE
        } else if (.eta$neta1[.i + 1] == .eta$neta2[.i + 1]){
            .lastBlock <- TRUE
        }
        if (.eta$neta1[.i] == .eta$neta2[.i]){
            .curLhs <- c(.curLhs, sprintf("ETA[%d]", .eta$neta1[.i]));
            .curRhs <- c(.curRhs, .curOme[.eta$neta1[.i], .eta$neta2[.i]]);
            if (.lastBlock){
                .ome[length(.ome) + 1] <- sprintf("%s ~ %s", paste(.curLhs, collapse=" + "),
                                                  paste(deparse(.curRhs), collapse=" "));
                .curLhs <- character();
                .curRhs <- numeric()
            }
        } else {
            .curRhs <- c(.curRhs, .curOme[.eta$neta1[.i], .eta$neta2[.i]]);
        }
    }
    .ome <- eval(parse(text=sprintf("list(%s)", paste(.ome, collapse=","))))
    return(.ome)
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
