.regFloat1 <- rex::rex(
  or(
    group(some_of("0":"9"), ".", any_of("0":"9")),
    group(any_of("0":"9"), ".", some_of("0":"9"))
  ),
  maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9")))
)
.regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"))
.regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
.regNum <- rex::rex(maybe("-"), or(.regDecimalint, .regFloat1, .regFloat2))

use.utf <- function() {
  opt <- getOption("cli.unicode", NULL)
  if (!is.null(opt)) {
    isTRUE(opt)
  } else {
    l10n_info()$`UTF-8` && !is.latex()
  }
}

is.latex <- function() {
  if (!("knitr" %in% loadedNamespaces())) {
    return(FALSE)
  }
  get("is_latex_output", asNamespace("knitr"))()
}

##' Control Options for FOCEi
##'
##' @param sigdig Optimization significant digits. This controls:
##'
##' \itemize{
##'
##'  \item The tolerance of the inner and outer optimization is \code{10^-sigdig}
##'
##'  \item The tolerance of the ODE solvers is
##'  \code{0.5*10^(-sigdig-2)}; For the sensitivity equations and
##'  steady-state solutions the default is \code{0.5*10^(-sigdig-1.5)}
##'  (sensitivity changes only applicable for liblsoda)
##'
##'  \item The tolerance of the boundary check is \code{5 * 10 ^ (-sigdig + 1)}
##'
##'  \item The significant figures that some tables are rounded to.
##' }
##'
##' @param atolSens Sensitivity atol, can be different than atol with
##'     liblsoda.  This allows a less accurate solve for gradients (if desired)
##'
##' @param rtolSens Sensitivity rtol, can be different than rtol with
##'     liblsoda.  This allows a less accurate solve for gradients (if desired)
##'
##' @param ssAtol Steady state absolute tolerance (atol) for calculating if steady-state
##'     has been archived.
##'
##' @param ssRtol Steady state relative tolerance (rtol) for
##'     calculating if steady-state has been achieved.
##'
##' @param ssAtolSens Sensitivity absolute tolerance (atol) for
##'     calculating if steady state has been achieved for sensitivity compartments.
##'
##' @param ssRtolSens Sensitivity relative tolerance (rtol) for
##'     calculating if steady state has been achieved for sensitivity compartments.
##'
##' @param epsilon Precision of estimate for n1qn1 optimization.
##'
##' @param maxstepsOde Maximum number of steps for ODE solver.
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
##'     value.  By default this is 1.
##'
##' @param derivEps Forward difference tolerances, which is a
##'     vector of relative difference and absolute difference.  The
##'     central/forward difference step size h is calculated as:
##'
##'         \code{h = abs(x)*derivEps[1] + derivEps[2]}
##'
##' @param derivMethod indicates the method for calculating
##'     derivatives of the outer problem.  Currently supports
##'     "switch", "central" and "forward" difference methods.  Switch
##'     starts with forward differences.  This will switch to central
##'     differences when abs(delta(OFV)) <= derivSwitchTol and switch
##'     back to forward differences when abs(delta(OFV)) >
##'     derivSwitchTol.
##'
##' @param derivSwitchTol The tolerance to switch forward to central
##'     differences.
##'
##' @param covDerivMethod indicates the method for calculating the
##'     derivatives while calculating the covariance components
##'     (Hessian and S).
##'
##' @param covMethod Method for calculating covariance.  In this
##'     discussion, R is the Hessian matrix of the objective
##'     function. The S matrix is the sum of individual
##'     gradient cross-product (evaluated at the individual empirical
##'     Bayes estimates).
##'
##' \itemize{
##'
##'  \item "\code{r,s}" Uses the sandwich matrix to calculate the
##'  covariance, that is: \code{solve(R) \%*\% S \%*\% solve(R)}
##'
##'  \item "\code{r}" Uses the Hessian matrix to calculate the
##'  covariance as \code{2 \%*\% solve(R)}
##'
##'  \item "\code{s}" Uses the cross-product matrix to calculate the
##'  covariance as \code{4 \%*\% solve(S)}
##'
##'  \item "" Does not calculate the covariance step.
##' }
##'
##' @param covTryHarder If the R matrix is non-positive definite and
##'     cannot be corrected to be non-positive definite try estimating
##'     the Hessian on the unscaled parameter space.
##'
##' @param hessEps is a double value representing the epsilon for the Hessian calculation.
##'
##' @param centralDerivEps Central difference tolerances.  This is a
##'     numeric vector of relative difference and absolute difference.
##'     The central/forward difference step size h is calculated as:
##'
##'         \code{h = abs(x)*derivEps[1] + derivEps[2]}
##'
##' @param lbfgsLmm An integer giving the number of BFGS updates
##'     retained in the "L-BFGS-B" method, It defaults to 7.
##'
##' @param lbfgsPgtol is a double precision variable.
##'
##'     On entry pgtol >= 0 is specified by the user.  The iteration
##'     will stop when:
##'
##'        \code{max(\| proj g_i \| i = 1, ..., n) <= lbfgsPgtol}
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
##'     of the \code{chol(solve(omega))}. This matrix and values are the
##'     parameters estimated in FOCEi. The possibilities are:
##'
##' \itemize{
##'  \item \code{sqrt} Estimates the sqrt of the diagonal elements of \code{chol(solve(omega))}.  This is the default method.
##'
##'  \item \code{log} Estimates the log of the diagonal elements of \code{chol(solve(omega))}
##'
##'  \item \code{identity} Estimates the diagonal elements without any transformations
##' }
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
##'     0.95 or 95\% confidence.
##'
##' @param useColor Boolean indicating if focei can use ASCII color codes
##'
##' @param boundTol Tolerance for boundary issues.
##'
##' @param calcTables This boolean is to determine if the foceiFit
##'     will calculate tables. By default this is \code{TRUE}
##'
##' @param ... Ignored parameters
##'
##' @param maxInnerIterations Number of iterations for n1qn1
##'     optimization.
##'
##' @param maxOuterIterations Maximum number of L-BFGS-B optimization
##'     for outer problem.
##'
##' @param n1qn1nsim Number of function evaluations for n1qn1
##'     optimization.
##'
##' @param eigen A boolean indicating if eigenvectors are calculated
##'     to include a condition number calculation.
##'
##' @param addPosthoc Boolean indicating if posthoc parameters are
##'     added to the table output.
##'
##' @param printNcol Number of columns to printout before wrapping
##'     parameter estimates/gradient
##'
##' @param noAbort Boolean to indicate if you should abort the FOCEi
##'     evaluation if it runs into troubles.  (default TRUE)
##'
##' @param interaction Boolean indicate FOCEi should be used (TRUE)
##'     instead of FOCE (FALSE)
##'
##' @param cholSEOpt Boolean indicating if the generalized Cholesky
##'     should be used while optimizing.
##'
##' @param cholSECov Boolean indicating if the generalized Cholesky
##'     should be used while calculating the Covariance Matrix.
##'
##' @param fo is a boolean indicating if this is a FO approximation routine.
##'
##' @param cholSEtol tolerance for Generalized Cholesky
##'     Decomposition.  Defaults to suggested (.Machine$double.eps)^(1/3)
##'
##' @param cholAccept Tolerance to accept a Generalized Cholesky
##'     Decomposition for a R or S matrix.
##'
##' @param outerOpt optimization method for the outer problem
##'
##' @param innerOpt optimization method for the inner problem (not
##'     implemented yet.)
##'
##' @param stateTrim Trim state amounts/concentrations to this value.
##'
##' @param resetEtaP represents the p-value for reseting the
##'     individual ETA to 0 during optimization (instead of the saved
##'     value).  The two test statistics used in the z-test are either
##'     chol(omega^-1) \%*\% eta or eta/sd(allEtas).  A p-value of 0
##'     indicates the ETAs never reset.  A p-value of 1 indicates the
##'     ETAs always reset.
##'
##' @param resetThetaP represents the p-value for reseting the
##'     population mu-referenced THETA parameters based on ETA drift
##'     during optimization, and resetting the optimization.  A
##'     p-value of 0 indicates the THETAs never reset.  A p-value of 1
##'     indicates the THETAs always reset and is not allowed.  The
##'     theta reset is checked at the beginning and when nearing a
##'     local minima.  The percent change in objective function where
##'     a theta reset check is initiated is controlled in
##'     \code{resetThetaCheckPer}.
##'
##' @param resetThetaCheckPer represents objective function
##'     \% percentage below which resetThetaP is checked.
##'
##' @param resetThetaFinalP represents the p-value for reseting the
##'     population mu-referenced THETA parameters based on ETA drift
##'     during optimization, and resetting the optimization one final time.
##'
##' @param resetHessianAndEta is a boolean representing if the
##'     individual Hessian is reset when ETAs are reset using the
##'     option \code{resetEtaP}.
##'
##' @param diagOmegaBoundUpper This represents the upper bound of the
##'     diagonal omega matrix.  The upper bound is given by
##'     diag(omega)*diagOmegaBoundUpper.  If
##'     \code{diagOmegaBoundUpper} is 1, there is no upper bound on
##'     Omega.
##'
##' @param diagOmegaBoundLower This represents the lower bound of the
##'     diagonal omega matrix.  The lower bound is given by
##'     diag(omega)/diagOmegaBoundUpper.  If
##'     \code{diagOmegaBoundLower} is 1, there is no lower bound on
##'     Omega.
##'
##' @param rhobeg Beginning change in parameters for bobyqa algorithm
##'     (trust region).  By default this is 0.2 or 20% of the initial
##'     parameters when the parameters are scaled to 1. rhobeg and
##'     rhoend must be set to the initial and final values of a trust
##'     region radius, so both must be positive with 0 < rhoend <
##'     rhobeg. Typically rhobeg should be about one tenth of the
##'     greatest expected change to a variable.  Note also that
##'     smallest difference abs(upper-lower) should be greater than or
##'     equal to rhobeg*2. If this is not the case then rhobeg will be
##'     adjusted.
##'
##' @param rhoend The smallest value of the trust region radius that
##'     is allowed. If not defined, then 10^(-sigdig-1) will be used.
##'
##' @param npt The number of points used to approximate the objective
##'     function via a quadratic approximation for bobyqa. The value
##'     of npt must be in the interval [n+2,(n+1)(n+2)/2] where n is
##'     the number of parameters in par. Choices that exceed 2*n+1 are
##'     not recommended. If not defined, it will be set to 2*n + 1
##' @param eval.max Number of maximum evaluations of the objective function
##'
##' @param iter.max Maximum number of iterations allowed.
##'
##' @param rel.tol Relative tolerance before nlminb stops.
##'
##' @param x.tol X tolerance for nlmixr optimizers
##'
##' @param abstol Absolute tolerance for nlmixr optimizer
##'
##' @param reltol  tolerance for nlmixr
##'
##' @param gillK The total number of possible steps to determine the
##'     optimal forward/central difference step size per parameter (by
##'     the Gill 1983 method).  If 0, no optimal step size is
##'     determined.  Otherwise this is the optimal step size
##'     determined.
##'
##' @param gillRtol The relative tolerance used for Gill 1983
##'     determination of optimal step size.
##'
##' @param scaleType The scaling scheme for nlmixr.  The supported types are:
##'
##' \itemize{
##' \item \code{nlmixr}  In this approach the scaling is performed by the following equation:
##'
##'    v_{scaled} = (v_{current} - v_{init})/scaleC[i] + scaleTo
##'
##' The \code{scaleTo} parameter is specified by the \code{normType},
##' and the scales are specified by \code{scaleC}.
##'
##' \item \code{norm} This approach uses the simple scaling provided
##'     by the \code{normType} argument.
##'
##' \item \code{mult} This approach does not use the data
##' normalization provided by \code{normType}, but rather uses
##' multiplicative scaling to a constant provided by the \code{scaleTo}
##' argument.
##'
##'   In this case:
##'
##'   v_{scaled} = v_{current}/v_{init}*scaleTo
##'
##' \item \code{multAdd} This approach changes the scaling based on
##' the parameter being specified.  If a parameter is defined in an
##' exponential block (ie exp(theta)), then it is scaled on a
##' linearly, that is:
##'
##'   v_{scaled} = (v_{current}-v_{init}) + scaleTo
##'
##' Otherwise the parameter is scaled multiplicatively.
##'
##'    v_{scaled} = v_{current}/v_{init}*scaleTo
##'
##' }
##'
##' @param scaleC The scaling constant used with
##'     \code{scaleType=nlmixr}.  When not specified, it is based on
##'     the type of parameter that is estimated.  The idea is to keep
##'     the derivatives similar on a log scale to have similar
##'     gradient sizes.  Hence parameters like log(exp(theta)) would
##'     have a scaling factor of 1 and log(theta) would have a scaling
##'     factor of ini_value (to scale by 1/value; ie
##'     d/dt(log(ini_value)) = 1/ini_value or scaleC=ini_value)
##'
##'    \itemize{
##'
##'    \item For parameters in an exponential (ie exp(theta)) or
##'    parameters specifying powers, boxCox or yeoJohnson
##'    transformations , this is 1.
##'
##'    \item For additive, proportional, lognormal error structures,
##'    these are given by 0.5*abs(initial_estimate)
##'
##'    \item Factorials are scaled by abs(1/digamma(inital_estimate+1))
##'
##'    \item parameters in a log scale (ie log(theta)) are transformed
##'    by log(abs(initial_estimate))*abs(initial_estimate)
##'
##'    }
##'
##'    These parameter scaling coefficients are chose to try to keep
##'    similar slopes among parameters.  That is they all follow the
##'    slopes approximately on a log-scale.
##'
##'    While these are chosen in a logical manner, they may not always
##'    apply.  You can specify each parameters scaling factor by this
##'    parameter if you wish.
##'
##' @param scaleC0 Number to adjust the scaling factor by if the initial
##'     gradient is zero.
##'
##' @param scaleCmax Maximum value of the scaleC to prevent overflow.
##'
##' @param scaleCmin Minimum value of the scaleC to prevent underflow.
##'
##' @param normType This is the type of parameter
##'     normalization/scaling used to get the scaled initial values
##'     for nlmixr.  These are used with \code{scaleType} of.
##'
##'     With the exception of \code{rescale2}, these come
##'     from
##'     \href{https://en.wikipedia.org/wiki/Feature_scaling}{Feature
##'     Scaling}. The \code{rescale2} The rescaling is the same type
##'     described in the
##'     \href{http://apmonitor.com/me575/uploads/Main/optimization_book.pdf}{OptdesX}
##'     software manual.
##'
##'     In general, all all scaling formula can be described by:
##'
##'     v_{scaled} = (v_{unscaled}-C_{1})/C_{2}
##'
##'     Where
##'
##'
##'     The other data normalization approaches follow the following formula
##'
##'     v_{scaled} = (v_{unscaled}-C_{1})/C_{2};
##'
##' \itemize{
##'
##' \item \code{rescale2} This scales all parameters from (-1 to 1).
##'     The relative differences between the parameters are preserved
##'     with this approach and the constants are:
##'
##'     C_{1} = (max(all unscaled values)+min(all unscaled values))/2
##'
##'     C_{2} = (max(all unscaled values) - min(all unscaled values))/2
##'
##'
##' \item \code{rescale} or min-max normalization. This rescales all
##'     parameters from (0 to 1).  As in the \code{rescale2} the
##'     relative differences are preserved.  In this approach:
##'
##'     C_{1} = min(all unscaled values)
##'
##'     C_{2} = max(all unscaled values) - min(all unscaled values)
##'
##'
##' \item \code{mean} or mean normalization.  This rescales to center
##'     the parameters around the mean but the parameters are from 0
##'     to 1.  In this approach:
##'
##'     C_{1} = mean(all unscaled values)
##'
##'     C_{2} = max(all unscaled values) - min(all unscaled values)
##'
##' \item \code{std} or standardization.  This standardizes by the mean
##'      and standard deviation.  In this approach:
##'
##'     C_{1} = mean(all unscaled values)
##'
##'     C_{2} = sd(all unscaled values)
##'
##' \item \code{len} or unit length scaling.  This scales the
##'    parameters to the unit length.  For this approach we use the Euclidean length, that
##'    is:
##'
##'     C_{1} = 0
##'
##'     C_{2} = sqrt(v_1^2 + v_2^2 + ... + v_n^2)
##'
##'
##' \item \code{constant} which does not perform data normalization. That is
##'
##'     C_{1} = 0
##'
##'     C_{2} = 1
##'
##' }
##'
##' @param gillStep When looking for the optimal forward difference
##'     step size, this is This is the step size to increase the
##'     initial estimate by.  So each iteration the new step size =
##'     (prior step size)*gillStep
##'
##' @param gillFtol The gillFtol is the gradient error tolerance that
##'     is acceptable before issuing a warning/error about the gradient estimates.
##'
##' @param gillKcov The total number of possible steps to determine
##'     the optimal forward/central difference step size per parameter
##'     (by the Gill 1983 method) during the covariance step.  If 0,
##'     no optimal step size is determined.  Otherwise this is the
##'     optimal step size determined.
##'
##' @param gillStepCov When looking for the optimal forward difference
##'     step size, this is This is the step size to increase the
##'     initial estimate by.  So each iteration during the covariance
##'     step is equal to the new step size = (prior step size)*gillStepCov
##'
##' @param gillFtolCov The gillFtol is the gradient error tolerance
##'     that is acceptable before issuing a warning/error about the
##'     gradient estimates during the covariance step.
##'
##' @param rmatNorm A parameter to normalize gradient step size by the
##'     parameter value during the calculation of the R matrix
##'
##' @param smatNorm A parameter to normalize gradient step size by the
##'     parameter value during the calculation of the S matrix
##'
##' @param covGillF Use the Gill calculated optimal Forward difference
##'     step size for the instead of the central difference step size
##'     during the central difference gradient calculation.
##'
##' @param optGillF Use the Gill calculated optimal Forward difference
##'     step size for the instead of the central difference step size
##'     during the central differences for optimization.
##'
##' @param covSmall The covSmall is the small number to compare
##'     covariance numbers before rejecting an estimate of the
##'     covariance as the final estimate (when comparing sandwich vs
##'     R/S matrix estimates of the covariance).  This number controls
##'     how small the variance is before the covariance matrix is
##'     rejected.
##'
##' @param adjLik In nlmixr, the objective function matches NONMEM's
##'     objective function, which removes a 2*pi constant from the
##'     likelihood calculation. If this is TRUE, the likelihood
##'     function is adjusted by this 2*pi factor.  When adjusted this
##'     number more closely matches the likelihood approximations of
##'     nlme, and SAS approximations.  Regardless of if this is turned
##'     on or off the objective function matches NONMEM's objective
##'     function.
##'
##' @param gradTrim The parameter to adjust the gradient to if the
##'     |gradient| is very large.
##'
##' @param gradCalcCentralSmall A small number that represents the value
##'     where |grad| < gradCalcCentralSmall where forward differences
##'     switch to central differences.
##'
##' @param gradCalcCentralLarge A large number that represents the value
##'     where |grad| > gradCalcCentralLarge where forward differences
##'     switch to central differences.
##'
##' @param etaNudge By default initial ETA estimates start at zero;
##'     Sometimes this doesn't optimize appropriately.  If this value
##'     is non-zero, when the n1qn1 optimization didn't perform
##'     appropriately, reset the Hessian, and nudge the ETA up by this
##'     value; If the ETA still doesn't move, nudge the ETA down by
##'     this value.  Finally if it doesn't move, reset it to zero and
##'     do not perform the optimization again.  This ETA nudge is only
##'     done on the first ETA optimization.
##'
##' @param maxOdeRecalc Maximum number of times to reduce the ODE
##'     tolerances and try to resolve the system if there was a bad
##'     ODE solve.
##'
##' @param odeRecalcFactor The factor to increase the rtol/atol with
##'     bad ODE solving.
##'
##' @param repeatGillMax If the tolerances were reduced when
##'     calculating the initial Gill differences, the Gill difference
##'     is repeated up to a maximum number of times defined by this
##'     parameter.
##'
##' @param stickyRecalcN The number of bad ODE solves before reducing
##'     the atol/rtol for the rest of the problem.
##'
##' @param nRetries If FOCEi doesn't fit with the current parameter
##'     estimates, randomly sample new parameter estimates and restart
##'     the problem.  This is similar to 'PsN' resampling.
##'
##' @param eventFD Finite difference step for forward or central
##'     difference estimation of event-based gradients
##'
##' @param eventCentral Use the central difference approximation when
##'     estimating event-based gradients
##' @param gradProgressOfvTime This is the time for a single objective
##'     function evaluation (in seconds) to start progress bars on gradient evaluations
##'
##' @inheritParams RxODE::rxSolve
##' @inheritParams minqa::bobyqa
##' @inheritParams foceiFit
##'
##' @details
##'
##' Note this uses the R's L-BFGS-B in \code{\link{optim}} for the
##' outer problem and the BFGS \code{\link[n1qn1]{n1qn1}} with that
##' allows restoring the prior individual Hessian (for faster
##' optimization speed).
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
##' @seealso \code{\link[RxODE]{rxSolve}}
##' @export
foceiControl <- function(sigdig = 3, ...,
                         epsilon = NULL, # 1e-4,
                         maxInnerIterations = 1000,
                         maxOuterIterations = 5000,
                         n1qn1nsim = NULL,
                         method = c("liblsoda", "lsoda", "dop853"),
                         transitAbs = NULL, atol = NULL, rtol = NULL,
                         atolSens = NULL, rtolSens = NULL,
                         ssAtol = NULL, ssRtol = NULL, ssAtolSens = NULL, ssRtolSens = NULL,
                         minSS = 10L, maxSS = 1000L,
                         maxstepsOde = 50000L, hmin = 0L, hmax = NA_real_, hini = 0, maxordn = 12L, maxords = 5L, cores,
                         covsInterpolation = c("locf", "linear", "nocb", "midpoint"),
                         print = 1L,
                         printNcol = floor((getOption("width") - 23) / 12),
                         scaleTo = 1.0,
                         scaleObjective = 0,
                         normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
                         scaleType = c("nlmixr", "norm", "mult", "multAdd"),
                         scaleCmax = 1e5,
                         scaleCmin = 1e-5,
                         scaleC = NULL,
                         scaleC0 = 1e5,
                         derivEps = rep(20 * sqrt(.Machine$double.eps), 2),
                         derivMethod = c("switch", "forward", "central"),
                         derivSwitchTol = NULL,
                         covDerivMethod = c("central", "forward"),
                         covMethod = c("r,s", "r", "s", ""),
                         hessEps = (.Machine$double.eps)^(1 / 3),
                         eventFD = sqrt(.Machine$double.eps),
                         eventCentral = TRUE,
                         centralDerivEps = rep(20 * sqrt(.Machine$double.eps), 2),
                         lbfgsLmm = 7L,
                         lbfgsPgtol = 0,
                         lbfgsFactr = NULL,
                         eigen = TRUE,
                         addPosthoc = TRUE,
                         diagXform = c("sqrt", "log", "identity"),
                         sumProd = FALSE,
                         optExpression = TRUE,
                         ci = 0.95,
                         useColor = crayon::has_color(),
                         boundTol = NULL,
                         calcTables = TRUE,
                         noAbort = TRUE,
                         interaction = TRUE,
                         cholSEtol = (.Machine$double.eps)^(1 / 3),
                         cholAccept = 1e-3,
                         resetEtaP = 0.15,
                         resetThetaP = 0.05,
                         resetThetaFinalP = 0.15,
                         diagOmegaBoundUpper = 5, # diag(omega) = diag(omega)*diagOmegaBoundUpper; =1 no upper
                         diagOmegaBoundLower = 100, # diag(omega) = diag(omega)/diagOmegaBoundLower; = 1 no lower
                         cholSEOpt = FALSE,
                         cholSECov = FALSE,
                         fo = FALSE,
                         covTryHarder = FALSE,
                         ## Ranking based on run 025
                         ## L-BFGS-B: 20970.53 (2094.004    429.535)
                         ## bobyqa: 21082.34 (338.677    420.754)
                         ## lbfgsb3* (modified for tolerances):
                         ## nlminb: 20973.468 (755.821    458.343)
                         ## mma: 20974.20 (Time: Opt: 3000.501 Cov: 467.287)
                         ## slsqp: 21023.89 (Time: Opt: 460.099; Cov: 488.921)
                         ## lbfgsbLG: 20974.74 (Time: Opt: 946.463; Cov:397.537)
                         outerOpt = c("nlminb", "bobyqa", "lbfgsb3c", "L-BFGS-B", "mma", "lbfgsbLG", "slsqp", "Rvmmin"),
                         innerOpt = c("n1qn1", "BFGS"),
                         ##
                         rhobeg = .2,
                         rhoend = NULL,
                         npt = NULL,
                         ## nlminb
                         rel.tol = NULL,
                         x.tol = NULL,
                         eval.max = 4000,
                         iter.max = 2000,
                         abstol = NULL,
                         reltol = NULL,
                         resetHessianAndEta = FALSE,
                         stateTrim = Inf,
                         gillK = 10L,
                         gillStep = 4,
                         gillFtol = 0,
                         gillRtol = sqrt(.Machine$double.eps),
                         gillKcov = 10L,
                         gillStepCov = 2,
                         gillFtolCov = 0,
                         rmatNorm = TRUE,
                         smatNorm = TRUE,
                         covGillF = TRUE,
                         optGillF = TRUE,
                         covSmall = 1e-5,
                         adjLik = TRUE, ## Adjust likelihood by 2pi for FOCEi methods
                         gradTrim = Inf,
                         maxOdeRecalc = 5,
                         odeRecalcFactor = 10^(0.5),
                         gradCalcCentralSmall = 1e-4,
                         gradCalcCentralLarge = 1e4,
                         etaNudge = 0.01, stiff,
                         nRetries = 3,
                         seed = 42,
                         resetThetaCheckPer = 0.1,
                         etaMat = NULL,
                         repeatGillMax = 3,
                         stickyRecalcN = 5,
                         gradProgressOfvTime = 10) {
  if (is.null(boundTol)) {
    boundTol <- 5 * 10^(-sigdig + 1)
  }
  if (is.null(epsilon)) {
    epsilon <- 10^(-sigdig - 1)
  }
  if (is.null(abstol)) {
    abstol <- 10^(-sigdig - 1)
  }
  if (is.null(reltol)) {
    reltol <- 10^(-sigdig - 1)
  }
  if (is.null(rhoend)) {
    rhoend <- 10^(-sigdig - 1)
  }
  if (is.null(lbfgsFactr)) {
    lbfgsFactr <- 10^(-sigdig - 1) / .Machine$double.eps
  }
  if (is.null(atol)) {
    atol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(rtol)) {
    rtol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(atolSens)) {
    atolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(rtolSens)) {
    rtolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(ssAtol)) {
    ssAtol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(ssRtol)) {
    ssRtol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(ssAtolSens)) {
    ssAtolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(ssRtolSens)) {
    ssRtolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(rel.tol)) {
    rel.tol <- 10^(-sigdig - 1)
  }
  if (is.null(x.tol)) {
    x.tol <- 10^(-sigdig - 1)
  }
  if (is.null(derivSwitchTol)) {
    derivSwitchTol <- 2 * 10^(-sigdig - 1)
  }
  ## if (is.null(gillRtol)){
  ##     ## FIXME: there is a way to calculate this according to the
  ##     ## Gill paper but it is buried in their optimization book.
  ##     gillRtol <- 10 ^ (-sigdig - 1);
  ## }
  .xtra <- list(...)
  if (is.null(transitAbs) && !is.null(.xtra$transit_abs)) { # nolint
    transitAbs <- .xtra$transit_abs # nolint
  }
  if (missing(covsInterpolation) && !is.null(.xtra$covs_interpolation)) { # nolint
    covsInterpolation <- .xtra$covs_interpolation # nolint
  }
  if (missing(maxInnerIterations) && !is.null(.xtra$max_iterations)) { # nolint
    maxInnerIterations <- .xtra$max_iterations # nolint
  }
  if (!missing(stiff) && missing(method)) {
    if (RxODE::rxIs(stiff, "logical")) {
      if (stiff) {
        method <- "lsoda"
        warning("stiff=TRUE has been replaced with method = \"lsoda\".")
      } else {
        method <- "dop853"
        warning("stiff=FALSE has been replaced with method = \"dop853\".")
      }
    }
  } else {
    if (!RxODE::rxIs(method, "integer")) {
      method <- match.arg(method)
    }
  }
  ## .methodIdx <- c("lsoda"=1L, "dop853"=0L, "liblsoda"=2L);
  ## method <- as.integer(.methodIdx[method]);
  if (RxODE::rxIs(scaleType, "character")) {
    .scaleTypeIdx <- c("norm" = 1L, "nlmixr" = 2L, "mult" = 3L, "multAdd" = 4L)
    scaleType <- as.integer(.scaleTypeIdx[match.arg(scaleType)])
  }
  if (RxODE::rxIs(normType, "character")) {
    .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6)
    normType <- as.integer(.normTypeIdx[match.arg(normType)])
  }
  derivMethod <- match.arg(derivMethod)
  .methodIdx <- c("forward" = 0L, "central" = 1L, "switch" = 3L)
  derivMethod <- as.integer(.methodIdx[derivMethod])
  covDerivMethod <- .methodIdx[match.arg(covDerivMethod)]
  if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1]
  if (!RxODE::rxIs(covsInterpolation, "integer")) {
    covsInterpolation <- tolower(match.arg(
      covsInterpolation,
      c("linear", "locf", "LOCF", "constant", "nocb", "NOCB", "midpoint")
    ))
  }

  ## if (covsInterpolation == "constant") covsInterpolation <- "locf";
  ## covsInterpolation  <- as.integer(which(covsInterpolation == c("linear", "locf", "nocb", "midpoint")) - 1);
  if (missing(cores)) {
    cores <- RxODE::rxCores()
  }
  if (missing(n1qn1nsim)) {
    n1qn1nsim <- 10 * maxInnerIterations + 1
  }
  if (length(covMethod) == 1) {
    if (covMethod == "") {
      covMethod <- 0L
    }
  }
  if (RxODE::rxIs(covMethod, "character")) {
    covMethod <- match.arg(covMethod)
    .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
    covMethod <- .covMethodIdx[match.arg(covMethod)]
  }
  .outerOptTxt <- "custom"
  if (RxODE::rxIs(outerOpt, "character")) {
    outerOpt <- match.arg(outerOpt)
    .outerOptTxt <- outerOpt
    if (outerOpt == "bobyqa") {
      RxODE::rxReq("minqa")
      outerOptFun <- .bobyqa
      outerOpt <- -1L
    } else if (outerOpt == "nlminb") {
      outerOptFun <- .nlminb
      outerOpt <- -1L
    } else if (outerOpt == "mma") {
      outerOptFun <- .nloptr
      outerOpt <- -1L
    } else if (outerOpt == "slsqp") {
      outerOptFun <- .slsqp
      outerOpt <- -1L
    } else if (outerOpt == "lbfgsbLG") {
      outerOptFun <- .lbfgsbLG
      outerOpt <- -1L
    } else if (outerOpt == "Rvmmin") {
      outerOptFun <- .Rvmmin
      outerOpt <- -1L
    } else {
      .outerOptIdx <- c("L-BFGS-B" = 0L, "lbfgsb3c" = 1L)
      outerOpt <- .outerOptIdx[outerOpt]
      if (outerOpt == 1L) {
        RxODE::rxReq("lbfgsb3c")
      }
      outerOptFun <- NULL
    }
  } else if (is(outerOpt, "function")) {
    outerOptFun <- outerOpt
    outerOpt <- -1L
  }
  if (RxODE::rxIs(innerOpt, "character")) {
    .innerOptFun <- c("n1qn1" = 1L, "BFGS" = 2L)
    innerOpt <- setNames(.innerOptFun[match.arg(innerOpt)], NULL)
  }

  if (resetEtaP > 0 & resetEtaP < 1) {
    .resetEtaSize <- qnorm(1 - (resetEtaP / 2))
  } else if (resetEtaP <= 0) {
    .resetEtaSize <- Inf
  } else {
    .resetEtaSize <- 0
  }

  if (resetThetaP > 0 & resetThetaP < 1) {
    .resetThetaSize <- qnorm(1 - (resetThetaP / 2))
  } else if (resetThetaP <= 0) {
    .resetThetaSize <- Inf
  } else {
    stop("Cannot always reset THETAs")
  }
  if (resetThetaFinalP > 0 & resetThetaFinalP < 1) {
    .resetThetaFinalSize <- qnorm(1 - (resetThetaFinalP / 2))
  } else if (resetThetaP <= 0) {
    .resetThetaFinalSize <- Inf
  } else {
    stop("Cannot always reset THETAs")
  }
  .ret <- list(
    maxOuterIterations = as.integer(maxOuterIterations),
    maxInnerIterations = as.integer(maxInnerIterations),
    method = method,
    transitAbs = transitAbs,
    atol = atol,
    rtol = rtol,
    atolSens = atolSens,
    rtolSens = rtolSens,
    ssAtol = ssAtol,
    ssRtol = ssRtol,
    ssAtolSens = ssAtolSens,
    ssRtolSens = ssRtolSens,
    minSS = minSS, maxSS = maxSS,
    maxstepsOde = maxstepsOde,
    hmin = hmin,
    hmax = hmax,
    hini = hini,
    maxordn = maxordn,
    maxords = maxords,
    cores = cores,
    covsInterpolation = covsInterpolation,
    n1qn1nsim = as.integer(n1qn1nsim),
    print = as.integer(print),
    lbfgsLmm = as.integer(lbfgsLmm),
    lbfgsPgtol = as.double(lbfgsPgtol),
    lbfgsFactr = as.double(lbfgsFactr),
    scaleTo = scaleTo,
    epsilon = epsilon,
    derivEps = derivEps,
    derivMethod = derivMethod,
    covDerivMethod = covDerivMethod,
    covMethod = covMethod,
    centralDerivEps = centralDerivEps,
    eigen = as.integer(eigen),
    addPosthoc = as.integer(addPosthoc),
    diagXform = match.arg(diagXform),
    sumProd = sumProd,
    optExpression = optExpression,
    outerOpt = as.integer(outerOpt),
    ci = as.double(ci),
    sigdig = as.double(sigdig),
    scaleObjective = as.double(scaleObjective),
    useColor = as.integer(useColor),
    boundTol = as.double(boundTol),
    calcTables = calcTables,
    printNcol = as.integer(printNcol),
    noAbort = as.integer(noAbort),
    interaction = as.integer(interaction),
    cholSEtol = as.double(cholSEtol),
    hessEps = as.double(hessEps),
    cholAccept = as.double(cholAccept),
    resetEtaSize = as.double(.resetEtaSize),
    resetThetaSize = as.double(.resetThetaSize),
    resetThetaFinalSize = as.double(.resetThetaFinalSize),
    diagOmegaBoundUpper = diagOmegaBoundUpper,
    diagOmegaBoundLower = diagOmegaBoundLower,
    cholSEOpt = as.integer(cholSEOpt),
    cholSECov = as.integer(cholSECov),
    fo = as.integer(fo),
    covTryHarder = as.integer(covTryHarder),
    outerOptFun = outerOptFun,
    ## bobyqa
    rhobeg = as.double(rhobeg),
    rhoend = as.double(rhoend),
    npt = npt,
    ## nlminb
    rel.tol = as.double(rel.tol),
    x.tol = as.double(x.tol),
    eval.max = eval.max,
    iter.max = iter.max,
    innerOpt = innerOpt,
    ## BFGS
    abstol = abstol,
    reltol = reltol,
    derivSwitchTol = derivSwitchTol,
    resetHessianAndEta = as.integer(resetHessianAndEta),
    stateTrim = as.double(stateTrim),
    gillK = as.integer(gillK),
    gillKcov = as.integer(gillKcov),
    gillRtol = as.double(gillRtol),
    gillStep = as.double(gillStep),
    gillStepCov = as.double(gillStepCov),
    scaleType = scaleType,
    normType = normType,
    scaleC = scaleC,
    scaleCmin = as.double(scaleCmin),
    scaleCmax = as.double(scaleCmax),
    scaleC0 = as.double(scaleC0),
    outerOptTxt = .outerOptTxt,
    rmatNorm = as.integer(rmatNorm),
    smatNorm = as.integer(smatNorm),
    covGillF = as.integer(covGillF),
    optGillF = as.integer(optGillF),
    gillFtol = as.double(gillFtol),
    gillFtolCov = as.double(gillFtolCov),
    covSmall = as.double(covSmall),
    adjLik = adjLik,
    gradTrim = as.double(gradTrim),
    gradCalcCentralSmall = as.double(gradCalcCentralSmall),
    gradCalcCentralLarge = as.double(gradCalcCentralLarge),
    etaNudge = as.double(etaNudge),
    maxOdeRecalc = as.integer(maxOdeRecalc),
    odeRecalcFactor = as.double(odeRecalcFactor),
    nRetries = nRetries,
    seed = seed,
    resetThetaCheckPer = resetThetaCheckPer,
    etaMat = etaMat,
    repeatGillMax = as.integer(repeatGillMax),
    stickyRecalcN = as.integer(max(1, abs(stickyRecalcN))),
    eventFD = eventFD,
    eventCentral = as.integer(eventCentral),
    gradProgressOfvTime = gradProgressOfvTime,
    ...
  )
  if (!missing(etaMat) && missing(maxInnerIterations)) {
    warning("By supplying etaMat, assume you wish to evaluate at ETAs, so setting maxInnerIterations=0")
    .ret$maxInnerIterations <- 0L
    .ret$etaMat
  }
  .tmp <- .ret
  .tmp$maxsteps <- maxstepsOde
  .tmp <- do.call(RxODE::rxControl, .tmp)
  .ret$rxControl <- .tmp
  class(.ret) <- "foceiControl"
  return(.ret)
}

.ucminf <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  RxODE::rxReq("ucminf")
  .ctl <- control
  .ctl$stepmax <- control$rhobeg
  .ctl$maxeval <- control$maxOuterIterations
  .ctl <- .ctl[names(.ctl) %in% c("stepmax", "maxeval")]
  .ret <- ucminf::ucminf(par, fn, gr = NULL, ..., control = list(), hessian = 2)
  .ret$x <- .ret$par
  return(.ret)
}

.bobyqa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- control
  if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
  .ctl$iprint <- 0L
  .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]
  .ret <- minqa::bobyqa(par, fn,
    control = .ctl,
    lower = lower,
    upper = upper
  )
  .ret$x <- .ret$par
  .ret$message <- .ret$msg
  .ret$convergence <- .ret$ierr
  .ret$value <- .ret$fval
  return(.ret)
}

.lbfgsb3c <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .w <- which(names(control) %in% c("trace", "factr", "pgtol", "abstol", "reltol", "lmm", "maxit", "iprint"))
  .control <- control[.w]
  .ret <- lbfgsb3c::lbfgsb3c(par = as.vector(par), fn = fn, gr = gr, lower = lower, upper = upper, control = .control)
  .ret$x <- .ret$par
  return(.ret)
}

.lbfgsbO <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .control <- control[names(control) %in% c("trace", "factr", "pgtol", "abstol", "reltol", "lmm", "maxit", "iprint")]
  .w <- which(sapply(.control, is.null))
  .control <- .control[-.w]
  .ret <- optim(
    par = par, fn = fn, gr = gr, method = "L-BFGS-B",
    lower = lower, upper = upper,
    control = .control, hessian = FALSE
  )
  .ret$x <- .ret$par
  return(.ret)
}

.mymin <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .control <- control[names(control) %in% c(
    "eval.max", "iter.max", "trace", "abs.tol",
    "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol", "scale.init", "diff.g"
  )]

  if (all(lower != -Inf) | all(upper != Inf)) {
    warning("Optimization: Boundaries not used in Nelder-Mead")
  }
  fit <- mymin(par, fn, control = .control)
  fit$message <- c("NON-CONVERGENCE", "NELDER_FTOL_REACHED")[1 + fit$convergence]
  fit$x <- fit$par
  return(fit)
}

.nlminb <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- control
  .ctl <- .ctl[names(.ctl) %in% c(
    "eval.max", "iter.max", "trace", "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol",
    "scale.inti", "diff.g"
  )]
  .ctl$trace <- 0
  .ret <- stats::nlminb(
    start = par, objective = fn, gradient = gr, hessian = NULL, control = .ctl,
    lower = lower, upper = upper
  )
  .ret$x <- .ret$par
  ## .ret$message   already there.
  ## .ret$convergence already there.
  return(.ret)
}

.Rvmmin <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  ## Also gives unreasonable estimates
  RxODE::rxReq("Rvmmin")
  .masked <- rep_len(1, length(par))
  .ctl <- list(
    maxit = control$maxOuterIterations,
    ## maxfevals
    trace = 0, dowarn = FALSE, checkgrad = FALSE, checkbounds = FALSE,
    keepinputpar = FALSE, eps = control$abstol
  )
  .ret <- Rvmmin::Rvmmin(par = par, fn = fn, gr = gr, lower = lower, upper = upper, bdmsk = .masked, control = list(), ...)
  .ret$x <- .ret$par
  .ret$message <- .ret$message
  return(.ret)
}

.nloptr <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ..., nloptrAlgoritm = "NLOPT_LD_MMA") {
  RxODE::rxReq("nloptr")
  .ctl <- list(
    algorithm = nloptrAlgoritm,
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn, eval_grad_f = gr,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  return(.ret)
}

.bobyqaNLopt <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- list(
    algorithm = "NLOPT_LN_BOBYQA",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  return(.ret)
}

.slsqp <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  return(.nloptr(par, fn, gr, lower, upper, control, ..., nloptrAlgoritm = "NLOPT_LD_SLSQP"))
}

.lbfgsbLG <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctlLocal <- list(
    algorithm = "NLOPT_LD_LBFGS",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ctl <- opts <- list(
    "algorithm" = "NLOPT_LD_AUGLAG",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    maxeval = control$maxOuterIterations,
    "local_opts" = .ctlLocal,
    "print_level" = 0
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn, eval_grad_f = gr,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  return(.ret)
}

.parseOM <- function(OMGA) {
  .re <- "\\bETA\\[(\\d+)\\]\\b"
  .offset <- as.integer(0)
  lapply(1:length(OMGA), function(.k) {
    .s <- OMGA[[.k]]
    .f <- eval(parse(text = (sprintf("y~%s", deparse(.s[[2]])))))
    .r <- unlist(lapply(attr(terms(.f), "variables"), deparse))[-(1:2)]
    .nr <- length(.r)

    .ix <- grep(.re, .r)
    if (.nr - length(.ix)) stop("invalid OMGA specs")

    .ix <- as.integer(sub(.re, "\\1", .r))
    if (any(.ix - (.offset + 1:.nr))) stop("invalid OMGA specs")
    .offset <<- .offset + .nr
    eval(.s[[3]])
  })
}

.genOM <- function(s) {
  .getNR <- function(.a) round(sqrt(2 * length(.a) + 0.25) - 0.1)
  .nr <- sum(sapply(s, .getNR))
  .mat <- matrix(0, .nr, .nr)
  .offset <- as.integer(0)
  j <- lapply(1:length(s), function(k) {
    .a <- s[[k]]
    .p <- .getNR(.a)
    .starts <- row(.mat) > .offset & col(.mat) > .offset
    .mat[col(.mat) >= row(.mat) & col(.mat) <= .offset + .p & .starts] <<- .a
    .offset <<- .offset + .p
  })
  .a <- .mat[col(.mat) >= row(.mat)]
  .mat <- t(.mat)
  .mat[col(.mat) >= row(.mat)] <- .a
  .mat
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
##' @param env An environment used to build the FOCEi or nlmixr object.
##' @return A focei fit or nlmixr fit object
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
##' ## Wang2007 prop error OBF 39.458 for NONMEM FOCEi, nlmixr matches.
##' fitPi <- foceiFit(w7, inits, mypar2,mod,pred,errProp,
##'      control=foceiControl(maxOuterIterations=0,covMethod=""))
##'
##' print(fitPi$objective)
##'
##' ## Wang2007 prop error OBF 39.207 for NONMEM FOCE; nlmixr matches.
##' fitP <- foceiFit(w7, inits, mypar2,mod,pred,errProp,
##'      control=foceiControl(maxOuterIterations=0,covMethod="",
##'      interaction=FALSE))
##'
##' print(fitP$objective)
##'
##' ## Wang 2007 prop error OBF 39.213 for NONMEM FO; nlmixr matches
##' fitPfo <- foceiFit(w7, inits, mypar2,mod,pred,errProp,
##'      control=foceiControl(maxOuterIterations=0,covMethod="",
##'      fo=TRUE))
##'
##' print(fitPfo$objective)
##'
##' ## Note if you have the etas you can evaluate the likelihood
##' ## of an arbitrary model.  It doesn't have to be solved by
##' ## FOCEi
##'
##' etaMat <- matrix(fitPi$eta[,-1])
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
##' ## Wang2007 add error of 0.026 for NONMEM FO; nlmixr matches
##'
##' fitAfo <- foceiFit(w7, inits, mypar2,mod,pred,errAdd,
##'      control=foceiControl(maxOuterIterations=0,fo=TRUE,covMethod=""))
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
##' ## You can also use the Box-Cox Transform of both sides with
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
##'     add.err <- c(0, 0.001)
##'     lambda <- c(-2, 1, 2)
##'     zeta <- c(0.1, 1, 10)
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
##'     cp ~ pow(prop.err, zeta) + add(add.err) + boxCox(lambda)
##' })}
##'
##' fitIVtbs <- nlmixr(one.compartment.IV.model, datr, "focei")
##'
##' ## If you want to use a variance normalizing distribution with
##' ## negative/positive data you can use the Yeo-Johnson transformation
##' ## as well.  This is implemented by the yeoJohnson(lambda) function.
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
##'     cp ~ pow(prop.err, delta) + add(add.err) + yeoJohnson(lambda)
##' })}
##'
##' fitIVyj <- nlmixr(one.compartment.IV.model, datr, "focei")
##'
##' ## In addition to using L-BFGS-B for FOCEi (outer problem) you may
##' ## use other optimizers.  An example is below
##'
##' one.cmt <- function() {
##'   ini({
##'       tka <- .5   # log Ka
##'       tcl <- -3.2 # log Cl
##'       tv <- -1    # log V
##'       eta.ka ~ 1
##'       eta.cl ~ 2
##'       eta.v ~ 1
##'       add.err <- 0.1
##'   })
##'   model({
##'       ka <- exp(tka + eta.ka)
##'       cl <- exp(tcl + eta.cl)
##'       v <- exp(tv + eta.v)
##'       linCmt() ~ add(add.err)
##'   })
##' }
##'
##' fit <- nlmixr(one.cmt, theo_sd, "focei", foceiControl(outerOpt="bobyqa"))
##'
##' ## You may also make an arbitrary optimizer work by adding a wrapper function:
##'
##' newuoa0 <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...){
##'   ## The function requires par, fn, gr, lower, upper and control
##'   ##
##'   ## The par, fn, gr, lower and upper and sent to the function from nlmixr's focei.
##'   ## The  control is the foceiControl list
##'   ##
##'   ##  The following code modifies the list control list for no warnings.
##'   .ctl <- control;
##'   if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
##'   ## nlmixr will print information this is to suppress the printing from the
##'   ## optimizer
##'   .ctl$iprint <- 0L;
##'   .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")];
##'   ## This does not require gradient and is an unbounded optimization:
##'   .ret <- minqa::newuoa(par, fn, control=.ctl);
##'   ## The return statement must be a list with:
##'   ##    - x for the final parameter message
##'   ##    - message for a minimization message
##'   ##    - convergence for a convergence code
##'   .ret$x <- .ret$par;
##'   .ret$message <- .ret$msg;
##'   .ret$convergence <- .ret$ierr
##'   ## you can access the final list from the optimization by fit$optReturn
##'   return(.ret);
##' }
##'
##' fit <- nlmixr(one.cmt, theo_sd, "focei", foceiControl(outerOpt=newuoa0))
##'
##' }
foceiFit <- function(data, ...) {
  UseMethod("foceiFit")
}

##' @rdname foceiFit
##' @export
focei.fit <- foceiFit
##' @rdname foceiFit
##' @export
foceiFit.data.frame <- function(data, ...) {
  call <- as.list(match.call(expand.dots = TRUE))[-1]
  return(.collectWarnings(do.call(foceiFit.data.frame0, call, envir = parent.frame(1))))
}

.updateParFixed <- function(.ret) {
  if (exists("uif", envir = .ret) & exists("shrink", envir = .ret)) {
    .uif <- .ret$uif
    .lab <- paste(.uif$ini$label[!is.na(.uif$ini$ntheta)])
    .lab[.lab == "NA"] <- ""
    .lab <- gsub(" *$", "", gsub("^ *", "", .lab))
    .muRef <- unlist(.uif$mu.ref)
    if (length(.muRef) > 0) {
      .nMuRef <- names(.muRef)
      .ome <- .ret$omega
      .omegaFix <- as.data.frame(.ret$uif$ini)
      .omegaFix <- .omegaFix[is.na(.omegaFix$ntheta), ]
      .omegaFix <- setNames(.omegaFix$fix, paste(.omegaFix$name))
      .muRef <- structure(as.vector(.nMuRef), .Names = as.vector(.muRef))
      .logEta <- .uif$log.eta
      .digs <- .ret$control$sigdig
      .cvOnly <- TRUE
      .sdOnly <- TRUE
      .cvp <- lapply(row.names(.ret$popDfSig), function(x) {
        .y <- .muRef[x]
        if (is.na(.y)) {
          return(data.frame(ch = " ", v = NA_real_))
        }
        .v <- .ome[.y, .y]
        if (any(.y == .logEta)) {
          .sdOnly <<- FALSE
          return(data.frame(
            ch = paste0(
              ifelse(.omegaFix[.y], "fix(", ""),
              formatC(signif(sqrt(exp(.v) - 1) * 100, digits = .digs),
                digits = .digs, format = "fg", flag = "#"
              ),
              ifelse(.omegaFix[.y], ")", "")
            ),
            v = sqrt(exp(.v) - 1) * 100
          ))
        } else {
          .cvOnly <<- FALSE
          return(data.frame(
            ch = paste0(
              ifelse(.omegaFix[.y], "fix(", ""),
              formatC(signif(sqrt(.v), digits = .digs),
                digits = .digs, format = "fg", flag = "#"
              ),
              ifelse(.omegaFix[.y], ")", "")
            ),
            v = .v
          ))
        }
      })
      .cvp <- do.call("rbind", .cvp)
      .shrink <- .ret$shrink
      .errs <- as.data.frame(.uif$ini)
      .errs <- paste(.errs[which(!is.na(.errs$err)), "name"])
      .sh <- lapply(row.names(.ret$popDfSig), function(x) {
        .y <- .muRef[x]
        if (is.na(.y)) {
          if (any(x == .errs)) {
            .v <- .shrink[7, "IWRES"]
            if (length(.v) != 0) {
              return(data.frame(ch = " ", v = NA_real_))
            }
            if (is.null(.v)) {
              return(data.frame(ch = " ", v = NA_real_))
            }
            if (is.na(.v)) {
              return(data.frame(ch = " ", v = NA_real_))
            }
          } else {
            return(" ")
          }
        } else {
          .v <- .shrink[7, .y]
        }
        if (length(.v) != 1) {
          return(data.frame(ch = " ", v = NA_real_))
        }
        if (is.na(.v)) {
          return(data.frame(ch = " ", v = NA_real_))
        }
        .t <- ">"
        if (.v < 0) {
        } else if (.v < 20) {
          .t <- "<"
        } else if (.v < 30) {
          .t <- "="
        }
        return(data.frame(
          ch = sprintf("%s%%%s", formatC(signif(.v, digits = .digs),
            digits = .digs, format = "fg", flag = "#"
          ), .t),
          v = .v
        ))
      })
      .sh <- do.call("rbind", .sh)
      .ret$parFixed <- data.frame(.ret$popDfSig, "BSD" = .cvp$ch, "Shrink(SD)%" = .sh$ch, check.names = FALSE)
      .ret$parFixedDf <- data.frame(.ret$popDf, "BSD" = .cvp$v, "Shrink(SD)%" = .sh$v, check.names = FALSE)
      .w <- which(names(.ret$parFixed) == "BSD")
      if (length(.w) >= 1) {
        names(.ret$parFixed)[.w] <- ifelse(.sdOnly, "BSV(SD)", ifelse(.cvOnly, "BSV(CV%)", "BSV(CV% or SD)"))
      }
      .w <- which(names(.ret$parFixedDf) == "BSD")
      if (length(.w) >= 1) {
        names(.ret$parFixedDf)[.w] <- ifelse(.sdOnly, "BSV(SD)", ifelse(.cvOnly, "BSV(CV%)", "BSV(CV% or SD)"))
      }
      if (!all(.lab == "")) {
        .ret$parFixed <- data.frame(Parameter = .lab, .ret$parFixed, check.names = FALSE)
      }
    } else {
      .ret$parFixed <- .ret$popDfSig
      .ret$parFixedDf <- .ret$popDfSig
    }
  } else {
    .ret$parFixed <- .ret$popDfSig
    .ret$parFixedDf <- .ret$popDf
  }
  class(.ret$parFixed) <- c("nlmixrParFixed", "data.frame")
}

.thetaReset <- new.env(parent = emptyenv())
##' @rdname foceiFit
##' @export
foceiFit.data.frame0 <- function(data,
                                 inits,
                                 PKpars,
                                 model = NULL,
                                 pred = NULL,
                                 err = NULL,
                                 lower = -Inf,
                                 upper = Inf,
                                 fixed = NULL,
                                 skipCov = NULL,
                                 control = foceiControl(),
                                 thetaNames = NULL,
                                 etaNames = NULL,
                                 etaMat = NULL,
                                 ...,
                                 env = NULL) {
  set.seed(control$seed)
  .pt <- proc.time()
  loadNamespace("n1qn1")
  if (!RxODE::rxIs(control, "foceiControl")) {
    control <- do.call(foceiControl, control)
  }
  if (is.null(env)) {
    .ret <- new.env(parent = emptyenv())
  } else {
    .ret <- env
  }
  .ret$origData <- data
  .ret$etaNames <- etaNames
  .ret$thetaFixed <- fixed
  .ret$control <- control
  .ret$control$focei.mu.ref <- integer(0)
  if (is(model, "RxODE") || is(model, "character")) {
    .ret$ODEmodel <- TRUE
    if (class(pred) != "function") {
      stop("pred must be a function specifying the prediction variables in this model.")
    }
  } else {
    ## Fixme
    .ret$ODEmodel <- TRUE
    model <- RxODE::rxGetLin(PKpars)
    pred <- eval(parse(text = "function(){return(Central);}"))
  }
  .square <- function(x) x * x
  .ret$diagXformInv <- c("sqrt" = ".square", "log" = "exp", "identity" = "identity")[control$diagXform]
  if (is.null(err)) {
    err <- eval(parse(text = paste0("function(){err", paste(inits$ERROR[[1]], collapse = ""), "}")))
  }
  .covNames <- .parNames <- c()
  .ret$adjLik <- control$adjLik
  if (!exists("noLik", envir = .ret)) {
    .atol <- rep(control$atol, length(RxODE::rxModelVars(model)$state))
    .rtol <- rep(control$rtol, length(RxODE::rxModelVars(model)$state))
    .ssAtol <- rep(control$ssAtol, length(RxODE::rxModelVars(model)$state))
    .ssRtol <- rep(control$ssRtol, length(RxODE::rxModelVars(model)$state))
    .ret$model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err,
      grad = (control$derivMethod == 2L),
      pred.minus.dv = TRUE, sum.prod = control$sumProd,
      theta.derivs = FALSE, optExpression = control$optExpression,
      interaction = (control$interaction == 1L),
      run.internal = TRUE
    )
    if (!is.null(.ret$model$inner)) {
      .atol <- c(.atol, rep(
        control$atolSens,
        length(RxODE::rxModelVars(.ret$model$inner)$state) -
          length(.atol)
      ))
      .rtol <- c(.rtol, rep(
        control$rtolSens,
        length(RxODE::rxModelVars(.ret$model$inner)$state) -
          length(.rtol)
      ))
      .ret$control$rxControl$atol <- .atol
      .ret$control$rxControl$rtol <- .rtol
      .ssAtol <- c(.ssAtol, rep(
        control$ssAtolSens,
        length(RxODE::rxModelVars(.ret$model$inner)$state) -
          length(.ssAtol)
      ))
      .ssRtol <- c(.ssRtol, rep(
        control$ssRtolSens,
        length(RxODE::rxModelVars(.ret$model$inner)$state) -
          length(.ssRtol)
      ))
      .ret$control$rxControl$ssAtol <- .ssAtol
      .ret$control$rxControl$ssRtol <- .ssRtol
    }
    .covNames <- .parNames <- RxODE::rxParams(.ret$model$pred.only)
    .covNames <- .covNames[regexpr(rex::rex(start, or("THETA", "ETA"), "[", numbers, "]", end), .covNames) == -1]
    colnames(data) <- sapply(names(data), function(x) {
      if (any(x == .covNames)) {
        return(x)
      } else {
        return(toupper(x))
      }
    })
    .lhs <- c(names(RxODE::rxInits(.ret$model$pred.only)), RxODE::rxLhs(.ret$model$pred.only))
    if (length(.lhs) > 0) {
      .covNames <- .covNames[regexpr(rex::rex(start, or(.lhs), end), .covNames) == -1]
    }
    if (length(.covNames) > 0) {
      if (!all(.covNames %in% names(data))) {
        message("Model:")
        RxODE::rxCat(.ret$model$pred.only)
        message("Needed Covariates:")
        nlmixrPrint(.covNames)
        stop("Not all the covariates are in the dataset.")
      }
      message("Needed Covariates:")
      print(.covNames)
    }
    .extraPars <- .ret$model$extra.pars
  } else {
    if (.ret$noLik) {
      .atol <- rep(control$atol, length(RxODE::rxModelVars(model)$state))
      .rtol <- rep(control$rtol, length(RxODE::rxModelVars(model)$state))
      .ret$model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err,
        grad = (control$derivMethod == 2L),
        pred.minus.dv = TRUE, sum.prod = control$sumProd,
        theta.derivs = FALSE, optExpression = control$optExpression, run.internal = TRUE,
        only.numeric = TRUE
      )
      if (!is.null(.ret$model$inner)) {
        .atol <- c(.atol, rep(
          control$atolSens,
          length(RxODE::rxModelVars(.ret$model$inner)$state) -
            length(.atol)
        ))
        .rtol <- c(.rtol, rep(
          control$rtolSens,
          length(RxODE::rxModelVars(.ret$model$inner)$state) -
            length(.rtol)
        ))
        .ret$control$rxControl$atol <- .atol
        .ret$control$rxControl$rtol <- .rtol
      }
      .covNames <- .parNames <- RxODE::rxParams(.ret$model$pred.only)
      .covNames <- .covNames[regexpr(rex::rex(start, or("THETA", "ETA"), "[", numbers, "]", end), .covNames) == -1]
      colnames(data) <- sapply(names(data), function(x) {
        if (any(x == .covNames)) {
          return(x)
        } else {
          return(toupper(x))
        }
      })
      .lhs <- c(names(RxODE::rxInits(.ret$model$pred.only)), RxODE::rxLhs(.ret$model$pred.only))
      if (length(.lhs) > 0) {
        .covNames <- .covNames[regexpr(rex::rex(start, or(.lhs), end), .covNames) == -1]
      }
      if (length(.covNames) > 0) {
        if (!all(.covNames %in% names(data))) {
          message("Model:")
          RxODE::rxCat(.ret$model$pred.only)
          message("Needed Covariates:")
          nlmixrPrint(.covNames)
          stop("Not all the covariates are in the dataset.")
        }
        message("Needed Covariates:")
        print(.covNames)
      }
      .extraPars <- .ret$model$extra.pars
    } else {
      .extraPars <- NULL
    }
  }
  .ret$skipCov <- skipCov
  if (is.null(skipCov)) {
    if (is.null(fixed)) {
      .tmp <- rep(FALSE, length(inits$THTA))
    } else {
      if (length(fixed) < length(inits$THTA)) {
        .tmp <- c(fixed, rep(FALSE, length(inits$THTA) - length(fixed)))
      } else {
        .tmp <- fixed[1:length(inits$THTA)]
      }
    }
    if (exists("uif", envir = .ret)) {
      .uifErr <- .ret$uif$ini$err[!is.na(.ret$uif$ini$ntheta)]
      .uifErr <- sapply(.uifErr, function(x) {
        if (is.na(x)) {
          return(FALSE)
        }
        return(!any(x == c("pow2", "tbs", "tbsYj")))
      })
      .tmp <- (.tmp | .uifErr)
    }
    .ret$skipCov <- c(
      .tmp,
      rep(TRUE, length(.extraPars))
    )
    .ret$control$focei.mu.ref <- .ret$uif$focei.mu.ref
  }
  if (is.null(.extraPars)) {
    .nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)))
  } else {
    .nms <- c(
      sprintf("THETA[%s]", seq_along(inits$THTA)),
      sprintf("ERR[%s]", seq_along(.extraPars))
    )
  }
  if (!is.null(thetaNames) && (length(inits$THTA) + length(.extraPars)) == length(thetaNames)) {
    .nms <- thetaNames
  }
  .ret$thetaNames <- .nms
  if (length(lower) == 1) {
    lower <- rep(lower, length(inits$THTA))
  } else if (length(lower) != length(inits$THTA)) {
    print(inits$THTA)
    print(lower)
    stop("Lower must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
  }
  if (length(upper) == 1) {
    upper <- rep(upper, length(inits$THTA))
  } else if (length(lower) != length(inits$THTA)) {
    stop("Upper must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
  }

  if (!is.null(.extraPars)) {
    .ret$model$extra.pars <- eval(call(control$diagXform, .ret$model$extra.pars))
    if (length(.ret$model$extra.pars) > 0) {
      inits$THTA <- c(inits$THTA, .ret$model$extra.pars)
      .lowerErr <- rep(control$atol[1] * 10, length(.ret$model$extra.pars))
      .upperErr <- rep(Inf, length(.ret$model$extra.pars))
      lower <- c(lower, .lowerErr)
      upper <- c(upper, .upperErr)
    }
  }
  if (is.null(data$ID)) stop('"ID" not found in data')
  if (is.null(data$DV)) stop('"DV" not found in data')
  if (is.null(data$EVID)) data$EVID <- 0
  if (is.null(data$AMT)) data$AMT <- 0
  ## Make sure they are all double amounts.
  for (.v in c("TIME", "AMT", "DV", .covNames)) {
    data[[.v]] <- as.double(data[[.v]])
  }
  .ret$dataSav <- data
  .ds <- data[data$EVID != 0 & data$EVID != 2, c("ID", "TIME", "AMT", "EVID", .covNames)]
  .w <- which(tolower(names(data)) == "limit")
  .limitName <- NULL
  if (length(.w) == 1L) {
    .limitName <- names(data)[.w]
  }
  .censName <- NULL
  .w <- which(tolower(names(data)) == "cens")
  if (length(.w) == 1L) {
    .censName <- names(data[.w])
  }
  data <- data[
    data$EVID == 0 | data$EVID == 2,
    c("ID", "TIME", "DV", "EVID", .covNames, .limitName, .censName)
  ]
  ## keep the covariate names the same as in the model
  .w <- which(!(names(.ret$dataSav) %in% .covNames))
  names(.ret$dataSav)[.w] <- tolower(names(.ret$dataSav[.w])) # needed in ev

  .lh <- .parseOM(inits$OMGA)
  .nlh <- sapply(.lh, length)
  .osplt <- rep(1:length(.lh), .nlh)
  .lini <- list(inits$THTA, unlist(.lh))
  .nlini <- sapply(.lini, length)
  .nsplt <- rep(1:length(.lini), .nlini)

  .om0 <- .genOM(.lh)
  if (length(etaNames) == dim(.om0)[1]) {
    .ret$etaNames <- .ret$etaNames
  } else {
    .ret$etaNames <- sprintf("ETA[%d]", seq(1, dim(.om0)[1]))
  }
  .ret$rxInv <- RxODE::rxSymInvCholCreate(mat = .om0, diag.xform = control$diagXform)
  .ret$xType <- .ret$rxInv$xType
  .om0a <- .om0
  .om0a <- .om0a / control$diagOmegaBoundLower
  .om0b <- .om0
  .om0b <- .om0b * control$diagOmegaBoundUpper
  .om0a <- RxODE::rxSymInvCholCreate(mat = .om0a, diag.xform = control$diagXform)
  .om0b <- RxODE::rxSymInvCholCreate(mat = .om0b, diag.xform = control$diagXform)
  .omdf <- data.frame(a = .om0a$theta, m = .ret$rxInv$theta, b = .om0b$theta, diag = .om0a$theta.diag)
  .omdf$lower <- with(.omdf, ifelse(a > b, b, a))
  .omdf$lower <- with(.omdf, ifelse(lower == m, -Inf, lower))
  .omdf$lower <- with(.omdf, ifelse(!diag, -Inf, lower))
  .omdf$upper <- with(.omdf, ifelse(a < b, b, a))
  .omdf$upper <- with(.omdf, ifelse(upper == m, Inf, upper))
  .omdf$upper <- with(.omdf, ifelse(!diag, Inf, upper))
  .ret$control$ntheta <- length(lower)
  .ret$control$nomega <- length(.omdf$lower)
  .ret$control$neta <- sum(.omdf$diag)
  .ret$control$nfixed <- sum(fixed)
  lower <- c(lower, .omdf$lower)
  upper <- c(upper, .omdf$upper)

  .ret$lower <- lower
  .ret$upper <- upper

  .ret$thetaIni <- inits$THTA

  .scaleC <- double(length(lower))
  if (is.null(control$scaleC)) {
    .scaleC <- rep(NA_real_, length(lower))
  } else {
    .scaleC <- as.double(control$scaleC)
    if (length(lower) > length(.scaleC)) {
      .scaleC <- c(.scaleC, rep(NA_real_, length(lower) - length(.scaleC)))
    } else if (length(lower) < length(.scaleC)) {
      .scaleC <- .scaleC[seq(1, length(lower))]
      warning("scaleC control option has more options than estimated population parameters, please check.")
    }
  }

  .ret$scaleC <- .scaleC
  if (exists("uif", envir = .ret)) {
    .ini <- as.data.frame(.ret$uif$ini)[!is.na(.ret$uif$ini$err), c("est", "err", "ntheta")]
    for (.i in seq_along(.ini$err)) {
      if (is.na(.ret$scaleC[.ini$ntheta[.i]])) {
        if (any(.ini$err[.i] == c("boxCox", "yeoJohnson", "pow2", "tbs", "tbsYj"))) {
          .ret$scaleC[.ini$ntheta[.i]] <- 1
        } else if (any(.ini$err[.i] == c("prop", "add", "norm", "dnorm", "logn", "dlogn", "lnorm", "dlnorm"))) {
          .ret$scaleC[.ini$ntheta[.i]] <- 0.5 * abs(.ini$est[.i])
        }
      }
    }

    for (.i in .ini$model$extraProps$powTheta) {
      if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- 1 ## Powers are log-scaled
    }
    .ini <- as.data.frame(.ret$uif$ini)
    for (.i in .ini$model$extraProps$factorial) {
      if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- abs(1 / digamma(.ini$est[.i] + 1))
    }
    for (.i in .ini$model$extraProps$gamma) {
      if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- abs(1 / digamma(.ini$est[.i]))
    }
    for (.i in .ini$model$extraProps$log) {
      if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- log(abs(.ini$est[.i])) * abs(.ini$est[.i])
    }
    ## FIXME: needs to be based on actual initial values in sin because typically change to correct scale
    ## Ctime is also usually used for circadian rhythm models
    ## for (.i in .ini$model$extraProps$sin){
    ##     if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- fabs(tan(.ini$est[.i]));
    ## }
    ## for (.i in .ini$model$extraProps$cos){
    ##     if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- fabs(1 / tan(.ini$est[.i]));
    ## }
    ## for (.i in .ini$model$extraProps$tan){
    ##     if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- fabs(2 * sin(2 * .ini$est[.i]));
    ## }
  }
  names(.ret$thetaIni) <- sprintf("THETA[%d]", seq_along(.ret$thetaIni))
  if (is.null(etaMat) & !is.null(control$etaMat)) {
    .ret$etaMat <- control$etaMat
  } else {
    .ret$etaMat <- etaMat
  }
  .ret$setupTime <- (proc.time() - .pt)["elapsed"]
  if (exists("uif", envir = .ret)) {
    .tmp <- .ret$uif$logThetasList
    .ret$logThetas <- .tmp[[1]]
    .ret$logThetasF <- .tmp[[2]]
  } else {
    .ret$logThetasF <- integer(0)
  }
  if (exists("noLik", envir = .ret)) {
    if (!.ret$noLik) {
      .ret$.params <- c(
        sprintf("THETA[%d]", seq_along(.ret$thetaIni)),
        sprintf("ETA[%d]", seq(1, dim(.om0)[1]))
      )
      .ret$.thetan <- length(.ret$thetaIni)
      .ret$nobs <- sum(data$EVID == 0)
    }
  }
  .ret$control$printTop <- TRUE
  .ret$control$nF <- 0
  .est0 <- .ret$thetaIni
  if (!is.null(.ret$model$pred.nolhs)) {
    .ret$control$predNeq <- length(.ret$model$pred.nolhs$state)
  } else {
    .ret$control$predNeq <- 0L
  }
  .fitFun <- function(.ret) {
    this.env <- environment()
    assign("err", "theta reset", this.env)
    while (this.env$err == "theta reset") {
      assign("err", "", this.env)
      .ret0 <- tryCatch(
        {
          foceiFitCpp_(.ret)
        },
        error = function(e) {
          if (regexpr("theta reset", e$message) != -1) {
            assign("err", "theta reset", this.env)
          } else {
            assign("err", e$message, this.env)
          }
        }
      )
      if (this.env$err == "theta reset") {
        .nm <- names(.ret$thetaIni)
        .ret$thetaIni <- setNames(.thetaReset$thetaIni + 0.0, .nm)
        .ret$rxInv$theta <- .thetaReset$omegaTheta
        .ret$control$printTop <- FALSE
        .ret$etaMat <- .thetaReset$etaMat
        .ret$control$etaMat <- .thetaReset$etaMat
        .ret$control$maxInnerIterations <- .thetaReset$maxInnerIterations
        .ret$control$nF <- .thetaReset$nF
        .ret$control$gillRetC <- .thetaReset$gillRetC
        .ret$control$gillRet <- .thetaReset$gillRet
        .ret$control$gillRet <- .thetaReset$gillRet
        .ret$control$gillDf <- .thetaReset$gillDf
        .ret$control$gillDf2 <- .thetaReset$gillDf2
        .ret$control$gillErr <- .thetaReset$gillErr
        .ret$control$rEps <- .thetaReset$rEps
        .ret$control$aEps <- .thetaReset$aEps
        .ret$control$rEpsC <- .thetaReset$rEpsC
        .ret$control$aEpsC <- .thetaReset$aEpsC
        .ret$control$c1 <- .thetaReset$c1
        .ret$control$c2 <- .thetaReset$c2
        message("Theta reset")
      }
    }
    if (this.env$err != "") {
      stop(this.env$err)
    } else {
      return(.ret0)
    }
  }
  .ret0 <- try(.fitFun(.ret))
  .n <- 1
  while (inherits(.ret0, "try-error") && control$maxOuterIterations != 0 && .n <= control$nRetries) {
    ## Maybe change scale?
    message(sprintf("Restart %s", .n))
    .ret$control$nF <- 0
    .estNew <- .est0 + 0.2 * .n * abs(.est0) * stats::runif(length(.est0)) - 0.1 * .n
    .estNew <- sapply(
      seq_along(.est0),
      function(.i) {
        if (.ret$thetaFixed[.i]) {
          return(.est0[.i])
        } else if (.estNew[.i] < lower[.i]) {
          return(lower + (.Machine$double.eps)^(1 / 7))
        } else if (.estNew[.i] > upper[.i]) {
          return(upper - (.Machine$double.eps)^(1 / 7))
        } else {
          return(.estNew[.i])
        }
      }
    )
    .ret$thetaIni <- .estNew
    .ret0 <- try(.fitFun(.ret))
    .n <- .n + 1
  }
  if (inherits(.ret0, "try-error")) stop("Could not fit data.")
  .ret <- .ret0
  if (!control$calcTables) {
    .etas <- .ret$ranef
    .thetas <- .ret$fixef
    .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas)
    .ret$shrink <- .Call(`_nlmixr_nlmixrShrink`, .ret$omega, .etas, .pars$eta.lst[-(dim(.ret$omega)[1] + 1)])
    .updateParFixed(.ret)
    return(.ret)
  }
  if (exists("parHistData", .ret)) {
    .tmp <- .ret$parHistData
    .tmp <- .tmp[.tmp$type == "Unscaled", names(.tmp) != "type"]
    .iter <- .tmp$iter
    .tmp <- .tmp[, names(.tmp) != "iter"]
    .ret$parHistStacked <- data.frame(stack(.tmp), iter = .iter)
    names(.ret$parHistStacked) <- c("val", "par", "iter")
    .ret$parHist <- data.frame(iter = .iter, .tmp)
  }
  .solve <- function(...) {
    .ret <- RxODE::rxSolve(..., warnIdSort = FALSE)
    if (names(.ret)[1] == "time") {
      ## For single subject ID is dropped.
      .ret <- data.frame(ID = 1, .ret)
    }
    return(.ret)
  }
  .solvePred <- function() {
    .res <- .solve(.ret$model$pred.only, .pars$pred, .ret$dataSav,
      returnType = "data.frame",
      atol = .ret$control$atol[1], rtol = .ret$control$rtol[1],
      maxsteps = .ret$control$maxstepsOde,
      hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini,
      transitAbs = .ret$control$transitAbs, maxordn = .ret$control$maxordn,
      maxords = .ret$control$maxords, method = .ret$control$method
    )
    if (any(is.na(.res$rx_pred_)) && .ret$control$method == 2L) {
      .res <- .solve(.ret$model$pred.only, .pars$pred, .ret$dataSav,
        returnType = "data.frame",
        atol = .ret$control$atol[1], rtol = .ret$control$rtol[1],
        maxsteps = .ret$control$maxstepsOde * 2,
        hmin = .ret$control$hmin, hmax = .ret$control$hmax / 2, hini = .ret$control$hini,
        transitAbs = .ret$control$transitAbs, maxordn = .ret$control$maxordn,
        maxords = .ret$control$maxords, method = "lsoda"
      )
      if (any(is.na(.res$rx_pred_))) {
        .res <- .solve(.ret$model$pred.only, .pars$pred, .ret$dataSav,
          returnType = "data.frame",
          atol = .ret$control$atol[1], rtol = .ret$control$rtol[1],
          maxsteps = .ret$control$maxstepsOde * 2,
          hmin = .ret$control$hmin, hmax = .ret$control$hmax / 2, hini = .ret$control$hini,
          transitAbs = .ret$control$transitAbs, maxordn = .ret$control$maxordn,
          maxords = .ret$control$maxords, method = "dop853"
        )
        if (any(is.na(.res$rx_pred_))) {
          warning("Problems solving pred/wres liblsoda, lsoda and dop853")
        } else {
          warning("Problems solving pred/wres liblsoda and lsoda switched to dop853")
        }
      } else {
        warning("Problems solving pred/wres liblsoda switched to lsoda")
      }
    }
    return(.res)
  }
  if (exists("noLik", envir = .ret)) {
    if (.ret$noLik) {
      message("Calculating residuals/tables")
      .pt <- proc.time()
      .etas <- .ret$ranef
      .thetas <- .ret$fixef
      .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas)
      .preds <- list(
        ipred = .solve(.ret$model$pred.only, .pars$ipred, .ret$dataSav,
          returnType = "data.frame.TBS",
          atol = .ret$control$atol[1], rtol = .ret$control$rtol[1], maxsteps = .ret$control$maxstepsOde,
          hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini,
          transitAbs = .ret$control$TransitAbs,
          maxordn = .ret$control$maxordn, maxords = .ret$control$maxords,
          method = .ret$control$method
        ),
        pred = .solvePred(),
        cwres = FALSE
      )
    } else {
      .pt <- proc.time()
      .etas <- .ret$ranef
      .thetas <- .ret$fixef
      .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas)
      .ret$shrink <- .Call(`_nlmixr_nlmixrShrink`, .ret$omega, .etas, .pars$eta.lst[-(dim(.ret$omega)[1] + 1)])
      .updateParFixed(.ret)
      return(.ret)
    }
  } else {
    if (exists("skipTable", envir = .ret)) {
      .etas <- .ret$ranef
      .thetas <- .ret$fixef
      .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas)
      .ret$shrink <- .Call(`_nlmixr_nlmixrShrink`, .ret$omega, .etas, .pars$eta.lst[-(dim(.ret$omega)[1] + 1)])
      .updateParFixed(.ret)
      if (.ret$skipTable) {
        return(.ret)
      }
    }
    message("Calculating residuals/tables")
    .pt <- proc.time()
    .etas <- .ret$ranef
    .thetas <- .ret$fixef
    .pars <- .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas)
    .preds <- list(
      ipred = .solve(.ret$model$inner, .pars$ipred, .ret$dataSav,
        returnType = "data.frame.TBS",
        atol = .ret$control$atol[1], rtol = .ret$control$rtol[1], maxsteps = .ret$control$maxstepsOde,
        hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini,
        transitAbs = .ret$control$TransitAbs,
        maxordn = .ret$control$maxordn, maxords = .ret$control$maxords,
        method = .ret$control$method
      ),
      pred = .solvePred()
    )
  }
  if (!is.null(.censName)) {
    .cens <- data[, .censName]
  } else {
    .cens <- rep(0L, length(data$DV))
  }
  if (!is.null(.limitName)) {
    .limit <- data[, .limitName]
  } else {
    .limit <- rep(-Inf, length(data$DV))
  }
  .lst <- .Call(
    `_nlmixr_nlmixrResid`, .preds, .ret$omega, data$DV, data$EVID, .preds$ipred$rxLambda, .preds$ipred$rxYj,
    .cens, .limit, .etas, .pars$eta.lst
  )
  if (is.null(.preds$cwres)) {
    .df <- RxODE::rxSolve(.ret$model$pred.only, .pars$ipred, .ret$dataSav,
      returnType = "data.frame",
      hmin = .ret$control$hmin, hmax = .ret$control$hmax, hini = .ret$control$hini, transitAbs = .ret$control$transitAbs,
      maxordn = .ret$control$maxordn, maxords = .ret$control$maxords,
      method = .ret$control$method,
      warnIdSort = FALSE
    )[, -(1:4)]
  } else {
    .df <- .preds$ipred[, -c(1:4, length(names(.preds$ipred)) - 0:1), drop = FALSE]
  }
  .df <- .df[, !(names(.df) %in% c("nlmixr_pred", names(.etas), names(.thetas)))]
  if (any(names(.df) %in% names(.lst[[1]]))) {
    warning("Calculated residuals like IPRED are masked by nlmixr calculated values")
    .df <- .df[, !(names(.df) %in% names(.lst[[1]]))]
  }
  if (any(names(.df) %in% names(.lst[[3]]))) {
    .df <- .df[, !(names(.df) %in% names(.lst[[3]]))]
  }
  .lst[[5]] <- .df
  .ret$shrink <- .lst[[2]]
  .updateParFixed(.ret)
  .df <- cbind(as.data.frame(data), .lst[[1]], .lst[[3]], .lst[[5]])
  .df$DV <- .lst[[4]]
  .ret$tableTime <- (proc.time() - .pt)["elapsed"]
  .ret$time <- data.frame(.ret$time, table = .ret$tableTime)
  .isDplyr <- requireNamespace("dplyr", quietly = TRUE)
  if (!.isDplyr) {
    .isDataTable <- requireNamespace("data.table", quietly = TRUE)
    if (.isDataTable) {
      .df <- data.table::data.table(.df)
    }
  } else {
    .df <- tibble::as_tibble(.df)
  }
  .cls <- class(.df)
  if (control$interaction) {
    .cls <- c(paste0("nlmixr", .ret$method, "i"), "nlmixrFitData", "nlmixrFitCore", .cls)
  } else {
    .cls <- c(paste0("nlmixr", .ret$method), "nlmixrFitData", "nlmixrFitCore", .cls)
  }
  class(.ret) <- "nlmixrFitCoreSilent"
  attr(.cls, ".foceiEnv") <- .ret
  class(.df) <- .cls
  message("done.")
  return(.df)
}

##' @export
`$.nlmixrFitCore` <- function(obj, arg, exact = FALSE) {
  .env <- obj
  if (arg == "posthoc") {
    return(nlmixrPosthoc(obj))
  } else if (arg == "notes") {
    return(.notesFit(obj))
  } else if (any(arg == c(
    "logLik", "value", "obf", "ofv",
    "objf", "OBJF", "objective", "AIC",
    "BIC"
  ))) {
    if (!is.null(obj$saem)) {
      .tmp <- obj$saem
      .curObj <- get("objective", .env)
      if (is.na(.curObj)) {
        .nnodes <- 3
        if (exists("nnodes.gq", .env)) {
          .nnodes <- .env$nnodes.gq
        }
        .nsd <- 1.6
        if (exists("nsd.gq", .env)) {
          .nsd <- .env$nsd.gq
        }
        if (.nnodes == 1) {
          setOfv(obj, paste0("laplace", .nsd))
        } else {
          setOfv(obj, paste0("gauss", .nnodes, "_", .nsd))
        }
      }
    }
  }
  if (any(arg == c("value", "obf", "ofv"))) arg <- "objf"
  if (arg == "sigma") {
    return(.sigma(obj))
  }
  if (arg == "coefficients") {
    return(list(
      fixed = fixef(obj),
      random = ranef(obj)
    ))
  }
  if (arg == "par.hist") arg <- "parHist"
  if (arg == "par.hist.stacked") arg <- "parHistStacked"
  if (arg == "omega.R") arg <- "omegaR"
  if (arg == "par.fixed") arg <- "parFixed"
  if (arg == "eta") arg <- "ranef"
  if (arg == "theta") arg <- "fixef"
  if (arg == "varFix") arg <- "cov"
  if (arg == "thetaMat") arg <- "cov"
  if (arg == "seed" && exists("saem", .env)) {
    return(attr(.env$saem, "saem.cfg")$seed)
  }
  if (arg == "saem.cfg" && exists("saem", .env)) {
    return(attr(.env$saem, "saem.cfg"))
  }
  if (exists(arg, envir = .env)) {
    return(get(arg, envir = .env))
  }
  if (arg == "env") {
    return(.env)
  }
  if (exists("uif", .env)) {
    .uif <- .env$uif
    if (arg == "modelName") arg <- "model.name"
    if (arg == "dataName") arg <- "data.name"
    .ret <- `$.nlmixrUI`(.uif, arg)
    if (!is.null(.ret)) {
      return(.ret)
    }
    .env2 <- `$.nlmixrUI`(.uif, "env")
    if (exists(arg, envir = .env2)) {
      return(get(arg, envir = .env2))
    }
  }
}

##' @export
`$.nlmixrFitCoreSilent` <- `$.nlmixrFitCore`

##' @export
`$.nlmixrFitData` <- function(obj, arg, exact = FALSE) {
  .ret <- obj[[arg]]
  if (is.null(.ret)) {
    if (arg == "posthoc") {
      return(nlmixrPosthoc(obj))
    }
    .cls <- class(obj)
    .env <- attr(.cls, ".foceiEnv")
    .ret <- `$.nlmixrFitCore`(.env, arg, exact)
    if (is.null(.ret)) {
      if (arg == "simInfo") {
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
as.nlme <- function(object, ...) {
  .ret <- object$nlme
  if (is.null(.ret)) stop("Cannot convert to nlme.")
  return(.ret)
}
##' Return composite saem/focei to saem
##'
##' @param x saem/focei object from common UI.
##' @return saem object or NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
as.saem <- function(x) {
  .ret <- x$saem
  if (is.null(.ret)) stop("Cannot convert to saem.")
  return(.ret)
}

##' @importFrom nlme VarCorr
##' @export
VarCorr.nlmixrFitCore <- function(x, sigma = NULL, ...) {
  .ret <- x$nlme
  if (is.null(.ret)) {
    .var <- diag(x$omega)
    .ret <- data.frame(
      Variance = .var, StdDev = sqrt(.var),
      row.names = names(.var)
    )
    .ret <- .ret[!is.na(.ret[, 1]), ]
    return(.ret)
  } else {
    VarCorr(.ret, ...)
  }
}

##' @export
VarCorr.nlmixrFitCoreSilent <- VarCorr.nlmixrFitCore

.sigma <- function(x) {
  .ret <- x$nlme
  if (is.null(.ret)) {
    if (exists("uif", envir = x$env)) {
      .df <- as.data.frame(x$uif$ini)
      .errs <- paste(.df[which(!is.na(.df$err)), "name"])
      return(fixef(x)[.errs])
    }
  } else {
    return(.ret$sigma)
  }
}

##' @export
str.nlmixrFitData <- function(object, ...) {
  NextMethod(object)
  .env <- object$env
  ## cat(" $ par.hist         : Parameter history (if available)\n")
  ## cat(" $ par.hist.stacked : Parameter history in stacked form for easy plotting (if available)\n")
  cat(" $ omega            : Omega matrix\n")
  cat(" $ omegaR           : Omega Correlation matrix\n")
  cat(" $ shrink           : Shrinkage table, includes skewness, kurtosis, and eta p-values\n")
  cat(" $ parFixed         : Fixed Effect Parameter Table\n")
  cat(" $ theta            : Fixed Parameter Estimates\n")
  cat(" $ eta              : Individual Parameter Estimates\n")
  cat(" $ seed             : Seed (if applicable)\n")
  cat(" $ coefficients     : Fixed and random coefficients\n")
  if (exists("uif", envir = object$env)) {
    cat(" $ meta             : Model meta information environment\n")
    cat(" $ modelName        : Model name (from R function)\n")
    cat(" $ dataName         : Name of R data input\n")
    cat(" $ simInfo          : RxODE list for simulation\n")
    cat(" $ sigma            : List of sigma components and their values\n")
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
residuals.nlmixrFitData <- function(object, ..., type = c("ires", "res", "iwres", "wres", "cwres", "cpred", "cres")) {
  return(object[, toupper(match.arg(type))])
}

##' Return the objective function
##'
##' @param x object to return objective function value
##' @param type Objective function type value to retrieve or add.
##'
##' \itemize{
##'
##' \item{focei} For most models you can specify "focei" and it will
##' add the focei objective function.
##'
##' \item{nlme} This switches/chooses the nlme objective function if
##'    applicable.  This objective function cannot be added if it
##'    isn't present.
##'
##' \item{fo} FO objective function value. Cannot be generated
##'
##' \item{foce} FOCE object function value. Cannot be generated
##'
##' \item{laplace#} This adds/retrieves  the Laplace objective function value.
##' The \code{#} represents the number of standard deviations
##' requested when expanding the Gaussian Quadrature.  This can
##' currently only be used with saem fits.
##'
##' \item{gauss#.#} This adds/retrieves the Gaussian Quadrature
##' approximation of the objective function.  The first number is the
##' number of nodes to use in the approximation. The second number is
##' the number of standard deviations to expand upon.
##'
##' }
##'
##' @param ... Other arguments sent to ofv for other methods.
##'
##' @return Objective function value
##'
##' @author Matthew Fidler
##'
##' @export
ofv <- function(x, type, ...) {
  UseMethod("ofv")
}

##' @export
ofv.nlmixrFitData <- function(x, type, ...) {
  if (!missing(type)) setOfv(x, type)
  return(x$ofv)
}

##' @export
logLik.nlmixrFitData <- function(object, ...) {
  .objName <- substitute(object)
  .lst <- list(...)
  if (!is.null(.lst$type)) {
    .new <- setOfv(object, .lst$type)
    .parent <- globalenv()
    .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
      if (.cur == .objName && identical(.parent[[.cur]]$env, object$env)) {
        return(.cur)
      }
      return(NULL)
    }))
    if (length(.bound) == 1) {
      if (exists(.bound, envir = .parent)) {
        assign(.bound, .new, envir = .parent)
      }
    }
    return(get("logLik", .new$env))
  } else {
    return(object$logLik)
  }
}

##' @export
logLik.nlmixrFitCore <- function(object, ...) {
  object$logLik
}

##' @export
nobs.nlmixrFitCore <- function(object, ...) {
  object$nobs
}

##' @export
vcov.nlmixrFitCore <- function(object, ...) {
  object$cov
}
##' This gets the parsed data in the lower-level manner that nlmixr expects.
##'
##' @param object nlmixr Object
##'
##' @export
##'
##' @author Matthew L. Fidler
##' @keywords internal
.nmGetData <- function(object) {
  .uif <- object$uif
  .tmp <- deparse(body(.uif$theta.pars))[-1]
  .tmp <- .tmp[-length(.tmp)]
  return(RxODE::etTrans(object$origData, paste(paste(.tmp, collapse = "\n"), "\n", .uif$rxode), TRUE, TRUE, TRUE))
}

##' @export
getData.nlmixrFitCore <- function(object) {
  object$origData
}

##' @export
ranef.nlmixrFitCore <- function(object, ...) {
  object$ranef
}

##' @export
fixef.nlmixrFitCore <- function(object, ...) {
  object$fixef
}
##' @export
fixef.nlmixrFitCoreSilent <- fixef.nlmixrFitCore

##' @export
ranef.nlmixrFitCoreSilent <- ranef.nlmixrFitCore

##' @export
getData.nlmixrFitCoreSilent <- getData.nlmixrFitCore

##' @export
logLik.nlmixrFitCoreSilent <- logLik.nlmixrFitCore

##' @export
nobs.nlmixrFitCoreSilent <- nobs.nlmixrFitCore

##' @export
vcov.nlmixrFitCoreSilent <- vcov.nlmixrFitCore

##' @export
`$.nlmixrGill83` <- function(obj, arg, exact = FALSE) {
  .ret <- obj[[arg]]
  if (is.null(.ret)) {
    .cls <- class(obj)
    .lst <- attr(.cls, ".nlmixrGill")
    return(.lst[[arg]])
  }
  return(.ret)
}

##' @export
nlmixrPosthoc <- function(x, ...) {
  UseMethod("nlmixrPosthoc")
}

##' @export
nlmixrPosthoc.default <- function(x, ...) {
  .posthoc <- (x$control$maxOuterIterations == 0L & x$control$maxInnerIterations > 0L)
  .posthoc <- ifelse(.posthoc, paste0(ifelse(x$method == "FO",
    ifelse(RxODE::rxIs(x, "nlmixrFitData"),
      paste0(
        " estimation with ", crayon::bold$yellow("FOCE"),
        gsub(rex::rex(any_spaces, "(", anything, ")"), "", x$extra),
        crayon::bold(" posthoc")
      ),
      ""
    ),
    crayon::bold(" posthoc")
  ), " estimation"), " fit")
  return(.posthoc)
}

.fmt3 <- function(name, bound, access, extra = "",
                  on = c(TRUE, TRUE)) {
  if (length(access) == 1) {
    on <- on[1]
  }
  cat(cli::cli_format_method({
    cli::cli_rule(paste0(
      crayon::bold(name), " (", extra,
      paste(crayon::bold$blue(paste0(ifelse(on, crayon::yellow(bound), ""), "$", access)), collapse = " or "), "):"
    ))
  }), sep = "\n")
}

.notesFit <- function(x) {
  .parent <- globalenv()
  .bound2 <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], x)) {
      return(.cur)
    }
    return(NULL)
  }))
  if (length(.bound2) > 0) {
    .bound <- .bound2[order(sapply(.bound2, nchar))]
    .bound <- .bound[1]
  } else {
    .bound <- ""
  }
  .c <- c(paste0(
    "  Covariance Type (", .bound, "$covMethod): ",
    x$covMethod
  ))
  if (is.na(get("objective", x$env))) {
    .c <- c(
      .c,
      "Missing Objective function; Can add by:",
      sprintf(
        " Gaussian/Laplacian Likelihoods: AIC(%s) or %s etc.",
        .bound, .bound, "$objf"
      ),
      sprintf(
        " FOCEi CWRES & Likelihoods: addCwres(%s)",
        .bound
      )
    )
  }
  if (x$message != "") {
    .c <- c(
      .c, paste0("  Minimization message (", .bound, "$message): "),
      paste0("    ", x$message)
    )
    if (x$message == "false convergence (8)") {
      .c <- c(
        .c, "  In an ODE system, false convergence may mean \"useless\" evaluations were performed.",
        "  See https://tinyurl.com/yyrrwkce",
        "  It could also mean the convergence is poor, check results before accepting fit",
        "  You may also try a good derivative free optimization:",
        "    nlmixr(...,control=list(outerOpt=\"bobyqa\"))"
      )
    }
  }
  .c
}

## FIXME fitted?

##' @export
getVarCov.nlmixrFitCore <- function(obj, ...) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .env <- obj
  if (RxODE::rxIs(obj, "nlmixrFitData")) {
    .env <- obj$env
  }
  .force <- FALSE
  .args <- list(...)
  if (!is.null(.args$force)) {
    .force <- .args$force
  }
  if (exists("cov", envir = .env) && !.force) {
    if (RxODE::rxIs(.env$cov, "matrix")) {
      return(.env$cov)
    }
  }
  .pt <- proc.time()
  .control <- .env$control
  ## .control$maxInnerIterations <- 0L;
  .control$maxOuterIterations <- 0L
  .control$boundTol <- 0
  .control$calcTables <- FALSE
  .lst <- list(...)
  if (!is.null(.lst$covMethod)) {
    if (length(.lst$covMethod) == 1) {
      if (.lst$covMethod == "") {
        .control$covMethod <- 0L
      }
    }
    if (RxODE::rxIs(.lst$covMethod, "character")) {
      .lst$covMethod <- match.arg(.lst$covMethod, c("r,s", "r", "s"))
      .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
      .control$covMethod <- .covMethodIdx[.lst$covMethod]
    }
  } else if (.control$covMethod == 0L) {
    .control$covMethod <- 1L
  }
  .lst$covMethod <- NULL
  ## covDerivMethod=c("central", "forward"),
  if (!is.null(.lst$hessEps)) {
    .control$hessEps <- .lst$hessEps
    .lst$hessEps <- NULL
  }
  if (!is.null(.lst$gillKcov)) {
    .control$gillKcov <- .lst$gillKcov
    .lst$gillKcov <- NULL
  }
  if (!is.null(.lst$gillStepCov)) {
    .control$gillStepCov <- .lst$gillStepCov
    .lst$gillStepCov <- NULL
  }
  if (!is.null(.lst$gillFtolCov)) {
    .control$gillFtolCov <- .lst$gillFtolCov
    .lst$gillFtolCov <- NULL
  }
  if (!is.null(.lst$rmatNorm)) {
    .control$rmatNorm <- .lst$rmatNorm
    .lst$rmatNorm <- NULL
  }
  if (!is.null(.lst$smatNorm)) {
    .control$smatNorm <- .lst$smatNorm
    .lst$smatNorm <- NULL
  }
  if (!is.null(.lst$covGillF)) {
    .control$covGillF <- .lst$covGillF
    .lst$covGillF <- NULL
  }
  if (!is.null(.lst$covSmall)) {
    .control$covSmall <- .lst$covSmall
    .lst$covSmall <- NULL
  }
  for (.n in names(.lst)) {
    .control[[.n]] <- .lst[[.n]]
  }
  .dat <- getData(obj)
  .uif <- obj$uif
  .mat <- as.matrix(nlme::random.effects(obj)[, -1])
  .skipCov <- obj$skipCov
  .inits <- list(
    THTA = as.vector(nlme::fixed.effects(obj)),
    OMGA = focei.eta.nlmixrFitCore(obj)
  )
  .fit2 <- foceiFit.data.frame0(
    data = .dat,
    inits = .inits,
    PKpars = .uif$theta.pars,
    ## par_trans=fun,
    model = .uif$rxode.pred,
    pred = function() {
      return(nlmixr_pred)
    },
    err = .uif$error,
    lower = .uif$focei.lower,
    upper = .uif$focei.upper,
    thetaNames = .uif$focei.names,
    etaNames = .uif$eta.names,
    etaMat = .mat,
    skipCov = .skipCov,
    control = .control
  )
  .env$cov <- .fit2$cov
  .env$popDf <- .fit2$popDf
  .env$popDfSig <- .fit2$popDfSig
  .env$covMethod <- .fit2$covMethod
  .updateParFixed(.env)
  .parent <- parent.frame(2)
  .bound <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], obj)) {
      return(.cur)
    }
    return(NULL)
  }))
  message(paste0("Updated original fit object ", ifelse(is.null(.bound), "", crayon::yellow(.bound))))
  .env$time$covariance <- (proc.time() - .pt)["elapsed"]
  return(.env$cov)
}

##' @export
getVarCov.nlmixrFitCoreSilent <- getVarCov.nlmixrFitCore

focei.eta.nlmixrFitCore <- function(object, ...) {
  .uif <- object$uif
  ## Reorder based on translation
  .df <- as.data.frame(.uif$ini)
  .eta <- .df[!is.na(.df$neta1), ]
  .len <- length(.eta$name)
  .curOme <- object$omega
  .curLhs <- character()
  .curRhs <- numeric()
  .ome <- character()
  for (.i in seq_along(.eta$name)) {
    .lastBlock <- FALSE
    if (.i == .len) {
      .lastBlock <- TRUE
    } else if (.eta$neta1[.i + 1] == .eta$neta2[.i + 1]) {
      .lastBlock <- TRUE
    }
    if (.eta$neta1[.i] == .eta$neta2[.i]) {
      .curLhs <- c(.curLhs, sprintf("ETA[%d]", .eta$neta1[.i]))
      .curRhs <- c(.curRhs, .curOme[.eta$neta1[.i], .eta$neta2[.i]])
      if (.lastBlock) {
        .ome[length(.ome) + 1] <- sprintf(
          "%s ~ %s", paste(.curLhs, collapse = " + "),
          paste(deparse(.curRhs), collapse = " ")
        )
        .curLhs <- character()
        .curRhs <- numeric()
      }
    } else {
      .curRhs <- c(.curRhs, .curOme[.eta$neta1[.i], .eta$neta2[.i]])
    }
  }
  .ome <- eval(parse(text = sprintf("list(%s)", paste(.ome, collapse = ","))))
  return(.ome)
}

##' @export
focei.eta.nlmixrFitCoreSilent <- focei.eta.nlmixrFitCore

##' Convert fit to FOCEi style fit
##'
##' @param object Fit object to convert to FOCEi-style fit.
##' @param uif Unified Interface Function
##' @param pt Proc time object
##' @param ... Other Parameters
##' @param data The data to pass to the FOCEi translation.
##' @param calcResid A boolean to indicate if the CWRES residuals
##'     should be calculated
##' @param nobs2 Number of observations without EVID=2
##' @return A FOCEi fit style object.
##' @author Matthew L. Fidler
as.focei <- function(object, uif, pt = proc.time(), ..., data, calcResid = TRUE) {
  UseMethod("as.focei")
}


##' Get the FOCEi theta or eta specification for model.
##'
##' @param object Fit object
##' @param uif User interface function or object
##' @param ... Other parameters
##' @return List for the OMGA list in FOCEi
##' @author Matthew L. Fidler
focei.eta <- function(object, uif, ...) {
  UseMethod("focei.eta")
}

##' Get the FOCEi theta specification for the model
##'
##' @inheritParams focei.eta
##' @return Parameter estimates for Theta
focei.theta <- function(object, uif, ...) {
  UseMethod("focei.theta")
}

##' Cox Box, Yeo Johnson and inverse transformation
##'
##' @param x data to transform
##' @param lambda Cox-box lambda parameter
##' @return Cox-Box Transformed Data
##' @author Matthew L. Fidler
##' @examples
##'
##' boxCox(1:3,1) ## Normal
##' iBoxCox(boxCox(1:3,1))
##'
##' boxCox(1:3,0) ## Log-Normal
##' iBoxCox(boxCox(1:3,0),0)
##'
##' boxCox(1:3,0.5) ## lambda=0.5
##' iBoxCox(boxCox(1:3,0.5),0.5)
##'
##' yeoJohnson(seq(-3,3),1) ## Normal
##' iYeoJohnson(yeoJohnson(seq(-3,3),1))
##'
##' yeoJohnson(seq(-3,3),0)
##' iYeoJohnson(yeoJohnson(seq(-3,3),0),0)
##' @export
boxCox <- function(x, lambda = 1) {
  .Call(`_nlmixr_boxCox_`, x, lambda, 0L)
}

##' @rdname boxCox
##' @export
iBoxCox <- function(x, lambda = 1) {
  .Call(`_nlmixr_iBoxCox_`, x, lambda, 0L)
}

##' @rdname boxCox
##' @export
yeoJohnson <- function(x, lambda = 1) {
  .Call(`_nlmixr_boxCox_`, x, lambda, 1L)
}

##' @rdname boxCox
##' @export
iYeoJohnson <- function(x, lambda = 1) {
  .Call(`_nlmixr_iBoxCox_`, x, lambda, 1L)
}

.setSaemExtra <- function(.env, type) {
  .uif <- .env$uif
  .txt <- paste0("(", crayon::italic(ifelse(is.null(.uif$nmodel$lin.solved), ifelse(.uif$predSys, "PRED", "ODE"), "Solved")), "); ")
  if (type == "FOCEi") {
    .txt <- paste0(.txt, crayon::blurred$italic("OBJF by FOCEi approximation"))
  } else if (type == "") {
    .txt <- paste0(.txt, crayon::blurred$italic("OBJF not calculated"))
  } else {
    .reg <- rex::rex(start, "laplace", capture(.regNum), end)
    .regG <- rex::rex(start, "gauss", capture(.regNum), "_", capture(.regNum), end)
    if (regexpr(.reg, type, perl = TRUE) != -1) {
      .nnode <- 1
      .nsd <- as.numeric(sub(.reg, "\\1", type, perl = TRUE))
    } else if (regexpr(.regG, type, perl = TRUE) != -1) {
      .nnode <- as.numeric(sub(.regG, "\\1", type, perl = TRUE))
      .nsd <- as.numeric(sub(.regG, "\\2", type, perl = TRUE))
    } else {
      stop("unknown error")
    }
    .txt <- paste0(.txt, crayon::blurred$italic(sprintf("OBJF by %s", paste0(ifelse(.nnode == 1, "Lapalcian (n.sd=", sprintf("Gaussian Quadrature (n.nodes=%s, n.sd=", .nnode)), .nsd, ")"))))
  }
  .env$extra <- .txt
  return(invisible(.txt))
}

##' Set Objective function type for a nlmixr object
##'
##' @param x nlmixr fit object
##' @param type Type of objective function to use for AIC, BIC, and
##'     $objective
##' @return Nothing
##' @author Matthew L. Fidler
##' @export
setOfv <- function(x, type) {
  if (inherits(x, "nlmixrFitCore") || inherits(x, "nlmixrFitCoreSilent")) {
    .objDf <- x$objDf
    .w <- which(tolower(row.names(.objDf)) == tolower(type))
    if (length(.w) == 1) {
      .env <- x$env
      .objf <- .objDf[.w, "OBJF"]
      .lik <- .objDf[.w, "Log-likelihood"]
      attr(.lik, "df") <- attr(get("logLik", .env), "df")
      attr(.lik, "nobs") <- attr(get("logLik", .env), "nobs")
      class(.lik) <- "logLik"
      .bic <- .objDf[.w, "BIC"]
      .aic <- .objDf[.w, "AIC"]
      assign("OBJF", .objf, .env)
      assign("objf", .objf, .env)
      assign("objective", .objf, .env)
      assign("logLik", .lik, .env)
      assign("AIC", .aic, .env)
      assign("BIC", .bic, .env)
      if (!is.null(x$saem)) {
        .setSaemExtra(.env, type)
      }
      invisible(x)
    } else {
      if (tolower(type) == "focei") {
        .tmp <- addCwres(x, TRUE, globalenv())
        return(setOfv(.tmp, "FOCEi"))
      } else if (!is.null(x$saem)) {
        .ret <- x$saem
        .reg <- rex::rex(start, "laplace", capture(.regNum), end)
        .regG <- rex::rex(start, "gauss", capture(.regNum), "_", capture(.regNum), end)
        if (regexpr(.reg, type, perl = TRUE) != -1) {
          .nnode <- 1
          .nsd <- as.numeric(sub(.reg, "\\1", type, perl = TRUE))
        } else if (regexpr(.regG, type, perl = TRUE) != -1) {
          .nnode <- as.numeric(sub(.regG, "\\1", type, perl = TRUE))
          .nsd <- as.numeric(sub(.regG, "\\2", type, perl = TRUE))
        } else {
          stop(sprintf("Cannot switch objective function to '%s' type.", type))
        }
        .likTime <- proc.time()
        .saemObf <- calc.2LL(x$saem, nnodes.gq = .nnode, nsd.gq = .nsd)
        .likTime <- proc.time() - .likTime
        .likTime <- .likTime["elapsed"]
        .env <- x$env
        .time <- .env$time
        if (any(names(.time) == "logLik")) {
          .time$logLik <- .time$logLik + .likTime
        } else {
          .time <- data.frame(.time, logLik = .likTime, check.names = FALSE)
        }
        .env$time <- .time
        .llik <- -.saemObf / 2
        .nobs <- .env$nobs
        attr(.llik, "df") <- attr(get("logLik", .env), "df")
        .objf <- ifelse(.env$adjObj, .saemObf - .nobs * log(2 * pi), .saemObf)
        .tmp <- data.frame(
          OBJF = .objf, AIC = .saemObf + 2 * attr(get("logLik", .env), "df"),
          BIC = .saemObf + log(.env$nobs) * attr(get("logLik", .env), "df"),
          "Log-likelihood" = as.numeric(.llik), check.names = FALSE
        )
        if (any(names(.env$objDf) == "Condition Number")) {
          .cn <- unique(.env$objDf[["Condition Number"]])
          .cn <- .cn[!is.na(.cn)]
          if (length(.cn) == 1) {
            .tmp <- data.frame(.tmp, "Condition Number" = .cn, check.names = FALSE)
          } else {
            .tmp <- data.frame(.tmp, "Condition Number" = NA, check.names = FALSE)
          }
        }
        .rn <- row.names(.env$objDf)
        .env$objDf <- rbind(
          .env$objDf,
          .tmp
        )
        row.names(.env$objDf) <- c(.rn, type)
        .env$objDf <- .env$objDf[order(row.names(.env$objDf)), ]
        .env$objDf <- .env$objDf[!is.na(.env$objDf$OBJF), ]
        return(setOfv(x, type))
      }
      stop(sprintf("Cannot switch objective function to '%s' type.", type))
    }
  } else {
    stop("Wrong type of object.")
  }
}

##' @importFrom utils capture.output
.captureOutput <- function(expr, envir = parent.frame()) {
  eval(
    {
      .file <- rawConnection(raw(0L), open = "w")
      on.exit({
        if (!is.null(.file)) close(.file)
      })
      capture.output(expr, file = .file)
      .ret <- rawConnectionValue(.file)
      close(.file)
      .file <- NULL
      .ret <- rawToChar(.ret)
      return(.ret)
    },
    envir = envir,
    enclos = envir
  )
}

##  LocalWords:  covariance
