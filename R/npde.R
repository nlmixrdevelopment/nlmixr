##' NPDE calculation for nlmixr
##'
##' @param object nlmixr fit object
##' @param nsim Number of simulations.  By default this is 300
##' @param ties When TRUE, the npde distribution can have ties.  When
##'     FALSE, the npde distribution uses uniform random numbers to
##'     prevent ties.
##' @param seed Seed for running nlmixr simulation.  By default 1009
##' @param updateObject Boolean indicating if original object should be updated.  By default this is TRUE.
##' @param ... Other ignored parameters.
##' @return New nlmixr fit object
##' @author Matthew L. Fidler
##' @export
addNpde <- function(object, nsim=300, ties=TRUE, seed=1009, updateObject=TRUE, ...){
    .objName <- substitute(object);
    if (any(names(object) == "NPDE")){
        warning("Already contains NPDE")
        return(object)
    }
    set.seed(seed);
    .si <- object$simInfo
    .rx <- .si$rx;
    .rx <- gsub(rex::rex(capture("ipred"), or("=", "~"),  except_any_of("\n;"), any_of("\n;")), "", .rx)
    .rx <- gsub(rex::rex("sim", or("=", "~"), "rxTBSi(", capture(except_any_of(",)")), ",", anything, any_of("\n;")),
                "sim=\\1", .rx)
    .si$rx <- .rx
    .dat <- nlmixrData(getData(fit))
    .dat <- .dat[.dat$EVID == 0, ]
    .si$object <- object;
    .si$returnType <- "data.frame.TBS";
    .si$nsim <- nsim;
    .si <- c(.si, list(...))
    .si$modelName <- "NPDE"
    .pt <- proc.time();
    .sim <- do.call("nlmixrSim", .si);
    .dv <- object$DV;
    .dvl <- length(.dv)
    .cls <- class(object)
    .new <- cbind(object, .Call(`_nlmixr_npde`, object$ID, .dv, .sim$sim, .sim$rxLambda, .sim$rxYj, ties=ties))
    class(.new) <- .cls;
    if (updateObject){
        .parent <- parent.frame(2);
        .bound <- do.call("c", lapply(ls(.parent), function(.cur){
                                   if (.cur == .objName && identical(.parent[[.cur]], object)){
                                       return(.cur)
                                   }
                                   return(NULL);
                               }))
        assign(.bound, .new, envir=.parent)
    }
    return(.new)
}
