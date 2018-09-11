##' VPC based on ui model
##'
##' @param fit nlmixr fit object
##' @param data this is the data to use to augment the VPC fit.  By
##'     default is the fitted data, (can be retrieved by
##'     \code{\link[nlme]{getData}}), but it can be changed by specifying
##'     this argument.
##' @param n Number of VPC simulations.  By default 100
##' @inheritParams vpc::vpc
##' @inheritParams RxODE::rxSolve
##' @param ... Args sent to \code{\link[RxODE]{rxSolve}}
##' @return Simulated dataset (invisibly)
##' @author Matthew L. Fidler
##' @export
vpc_ui <- function(fit, data=NULL, n=100, bins = "jenks",
                   n_bins = "auto", bin_mid = "mean",
                   show = NULL, stratify = NULL, pred_corr = FALSE,
                   pred_corr_lower_bnd = 0, pi = c(0.05, 0.95), ci = c(0.05, 0.95),
                   uloq = NULL, lloq = NULL, log_y = FALSE, log_y_min = 0.001,
                   xlab = NULL, ylab = NULL, title = NULL, smooth = TRUE, vpc_theme = NULL,
                   facet = "wrap", labeller = NULL, vpcdb = FALSE, verbose = FALSE, ...){
    if (is(data, "numeric") | is(data, "integer")){
        if (missing(n)){
            n <- data;
            data <- NULL;
        } else {
            stop("Data needs to be a data frame instead of a numeric value.")
        }
    }
    .xtra <- list(...);
    if (is.null(.xtra$nsim)){
        .xtra$nsim <- n;
    }
    if (inherits(fit, "nlmixrVpc")){
        sim <- fit;
    } else {
        .xtra$object <- fit;
        ## .xtra$returnType <- "data.frame";
        .xtra$returnType <- "rxSolve";
        pt <- proc.time();
        sim <- do.call("nlmixrSim", .xtra);
        sim0 <- sim;
        sim <- sim[, c("sim.id", "id", "time", "sim")]
        names(sim)[4] <- "dv";
        if (is.null(data)){
            dat <- nlmixrData(getData(fit));
        } else {
            dat <- data
        }
        diff <- proc.time() - pt;
        message(sprintf("done (%.2f sec)", diff["elapsed"]));
        onames <- names(dat)
        names(dat) <- tolower(onames)
        w <- which(duplicated(names(dat)));
        if (length(w) > 0){
            warning(sprintf("Dropping duplicate columns (case insensitive): %s", paste(onames, collapse=", ")))
            dat <- dat[, -w];
        }
        if (!is.null(stratify)){
            cols <- c(tolower(stratify), "dv")
            stratify <- tolower(stratify);
        }  else {
            cols <- c("dv");
        }

        dat <- dat[dat$evid == 0, ];
        ## Assume this is in the observed dataset. Add it to the current dataset
        if(!all(names(sim) %in% cols)){
            w <- cols[!(cols %in% names(sim))]
            if (length(w) >= 1){
                n <- names(sim)
                sim <- cbind(sim, dat[, w, drop = FALSE]);
                names(sim) <- c(n, w);
            }
        }
        sim <- list(rxsim=sim0, sim=sim, obs=dat)
        attr(sim, "nsim") <- .xtra$nsim;
        class(sim) <- "nlmixrVpc";
    }
    ns <- loadNamespace("vpc");
    if (exists("vpc_vpc",ns)){
        vpcn <- "vpc_vpc"
    } else {
        vpcn <- "vpc"
    }
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    call <- call[names(call) %in% methods::formalArgs(getFromNamespace(vpcn,"vpc"))]
    call$obs_cols = list(id="id", dv="dv", idv="time")
    call$sim_cols = list(id="id", dv="dv", idv="time")
    call$stratify = stratify
    p = do.call(getFromNamespace(vpcn,"vpc"), c(sim, call), envir = parent.frame(1))
    print(p);
    sim$gg <- p;

    return(invisible(sim));
}

##'@export
print.nlmixrVpc <- function(x, ...){
    cat(sprintf("nlmixr vpc object of %d simulations.\n", attr(x, "nsim")))
    cat("  $sim = simulated data\n")
    cat("  $obs = observed data\n")
    cat("  $gg = vpc ggplot\n")
    cat("use vpc(...) to change plot options\n")
}

##' @rdname vpc_ui
##' @export
vpc.nlmixrFitData <- function(sim, ...){
    vpc_ui(fit=sim, ...);
}

##' @rdname vpc_ui
##' @export
vpc.nlmixrVpc <- function(sim, ...){
    vpc_ui(fit=sim, ...);
}

##' @rdname vpc_ui
##' @S3method vpc ui
##' @export vpc.ui
vpc.ui <- function(sim, ...){
    vpc_ui(fit=sim, ...);
}
