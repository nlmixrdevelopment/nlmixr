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
    if (inherits(fit, "nlmixrVpc")){
        sim <- fit;
    } else {
        .xtra$object <- fit;
        ## .xtra$returnType <- "data.frame";
        .xtra$returnType <- "rxSolve";
        .xtra$modelName <- "VPC"
        pt <- proc.time();
        if (is.null(data)){
            dat <- nlmixrData(getData(fit));
        } else {
            dat <- data
        }
        .xtra$nStud <- n;
        if (!is.null(.xtra$nsim)){
            .xtra$nStud <- .xtra$nsim
            .xtra$nsim <- NULL;
        }
        .xtra$dfObs <- 0
        .xtra$dfSub <- 0
        .xtra$thetaMat <- NA
        sim <- do.call("nlmixrSim", .xtra);
        sim0 <- sim;
        sim <- sim[, c("id", "time", "sim")]
        names(sim)[3] <- "dv";
        .si <- fit$simInfo;
        if (pred_corr){
            .xtra.prd <- .xtra;
            .xtra.prd$modelName <- "Pred (for pcVpc)"
            .xtra.prd$params <- c(.si$params, setNames(rep(0, dim(.si$omega)[1]), dimnames(.si$omega)[[2]]),
                                  setNames(rep(0, dim(.si$sigma)[1]), dimnames(.si$sigma)[[2]]))
            .xtra.prd$omega <- NA
            .xtra.prd$sigma <- NA
            .xtra.prd$returnType <- "data.frame";
            .xtra.prd$nStud <- 1;
            .xtra.prd$nsim <- NULL;
            sim2 <- do.call("nlmixrSim", .xtra.prd);
            sim$pred <- sim2$sim
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
        if (any(names(sim)=="cmt") && any(names(fit)=="CMT")){
            if (is(fit$CMT, "factor")){
                sim$cmt  <- factor(sim$cmt,sort(unique(sim$cmt)), labels=levels(fit$CMT))
            }
        }
        if (any(names(dat)=="cmt") && any(names(fit)=="CMT")){
            if (is(fit$CMT, "factor")){
                dat$cmt  <- factor(dat$cmt,sort(unique(dat$cmt)), labels=levels(fit$CMT))
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
    p = do.call(getFromNamespace(vpcn, "vpc"), c(sim, call), envir = parent.frame(1))
    print(p);
    sim$gg <- p;
    return(invisible(sim));
}

##'@export
print.nlmixrVpc <- function(x, ...){
    cat(sprintf("nlmixr vpc object of %d simulations.\n", attr(x, "nsim")))
    cat("  $rxsim = original simulated data\n")
    cat("  $sim = merge simulated data\n")
    cat("  $obs = observed data\n")
    cat("  $gg = vpc ggplot\n")
    cat("use vpc(...) to change plot options\n")
}

##'@export
plot.nlmixrVpc <- function(x, ...){
    return(x$gg)
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
