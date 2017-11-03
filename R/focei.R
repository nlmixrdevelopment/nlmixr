regEta <- rex::rex(start, "ETA[", capture("1":"9", any_of("0":"9")), "]")
rxFoceiEtaSetup <- function(object, ..., dv, eta, theta, nonmem=FALSE, table=TRUE, inv.env=parent.frame(1), id= -1,
                            inits.vec=NULL, atol.outer=1e-8, rtol.outer=1e-6,
                            switch.solver=FALSE, pred.minus.dv=TRUE, scale.to=1,
                            numeric=FALSE){
    args <- list(object=object, ..., eta=eta, theta=theta);
    args$do.solve <- FALSE;
    setup <- c(do.call(object$solve, args), as.list(inv.env));
    return(with(setup, {
        tmp <- environment(object$.call);
        if (any(names(tmp) == "eta.trans")){
            eta.trans <- tmp$eta.trans;
        } else {
            n <- names(params)
            max.eta <- n[regexpr(regEta, n) != -1]
            max.eta <- max(as.numeric(gsub(regEta, "\\1", max.eta)));
            eta.trans <- sapply(seq(1, max.eta), function(x){return(which(n == sprintf("ETA[%s]", x)) - 1)});
            tmp$eta.trans <- eta.trans;
        }
        setup$nonmem <- as.integer(nonmem)
        setup$table <- as.integer(table)
        setup$DV <- as.double(dv);
        setup$eta.trans <- as.integer(eta.trans);
        setup$object <- object;
        setup$eta <- eta;
        setup$eta.mat <- matrix(eta, ncol=1);
        setup$neta <- as.integer(length(eta))
        setup$theta <- theta;
        setup$ntheta <- as.integer(length(theta))
        setup$id <- as.integer(id);
        setup$atol.outer <- atol.outer;
        setup$rtol.outer <- rtol.outer;
        setup$nearPD <- rxNearPd;
        setup$rc <- 0L;
        setup$switch.solver <- as.integer(switch.solver);
        setup$pred.minus.dv <- as.integer(pred.minus.dv);
        setup$numeric <- as.integer(numeric)
        if (!is.null(inits.vec)){
            setup$inits.vec <- inits.vec;
        }
        setup$scale.to <- scale.to;
        return(list2env(setup));
    }))
}

rxNearPd <- function(mat, env){
    if (any(is.na(mat))){
        cat("Bad matrix:\n");
        print(mat);
        env$reset <- 1;
        return(mat)
    } else {
        return(as.matrix(Matrix::nearPD(mat)$mat));
    }
}

##' FOCEI ETA setup,
##'
##' This is basically for testing
##'
##' @param object RxODE object
##' @param ... values sent to rxFoceiEtaSetup
##' @param dv dependent variable
##' @param eta eta values.
##' @param env Environment to use, instead of rxFoceiEtaSetup environment
##' @param nonmem Match NONMEMs approximation of the hessian.
##' @return environment of solved information.
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxFoceiEta <- function(object, ..., dv, eta,  env, nonmem=FALSE){
    if (class(object) == "rxFocei"){
        object <- object$inner;
    }
    if (missing(env)){
        args <- c(as.list(match.call(expand.dots=TRUE))[-1], list(do.solve=FALSE));
        args[[1]] <- object;
        env <- do.call(getFromNamespace("rxFoceiEtaSetup", "nlmixr"), args, envir = parent.frame(1));
    }
    object$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
    rxInner(eta, env);
    return(env);
}

##' Get -LLIK for individual
##'
##' Used for testing.
##'
##' @inheritParams rxFoceiEta
##' @return -llik
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiLik <- function(object, ..., dv, eta){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    if (class(object) == "rxFocei"){
        object <- object$inner;
        args[[1]] <- object;
    }
    object$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "nlmixr"), args, envir = parent.frame(1));
    return(RxODE_focei_eta_lik(eta, env))
}

##' Get -Grad for individual
##'
##' For testing
##' @inheritParams rxFoceiEta
##' @return -grad
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiLp <- function(object, ..., dv, eta){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    if (class(object) == "rxFocei"){
        object <- object$inner;
        args[[1]] <- object;
    }
    object$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "nlmixr"), args, envir = parent.frame(1));
    return(RxODE_focei_eta_lp(eta, env));
}

##' Calculate the Gradient for an individual
##'
##' This can be based on a combination of numerical differention and
##' symbolic differentiation or purely numerical differentiation.
##'
##' @param object focei object
##' @param ret individual return object to add "grad" attribute to...
##' @inheritParams rxFoceiEta
##' @keywords internal
##' @export
rxFoceiGrad <- function(object, ret, ..., theta, eta=NULL, dv,
                        numDeriv.method="Richardson"){
    grad <- attr(ret,"grad");
    if (!is.null(grad) && length(grad) > 0){
        return(ret);
    } else if (!is.null(object$theta)) {
        ## Use Almquist's method BUT the dH/dTheta is numerically caluclated...
        theta.rxode <- object$theta;
        inner.rxode <- object$inner;
        if (is.null(eta)){
            eta <- as.vector(attr(ret, "posthoc"))
        }
        ome.28 <- as.vector(attr(ret,"omega.28"))
        c.hess <- attr(ret,"c.hess")
        ## if (length(theta) == length(rxParams(theta.rxode)) - length(eta)){
        ##     new.theta <- theta;
        ##     rest.theta <- c();
        ## } else {
        ##     new.theta <- theta[seq(1, length(theta) - length(ome.28))]
        ##     rest.theta <- theta[-seq(1, length(theta) - length(ome.28))]
        ## }
        ## print(new.theta)
        ## print(rest.theta)
        args <- as.list(match.call(expand.dots=TRUE))[-1];
        args$dv <- dv
        args$theta <- theta
        args$eta <- eta
        args$object <- theta.rxode;
        env <- do.call(getFromNamespace("rxFoceiEtaSetup", "nlmixr"), args, envir = parent.frame(1));
        theta.rxode$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
        rxGrad(env);
        args$object <- object
        ## Get dH/dTheta by numeric differentiation...
        ##inner.rxode$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
        args$return.env <- TRUE
        args$c.hess <- c.hess
        fun <- function(theta){
            args$theta <- theta## c(theta, rest.theta) ;
            args$estimate <- FALSE
            env2 <- do.call(getFromNamespace("rxFoceiInner", "nlmixr"), args, envir = parent.frame(1));
            rxHessian(env2);
            return(as.vector(env2$H));
        }
        H <- fun(theta)
        Hinv <- RxODE::rxInv(matrix(H, length(args$eta)));
        jac <- numDeriv::jacobian(fun, theta, method=numDeriv.method);
        gr <- sapply(seq_along(env$lp), function(x){
            return(-env$lp[x] - 0.5 * sum(diag(Hinv %*% matrix(jac[, x], length(args$eta)))))
        })
        gr <- c(gr, ome.28)
        if (any(ls(env) == "inits.vec") && !is.null(args$scale.to)){
            gr <- gr / (env$inits.vec / env$scale.to);
        }
        attr(ret,"grad") <- gr;
        return(ret);
    } else {
        ## Use brute force
        inner.rxode <- object$inner;
        if (is.null(eta)){
            eta <- as.vector(attr(ret, "posthoc"))
        }
        ome.28 <- as.vector(attr(ret,"omega.28"))
        c.hess <- attr(ret,"c.hess")
        args <- as.list(match.call(expand.dots=TRUE))[-1];
        args$object <- object;
        args$dv <- dv
        args$theta <- theta
        args$eta <- eta
        args$c.hess <- c.hess;
        args$add.grad <- FALSE
        func <- function(theta){
            args$theta <- theta;
            ret <- do.call(getFromNamespace("rxFoceiInner","nlmixr"), args)
        }
        gr <- c(numDeriv::grad(func, theta, method=numDeriv.method), attr(ret, "omega.28"))
        if (any(names(args) == "inits.vec") && !is.null(args$scale.to)){
            gr <- gr / (args$inits.vec / args$scale.to);
        }
        attr(ret,"grad") <- gr;
        return(ret);
    }
}

##' Solve the FOCEI inner problem
##'
##' @param object RxODE object
##' @param ... values sent to rxFoceiEtaSetup and lbfgs
##' @param dv dependent variable
##' @param eta eta value
##' @param c.hess Hessian curvature matrix in compressed form (Used for n1qn1)
##' @param eta.bak backup eta value in case eta estimation fails
##' @param estimate Boolean indicating if the optimization should be
##'     perfomed(TRUE), or just keep the eta, and calculate the Loglik
##'     with fitted/posthoc attributes(FALSE).
##' @param inner.opt Inner optimization method; Either n1qn1 or lbfgs
##' @param return.env Return the environment instead of the llik.
##' @param add.grad Add A gradient attribute, if not present.
##' @return Loglik with fitted and posthoc attributes
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiInner <- function(object, ..., dv, eta, c.hess=NULL, eta.bak=NULL,
                         estimate=TRUE, inner.opt=c("n1qn1", "lbfgs"), return.env=FALSE,
                         add.grad=FALSE){
    inner.opt <- match.arg(inner.opt);
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$dv <- dv
    args$eta <- eta
    args$eta.bak <- eta.bak
    args$estimate <- estimate
    if (!any(names(args) == "numeric")){
        args$numeric <- FALSE;
    }
    if (args$numeric){
        inner.rxode <- object$pred.only;
        inner.rxode$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
    } else {
        inner.rxode <- object$inner;
        inner.rxode$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
    }
    args$object <- inner.rxode;
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "nlmixr"), args, envir = parent.frame(1));

    lik <- RxODE_focei_eta("lik");
    lp <- RxODE_focei_eta("lp")
    c.hess <- NULL
    est <- function(){
        if (estimate){
            args$call_eval <- lik;
            args$call_grad <- lp;
            args$vars <- args$eta;
            args$environment <- env
            env$eta <- args$eta;
            est0 <- function(){
                if (inner.opt == "n1qn1"){
                    if (is.null(c.hess)){
                        args$c.hess <- c.hess;
                    }
                    ## output <- do.call(getFromNamespace("n1qn1","n1qn1"), args, envir = parent.frame(1))
                    output <- try(do.call(getFromNamespace("n1qn1","n1qn1"), args, envir = parent.frame(1)), silent=TRUE)
                    if (!inherits(output, "try-error")){
                        env$c.hess <- output$c.hess
                        c.hess <<- output$c.hess;
                    }
                } else {
                    output <- try(do.call(getFromNamespace("lbfgs","lbfgs"), args, envir = parent.frame(1)), silent=TRUE)
                }
                return(output);
            }
            output <- est0();
            if (inherits(output, "try-error")){
                args$vars <- rep(0, length(args$eta));
                output <- est0();
                if (inherits(output, "try-error")){
                    if (env$numeric == 1){
                        output <- rxInnerNum(rep(0, length(args$eta)), env);
                    } else {
                        output <- rxInner(rep(0, length(args$eta)), env);
                    }
                }
            }
            return(output);
        } else {
            ret <- try(R.utils::captureOutput(rxInner(args$eta, env)));
            ## ret <- rxInner(args$eta, env);
            ## print(as.list(env))
            if (inherits(ret, "try-error")){
                pred.only <- object$pred.only;
                pred.only$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
                ret <- R.utils::captureOutput(rxInnerNum(args$eta, env)); ## Use finite difference instead.
            }
            return(ret);
        }
    }
    est();
    if (estimate){
        if (any(is.na(env$eta))){
            if (!is.null(args$eta.bak)){
                args$eta <- args$eta.bak
                est()
                if (any(is.na(env$eta))){
                    args$eta <- rep(0, length(args$eta));
                    env <- do.call(getFromNamespace("rxFoceiEta", "nlmixr"), args, envir = parent.frame(1));
                }
            } else {
                args$eta <- rep(0, length(args$eta))
                est();
                if (any(is.na(env$eta))){
                    args$eta <- rep(0, length(args$eta))
                    env <- do.call(getFromNamespace("rxFoceiEta", "nlmixr"), args, envir = parent.frame(1));
                }
            }
        }
        if (any(abs(env$eta) > 1e4)){
            env <- do.call(getFromNamespace("rxFoceiEta", "nlmixr"), args, envir = parent.frame(1));
            if (any(abs(env$eta) > 1e4)){
                args$eta <- rep(0, length(args$eta))
                est();
                if (any(is.na(env$eta))){
                    args$eta <- rep(0, length(args$eta))
                    env <- do.call(getFromNamespace("rxFoceiEta", "nlmixr"), args, envir = parent.frame(1));
                }
            }
        }
    }
    if (return.env){
        return(env);
    }
    if (is.null(object$outer)){
        ret <- try(RxODE_focei_finalize_llik(env), silent=TRUE);
        if ((attr(ret, "corrected") == 1) || inherits(ret, "try-error")){
            if (estimate){
                RxODE::rxCat(sprintf("Warning: Problem with Hessian or ETA estimate, resetting ETAs to 0 (ID=%s).\n", env$id));
                args$eta <- rep(0, length(env$eta));
                env$eta <- args$eta;
                args$orthantwise_end <- length(args$eta);
                est();
                if (any(is.na(env$eta))){
                    warning("ETA estimate failed; Assume ETA=0");
                    args$eta <- rep(0, length(args$eta))
                    env <- do.call(getFromNamespace("rxFoceiEta", "nlmixr"), args, envir = parent.frame(1));
                }
                if (any(env$eta > 1e4)){
                    warning("ETA estimate overflow; Assume ETA");
                    env <- do.call(getFromNamespace("rxFoceiEta", "nlmixr"), args, envir = parent.frame(1));
                }
                ret <- try(RxODE_focei_finalize_llik(env), silent=TRUE)
                if (inherits(ret, "try-error")){
                    warning("ETA estimate failed; Assume ETA=0");
                    args$eta <- rep(0, length(args$eta))
                    env <- do.call(getFromNamespace("rxFoceiEta", "nlmixr"), args, envir = parent.frame(1));
                    ret <- RxODE_focei_finalize_llik(env)
                }
                attr(ret, "corrected") <- 1L;
            } else {
                if (inherits(ret, "try-error")){
                    cat("Hessian:\n");
                    print(env$H)
                    cat("log.det.H.neg.5:\n");
                    print(env$log.det.H.neg.5)
                    cat("log.det.OMGAinv.5:\n");
                    print(env$log.det.OMGAinv.5)
                    ## Try with numeric differences instead.
                    pred.only <- object$pred.only;
                    pred.only$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
                    rxInnerNum(args$eta, env); ## Use finite difference instead.
                    ret <- RxODE_focei_finalize_llik(env)
                }
            }
        }
        if (add.grad){
            args <- as.list(match.call(expand.dots=TRUE))[-1];
            args$eta <- NULL;
            args$ret <- ret;
            ret <- do.call(getFromNamespace("rxFoceiGrad", "nlmixr"), args, envir = parent.frame(1))
        }
        return(ret);
    } else {
        args$object <- object;
        args$eta <- env$eta;
        args$theta <- env$theta;
        args$orthantwise_end <- length(args$eta);
        ## print(env$id)
        ## print(args$eta);
        ## print(eta);
        ## print(output);
        env <- do.call(getFromNamespace("rxFoceiTheta", "nlmixr"), args, envir = parent.frame(1));
        if (env$reset == 1){
            cat(sprintf("Warning: Problem with Hessian or ETA estimate, resetting ETAs to 0 (ID=%s).\n", env$id));
            inner.rxode$assignPtr();
            args$eta <- rep(0, length(env$eta));
            args$orthantwise_end <- length(args$eta);
            args$object <- inner.rxode;
            env <- do.call(getFromNamespace("rxFoceiEtaSetup", "nlmixr"), args, envir = parent.frame(1));
            est();
            args$object <- object;
            args$eta <- env$eta;
            args$orthantwise_end <- length(args$eta);
            env <- do.call(getFromNamespace("rxFoceiTheta", "nlmixr"), args, envir = parent.frame(1));
            if (env$reset == 1){
                print(env$par)
                stop("Cannot reset problem.");
            }
            ret <- env$ret;
            attr(ret, "corrected") <- 1L;
        }
        return(env$ret);
    }
}


##' FOCEI THETA setup,
##'
##' This is basically for testing
##'
##' @param object RxODE object
##' @param ... values sent to rxFoceiEtaSetup
##' @param dv dependent variable
##' @param eta eta values.
##' @param omegaInv Inverse omega
##' @param env Environment to use, instead of rxFoceiEtaSetup environment
##' @param nonmem Match NONMEMs approximation of the hessian.
##' @return environment of solved information.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiTheta <- function(object, ..., dv, eta, env, nonmem=FALSE){
    if (class(object) == "rxFocei"){
        object <- object$outer;

    }
    if (missing(env)){
        args <- c(as.list(match.call(expand.dots=TRUE))[-1], list(do.solve=FALSE));
        args[[1]] <- object;
        env <- do.call(getFromNamespace("rxFoceiEtaSetup", "nlmixr"), args, envir = parent.frame(1));
    }
    object$assignPtr(); ## Assign the ODE pointers (and Jacobian Type)
    env$atol.inner <- env$atol;
    env$rtol.inner <- env$rtol;
    env$atol <- env$atol.outer;
    env$rtol <- env$rtol.outer;
    on.exit({
        env$atol <- env$atol.inner;
        env$rtol <- env$rtol.inner;
    })
    rxOuter(env)
    return(env);
}
