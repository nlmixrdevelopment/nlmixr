## TODO:
## Lincmt infusion checks...
## High Covariance in Omgea for nlme
## Multiple endpoints check.
## Initial conditions between methods -- OK
## Dots in variable names (especially THETAs for SAEM)
## Fixing components?

## Unified UI observations
## - nlme
## - Initial conditions for residuls.
## - Error structure
## -
## Check on unified SD for SAEM
## Check on output for OMEGA (is is Var or SD)

## Also need to add check for types of omega blocks

## SAEM:
## - Covariate part

##' Prepares the UI function and returns a list.
##'
##' @param fun UI function
##' @return nlmixr UI function
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrUI <- function(fun){
    lhs0 <- nlmixrfindLhs(body(fun))
    dum.fun <- function(){return(TRUE)}
    env.here <- environment(dum.fun)
    env <- new.env(parent=.GlobalEnv)
    assign("fun", fun, env)
    fun2 <- attr(fun,"srcref");
    if (is.null(fun2)){
        stop("option \"keep.source\" must be TRUE for nlmixr models.")
    }
    fun2 <- as.character(fun2, useSource=TRUE)
    rg <- rex::rex("function", any_spaces, "(", anything, ")")
    w <- which(regexpr(rg, fun2) != -1);
    if (length(w) > 0){
        w <- w[1];
        fun2[w] <- sub(rg, "", fun2[w]);
    }
    fun2 <- gsub(rex::rex(boundary, "ini(",any_spaces,"{"), "ini <- function()({", fun2)
    fun2 <- gsub(rex::rex(boundary, "model(",any_spaces,"{"), "model <- function()({", fun2)
    if (fun2[length(fun2)] != "}"){
        fun2[length(fun2)] <- sub(rex::rex("}", end), "", fun2[length(fun2)]);
        fun2[length(fun2) + 1] <- "}"
    }
    w <- which(regexpr(rex::rex(start, any_spaces, "#", anything), fun2) != -1);
    if (length(w) > 0 && all(lhs0 != "desc")){
        w2 <- w[1];
        if (length(w) > 1){
            for (i in 2:length(w)){
                if (w[i] - 1 == w[i - 1]){
                    w2[i] <- w[i];
                } else {
                    break;
                }
            }
        }
        desc <- paste(gsub(rex::rex(any_spaces, end), "", gsub(rex::rex(start, any_spaces, any_of("#"), any_spaces), "", fun2[w2])), collapse=" ");
        lhs0 <- c(lhs0, "desc");
    }
    model <- ini <- NULL; ## Rcheck hack
    eval(parse(text=paste(fun2, collapse="\n"), keep.source=TRUE))
    ini <- nlmixrBounds(ini);
    meta <- list();
    for (var in lhs0){
        if (!any(var == c("ini", "model")) &&
            exists(var, envir=env.here)){
            meta[[var]] <- get(var, envir=env.here);
        }
    }
    if (inherits(ini, "try-error")){
        stop("Error parsing initial estimates.")
    }
    fun2 <- try(nlmixrUIModel(model,ini,fun));
    if (inherits(fun2, "try-error")){
        stop("Error parsing model.")
    }
    class(fun2) <- "nlmixrUI";
    var.ini <- c(fun2$focei.names, fun2$eta.names)
    var.def <- fun2$all.vars;
    diff <- setdiff(var.ini, var.def);
    if (length(diff) > 0){
        stop(sprintf("Model error: initial estimates provided without variables being used: %s", paste(diff, collapse=", ")))
    }
    ns <- fun2$name[which(!is.na(fun2$neta1) & !is.na(fun2$err))]
    if (length(ns) > 0){
        stop(sprintf("Residual error component(s) need to be defined with assignment ('=' or '<-') in ini block (not '~'): %s", paste(ns, collapse=", ")))
    }
    ns <- fun2$name[is.na(fun2$est)];
    if (length(ns) > 0){
        stop(sprintf("The following parameters initial estimates are NA: %s", paste(ns, collapse=", ")))
    }
    fun2$meta <- list2env(meta, parent=emptyenv());
    return(fun2)
}
##' Print UI function
##'
##' @param x  UI function
##' @param ... other arguments
##' @author Matthew L. Fidler
##' @export
print.nlmixrUI <- function(x, ...){
    message(cli::rule(x$model.desc, line="bar2"))
    message(cli::rule(crayon::bold("Initialization:")))
    print(x$ini)
    if (length(x$all.covs) > 0){
        message("\n Covariates or Uninitialized Parameters ($all.covs)")
        print(x$all.covs);
    }
    message(cli::rule(crayon::bold(sprintf("Model%s:", ifelse(class(x$rxode) == "RxODE", " (RxODE)", "")))))
    message(x$fun.txt)
    message(cli::rule(line="bar2"))
}

## This is a list of supported distributions with the number of arguments they currently support.
dists <- list("dpois"=1,
              "dbinom"=2,
              "dbeta"=2:3,
              ##
              ## "dnbinom"=2:3,  ## dnbinom is in R; FIXME: how does ot compare to dneg_binomial
              ## "dneg_binomial", ## not in base R (but in glnmm2)
              ##
              ## Available as external package http://ugrad.stat.ubc.ca/R/library/rmutil/html/BetaBinom.html
              ## "dbetabinomial", ## not in base R (but in glnmm2)
              "dt"=1:2,
              "pois"=1,
              "binom"=2,
              "beta"=2:3,
              "t"=1:2,
              "add"=1,
              "prop"=1,
              "norm"=1,
              "dnorm"=1
              )

allVars <- function(x){
    defined <- character()
    this.env <- environment()
    f <- function(x){
        if (is.atomic(x)) {
            character()
        } else if (is.name(x)) {
            as.character(x)
        } else if (is.call(x) || is.pairlist(x)) {
            if (identical(x[[1]], quote(`~`)) ||
                identical(x[[1]], quote(`=`)) ||
                identical(x[[1]], quote(`<-`))){
                if (is.call(x[[3]])){
                    ret <- unique(unlist(lapply(x[[3]][-1], f)));
                } else {
                    ret <- unique(unlist(lapply(x[[3]], f)));
                }
                ret <- ret[!(ret %in% defined)]
                assign("defined", unique(c(defined, x[[2]])), this.env)
                return(ret)
            } else {
                children <- lapply(x[-1], f)
                unique(unlist(children))
            }
        } else {
            stop("Don't know how to handle type ", typeof(x),
                 call. = FALSE)
        }
    }
    f(x);
}

allNames <- function(x) {
    if (is.atomic(x)) {
        character()
    } else if (is.name(x)) {
        as.character(x)
    } else if (is.call(x) || is.pairlist(x)) {
        children <- lapply(x[-1], allNames)
        unique(unlist(children))
    } else {
        stop("Don't know how to handle type ", typeof(x),
             call. = FALSE)
    }
}

allCalls <- function(x) {
    if (is.atomic(x) || is.name(x)) {
        character()
    } else if (is.call(x)) {
        fname <- as.character(x[[1]])
        children <- lapply(x[-1], allCalls)
        unique(c(fname, unlist(children)))
    } else if (is.pairlist(x)) {
        unique(unlist(lapply(x[-1], allCalls), use.names = FALSE))
    } else {
        stop("Don't know how to handle type ", typeof(x), call. = FALSE)
    }
}

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

## cauchy = t w/df=1; Support?
## dgeom is a special negative binomial; Support?
## fixme
unsupported.dists <- c("dchisq", "chisq", "dexp", "df", "f", "dgeom", "geom",
                       "dhyper", "hyper", "dlnorm", "lnorm", "dunif", "unif",
                       "dweibull", "weibull",
                       ## for testing...
                       "nlmixrDist")
add.dists <- c("add", "prop");

##' Build linear solved information based on defined parameters.
##'
##' @param lhs List of defined parameters
##' @return List containing the translation parmateres (extra lines)
##'     to give what nlmixr expects, and model properties
##'     (parameterization, ncmt, oral, and tlag).  The infusion needs
##'     to be guessed from the data and is not included in this function's output
##' @author Matthew L. Fidler
nlmixrUILinCmt <- function(lhs){
    reg <- rex::rex(start, one_of("V", "v"), capture(number), end);
    w <- which(regexpr(reg, lhs) != -1);
    if (any(lhs == "V") || any(lhs == "v")){
        if (length(w) > 0){
            min.v <- min(as.numeric(gsub(reg, "\\1", lhs[w])));
            vs <- c("V", paste0("V", seq(min.v, min.v + 1)))
        } else {
            vs <- c("V", paste0("V", 1:2));
        }
    } else {
        if (length(w) > 0){
            min.v <- min(as.numeric(gsub(reg, "\\1", lhs[w])));
            vs <- paste0("V", seq(min.v, min.v + 2))
        } else {
            vs <- paste0("V", 1:3);
        }
    }

    par1 <- list(CL  =c("CL"),
                 V   =c("V", "VC", vs[1]),
                 CLD =c("Q", "CLD"),
                 VT  =c("VT", "VP", vs[2]),
                 CLD2=c("CLD2", "Q2"),
                 VT2 =c("VT2", "VP2", vs[3]));
    ## If Cmt #1 is the central depot compartment then
    ## K12, K21 = 2 cmt model
    ## K13, K31 = 3 compartment model
    ## K10 = elimination
    ## K01 = absoprtion
    ##
    ## If Cmt #2 is the central compartment then
    ## K23, K32 = 2 cmt model
    ## K24, K42 = 3 cmt model
    ## K20 = elimination
    ## K12 = absorption
    ##
    ## Both ways completely identify the 2/3 compartment model and are unique.
    par2 <- list(KE=c("KE", "KEL", "K", "K10", "K20"),
                 K12=c("K12", "K23"),
                 K21=c("K21", "K32"),
                 K13=c("K13", "K24"),
                 K31=c("K31", "K42"));
    oral.pars <- c("KA", "K12", "K01");
    tlag.pars <- c("TLAG");
    lhs.up <- toupper(lhs);
    npars <- 0;
    extra.lines <- c();
    par.trans <- list();
    for (i in seq_along(par1)){
        possible <- par1[[i]];
        val <- intersect(possible, lhs.up);
        if (length(val) == 1){
            npars <- npars + 1;
            if (i == 2 && npars == 1){
                stop(sprintf("Clearance (%s) and Volume (%s) both have to be specified in this 1-cmt solved paramterization.",
                             paste(par1[["CL"]], collapse=", "), paste(par1[["V"]], collapse=", ")));
            }
            if (i == 4 && npars < 4){
                stop(sprintf("Distribtuional Clearance (%s) and Peripheral Volume (%s) both have to be specified in this 2-cmt solved paramterization.",
                             paste(par1[["CLD"]], collapse=", "), paste(par1[["VT"]], collapse=", ")))
            }
            if (i == 6 && npars < 6){
                stop(sprintf("Distribtuional Clearance #2 (%s) and Peripheral Volume #2 (%s) both have to be specified in 3-cmt this solved paramterization.",
                             paste(par1[["CLD"]], collapse=", "), paste(par1[["VT"]], collapse=", ")))
            }
            cur <- names(par1)[i]
            w <- which(lhs.up == val);
            cur.lhs <- lhs[w];
            par.trans[[length(par.trans) + 1]] <- c(cur, cur.lhs);
            if (cur != cur.lhs){
                extra.lines[length(extra.lines) + 1] <- sprintf("%s = %s;", cur, cur.lhs);
            }
        } else if (length(val) > 2){
            stop(sprintf("Need clearer paramterization for %s; Currently defined these similiar parameters: %s",
                         names(par1)[i], paste(val, collapse=", ")))
        }
    }
    param <- 0;
    if (npars > 1){
        param <- 1;
        ncmt <- npars / 2;
    } else {
        for (i in seq_along(par2)){
            possible <- par2[i];
            val <- intersect(possible, lhs.up);
            if (length(val) == 1){
                npars <- npars + 1;
                if ((i == 2 && npars == 1) ||
                    (i == 4 && npars < 4) ||
                    (i == 6 && npars < 6)){
                    stop("Not all the appropriate micro-constants have been specified this solved paramterization.");
                }
                cur <- names(par2)[i]
                w <- which(lhs.up == val);
                cur.lhs <- lhs[w];
                par.trans[[length(par.trans) + 1]] <- c(cur, cur.lhs);
                if (cur != cur.lhs){
                    extra.lines[length(extra.lines) +1] <- sprintf("%s <- %s;", cur, cur.lhs);
                }
            } else if (length(val) > 2){
                stop(sprintf("Need clearer paramterization for %s; Currently defined these similiar parameters: %s",
                             names(par2)[i], paste(val, collapse=", ")))
            }
        }
        if (npars > 1){
            param <- 2;
            ncmt <- npars / 2;
        }
    }
    if (param == 0){
        return(NULL);
    }
    ## Now get oral / tlag
    oral <- FALSE;
    val <- intersect(oral.pars, lhs.up);
    if (length(val) == 1){
        w <- which(val == lhs.up);
        if (length(w) > 0){
            cur.lhs <- lhs[w]
            par.trans[[length(par.trans) + 1]] <- c("KA", cur.lhs);
            if (cur.lhs != "KA"){
                extra.lines[length(extra.lines) +1] <- sprintf("KA <- %s;", cur.lhs)
            }
        }
        oral <- TRUE
    } else if (length(val) == 2){
        stop(sprintf("Ambiguous Absorption constant; Could be: %s", paste(val, collapse=", ")));
    }
    tlag <- FALSE
    val <- intersect(tlag.pars, lhs.up);
    if (length(val) == 1){
        w <- which(val == lhs.up);
        cur.lhs <- lhs[w]
        par.trans[[length(par.trans) + 1]] <- c("TLAG", cur.lhs);
        if (cur.lhs != "TLAG"){
            extra.lines[length(extra.lines) +1] <- sprintf("TLAG <- %s;", cur.lhs)
        }
        if (!oral){
            stop("Absorpation lag time requires an solved oral compartmental model.");
        }
        tlag <- TRUE
    } else if (length(val) == 2){
        stop(sprintf("Ambiguous Lag-time constant; Could be: %s", paste(val, collapse=", ")));
    }
    extra.lines <- paste(extra.lines, collapse="\n");
    return(list(extra.lines=extra.lines,
                ncmt=ncmt, parameterization=param, oral=oral, tlag=tlag, par.trans=par.trans))
}

nlmixrUIModel <- function(fun, ini=NULL, bigmodel=NULL){
    ## Parses the UI function to extract predictions and errors, and the other model specification.
    rxode <- FALSE
    all.names <- allNames(body(fun));
    if (any(regexpr(rex::rex(start, or("rx_", "nlmixr_")), paste(all.names)) != -1)){
        stop("Parameters/States/Variables cannot start with `rx_` or `nlmixr_`");
    }
    all.vars <- allVars(body(fun));
    all.funs <- allCalls(body(fun));
    all.lhs <- nlmixrfindLhs(body(fun));
    errs.specified <- c()
    add.prop.errs <- data.frame(y=character(), add=logical(), prop=logical());
    bounds <- ini;
    theta.names <-  c()
    theta.ord <- c();
    eta.names <- c();
    .mu.ref <- list();
    cov.ref <- list();
    log.theta <- c();
    log.eta <- c();
    this.env <- environment();
    if (!is.null(ini)){
        unnamed.thetas <- ini$ntheta[(!is.na(ini$ntheta) & is.na(ini$name))];
        if (length(unnamed.thetas) > 0){
            stop(sprintf("The following THETAs are unnamed: %s", paste(sprintf("THETA[%d]", unnamed.thetas), collapse=", ")))
        }
        unnamed.etas <- ini$neta1[!is.na(ini$neta1) & (ini$neta1 == ini$neta2) & is.na(ini$name)];
        if (length(unnamed.etas) > 0){
            stop(sprintf("The following ETAs are unnamed: %s", paste(sprintf("ETA[%d]", unnamed.etas), collapse=", ")))
        }
        theta.names <- ini$theta.names;
        eta.names <- ini$eta.names;
    }
    errn <- 0

    any.theta.names <- function(what, th.names){
        ## for (i in th.names){
        ##     if (any(what == i)){
        ##         return(TRUE)
        ##     }
        ## }
        ## return(FALSE)
        return(any(what[1] == th.names))
    }

    find.theta <- function(x){
        if (length(x) == 1 && any.theta.names(as.character(x), theta.names)){
            return(as.character(x));
        } else if (identical(x[[1]], quote(`+`))){
            th <- c();
            if (length(x) >= 3){
                if (any.theta.names(as.character(x[[3]]), theta.names)){
                    th <- as.character(x[[3]]);
                }
            }
            if (length(x) >= 2){
                if (length(x[[2]]) > 1){
                    return(c(th, find.theta(x[[2]])))
                } else {
                    if (any.theta.names(as.character(x[[2]]), theta.names)){
                        th <- c(th, as.character(x[[2]]));
                    }
                    return(th)
                }
            }
        }
    }

    f <- function(x) {
        if (is.name(x)) {
            if (any.theta.names(as.character(x), theta.names)){
                assign("theta.ord", unique(c(theta.ord, as.character(x))), this.env)
            }
            return(x)
        } else if (is.call(x)) {
            if (identical(x[[1]], quote(`~`)) &&
                any(as.character(x[[3]][[1]]) == c(names(dists), unsupported.dists))){
                ch.dist <- as.character(x[[3]]);
                if (any(ch.dist[1] == unsupported.dists)){
                    stop(sprintf("The %s distribution is currently unsupported.", ch.dist[1]))
                }
                nargs <- dists[[ch.dist[1]]]
                if (any((length(ch.dist) - 1) == nargs)){
                    assign("errs.specified", unique(errs.specified, as.character(x[[3]][[1]])))
                    if (do.pred == 1){
                        return(bquote(nlmixr_pred <- .(x[[2]])));
                    } else if (do.pred == 0){
                        dist.name <- ch.dist[1];
                        dist.args <- ch.dist[-1];
                        for (i in seq_along(dist.args)){
                            tmp <- suppressWarnings(as.numeric(dist.args[i]))
                            errn <- errn + 1;
                            if (!is.na(tmp)){
                                ## FIXME: allow numeric estimates...?
                                stop("Distribution parameters cannot be numeric, but need to be estimated.")
                            }
                            w <- which(bounds$name == dist.args[i]);
                            tmp <- bounds;
                            tmp$err[w] <- dist.name;
                            assign("bounds", tmp, this.env)
                        }
                        if (any(ch.dist[1] == c("add", "norm"))){
                            assign("errn", errn + 1, this.env);
                            assign("add.prop.errs", rbind(add.prop.errs,
                                                          data.frame(y=sprintf("Y%02d", errn), add=TRUE, prop=FALSE)), this.env)
                        } else if (ch.dist[1] == "prop"){
                            assign("errn", errn + 1, this.env);
                            assign("add.prop.errs", rbind(add.prop.errs,
                                                          data.frame(y=sprintf("Y%02d", errn), add=FALSE, prop=TRUE)), this.env);
                        }
                        return(bquote(return(.(eval(parse(text=sprintf("quote(%s(%s))", dist.name, paste(dist.args, collapse=", "))))))));
                    } else if (do.pred == 3){ ## Dataset preparation function for nlme
                        return(bquote(return(.(sprintf("Y%02d", errn)))))
                    } else {
                        return(quote(nlmixrIgnore()))
                    }
                } else {
                    if (any(do.pred == c(2, 4, 5))){
                        return(x);
                    } else {
                        if (length(nargs) == 1){
                            stop(sprintf("The %s distribution requires %s arguments.", ch.dist[1], nargs))
                        } else {
                            stop(sprintf("The %s distribution requires %s-%s arguments.", ch.dist[1], min(nargs), max(nargs)))
                        }
                    }
                }
            } else if (identical(x[[1]], quote(`~`)) &&
                       identical(x[[3]][[1]], quote(`+`)) &&
                       identical(x[[3]][[2]][[1]], quote(`+`)) &&
                       any(as.character(x[[3]][[2]][[2]][[1]]) == c(names(dists), unsupported.dists))){
                stop(sprintf("Only 2 distributions can be combined.\nCurrently can combine: %s",
                             paste(add.dists, collapse=", ")))
            } else if (identical(x[[1]], quote(`~`)) &&
                       identical(x[[3]][[1]], quote(`+`)) &&
                       length(as.character(x[[3]])) == 3  &&
                       any(as.character(x[[3]][[2]][[1]]) == c(names(dists), unsupported.dists)) &&
                       any(as.character(x[[3]][[3]][[1]]) == c(names(dists), unsupported.dists))){
                err1 <- as.character(x[[3]][[2]][[1]]);
                err1.v <- as.character(x[[3]][[2]][[2]])
                err2 <- as.character(x[[3]][[3]][[1]]);
                err2.v <- as.character(x[[3]][[3]][[2]])
                if (!is.na(suppressWarnings(as.numeric(err1.v)))){
                    stop("Distribution parameters cannot be numeric, but need to be estimated.")
                }
                if (!is.na(suppressWarnings(as.numeric(err2.v)))){
                    stop("Distribution parameters cannot be numeric, but need to be estimated.")
                }
                if (any(err1 == add.dists) &&
                    any(err2 == add.dists)){
                    tmp <- paste(sort(c(err1, err2)), collapse="+");
                    assign("errs.specified", unique(errs.specified, tmp), this.env)
                    if (any(do.pred == c(2, 4, 5))){
                        return(quote(nlmixrIgnore()));
                    }
                    else if (any(do.pred == c(1, 4, 5))){
                        return(bquote(nlmixr_pred <- .(x[[2]]))) ;
                    } else if (do.pred == 3){
                        w <- which(bounds$name == err1.v);
                        tmp <- bounds;
                        tmp$err[w] <- err1;
                        w <- which(bounds$name == err2.v);
                        tmp$err[w] <- err2;
                        assign("bounds", tmp, this.env);
                        if ((any(paste(tmp$err) == "add") || any(paste(tmp$err) == "norm")) && any(paste(tmp$err) == "prop")){
                            assign("errn", errn + 1, this.env);
                            assign("add.prop.errs", rbind(add.prop.errs,
                                                          data.frame(y=sprintf("Y%02d", errn), add=TRUE, prop=TRUE)), this.env);
                        }
                        return(bquote(return(.(sprintf("Y%02d", errn)))));
                    } else {
                        return(bquote(return(.(x[[3]]))));
                    }
                } else {
                    stop(sprintf("The %s and %s distributions cannot be combined\nCurrently can combine: %s",
                                 as.character(x[[3]][[2]][[1]]), as.character(x[[3]][[3]][[1]]),
                                 paste(add.dists, collapse=", ")))
                }
            } else if (identical(x[[1]], quote(`~`)) && (do.pred != 2)){
                return(quote(nlmixrIgnore()))
            } else if (identical(x[[1]], quote(`<-`)) && !any(do.pred == c(2, 4, 5) )){
                return(quote(nlmixrIgnore()))
            } else if (identical(x[[1]], quote(`=`)) && !any(do.pred == c(2, 4, 5))){
                return(quote(nlmixrIgnore()))
            } else if (identical(x[[1]], quote(`<-`)) && do.pred == 4){
                ## SAEM requires = instead of <-
                x[[1]] <- quote(`=`);
                return(as.call(lapply(x, f)))
            } else if (identical(x[[1]], quote(`exp`)) && any(do.pred == c(4, 5))){
                ## Need traverse the parsing tree to get log theta/eta
                ## parameters.
                find.log <- function(x){
                    if (is.atomic(x) || is.name(x)) {
                        if (any.theta.names(as.character(x), theta.names)){
                            assign("log.theta", unique(c(log.theta, as.character(x))), this.env)
                        } else if (any.theta.names(as.character(x), eta.names)){
                            assign("log.eta", unique(c(log.eta, as.character(x))), this.env);
                        }
                        return(x)
                    } else if (is.pairlist(x)) {
                        return(lapply(x, find.log));
                    } else if (is.call(x)) {
                        return(lapply(x, find.log));
                    } else {
                        stop("Don't know how to handle type ", typeof(x),
                             call. = FALSE)
                    }
                }
                find.log(x[[2]])
                return(as.call(lapply(x, f)))
            } else if (identical(x[[1]], quote(`+`))){
                ## print(as.character(x))
                if (any(do.pred == c(4, 5))){
                    if (length(x) >= 3){
                        ## message("---")
                        ## print(x[[1]])
                        ## print(x[[2]]);
                        ## print(x[[3]]);
                        wm <- NULL
                        if (length(x[[3]]) == 3){
                            if (identical(x[[3]][[1]], quote(`*`))){
                                wm <- 3;
                                wm0 <- 2;
                            }
                        }
                        if (length(x[[2]]) == 3){
                            if (identical(x[[2]][[1]], quote(`*`))){
                                wm <- 2;
                                wm0 <- 3;
                            }
                        }
                        if (!is.null(wm)){
                            cur <- 2;
                            th <- 3;
                            w <- which(all.covs == as.character(x[[wm]][[2]])[1])
                            if (length(w) == 0){
                                cur <- 3;
                                th <- 2;
                                w <- which(all.covs == as.character(x[[wm]][[3]])[1])
                            }
                            if (length(w) == 1){
                                cov <- all.covs[w];
                                th <- as.character(x[[wm]][[th]]);
                                th0 <- find.theta(x[[wm0]]);
                                if (length(th0) == 1 &&  do.pred == 4){
                                    tmp <- get("cov.ref", this.env)
                                    tmp[[cov]] <- c(tmp[[cov]], structure(th0, .Names=th));
                                    assign("cov.ref", tmp, this.env);
                                    return(f(x[[wm0]]))
                                }
                            }
                        }
                        if (any.theta.names(as.character(x[[2]]), eta.names) &&
                            any.theta.names(as.character(x[[3]]), theta.names)){
                            ## Found ETA+THETA
                            tmp <- .mu.ref;
                            tmp[[as.character(x[[2]])]] <- as.character(x[[3]]);
                            ## assign("mu.ref", tmp, this.env);
                            .mu.ref <<- tmp
                            ## Collapse to THETA
                            return(x[[3]])
                        } else if (any.theta.names(as.character(x[[3]]), eta.names) &&
                                   any.theta.names(as.character(x[[2]]), theta.names)){
                            ## Found THETA+ETA
                            tmp <- .mu.ref
                            tmp[[as.character(x[[3]])]] <- as.character(x[[2]]);
                            ## assign(".mu.ref", tmp, this.env)
                            .mu.ref <<- tmp
                            ## Collapse to THETA
                            ## model$omega=diag(c(1,1,0))
                            ## 0 is not estimated.
                            ## inits$omega has the initial estimate
                            ## mod$res.mod = 1 = additive
                            ## mod$res.mod = 2 = proportional
                            ## mod$res.mod = 3 = additive + proportional
                            ## a+b*f
                            ## mod$ares = initial estimate of res
                            ## mod$bres = initial estimate of
                            return(x[[2]]);
                        } else if (any.theta.names(as.character(x[[3]]), eta.names) &&
                                   length(x[[2]]) > 1){
                            ## This allows 123 + Cl + 123 + eta.Cl + 123
                            ## And collapses to 123 + Cl + 123 + 123
                            ## Useful for covariates...
                            eta <- as.character(x[[3]]);
                            th <- find.theta(x[[2]]);
                            if (length(th) == 1){
                                tmp <- .mu.ref
                                tmp[[eta]] <- th;
                                ## assign("tmp", .mu.ref, this.env)
                                .mu.ref <<- tmp
                                return(f(as.call(x[[2]])));
                            }
                        } else if (length(x) < 3) {
                        } else if (any.theta.names(as.character(x[[3]]), theta.names) &&
                                   length(x[[2]]) > 1){
                            ## This allows 123 + eta.Cl + 123 + Cl + 123
                            ## And collapses to 123  + 123 + Cl + 123
                            ## Useful for covariates...
                            theta <- as.character(x[[3]]);
                            .etas <- c();
                            find.etas <- function(x){
                                if (is.atomic(x) || is.name(x)) {
                                    return(x)
                                } else if (is.pairlist(x)) {
                                    return(lapply(x, find.etas));
                                } else if (is.call(x)) {
                                    if (identical(x[[1]], quote(`+`)) &&
                                        any.theta.names(as.character(x[[3]]), eta.names)){
                                        .etas <<- c(.etas,as.character(x[[3]]));
                                        return(x[[2]]);
                                    }
                                    return(as.call(lapply(x, find.etas)));
                                } else {
                                    stop("Don't know how to handle type ", typeof(x),
                                         call. = FALSE)
                                }
                            }
                            new <- find.etas(x[[2]]);
                            if (length(.etas) == 1){
                                tmp <- .mu.ref;
                                tmp[[.etas]] <- theta;
                                ## assign("tmp", .mu.ref, this.env);
                                .mu.ref <<- tmp
                                x[[2]] <- new;
                            }
                        }
                    }
                }
                return(as.call(lapply(x, f)))
            } else {
                return(as.call(lapply(x, f)));
            }
        } else if (is.pairlist(x)) {
            as.pairlist(lapply(x, f));
        } else if (is.atomic(x)){
            return(x)
        } else {
            stop("Don't know how to handle type ", typeof(x),
                 call. = FALSE)
        }
    }
    rm.empty <- function(x){
        ## empty if/else
        ## First remove if () followed by nlmixrIgnore()
        ignoreit <- rex::rex(any_spaces, "nlmixrIgnore()", any_spaces)
        w1 <- which(regexpr(rex::rex(start, any_spaces, or(group("if", any_spaces, "(", anything, ")"), "else"),any_spaces, end), x) != -1);
        if (length(w1) > 0){
            w2 <- which(regexpr(ignoreit, x[w1 + 1]) != -1)
            if (length(w2) > 0){
                x <- x[-w1[w2]];
            }
        }
        x <- x[regexpr(ignoreit, x, perl=TRUE) == -1];
        w1 <- which(regexpr(rex::rex(start, any_spaces, or("if", "else"), anything, "{", end), x) != -1)
        if (length(w1) > 0){
            w2 <- w1 + 1;
            w3 <- which(regexpr(rex::rex(start, any_spaces, "}", end), x[w2]) != -1);
            if (length(w3) > 0){
                w1 <- w1[w3];
                w2 <- w2[w3];
                return(x[-c(w1, w2)])
            } else {
                return(x);
            }
        } else {
            return(x)
        }
    }
    new.fn <- function(x){
        x <- rm.empty(x);
        if (do.pred == 2){
            assign("rxode", any(regexpr(rex::rex(start, any_spaces, "d/dt(", anything, ")", any_spaces, or("=", "<-")), x) != -1), this.env)
        }
        x[1] <- paste0("function(){");
        x[length(x)] <- "}"
        x <- eval(parse(text=paste(x, collapse="\n")))
        return(x)
    }
    all.covs <- character();
    do.pred <- 1;
    pred.txt <- deparse(f(body(fun)))
    pred <- new.fn(pred.txt);
    do.pred <- 0;
    err <- new.fn(deparse(f(body(fun))));
    do.pred <- 2;
    rest.txt <- deparse(f(body(fun)))
    rest <- new.fn(rest.txt);
    rest.funs <- allCalls(body(rest));
    rest.vars <- allVars(body(rest));
    all.covs <- setdiff(rest.vars,paste0(bounds$name))
    do.pred <- 3;
    grp.fn <- new.fn(deparse(f(body(fun))));
    do.pred <- 4;
    saem.pars <- try(deparse(f(body(fun))), silent=TRUE);
    nlme.mu.fun2 <- NULL
    if (inherits(saem.pars, "try-error")){
        saem.pars <- NULL
    }
    do.pred <- 5;
    nlme.mu.fun <- try(deparse(f(body(fun))), silent=TRUE);
    if (inherits(nlme.mu.fun, "try-error")){
        nlme.mu.fun <- NULL
    }
    if (rxode){
        rx.txt <- deparse(body(rest))[-1]
        rx.txt <- rx.txt[-length(rx.txt)];
        reg <- rex::rex(boundary, or(ini$name), boundary);
        w <- which(regexpr(reg, rx.txt, perl=TRUE) != -1);
        if (length(w) == 0){
            stop("Error parsing model -- no parameters found.")
        }
        ## Separate ode and pred
        w <- max(w);
        if (any(regexpr(rex::rex(or("d/dt(", group("(0)", any_spaces, or("=", "~")))), rx.txt[1:w]) != -1)){
            ## mixed PK parameters and ODEs
            stop("Mixed PK/ODEs")
        } else {
            rx.ode <- rx.txt[-(1:w)];
            rx.pred <- eval(parse(text=paste(c("function() {", rx.txt[1:w], "}"), collapse="\n")))
            ## Now separate out parameters for SAEM.
            w <- max(which(regexpr(reg, saem.pars, perl=TRUE) != -1));
            saem.pars <- c(saem.pars[1:w], "");
            nlme.mu.fun2 <- saem.pars;
            w <- max(which(regexpr(reg, nlme.mu.fun, perl=TRUE) != -1));
            nlme.mu.fun <- c(nlme.mu.fun[1:w], "");
        }
        rxode <- paste(rx.ode, collapse="\n")
        rest <- rx.pred;
        all.vars <- all.vars[!(all.vars %in% RxODE::rxState(rxode))]
        rest.vars <- rest.vars[!(rest.vars %in% RxODE::rxState(rxode))]
        all.covs <- setdiff(rest.vars,paste0(bounds$name))
        all.covs <- all.covs[!(all.covs %in% RxODE::rxLhs(rxode))]
        all.covs <- all.covs[!(all.covs %in% RxODE::rxState(rxode))]
        all.covs <- setdiff(all.covs,c("t", "time", "podo", "M_E","M_LOG2E","M_LOG10E","M_LN2","M_LN10","M_PI","M_PI_2","M_PI_4","M_1_PI",
                                       "M_2_PI","M_2_SQRTPI","M_SQRT2","M_SQRT1_2","M_SQRT_3","M_SQRT_32","M_LOG10_2","M_2PI","M_SQRT_PI",
                                       "M_1_SQRT_2PI","M_SQRT_2dPI","M_LN_SQRT_PI","M_LN_SQRT_2PI","M_LN_SQRT_PId2","pi"))
    } else {
        all.covs <- setdiff(rest.vars,paste0(bounds$name))
        nlme.mu.fun2 <- saem.pars
        rxode <- NULL
    }
    fun2 <- as.character(attr(fun, "srcref"), useSource=TRUE)
    fun2[1] <- "function(){"
    fun2[length(fun2)] <- "}";
    fun2 <- try(eval(parse(text=paste0(fun2, collapse = "\n"), keep.source=TRUE)), silent=TRUE);
    if (inherits(fun2, "try-error")){
        stop("Error parsing model")
    }
    fun2 <- as.character(attr(fun2, "srcref"), useSource=TRUE);
    fun3 <- fun2
    fun3 <- fun3[-1];
    fun3 <- fun3[-length(fun3)]
    fun3 <- fun3[-length(fun3)]
    fun3 <- paste0(fun3, collapse="\n");

    misplaced.dists <- intersect(rest.funs, c(names(dists), unsupported.dists));
    if (length(misplaced.dists) == 1){
        if (misplaced.dists == "dt"){
            if (!any(regexpr("[^/]\\bdt[(]", deparse(rest), perl=TRUE) != -1)){
                misplaced.dists <- character();
            }
        }
    }
    if (length(misplaced.dists) > 0){
        stop(sprintf("Distributions need to be on residual model lines (like f ~ add(add.err)).\nMisplaced Distribution(s): %s", paste(misplaced.dists, collapse=", ")))
    }
    lin.solved <- NULL;
    tmp <- deparse(pred);
    tmp <- tmp[regexpr(rex::rex("nlmixr_pred <- "), tmp) != -1];
    if (length(tmp) > 0){
        if (all(regexpr(rex::rex("nlmixr_pred <- linCmt()"), tmp) != -1)){
            lin.solved <- nlmixrUILinCmt(all.lhs)
            par.ord <- as.list(lin.solved$par.trans);
            ## The theta.ord has to be adjusted...
            theta.ord2 <- c();
            log.theta2 <- c();
            log.eta2 <- c();
            saem.pars <- c("({", unlist(lapply(par.ord, function(x){
                x1 <- gsub(" *$", "", x[1]);
                x2 <- gsub(";", "", gsub("^ *", "", x[2]));
                reg <- rex::rex(start, any_spaces, x2, any_spaces, or("=", "<-"), any_spaces)
                w <- which(regexpr(reg, saem.pars) != -1);
                w2 <- which(regexpr(reg, rest.txt) != -1);
                if (length(w) == 1){
                    for (v in theta.ord){
                        if (regexpr(rex::rex(boundary, v, boundary), rest.txt[w2], perl=TRUE) != -1){
                            assign("theta.ord2", c(theta.ord2, v), this.env);
                        }
                    }
                    for (v in log.theta){
                        if (regexpr(rex::rex(boundary, v, boundary), rest.txt[w2], perl=TRUE) != -1){
                            assign("log.theta2", c(log.theta2, v), this.env);
                        }
                    }
                    for (v in log.eta){
                        if (regexpr(rex::rex(boundary, v, boundary), rest.txt[w2], perl=TRUE) != -1){
                            assign("log.eta2", c(log.eta2, v), this.env);
                        }
                    }
                    return(gsub(reg, paste0(x1, " = "), saem.pars[w]));
                } else {
                    return("error")
                }
                })), "})");
            theta.ord <- theta.ord2;
            log.theta <- log.theta2;
            log.eta <- log.eta2;
        }
    } else {
        add <- linCmt <- function(...) NULL
        pred <- function(){return(linCmt())}
        err <- function(){return(add(0.1))}
        grp.fn <- function(){return("Y01")};
        errs.specified <- c("add");
        add.prop.errs <- data.frame(y="Y1", add=TRUE, prop=FALSE);
        lin.solved <- nlmixrUILinCmt(all.lhs)
    }
    if (!is.null(saem.pars)){
        saem.pars <- new.fn(saem.pars)
        saem.theta.trans <- rep(NA, length(theta.names));
    } else {
        saem.pars <- NULL
        saem.theta.trans <- NULL
    }
    if (!is.null(nlme.mu.fun)){
        nlme.mu.fun <- new.fn(nlme.mu.fun)
    }
    if (!is.null(nlme.mu.fun2)){
        nlme.mu.fun2 <- new.fn(nlme.mu.fun2)
    }

    cov.theta.pars <- gsub(rex::rex(or(all.covs), "."), "", names(unlist(cov.ref)))
    for (i in seq_along(theta.names)){
        if (!any(theta.names[i] == cov.theta.pars)){
            w <- which(theta.names[i] == theta.ord);
            if (length(w) == 1){
                if (!any(theta.names[i] == cov.theta.pars)){
                    saem.theta.trans[i] <- w;
                }
            }
        }
    }
    cur <- 1;
    if (!all(is.na(saem.theta.trans))){
        while (cur <= max(saem.theta.trans, na.rm=TRUE)){
            while(!any(saem.theta.trans[!is.na(saem.theta.trans)] == cur)){
                w <- which(saem.theta.trans > cur);
                saem.theta.trans[w] <- saem.theta.trans[w] - 1;
            }
            cur <- cur + 1;
        }
    }
    env <- new.env(parent=emptyenv());
    env$infusion <- FALSE
    env$sum.prod <- FALSE
    ## Split out inPars
    saem.all.covs <- all.covs[all.covs %in% names(cov.ref)]
    saem.inPars <- all.covs[!(all.covs %in% names(cov.ref))]
    ret <- list(ini=bounds, model=bigmodel,
                nmodel=list(fun=fun2, fun.txt=fun3, pred=pred, error=err, rest=rest, rxode=rxode,
                            all.vars=all.vars, rest.vars=rest.vars, all.names=all.names, all.funs=all.funs, all.lhs=all.lhs,
                            all.covs=all.covs, saem.all.covs=saem.all.covs, saem.inPars=saem.inPars, lin.solved=lin.solved,
                            errs.specified=errs.specified, add.prop.errs=add.prop.errs, grp.fn=grp.fn, mu.ref=.mu.ref, cov.ref=cov.ref,
                            saem.pars=saem.pars, nlme.mu.fun=nlme.mu.fun, nlme.mu.fun2=nlme.mu.fun2, log.theta=log.theta,
                            log.eta=log.eta, theta.ord=theta.ord, saem.theta.trans=saem.theta.trans,
                            env=env))
    return(ret)
}
##' Create the nlme specs list for nlmixr nlme solving
##' @inheritParams nlmixrUI.nlmefun
##' @param mu.type is the mu-referencing type of model hat nlme will be using.
##' @return specs list for nlme
##' @author Matthew L. Fidler
nlmixrUI.nlme.specs <- function(object, mu.type=c("thetas", "covariates", "none")){
    mu.type <- match.arg(mu.type);
    if (mu.type == "thetas"){
        return(list(fixed=object$fixed.form,
                    random=object$random.mu,
                    start=object$theta))
    } else if (mu.type == "covariates") {
        theta <- names(object$theta);
        cov.ref <- object$cov.ref;
        cov.theta <- unique(as.vector(unlist(cov.ref)))
        cov.base <- theta[!(theta %in% cov.theta)]
        cov.base <- cov.base[!(cov.base %in% unlist(lapply(names(cov.ref), function(x){names(cov.ref[[x]])})))];
        cov.lst <- list();
        new.theta <- cov.base;
        for (n in names(cov.ref)){
            cov.base <- cov.base[!(cov.base %in% (names(cov.ref[[n]])))]
            cur <- cov.ref[[n]]
            for (i in seq_along(cur)){
                m <- cur[i]
                cov.lst[[m]] <- c(cov.lst[[m]], n);
                new.theta <- c(new.theta, as.vector(m), names(m))
            }
        }
        e1 <- paste(paste(cov.base, collapse="+"), "~ 1");
        fixed.form <- paste(c(e1, sapply(names(cov.lst), function(x){paste(x, "~", paste(cov.lst[[x]], collapse="+"))})), collapse=", ")
        fixed.form <- eval(parse(text=sprintf("list(%s)", fixed.form)))
        if (length(cov.base) == 0){
            fixed.form <- fixed.form[-1];
        }
        theta <- theta[new.theta]
        return(list(fixed=fixed.form,
                    random=object$random.mu,
                    start=object$theta))
    } else {
        return(list(fixed=object$fixed.form,
                    random=object$random,
                    start=object$theta))
    }
}
##' Create the nlme parameter transform function from the UI object.
##'
##' @param object UI object
##' @param mu Is the model mu referenced?
##' \itemize{
##'
##' \item With the "thetas" only the population parameters are
##' mu-referenced; All covariates are included in the model parameter
##' function.  The between subject variability pieces are specified in
##' the \code{random} specs parameter.
##'
##' \item With the "covariates" option, the population parameters are
##' mu referenced and covariates are removed from the model function.
##' The covariates will be specified used in the fixed effects
##' parameterization of nlme, like \code{list(lKA+lCL~1, lV~WT)}
##'
##' \item With the "none" option, the model function is given to nlme
##' without any modification.
##'
##' }
##' @return Parameter function for nlme
##' @author Matthew L. Fidler
##' @keywords internal
nlmixrUI.nlmefun <- function(object, mu.type=c("thetas", "covariates", "none")){
    ## create nlme function
    mu.type <- match.arg(mu.type);
    if (!is.null(object$lin.solved)){
        ## This is only a solved system.
        if (mu.type == "thetas"){
            bod <- deparse(body(object$nlme.mu.fun));
            bod[length(bod)] <- paste0(object$lin.solved$extra.lines, "\n}");
            bod <- eval(parse(text=sprintf("quote(%s)", paste0(bod, collapse="\n"))));
            fn <- eval(parse(text=sprintf("function(%s) NULL", paste(unique(c(names(object$ini$theta), object$all.covs)), collapse=", "))));
            body(fn) <- bod
            return(fn);
        } else if (mu.type == "covariates"){
            bod <- deparse(body(object$nlme.mu.fun2));
            bod[length(bod)] <- paste0(object$lin.solved$extra.lines, "\n}");
            bod <- eval(parse(text=sprintf("quote(%s)", paste0(bod, collapse="\n"))));
            vars <- unique(c(unlist(object$mu.ref), unlist(object$cov.ref)));
            fn <- eval(parse(text=sprintf("function(%s) NULL", paste(vars, collapse=", "))));
            vars2 <- allVars(bod);
            body(fn) <- bod;
            if (length(vars) != length(vars2)) return(NULL);
            return(fn);
        } else {
            bod <- deparse(body(object$rest));
            bod[length(bod)] <- paste0(object$lin.solved$extra.lines, "\n}");
            bod <- eval(parse(text=sprintf("quote(%s)", paste0(bod, collapse="\n"))));
            fn <- eval(parse(text=sprintf("function(%s) NULL", paste(object$rest.vars, collapse=", "))));
            body(fn) <- bod
            return(fn);
        }
    } else {
        if (mu.type == "thetas"){
            fn <- eval(parse(text=sprintf("function(%s) NULL", paste(unique(c(names(object$ini$theta), object$all.covs)), collapse=", "))))
            body(fn) <- body(object$nlme.mu.fun);
        } else if (mu.type == "covariates"){
            vars <- unique(c(unlist(object$mu.ref), unlist(object$cov.ref)));
            fn <- eval(parse(text=sprintf("function(%s) NULL", paste(vars, collapse=", "))))
            body(fn) <- body(object$nlme.mu.fun2);
            vars2 <- allVars(body(fn));
            if (length(vars) != length(vars2)){
                return(NULL);
            }
        } else {
            fn <- eval(parse(text=sprintf("function(%s) NULL", paste(object$rest.vars, collapse=", "))))
            body(fn) <- body(object$rest);
        }
    }
    return(fn)
}
##' Get the variance for the nlme fit process based on UI
##'
##' @param object UI object
##' @return nlme/lme variance object
##' @author Matthew L. Fidler
##' @keywords internal
nlmixrUI.nlme.var <- function(object){
    ## Get the variance for the nlme object
    add.prop.errs <- object$add.prop.errs;
    w.no.add <- which(!add.prop.errs$add);
    w.no.prop <- which(!add.prop.errs$prop);
    const <- grp <- ""
    power <- ", fixed=c(1)"
    powera <- ", fixed=list(power=1)"
    if (length(add.prop.errs$y) > 1){
        grp <- " | nlmixr.grp";
    }
    if (length(w.no.add) > 0){
        const <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y[w.no.add], "=0"), collapse=", "))
    }
    if (length(w.no.prop) > 0){
        power <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y, "=", ifelse(add.prop.errs$prop, 1, 0)), collapse=", "))
        powera <- sprintf(", fixed=list(power=list(%s))", paste(paste0(add.prop.errs$y, "=", ifelse(add.prop.errs$prop, 1, 0)), collapse=", "))
    }
    tmp <- sprintf("varConstPower(form=~fitted(.)%s%s)", grp, powera)
    if (all(!add.prop.errs$prop)){
        tmp <- sprintf("varIdent(form = ~ 1%s)", grp);
        if (tmp == "varIdent(form = ~ 1)"){
            warning("Initial condition for additive error ignored with nlme")
            return(NULL)
        }
    } else if (all(!add.prop.errs$add)){
        tmp <- sprintf("varPower(form = ~ fitted(.)%s%s)", grp, power);
    }
    return(eval(parse(text=tmp)))
}
##' Return RxODE model with predictions appended
##'
##' @param object UI object
##' @return String or NULL if RxODE is not specified by UI.
##' @author Matthew L. Fidler
nlmixrUI.rxode.pred <- function(object){
    if (is.null(object$rxode)){
        return(NULL)
    } else {
        tmp <- deparse(body(object$pred))[-1]
        tmp <- tmp[-length(tmp)]
        return(paste(c(object$rxode, tmp), collapse="\n"));
    }
}
##' Get the Parameter  function with THETA/ETAs defined
##'
##' @param obj UI object
##' @return parameters function defined in THETA[#] and ETA[#]s.
##' @author Matthew L. Fidler
nlmixrUI.theta.pars <- function(obj){
    df <- as.data.frame(obj$ini)
    dft <- df[!is.na(df$ntheta), ];
    dft.fixed <- dft[dft$fix, ];
    dft.unfixed <- dft[!dft$fix, ];
    fixed <- with(dft.fixed, sprintf("%s=%s", name, est))
    unfixed <- with(dft.unfixed, sprintf("%s=THETA[%d]", name, seq_along(dft.unfixed$name)))
    eta <- df[!is.na(df$neta1), ];
    eta <- eta[eta$neta1 == eta$neta2, ];
    eta <- with(eta, sprintf("%s=ETA[%d]", name, eta$neta1))
    f <- deparse(body(obj$rest))[-1]
    f <- eval(parse(text=paste(c("function(){", unfixed, eta, fixed, f[-length(f)], "}"), collapse="\n")))
    return(f)
}

##' Get the FOCEi initilizations
##'
##' @param obj UI object
##' @return list with FOCEi style initilizations
##' @author Matthew L. Fidler
nlmixrUI.focei.inits <- function(obj){
    df <- as.data.frame(obj$ini);
    dft <- df[!is.na(df$ntheta), ];
    dft.unfixed <- dft[!dft$fix, ];
    eta <- df[!is.na(df$neta1), ];
    len <- length(eta$name)
    cur.lhs <- character()
    cur.rhs <- numeric()
    ome <- character()
    for (i in seq_along(eta$name)){
        last.block <- FALSE;
        if (i == len){
            last.block <- TRUE
        } else if (eta$neta1[i + 1] == eta$neta2[i + 1]){
            last.block <- TRUE
        }
        if (eta$neta1[i] == eta$neta2[i]){
            cur.lhs <- c(cur.lhs, sprintf("ETA[%d]", eta$neta1[i]));
            cur.rhs <- c(cur.rhs, eta$est[i]);
            if (last.block){
                ome[length(ome) + 1] <- sprintf("%s ~ %s", paste(cur.lhs, collapse=" + "),
                                                paste(deparse(cur.rhs), collapse=" "));
                cur.lhs <- character();
                cur.rhs <- numeric()
            }
        } else {
            cur.rhs <- c(cur.rhs, eta$est[i]);
        }
    }
    ome <- eval(parse(text=sprintf("list(%s)", paste(ome, collapse=","))))
    return(list(THTA=dft.unfixed$est,
                OMGA=ome));
}
##' Get the eta->eta.trans for SAEM
##'
##' @param obj ui object
##' @return list of eta to eta.trans
##' @author Matthew L. Fidler
nlmixrUI.saem.eta.trans <- function(obj){
    eta.names <- obj$eta.names;
    theta.names <- obj$theta.names;
    theta.trans <- obj$saem.theta.trans;
    mu.ref <- obj$mu.ref
    trans <- rep(NA, length(eta.names))
    for (i in seq_along(eta.names)){
        ref <- mu.ref[[eta.names[i]]]
        if (!is.null(ref)){
            w <- which(ref == theta.names)
            if (length(w) == 1){
                trans[i] <- theta.trans[w];
            }
        }
    }
    if (any(is.na(trans))){
        stop("Could not figure out the mu-referencing for this model.")
    }
    return(trans)
}
##' Get the SAEM model Omega
##'
##' @param obj UI model
##' @return SAEM model$omega spec
##' @author Matthew L. Fidler
nlmixrUI.saem.model.omega <- function(obj){
    dm <- sum(!is.na(obj$saem.theta.trans));
    et <- obj$saem.eta.trans;
    mat <- matrix(rep(0, dm * dm), dm);
    etd <- which(!is.na(obj$neta1));
    for (i in etd){
        mat[et[obj$neta1[i]], et[obj$neta2[i]]] <- mat[et[obj$neta2[i]], et[obj$neta1[i]]] <- 1;
    }
    return(mat)
}
##' Get the SAEM model$res.mod code
##'
##' @param obj UI model
##' @return SAEM model$res.mod spec
##' @author Matthew L. Fidler
nlmixrUI.saem.res.mod <- function(obj){
    obj <- obj$add.prop.errs
    if (length(obj$add) == 1){
        if (obj$add && !obj$prop){
            return(1)
        } else if (obj$add && obj$prop){
            return(3)
        } else {
            return(2)
        }
    } else if (length(obj$add) == 0) {
        ## Use default...
        return(1)
    } else {
        stop("Currently SAEM only supports one type of residual error model.")
    }
}
##' Get error names for SAEM
##'
##' @param obj SAEM user interface function.
##' @return Names of error estimates for SAEM
##' @author Matthew L. Fidler
nlmixrUI.saem.res.name <- function(obj){
    w <- which(obj$err == "add");
    ret <- c();
    if (length(w) == 1){
        ret[length(ret) + 1] <- paste(obj$name[w])
    }
    w <- which(obj$err == "prop");
    if (length(w) == 1){
        ret[length(ret) + 1] <- paste(obj$name[w])
    }
    return(ret);
}

##' Get initial estimate for ares SAEM.
##'
##' @param obj UI model
##' @return SAEM model$ares spec
##' @author Matthew L. Fidler
nlmixrUI.saem.ares <- function(obj){
    w <- which(obj$err == "add");
    if (length(w) == 1){
        return(obj$est[w]);
    } else {
        return(10); ## SAME as SAEM
    }
}

##' Get initial estimate for bres SAEM.
##'
##' @param obj UI model
##' @return SAEM model$ares spec
##' @author Matthew L. Fidler
nlmixrUI.saem.bres <- function(obj){
    w <- which(obj$err == "prop");
    if (length(w) == 1){
        return(obj$est[w]);
    } else {
        return(1); ## SAME as SAEM
    }
}

##' Get model$log.eta for SAEM
##'
##' @param obj UI model
##' @return SAEM model$log.eta
##' @author Matthew L. Fidler
nlmixrUI.saem.log.eta <- function(obj){
    lt <- obj$log.theta;
    dm <- sum(!is.na(obj$saem.theta.trans));
    ret <- rep(FALSE, dm);
    theta.trans <- obj$saem.theta.trans;
    theta.names <- obj$theta.names;
    for (n in lt){
        w <- which(n == theta.names);
        if (length(w) == 1){
            ret[theta.trans[w]] <- TRUE
        }
    }
    return(ret)
}

##' Generate saem.fit user function.
##'
##' @param obj UI object
##' @return saem user function
##' @author Matthew L. Fidler
nlmixrUI.saem.fit <- function(obj){
    if (any(ls(envir=obj$env) == "saem.fit")){
        return(obj$env$saem.fit)
    } else if (!is.null(obj$rxode.pred)) {
        ## RxODE function
        message("Compiling RxODE differential equations...", appendLF=FALSE)
        if (obj$env$sum.prod){
            ode <- RxODE::RxODE(RxODE::rxSumProdModel(obj$rxode.pred));
        } else {
            ode <- RxODE::RxODE(obj$rxode.pred);
        }
        RxODE::rxLoad(ode);
        obj$env$saem.ode <- ode;
        RxODE::rxLoad(ode);
        message("done.")
        inPars <- obj$saem.inPars;
        if (length(inPars) == 0) inPars <- NULL
        saem.fit <- gen_saem_user_fn(model=ode, obj$saem.pars, pred=function() nlmixr_pred, inPars=inPars);
        message("done.")
        obj$env$saem.fit <- saem.fit;
        return(obj$env$saem.fit);
    } else if (!is.null(obj$lin.solved)) {
        message("Compiling SAEM user function...", appendLF=FALSE)
        saem.fit <- gen_saem_user_fn(model=lincmt(ncmt=obj$lin.solved$ncmt,
                                                  oral=obj$lin.solved$oral,
                                                  tlag=obj$lin.solved$tlag,
                                                  infusion = obj$env$infusion,
                                                  parameterization = obj$lin.solved$parameterization))
        message("done.")
        obj$env$saem.fit <- saem.fit;
        return(obj$env$saem.fit);
    }
}
##' Generate SAEM model list
##'
##' @param obj  nlmixr UI object
##' @return SAEM model list
##' @author Matthew L. Fidler
nlmixrUI.saem.model <- function(obj){
    mod <- list(saem_mod=obj$saem.fit);
    if (length(obj$saem.all.covs > 0)){
        mod$covars <- obj$saem.all.covs;
    }
    mod$res.mod <- obj$saem.res.mod;
    mod$log.eta <- obj$saem.log.eta;
    ## if (FALSE){
    ## FIXME option/warning
    mod$ares <- obj$saem.ares;
    mod$bres <- obj$saem.bres;
    ## }
    mod$omega <- obj$saem.model.omega;
    return(mod)
}
##' Get THETA names for nlmixr's SAEM
##'
##' @param uif nlmixr UI object
##' @return SAEM theta names
##' @author Matthew L. Fidler
nlmixrUI.saem.theta.name <- function(uif){
    trans <- uif$saem.theta.trans
    trans.name <- paste(uif$ini$name[which(!is.na(trans))]);
    trans <- trans[!is.na(trans)]
    theta.name <- trans.name[order(trans)]
    all.covs <- uif$saem.all.covs
    lc <- length(all.covs);
    if (lc > 0){
        m <- matrix(rep(NA, length(theta.name) * (lc + 1)), nrow=lc + 1);
        dimnames(m) <- list(c("_name",  all.covs), theta.name);
        m["_name", ] <- theta.name;
        for (cn in names(uif$cov.ref)){
            v <- uif$cov.ref[[cn]];
            for (var in names(v)){
                rn <- v[var];
                m[cn, rn] <- var
            }
        }
        ret <- unlist(m)
        ret <- ret[!is.na(ret)]
        return(ret)
    }
    return(theta.name);
}
##' Generate SAEM initial estimates for THETA.
##'
##' @param obj nlmixr UI object
##' @return SAEM theta initial estimates
##' @author Matthew L. Fidler
nlmixrUI.saem.init.theta <- function(obj){
    theta.name <- obj$saem.theta.name
    cov.names <- unique(names(unlist(structure(obj$cov.ref, .Names=NULL))));
    theta.name <- theta.name[!(theta.name %in% cov.names)];
    nm <- paste(obj$ini$name)
    lt <- obj$log.theta;
    i <- 0;
    this.env <- environment();
    theta.ini <- sapply(theta.name, function(x){
        w <- which(x == nm)
        assign("i", i + 1, this.env);
        if (any(lt == x)){
            return(exp(obj$ini$est[w]))
        } else {
            return(obj$ini$est[w])
        }
    })
    all.covs <- obj$saem.all.covs;
    lc <- length(all.covs);
    if (lc > 0){
        m <- matrix(rep(NA, lc * length(theta.name)), ncol=lc)
        dimnames(m) <- list(theta.name, all.covs);
        for (cn in names(obj$cov.ref)){
            v <- obj$cov.ref[[cn]];
            for (var in names(v)){
                rn <- v[var];
                w <- which(var == nm);
                m[rn, cn] <- obj$ini$est[w]
            }
        }
        return(as.vector(c(theta.ini, as.vector(m))))
    }
    return(as.vector(theta.ini))
}
##' SAEM's init$omega
##'
##' @param obj nlmixr UI object
##' @param names When \code{TRUE} return the omega names.  By default
##'     this is \code{FALSE}.
##' @return Return initial matrix
##' @author Matthew L. Fidler
nlmixrUI.saem.init.omega <- function(obj, names=FALSE){
    dm <- sum(!is.na(obj$saem.theta.trans));
    et <- obj$saem.eta.trans;
    ret <- rep(NA, dm);
    etd <- which(obj$neta1 == obj$neta2);
    for (i in etd){
        if (names){
            ret[et[obj$neta1[i]]] <- paste(obj$name[i])
        } else {
            ret[et[obj$neta1[i]]] <- obj$est[i]
        }
    }
    if (names){
        ret <- ret[!is.na(ret)];
        return(ret);
    } else {
        tmp <- unique(ret[!is.na(ret)])
        if (length(tmp) == 1){
            ret[is.na(ret)] <- tmp;
        } else {
            ret[is.na(ret)] <- 1;
        }
    }
    return(ret)
}
##' Get saem initilization list
##'
##' @param obj nlmixr UI object
##' @return Return SAEM inits list.
##' @author Matthew L. Fidler
nlmixrUI.saem.init <- function(obj){
    ret <- list();
    ret$theta <- obj$saem.init.theta;
    ## if (FALSE){
    ret$omega <-obj$saem.init.omega;
    ## }
    return(ret);
}

nlmixrUI.model.desc <- function(obj){
    if (!is.null(obj$rxode.pred)){
        return("RxODE-based model")
        ## n.cmt <- length(RxODE::rxState(RxODE::rxGetModel(obj$rxode.pred)));
        ## if (n.cmt == 0){
        ##     return("Compiled Model (with no ODEs)");
        ## } else {
        ##     return(sprintf("ODE(%d compartments)", n.cmt));
        ## }
    } else {
        return(sprintf("%s-compartment model%s%s", obj$lin.solved$ncmt,
                       ifelse(obj$lin.solved$oral, " with first-order absorption", ""),
                       ifelse(obj$lin.solved$parameterization == 1, " in terms of Cl", " in terms of micro-constants"),
                       ifelse(obj$lin.solved$tlag, " (with lag time)", "")));
    }
}

##' @export
`$.nlmixrUI` <- function(obj, arg, exact = TRUE){
    x <- obj;
    class(x) <- "list"
    if (arg == "ini"){
        return(x$ini);
    } else if (arg == "nmodel"){
        return(x$nmodel);
    } else if (arg == "model"){
        return(x$model);
    } else if (arg == "nlme.fun.mu"){
        return(nlmixrUI.nlmefun(obj, "thetas"))
    } else if (arg == "nlme.fun"){
        return(nlmixrUI.nlmefun(obj, "none"))
    } else if (arg == "nlme.fun.mu.cov"){
        return(nlmixrUI.nlmefun(obj, "covariates"))
    } else if (arg == "nlme.specs"){
        return(nlmixrUI.nlme.specs(obj, "none"))
    } else if (arg == "nlme.specs.mu"){
        return(nlmixrUI.nlme.specs(obj, "thetas"));
    } else if (arg == "nlme.specs.mu.cov"){
        return(nlmixrUI.nlme.specs(obj, "covariates"));
    } else if (arg == "nlme.var"){
        return(nlmixrUI.nlme.var(obj))
    } else if (arg == "rxode.pred"){
        return(nlmixrUI.rxode.pred(obj))
    } else if (arg == "theta.pars"){
        return(nlmixrUI.theta.pars(obj))
    } else if (arg == "focei.inits"){
        return(nlmixrUI.focei.inits(obj));
    } else if (arg == "saem.theta.name"){
        return(nlmixrUI.saem.theta.name(obj))
    } else if (arg == "saem.eta.trans"){
        return(nlmixrUI.saem.eta.trans(obj))
    } else if (arg == "saem.model.omega"){
        return(nlmixrUI.saem.model.omega(obj))
    } else if (arg == "saem.res.mod"){
        return(nlmixrUI.saem.res.mod(obj));
    } else if (arg == "saem.ares"){
        return(nlmixrUI.saem.ares(obj));
    } else if (arg == "saem.bres"){
        return(nlmixrUI.saem.bres(obj));
    } else if (arg == "saem.log.eta"){
        return(nlmixrUI.saem.log.eta(obj))
    } else if (arg == "saem.fit"){
        return(nlmixrUI.saem.fit(obj))
    } else if (arg == "saem.model"){
        return(nlmixrUI.saem.model(obj))
    } else if (arg == "saem.init.theta"){
        return(nlmixrUI.saem.init.theta(obj))
    } else if (arg == "saem.init.omega"){
        return(nlmixrUI.saem.init.omega(obj))
    } else if (arg == "saem.init"){
        return(nlmixrUI.saem.init(obj))
    } else if (arg == "saem.omega.name"){
        return(nlmixrUI.saem.init.omega(obj, TRUE))
    } else if (arg == "saem.res.name"){
        return(nlmixrUI.saem.res.name(obj));
    } else if (arg == "model.desc"){
        return(nlmixrUI.model.desc(obj))
    } else if (arg == "meta"){
        return(x$meta);
    } else if (arg == ".clean.dll"){
        if (exists(".clean.dll", envir=x$meta)){
            clean <- x$meta$.clean.dll;
            if (is(clean, "logical")){
                return(clean)
            }
        }
        return(TRUE);
    } else if (arg == "random.mu"){
        return(nlmixrBoundsOmega(x$ini,x$nmodel$mu.ref))
    }
    m <- x$ini;
    ret <- `$.nlmixrBounds`(m, arg, exact=exact)
    if (is.null(ret)){
        m <- x$nmodel;
        ret <- m[[arg, exact = exact]];
        if (is.null(ret)){
            if (exists(arg, envir=x$meta)){
                ret <- get(arg, envir=x$meta);
            }
        }
    }
    ret
}

##' @export
str.nlmixrUI <- function(object, ...){
    obj <- object;
    class(obj) <- "list";
    str(obj$ini);
    str(obj$nmodel)
    cat(" $ ini       : Model initilizations/bounds object\n");
    cat(" $ model     : Original Model\n");
    cat(" $ nmodel    : Parsed Model List\n");
    cat(" $ nlme.fun  : The nlme model function.\n");
    cat(" $ nlme.specs: The nlme model specs.\n");
    cat(" $ nlme.var  : The nlme model varaince.\n")
    cat(" $ rxode.pred: The RxODE block with pred attached (final pred is nlmixr_pred)\n")
    cat(" $ theta.pars: Parameters in terms of THETA[#] and ETA[#]\n")
    cat(" $ focei.inits: Initilization for FOCEi style blocks\n")
    cat(" $ saem.eta.trans: UI ETA -> SAEM ETA\n")
    cat(" $ saem.model.omega: model$omega for SAEM\n")
    cat(" $ saem.res.mod: model$res.mod for SAEM\n")
    cat(" $ saem.ares: model$ares for SAEM\n")
    cat(" $ saem.bres: model$bres for SAEM\n")
    cat(" $ saem.log.eta: model$log.eta for SAEM\n")
    cat(" $ saem.fit  : The SAEM fit user function\n")
    cat(" $ saem.model: The SAEM model list\n")
    cat(" $ saem.init.theta: The SAEM init$theta\n")
    cat(" $ saem.init.omega: The SAEM init$omega\n")
    cat(" $ saem.init : The SAEM inits list\n")
    cat(" $ saem.theta.name : The SAEM theta names\n")
    cat(" $ saem.omega.name : The SAEM theta names\n")
    cat(" $ saem.res.name : The SAEM omega names\n")
    cat(" $ model.desc : Model description\n")
    cat(" $ .clean.dll : boolean representing if dlls are cleaned after running.\n")
}
