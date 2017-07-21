##' Prepares the UI function and returns a list.
##'
##' @param fun UI function
##' @return nlmixr UI function
##' @author Matthew L. Fidler
nlmixrUI <- function(fun){
    temp <- tempfile();
    on.exit({while(sink.number() != 0){sink()};if (file.exists(temp)){unlink(temp)}});
    sink(temp);
    print(fun);
    sink();
    fun2 <- readLines(temp);
    unlink(temp);
    if (regexpr("^<.*>$", fun2[length(fun2)]) != -1){
        fun2 <- fun2[-length(fun2)];
    }
    fun2 <- gsub(rex::rex(boundary, "ini(",any_spaces,"{"), "ini <- function()({", fun2)
    fun2 <- gsub(rex::rex(boundary, "model(",any_spaces,"{"), "model <- function()({", fun2)
    fun2[length(fun2)] <- "ini <- nlmixrBounds(ini);return(nlmixrUIModel(model,ini,fun))\n}"
    fun2 <- try(eval(parse(text=paste0(fun2, collapse = "\n"))), silent=TRUE);
    if (inherits(fun2, "try-error")){
        stop("Error parsing model")
    }
    fun2 <- fun2();
    class(fun2) <- "nlmixrUI";
    return(fun2)
}
##' Print UI function
##'
##' @param x  UI function
##' @param ... other arguments
##' @author Matthew L. Fidler
##' @export
print.nlmixrUI <- function(x, ...){
    ## nlmixrLogo("", "Model")
    message("\n## Initialization:")
    message("################################################################################")
    print(x$ini)
    if (length(x$all.covs) > 0){
        message("\n Covariates or Uninitialized Parameters ($all.covs)")
        print(x$all.covs);
    }
    message(sprintf("\n## Model%s:", ifelse(class(x$rxode) == "RxODE", " (RxODE)", "")))
    message("################################################################################")
    message(x$fun.txt)
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
                defined <<- unique(c(defined, x[[2]]))
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
    par1 <- list(CL  =c("CL"),
                 V   =c("V", "VC"),
                 CLD =c("Q", "CLD"),
                 VT  =c("VT", "VP"),
                 CLD2=c("CLD2", "Q2"),
                 VT2 =c("VT2", "VP2"));
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
            if (cur != cur.lhs){
                extra.lines[length(extra.lines) + 1] <- sprintf("%s <- %s;", cur, cur.lhs);
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
        w <- which(oral.pars == lhs.up);
        cur.lhs <- lhs[w]
        if (cur.lhs != "KA"){
            extra.lines[length(extra.lines) +1] <- sprintf("KA <- %s;", cur.lhs)
        }
        oral <- TRUE
    } else if (length(val) == 2){
        stop(sprintf("Ambiguous Absorption constant; Could be: %s", paste(val, collapse=", ")));
    }
    tlag <- FALSE
    val <- intersect(tlag.pars, lhs.up);
    if (length(val) == 1){
        w <- which(tlag.pars == lhs.up);
        cur.lhs <- lhs[w]
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
                ncmt=ncmt, parameterization=param, oral=oral, tlag=tlag))
}

nlmixrUIModel <- function(fun, ini=NULL, bigmodel=NULL){
    ## Parses the UI function to extract predictions and errors, and the other model specification.
    rxode <- FALSE
    all.names <- allNames(body(fun));
    all.vars <- allVars(body(fun));
    all.funs <- allCalls(body(fun));
    all.lhs <- nlmixrfindLhs(body(fun));
    errs.specified <- c()
    add.prop.errs <- data.frame(y=character(), add=logical(), prop=logical());
    bounds <- ini;
    errn <- 0
    f <- function(x) {
        if (is.name(x)) {
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
                    errs.specified <<- unique(errs.specified,
                                              as.character(x[[3]][[1]]))
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
                            bounds$err[w] <<- dist.name;
                        }
                        if (any(ch.dist[1] == c("add", "norm"))){
                            errn <<- errn + 1;
                            add.prop.errs <<- rbind(add.prop.errs,
                                                    data.frame(y=sprintf("Y%02d", errn), add=TRUE, prop=FALSE))
                        } else if (ch.dist[1] == "prop"){
                            errn <<- errn + 1;
                            add.prop.errs <<- rbind(add.prop.errs,
                                                    data.frame(y=sprintf("Y%02d", errn), add=FALSE, prop=TRUE))
                        }
                        return(bquote(return(.(eval(parse(text=sprintf("quote(%s(%s))", dist.name, paste(dist.args, collapse=", "))))))));
                    } else if (do.pred == 3){ ## Dataset preparation function for nlme
                        return(bquote(return(.(sprintf("Y%02d", errn)))))
                    } else {
                        return(quote(nlmixrIgnore()))
                    }
                } else {
                    if (do.pred == 2){
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
                    errs.specified <<- unique(errs.specified, tmp);
                    if (do.pred == 2){
                        return(quote(nlmixrIgnore()));
                    }
                    else if (do.pred == 1){
                        return(bquote(nlmixr_pred <- .(x[[2]]))) ;
                    } else if (do.pred == 3){
                        w <- which(bounds$name == err1.v);
                        bounds$err[w] <<- err1;
                        w <- which(bounds$name == err2.v);
                        bounds$err[w] <<- err2;
                        if (any(tmp == c("norm+prop", "add+prop"))){
                            errn <<- errn + 1;
                            add.prop.errs <<- rbind(add.prop.errs,
                                                    data.frame(y=sprintf("Y%02d", errn), add=TRUE, prop=TRUE))
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
            } else if (identical(x[[1]], quote(`~`)) && do.pred != 2){
                return(quote(nlmixrIgnore()))
            } else if (identical(x[[1]], quote(`<-`)) && do.pred != 2){
                return(quote(nlmixrIgnore()))
            } else if (identical(x[[1]], quote(`=`)) && do.pred != 2){
                return(quote(nlmixrIgnore()))
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
        x <- x[regexpr(rex::rex(any_spaces, "nlmixrIgnore()", any_spaces), x, perl=TRUE) == -1];
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
            rxode <<- any(regexpr(rex::rex(start, any_spaces, "d/dt(", anything, ")", any_spaces, or("=", "<-")), x) != -1)
        }
        x[1] <- paste0("function(){");
        x[length(x)] <- "}"
        x <- eval(parse(text=paste(x, collapse="\n")))
        return(x)
    }
    do.pred <- 1;
    pred <- new.fn(deparse(f(body(fun))))
    do.pred <- 0;
    err <- new.fn(deparse(f(body(fun))));
    do.pred <- 2;
    rest <- new.fn(deparse(f(body(fun))));
    do.pred <- 3;
    grp.fn <- new.fn(deparse(f(body(fun))));
    rest.funs <- allCalls(body(rest));
    rest.vars <- allVars(body(rest));
    if (rxode){
        rx.txt <- deparse(body(rest))[-1]
        rx.txt <- rx.txt[-length(rx.txt)];
        ## now parse ode
        w <- seq(which(regexpr(rex::rex("d/dt("), rx.txt) != -1)[1], length(rx.txt));
        rx.ode <- rx.txt[w];
        rx.pred <- eval(parse(text=paste(c("function() {", rx.txt[-w], "}"), collapse="\n")))
        rxode <- paste(rx.ode, collapse="\n")
        rest <- rx.pred;
        all.vars <- all.vars[!(all.vars %in% RxODE::rxState(rxode))]
        rest.vars <- rest.vars[!(rest.vars %in% RxODE::rxState(rxode))]
    } else {
        rxode <- NULL
    }
    temp <- tempfile();
    on.exit({while(sink.number() != 0){sink()};if (file.exists(temp)){unlink(temp)}});
    sink(temp);
    print(fun);
    sink();
    fun2 <- readLines(temp);
    unlink(temp);
    if (regexpr("^<.*>$", fun2[length(fun2)]) != -1){
        fun2 <- fun2[-length(fun2)];
    }
    fun2[1] <- "function(){"
    fun2[length(fun2)] <- "}";
    fun2 <- try(eval(parse(text=paste0(fun2, collapse = "\n"))), silent=TRUE);
    if (inherits(fun2, "try-error")){
        stop("Error parsing model")
    }
    temp <- tempfile();
    on.exit({while(sink.number() != 0){sink()};if (file.exists(temp)){unlink(temp)}});
    sink(temp);
    print(fun2);
    sink();
    fun3 <- readLines(temp);
    unlink(temp);
    fun3 <- fun3[-1];
    fun3 <- fun3[-length(fun3)]
    fun3 <- fun3[-length(fun3)]
    fun3 <- paste0(fun3, collapse="\n");
    all.covs <- setdiff(rest.vars,paste0(bounds$name))
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
        }
    } else {
        pred <- function(){return(linCmt())}
        err <- function(){return(add(0.1))}
        grp.fn <- function(){return("Y01")};
        errs.specified <- c("add");
        add.prop.errs <- data.frame(y="Y1", add=TRUE, prop=FALSE);
        lin.solved <- nlmixrUILinCmt(all.lhs)
    }

    ret <- list(ini=bounds, model=bigmodel,
                nmodel=list(fun=fun2, fun.txt=fun3, pred=pred, error=err, rest=rest, rxode=rxode,
                            all.vars=all.vars, rest.vars=rest.vars, all.names=all.names, all.funs=all.funs, all.lhs=all.lhs,
                            all.covs=all.covs, lin.solved=lin.solved,
                            errs.specified=errs.specified,
                            add.prop.errs=add.prop.errs,
                            grp.fn=grp.fn))
    return(ret)
}
##' Create the nlme specs list for nlmixr nlme solving
##'
##' @param object UI object
##' @return specs list for nlme
##' @author Matthew L. Fidler
nlmixrUI.nlme.specs <- function(object){
    return(list(fixed=object$fixed.form,
                random=object$random,
                start=object$theta))
}
##' Create the nlme parameter transform function from the UI object.
##'
##' @param object UI object
##' @return parameter function for nlme
##' @author Matthew L. Fidler
##' @keywords internal
nlmixrUI.nlmefun <- function(object){
    ## create nlme function
    if (!is.null(object$lin.solved)){
        ## This is only a solved system.
        bod <- deparse(body(object$rest));
        bod[length(bod)] <- paste0(object$lin.solved$extra.lines, "\n}");
        bod <- eval(parse(text=sprintf("quote(%s)", paste0(bod, collapse="\n"))));
        fn <- eval(parse(text=sprintf("function(%s) NULL", paste(object$rest.vars, collapse=", "))));
        body(fn) <- bod
        return(fn);
    } else {
        fn <- eval(parse(text=sprintf("function(%s) NULL", paste(object$rest.vars, collapse=", "))))
        body(fn) <- body(object$rest);
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
    if (length(add.prop.errs$y) > 1){
        grp <- " | nlmixr.grp";
    }
    if (length(w.no.add) > 0){
        const <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y[w.no.add], "=0"), collapse=", "))
    }
    if (length(w.no.prop) > 0){
        power <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y, "=", ifelse(add.prop.errs$prop, 1, 0)), collapse=", "))
    }
    tmp <- sprintf("varComb(varIdent(form = ~ 1%s%s), varPower(form=~fitted(.)$s%s))", grp, const, grp, power)
    if (all(!add.prop.errs$prop)){
        tmp <- sprintf("varIdent(form = ~ 1%s)", grp);
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
    } else if (arg == "nlme.fun"){
        return(nlmixrUI.nlmefun(obj))
    } else if (arg == "nlme.specs"){
        return(nlmixrUI.nlme.specs(obj))
    } else if (arg == "nlme.var"){
        return(nlmixrUI.nlme.var(obj))
    } else if (arg == "rxode.pred"){
        return(nlmixrUI.rxode.pred(obj))
    }
    m <- x$ini;
    ret <- `$.nlmixrBounds`(m, arg, exact=exact)
    if (is.null(ret)){
        m <- x$nmodel;
        return(m[[arg, exact = exact]]);
    } else {
        return(ret)
    }
}

##' @export
str.nlmixrUI <- function(object, ...){
    obj <- object;
    class(obj) <- "list";
    str(obj$ini);
    message(" $ ini       : Model initilizations/bounds object");
    message(" $ model     : Original Model");
    message(" $ nmodel    : Parsed Model List");
    message(" $ nlme.fun  : The nlme model function.");
    message(" $ nlme.specs: The nlme model specs.");
    message(" $ nlme.var  : The nlme model varaince.")
    message(" $ rxode.pred: The RxODE block with pred attached (final pred is nlmixr_pred)")
    str(obj$nmodel)
}
