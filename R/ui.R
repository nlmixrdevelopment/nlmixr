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
    if (fun2[length(fun2)] == "}"){
        fun2[length(fun2)] <- "ini <- nlmixrBounds(ini);return(nlmixrUIModel(model,ini,fun))\n}"
    }
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
    message("\n##Initialization:")
    message("################################################################################")
    print(x$ini)
    if (length(x$all.covs) > 0){
        message("\n Covariates or Uninitialized Parameters ($all.covs)")
        print(x$all.covs);
    }
    message("\n## Model:")
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
    if (is.atomic(x)) {
        character()
    } else if (is.name(x)) {
        as.character(x)
    } else if (is.call(x) || is.pairlist(x)) {
        if (identical(x[[1]], quote(`~`)) ||
            identical(x[[1]], quote(`=`)) ||
            identical(x[[1]], quote(`<-`))){
            if (is.call(x[[3]])){
                unique(unlist(lapply(x[[3]][-1], allVars)));
            } else {
                unique(unlist(lapply(x[[3]], allVars)));
            }
        } else {
            children <- lapply(x[-1], allVars)
            unique(unlist(children))
        }
    } else {
        stop("Don't know how to handle type ", typeof(x),
             call. = FALSE)
    }
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
                         names(par)[i], paste(val, collapse=", ")))
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
                    if (do.pred == 1){
                        return(bquote(return(.(x[[2]]))));
                    } else if (do.pred == 0){
                        return(bquote(return(.(x[[3]]))));
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
                if (any(as.character(x[[3]][[2]][[1]]) == add.dists) &&
                    any(as.character(x[[3]][[3]][[1]]) == add.dists)){
                    if (do.pred == 2){
                        return(quote(nlmixrIgnore()));
                    }
                    else if (do.pred == 1){
                        return(bquote(return(.(x[[2]]))));
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
            rxode <- any(regexpr(rex::rex(start, any_spaces, "d/dt(", anything, ")", any_spaces, or("=", "<-")), x) != -1)
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
    all.covs <- setdiff(all.vars,paste0(ini$name))
    rest.funs <- allCalls(body(rest));
    misplaced.dists <- intersect(rest.funs, c(names(dists), unsupported.dists));
    if (length(misplaced.dists) > 0){
        stop(sprintf("Distributions need to be on residual model lines (like f ~ add(add.err)).\nMisplaced Distribution(s): %s", paste(misplaced.dists, collapse=", ")))
    }
    lin.solved <- NULL;
    tmp <- deparse(pred);
    tmp <- tmp[regexpr(rex::rex("return("), tmp) != -1];
    if (all(regexpr(rex::rex("return(linCmt())"), tmp) != -1)){
        lin.solved <- nlmixrUILinCmt(all.lhs)
    }
    ret <- list(ini=ini, model=bigmodel,
                nmodel=list(fun=fun2, fun.txt=fun3, pred=pred, err=err, rest=rest, rxode=rxode,
                            all.vars=all.vars, all.names=all.names, all.funs=all.funs, all.lhs=all.lhs,
                            all.covs=all.covs, lin.solved=lin.solved))
    return(ret)
}

nlmixrUI.nlme.specs <- function(object){
    return(list(fixed=object$fixed.form,
                random=object$random,
                start=object$theta))
}

nlmixrUI.nlmefun <- function(object){
    ## create nlme function
    if (!is.null(object$lin.solved)){
        ## This is only a solved system.
        bod <- deparse(body(object$rest));
        bod[length(bod)] <- paste0(object$lin.solved$extra.lines, "\n}");
        bod <- eval(parse(text=sprintf("quote(%s)", paste0(bod, collapse="\n"))));
        fn <- eval(parse(text=sprintf("function(%s) NULL", paste(object$all.vars, collapse=", "))));
        body(fn) <- bod
        return(fn);
    } else {
        fn <- eval(parse(text=sprintf("function(%s) NULL", paste(object$all.vars, collapse=", "))))
        body(fn) <- body(object$rest);
    }
    return(fn)
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
    str(obj$nmodel)
}
