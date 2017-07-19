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
        fun2[length(fun2)] <- "return(list(ini=nlmixrBounds(ini),model=nlmixrUIModel(model)))\n}"
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
    message("################################################################################")
    message("## Initial Conditions:")
    message("################################################################################")
    print(x$ini)
    message("################################################################################")
    message("## Model:")
    message("################################################################################")
    temp <- tempfile();
    on.exit({while(sink.number() != 0){sink()};if (file.exists(temp)){unlink(temp)}});
    sink(temp);
    print(x$model$fun);
    sink();
    fun2 <- readLines(temp);
    unlink(temp);
    fun2 <- fun2[-1];
    fun2 <- fun2[-length(fun2)]
    fun2 <- fun2[-length(fun2)]
    message(paste0(fun2, collapse="\n"));
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
## cauchy = t w/df=1; Support?
## dgeom is a special negative binomial; Support?
## fixme
unsupported.dists <- c("dchisq", "chisq", "dexp", "df", "f", "dgeom", "geom",
                       "dhyper", "hyper", "dlnorm", "lnorm", "dunif", "unif",
                       "dweibull", "weibull",
                       ## for testing...
                       "nlmixrDist")

add.dists <- c("add", "prop");

nlmixrUIModel <- function(fun){
    rxode <- FALSE
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
                       length(as.character(x[[3]])) == 3  &&
                       any(as.character(x[[3]][[2]][[1]]) == add.dists) &&
                       any(as.character(x[[3]][[2]][[1]]) == add.dists)){
                if (do.pred == 1){
                    return(bquote(return(.(x[[2]]))));
                } else {
                    return(bquote(return(.(x[[3]]))));
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
    ret <- list(fun=fun2, pred=pred, err=err, rest=rest, rxode=rxode)
    return(ret)
}

