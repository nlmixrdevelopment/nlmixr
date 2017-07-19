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
        fun2[length(fun2)] <- "return(list(ini=nlmixrBounds(ini),model=model))\n}"
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
    message("## Initial Conditions:")
    message("################################################################################")
    print(x$ini)
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

nlmixrUIModelPredErr <- function(fun){
    do.pred <- TRUE;
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
                    if (do.pred){
                        return(bquote(return(.(x[[2]]))));
                    } else {
                        return(bquote(return(.(x[[3]]))));
                    }
                } else {
                    if (length(nargs) == 1){
                        stop(sprintf("The %s distribution requires %s arguments.", ch.dist[1], nargs))
                    } else {
                        stop(sprintf("The %s distribution requires %s-%s arguments.", ch.dist[1], min(nargs), max(nargs)))
                    }
                }
            } else if (identical(x[[1]], quote(`~`)) &&
                       identical(x[[3]][[1]], quote(`+`)) &&
                       length(as.character(x[[3]])) == 3  &&
                       any(as.character(x[[3]][[2]][[1]]) == add.dists) &&
                       any(as.character(x[[3]][[2]][[1]]) == add.dists)){
                if (do.pred){
                    return(bquote(return(.(x[[2]]))));
                } else {
                    return(bquote(return(.(x[[3]]))));
                }
            } else if (identical(x[[1]], quote(`~`))){
                return(quote(nlmixrIgnore()))
            } else if (identical(x[[1]], quote(`<-`))){
                return(quote(nlmixrIgnore()))
            } else if (identical(x[[1]], quote(`=`))){
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
    pred <- deparse(f(body(fun)))
    pred <- pred[regexpr(rex::rex(any_spaces, "nlmixrIgnore()", any_spaces), pred, perl=TRUE) == -1];
    pred[1] <- paste0("function()", pred[1]);
    pred <- eval(parse(text=paste(pred, collapse="\n")))
    do.pred <- FALSE;
    err <- deparse(f(body(fun)));
    err <- err[regexpr(rex::rex(any_spaces, "nlmixrIgnore()", any_spaces), err, perl=TRUE) == -1];
    err[1] <- paste0("function()", err[1]);
    err <- eval(parse(text=paste(err, collapse="\n")))
    ret <- list(pred=pred, err=err)
    return(ret)
}

