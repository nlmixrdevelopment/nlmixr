## Start DSL tmp on http://adv-r.had.co.nz/dsl.tmp
## These operators are called to create the language and are not called in tests.

## These operators are called to create the language and are not called in tests.
## nocov start
unaryOp <- function(left, right) {
    force(left)
    force(right)
    function(e1) {
        paste0(left, e1, right)
    }
}

binaryOp <- function(sep) {
    force(sep)
    function(e1, e2) {
        if (missing(e2)){
            paste0(gsub(" ", "", sep), e1)
        } else {
            paste0(e1, sep, e2)
        }
    }
}

functionOp <- function(fn){
    force(fn);
    function(...){
        paste0(fn, "(", paste(unlist(list(...)), collapse=", "), ")")
    }
}

base.dsl.ui <- new.env(parent = emptyenv());

for (tmp in c("+", "-", "*", "/", "**", "^")){
    base.dsl.ui[[tmp]] <- functionOp(paste0(" ", tmp, " "));
}

for (tmp in c("abs", "acos", "acosh", "asin", "atan", "atan2", "atanh", "beta",
              "cos", "cosh", "digamma", "exp", "factorial","gamma", "log", "log10",
              "sin", "sinh", "sqrt", "tan", "tanh", "trigamma",
              ## Shold be supported?
              "choose", "lchoose", "psigamma")){
    base.dsl.ui[[tmp]] <- functionOp(tmp)
}

base.dsl.ui$"(" <- unaryOp("(", ")")

uiCurentBounds <- NULL;
uiCurrentTrans <- "nlme";
uiSupportUnkown <- TRUE

uiThetaEta <- function(name, val){
    if (is.null(uiCurrentBounds)){
        stop("The bounds are not set.")
    } else {
        n <- toupper(name);
        if (n == "THETA" && is.numeric(val)){
            if (uiCurentTrans == "nlme"){
                return(sprintf("theta.%s.", val));
            } else {
                stop("Unknown Translation")
            }
        } else if (n == "ETA" && is.numeric(val)){
            if (uiCurentTrans == "nlme"){
                return(sprintf("eta.%s.", val));
            } else {
                stop("Unknown Translation")
            }
        } else {
            stop("This only supports THETA[#] an ETA[#] vectors.")
        }
    }
}

base.dsl.ui$"[" <- uiThetaEta;

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

unknown.base.pars <- function(op){
    force(op)
    function(...){
        if (uiSupportUnkown){
            warning(sprintf("Other estimation methods may not support the function '%s'.", op))
            return(paste0(fn, "(", paste(unlist(list(...)), collapse=", "), ")"));
        } else {
            stop("This estimation method does not support the function '%s'.", op);
        }
    }
}

base.dsl.ui.env <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknown.base.pars), calls)
    callEnv <- list2env(callList);
    rxSymPyFEnv <- cloneEnv(base.dsl.ui, callEnv);
    names <- allNames(expr)
    n1 <- names;
    n2 <- names;
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSymPyFEnv);
    return(symbol.env)
}
