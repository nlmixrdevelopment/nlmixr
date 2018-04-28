##' Extract the Nlmixr bound information from a function.
##'
##' @param fun Function to extract bound information from.
##' @return a dataframe with bound information.
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
nlmixrBounds <- function(fun){
    fun2 <- as.character(attr(fun,"srcref"),useSource=TRUE);
    w <- which(regexpr("^ *#+.*", fun2) == 1);
    if (length(w) > 0){
        fun2 <- fun2[-w];
    }
    w <- which(regexpr("#+.*", fun2) != -1);
    if (length(w) > 0){
        labels <- gsub(".*#+ *(.*) *$", "\\1", fun2[w]);
                           labels <- sapply(labels,
                                            function(x){
                               return(sprintf("label(%s)", paste0(deparse(x))))});
        fun2[w] <- paste0(fun2[w], "\n", labels);
    }
    w <- which(regexpr("^ *[|].*", fun2) != -1);
    if (length(w) > 0){
        stop("A conditional statement cannot be on a line by itself.")
    }
    w <- which(regexpr("^.*[|]", fun2) != -1);
    if (length(w) > 0){
        condition <- gsub(".*[|] *(.*) *$", "\\1", fun2[w]);
        condition <- sapply(condition,
                            function(x){
            return(sprintf("condition(%s)", paste0(deparse(x))))});
        fun2[w] <- paste0(gsub("[|].*", "", fun2[w]), "\n", condition);
    }
    fun2 <- try(eval(parse(text=paste0(fun2, collapse = "\n"))), silent=TRUE);
    if (inherits(fun2, "try-error")){
        stop("Error parsing bounds; Perhaps there is an (unsupported) comment/condition inside the bounds themselves.")
    }
    theta <- 0;
    eta1 <- 0;
    df <- data.frame(ntheta=numeric(),
                     neta1=numeric(),
                     neta2=numeric(),
                     name=character(),
                     lower=numeric(),
                     est=numeric(),
                     upper=numeric(),
                     fix=logical(),
                     err=character(),
                     label=character(),
                     condition=character());
    netas <- 0;
    nerr <- 0;
    f <- function(x, env) {
        if (is.name(x)) {
            character()
        } else if (is.call(x)) {
            nenv <- new.env(parent=emptyenv());
            nenv$do.fixed <- FALSE;
            fix <- FIX <- fixed <- FIXED <- function(x){
                assign("do.fixed", TRUE, nenv)
                return(x)
            }
            if ((identical(x[[1]], quote(`<-`)) ||
                 identical(x[[1]], quote(`=`))) &&
                is.name(x[[2]])) {
                if (length(x[[3]]) > 5){
                    stop(sprintf("%s %s c(%s) syntax is not supported for thetas",
                                 as.character(x[[2]]), as.character(x[[1]]), paste(sapply(x[[3]][-1], as.character), collapse=", ")))
                } else if (length(x[[3]]) == 5){
                    if (as.character(x[[3]][[1]]) == "c" &&
                        any(tolower(as.character(x[[3]][[5]])) == c("fix", "fixed")) &&
                        length(x[[3]][[5]]) == 1){
                        ## a = c(1,2,3,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower=as.numeric(eval(x[[3]][[2]])),
                                                   est=as.numeric(eval(x[[3]][[3]])),
                                                   upper=as.numeric(eval(x[[3]][[4]])),
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA))
                    } else {
                        stop(sprintf("%s %s c(%s) syntax is not supported for thetas", as.character(x[[2]]),
                                     as.character(x[[1]]), paste(sapply(x[[3]][-1], as.character), collapse=", ")))
                    }
                } else if (length(x[[3]]) == 4 &&
                           any(tolower(as.character(x[[3]][[1]])) == c("c", "fix", "fixed"))){
                    nenv$do.fixed <- any(tolower(as.character(x[[3]][[1]])) == c("fix", "fixed"));
                    if (any(tolower(as.character(x[[3]][[4]])) == c("fix", "fixed")) &&
                        length(x[[3]][[4]]) == 1){
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower=as.numeric(eval(x[[3]][[2]])),
                                                   est=as.numeric(eval(x[[3]][[3]])),
                                                   upper=Inf,
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    } else {
                        ## a = c(1,2,3)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower=as.numeric(eval(x[[3]][[2]])),
                                                   est=as.numeric(eval(x[[3]][[3]])),
                                                   upper=as.numeric(eval(x[[3]][[4]])),
                                                   fix=nenv$do.fixed,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    }
                } else if (length(x[[3]]) == 3 &&
                           any(tolower(as.character(x[[3]][[1]])) == c("c", "fix", "fixed"))){
                    nenv$do.fixed <- any(tolower(as.character(x[[3]][[1]])) == c("fix", "fixed"));
                    if  (any(tolower(as.character(x[[3]][[3]])) == c("fix", "fixed"))
                         && length(x[[3]][[3]]) == 1) {
                        ## a = c(1,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower= -Inf,
                                                   est=as.numeric(eval(x[[3]][[2]])),
                                                   upper=Inf,
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    }  else {
                        ## a = c(1,2)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower=as.numeric(eval(x[[3]][[2]])),
                                                   est=as.numeric(eval(x[[3]][[3]])),
                                                   upper=Inf,
                                                   fix=nenv$do.fixed,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    }
                } else if (length(x[[3]]) == 2 &&
                           any(tolower(as.character(x[[3]][[1]])) == c("c", "fix", "fixed"))){
                    ## a = c(1)
                    nenv$do.fixed <- any(tolower(as.character(x[[3]][[1]])) == c("fix", "fixed"));
                    env$theta <- env$theta + 1;
                    env$df <- rbind(env$df,
                                    data.frame(ntheta=env$theta,
                                               neta1=NA,
                                               neta2=NA,
                                               name=as.character(x[[2]]),
                                               lower=-Inf,
                                               est=as.numeric(eval(x[[3]][[2]])),
                                               upper=Inf,
                                               fix=nenv$do.fixed,
                                               err=NA,
                                               label=NA,
                                               condition=NA))
                } else if (length(x[[3]]) == 1){
                    ## a = 1
                    env$theta <- env$theta + 1;
                    env$df <- rbind(env$df,
                                    data.frame(ntheta=env$theta,
                                               neta1=NA,
                                               neta2=NA,
                                               name=as.character(x[[2]]),
                                               lower=-Inf,
                                               est=as.numeric(eval(x[[3]][[1]])),
                                               upper=Inf,
                                               fix=nenv$do.fixed,
                                               err=NA,
                                               label=NA,
                                               condition=NA))
                } else {
                    num <- try(eval(x[[3]]), silent=TRUE)
                    if (is.numeric(num)){
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower= -Inf,
                                                   est=num,
                                                   upper=Inf,
                                                   fix=nenv$do.fixed,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA))
                    }
                }
            } else if (any(tolower(as.character(x[[1]])) == c("c", "fix", "fixed"))){
                nenv$do.fixed <- any(tolower(as.character(x[[1]])) == c("fix", "fixed"));
                if (length(x) > 5){
                    stop(sprintf("c(%s) syntax is not supported for thetas", paste(sapply(x[-1], as.character), collapse=", ")))
                } else if (length(x) == 5){
                    if (any(tolower(as.character(x[[5]])) == c("fix", "fixed")) &&
                        length(x[[5]]) == 1){
                        ## c(1,2,3,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=NA,
                                                   lower=as.numeric(eval(x[[2]])),
                                                   est=as.numeric(eval(x[[3]])),
                                                   upper=as.numeric(eval(x[[4]])),
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    } else {
                        stop(sprintf("c(%s) syntax is not supported for thetas", paste(sapply(x[-1], as.character), collapse=", ")))
                    }
                } else if (length(x) == 4){
                    ## No assignment...
                    if (any(tolower(as.character(x[[4]])) == c("fix", "fixed")) &&
                        length(x[[4]]) == 1){
                        ## c(1,2,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=NA,
                                                   lower=as.numeric(eval(x[[2]])),
                                                   est=as.numeric(eval(x[[3]])),
                                                   upper=Inf,
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    } else {
                        ## c(1,2,3)
                        env$theta <- env$theta + 1
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=NA,
                                                   lower=as.numeric(eval(x[[2]])),
                                                   est=as.numeric(eval(x[[3]])),
                                                   upper=as.numeric(eval(x[[4]])),
                                                   fix=nenv$do.fixed,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    }
                } else if (length(x) == 3){
                    if (any(tolower(as.character(x[[3]])) == c("fix", "fixed")) &&
                        length(x[[3]]) == 1){
                        ## c(1, fixed)
                        env$theta <- env$theta + 1
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=NA,
                                                   lower= -Inf,
                                                   est=as.numeric(eval(x[[2]])),
                                                   upper=Inf,
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    } else {
                        ## c(1,2)
                        env$theta <- env$theta + 1
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=NA,
                                                   lower=as.numeric(eval(x[[2]])),
                                                   est=as.numeric(eval(x[[3]])),
                                                   upper=Inf,
                                                   fix=nenv$do.fixed,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    }
                } else if (length(x) == 2){
                    ## c(1)
                    env$theta <- env$theta + 1
                    env$df <- rbind(env$df,
                                    data.frame(ntheta=env$theta,
                                               neta1=NA,
                                               neta2=NA,
                                               name=NA,
                                               lower=-Inf,
                                               est=as.numeric(eval(x[[2]])),
                                               upper=Inf,
                                               fix=FALSE,
                                               err=NA,
                                               label=NA,
                                               condition=NA));
                }
            } else if (identical(x[[1]], quote(`~`))){
                if (length(x) == 3){
                    if (length(x[[3]]) == 1){
                        ## et1 ~ 0.2
                        env$netas <- 1;
                        env$eta1 <- env$eta1 + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=NA,
                                                   neta1=env$eta1,
                                                   neta2=env$eta1,
                                                   name=as.character(x[[2]]),
                                                   lower=-Inf,
                                                   est=as.numeric(eval(x[[3]])),
                                                   upper=Inf,
                                                   fix=FALSE,
                                                   err=NA,
                                                   label=NA,
                                                   condition="ID"));
                    } else {
                        ## et1+et2+et3~c() lower triangular matrix
                        if (any(tolower(as.character(x[[3]][[1]])) == c("c", "fix", "fixed"))){
                            full.fixed <- any(tolower(as.character(x[[3]][[1]])) == c("fix", "fixed"));
                            env$netas <- length(x[[3]]) - 1;
                            num <- sqrt(1+env$netas*8)/2-1/2
                            if (round(num) == num){
                                n <- unlist(strsplit(as.character(x[[2]]), " +[+] +"));
                                n <- n[n != "+"];
                                if(length(n) == num){
                                    r <- x[[3]][-1];
                                    r <- t(sapply(r, function(x){
                                        nenv$do.fixed <- FALSE;
                                        ret <- as.numeric(eval(x));
                                        return(c(v=ret, do.fixed=nenv$do.fixed));
                                    }));
                                    r <- data.frame(r)
                                    r$do.fixed <- as.logical(r$do.fixed);
                                    i <- 0
                                    j <- 1;
                                    for (k in seq_along(r$do.fixed)){
                                        v <- r$v[k];
                                        do.fixed <- r$do.fixed[k]
                                        i <- i + 1;
                                        if (i == j){
                                            nm <- n[i];
                                        } else {
                                            nm <- sprintf("(%s,%s)", n[j], n[i]);
                                        }
                                        env$df <- rbind(env$df,
                                                        data.frame(ntheta=NA,
                                                                   neta1=env$eta1 + j,
                                                                   neta2=env$eta1 + i,
                                                                   name=nm,
                                                                   lower=-Inf,
                                                                   est=v,
                                                                   upper=Inf,
                                                                   fix=full.fixed | do.fixed,
                                                                   err=NA,
                                                                   label=NA,
                                                                   condition="ID"))
                                        if (i == j){
                                            j <- j + 1;
                                            i <- 0;
                                        }
                                    }
                                    env$eta1 <- env$eta1 + num;
                                }  else {
                                    stop("The left handed side of the expression must match the number of ETAs in the lower triangular matrix.");
                                }
                            } else {
                                n <- unlist(strsplit(as.character(x[[2]]), " +[+] +"));
                                n <- n[n != "+"];
                                stop(sprintf("%s ~ c(%s) does not have the right dimensions for a lower triangular matrix.", paste(n, collapse=" + "), paste(sapply(x[[3]][-1], as.character), collapse=", ")))
                            }
                        } else if (any(as.character(x[[3]][[1]]) == c("add", "prop"))){
                            env$theta <- env$theta + 1;
                            env$nerr <- 1;
                            env$df <- rbind(env$df,
                                            data.frame(ntheta=env$theta,
                                                       neta1=NA,
                                                       neta2=NA,
                                                       name=as.character(x[[2]]),
                                                       lower=0,
                                                       est=as.numeric(eval(x[[3]][[2]])),
                                                       upper=Inf,
                                                       fix=TRUE,
                                                       err=as.character(x[[3]][[1]]),
                                                       label=NA,
                                                       condition=NA))
                        } else if (as.character(x[[3]][[1]]) == "+" && length(x[[3]])  == 3){
                            env$nerr <- 2;
                            env$theta <- env$theta + 1;
                            env$df <- rbind(env$df,
                                            data.frame(ntheta=env$theta,
                                                       neta1=NA,
                                                       neta2=NA,
                                                       name=as.character(x[[2]]),
                                                       lower=0,
                                                       est=as.numeric(eval(x[[3]][[2]][[2]])),
                                                       upper=Inf,
                                                       fix=TRUE,
                                                       err=as.character(x[[3]][[2]][[1]]),
                                                       label=NA,
                                                       condition=NA))
                            env$theta <- env$theta + 1;
                            env$df <- rbind(env$df,
                                            data.frame(ntheta=env$theta,
                                                       neta1=NA,
                                                       neta2=NA,
                                                       name=as.character(x[[2]]),
                                                       lower=0,
                                                       est=as.numeric(eval(x[[3]][[3]][[2]])),
                                                       upper=Inf,
                                                       fix=TRUE,
                                                       err=as.character(x[[3]][[3]][[1]]),
                                                       label=NA,
                                                       condition=NA))
                        }
                    }
                } else if (length(x) == 2){
                    if (length(x[[2]]) == 1){
                        ## ~ 0.1
                        env$netas <- 1;
                        env$eta1 <- env$eta1 + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=NA,
                                                   neta1=env$eta1,
                                                   neta2=env$eta1,
                                                   name=NA,
                                                   lower=-Inf,
                                                   est=as.numeric(eval(x[[2]])),
                                                   upper=Inf,
                                                   fix=FALSE,
                                                   err=NA,
                                                   label=NA,
                                                   condition="ID"))
                    } else {
                        ## ~c() lower triangular matrix
                        if (any(tolower(as.character(x[[2]][[1]])) == c("c", "fix", "fixed"))){
                            full.fixed <- any(tolower(as.character(x[[2]][[1]])) == c("fix", "fixed"))
                            env$netas <- length(x[[2]]) - 1;
                            num <- sqrt(1+env$netas*8)/2-1/2
                            if (round(num) == num){
                                r <- x[[2]][-1];
                                r <- t(sapply(r, function(x){
                                    nenv$do.fixed <- FALSE;
                                    ret <- as.numeric(eval(x));
                                    return(c(v=ret, do.fixed=nenv$do.fixed));
                                }));
                                r <- data.frame(r)
                                r$do.fixed <- as.logical(r$do.fixed);
                                i <- 0
                                j <- 1;
                                for (k in seq_along(r$do.fixed)){
                                    v <- r$v[k];
                                    do.fixed <- r$do.fixed[k]
                                    i <- i + 1;
                                    env$df <- rbind(env$df,
                                                    data.frame(ntheta=NA,
                                                               neta1=eta1 + j,
                                                               neta2=eta1 + i,
                                                               name=NA,
                                                               lower=-Inf,
                                                               est=v,
                                                               upper=Inf,
                                                               fix=full.fixed | do.fixed,
                                                               err=NA,
                                                               label=NA,
                                                               condition="ID"))
                                    if (i == j){
                                        j <- j + 1;
                                        i <- 0;
                                    }
                                }
                                env$eta1 <- env$eta1 + num;
                            } else {
                                stop(sprintf("~c(%s) does not have the right dimensions for a lower triangular matrix.", paste(sapply(x[[2]][-1], as.character), collapse=", ")))
                            }
                        }
                    }
                }
            } else if (identical(x[[1]], quote(`label`))) {
                lab <- as.character(x[[2]]);
                lab <- gsub(";+", "", lab)
                len <- length(env$df$ntheta);
                if (len > 0){
                    if (!is.na(env$df$ntheta[len])){
                        tmp <- as.character(env$df$label);
                        tmp[len] <- lab
                        env$df$label <- tmp
                    } else {
                        ## stop("Currently only thetas can be labeled");
                        ## warning("Currently only thetas can be labeled")
                    }
                }
            } else if (identical(x[[1]], quote(`condition`))){
                lab <- as.character(x[[2]]);
                lab <- gsub(";+", "", lab)
                len <- length(env$df$ntheta);
                if (!is.na(env$df$ntheta[len]) && is.na(env$df$err[len])){
                    stop("Currently only etas/errs can be conditioned on...");
                } else if (is.na(env$df$ntheta[len])){
                    cnd <- as.character(env$df$condition);
                    cnd[seq(len - env$netas + 1, len)] <- lab;
                    env$df$condition <- cnd;
                } else {
                    cnd <- as.character(env$df$condition);
                    cnd[seq(len - env$nerr + 1, len)] <- lab;
                    env$df$condition <- cnd;
                }
            } else {
                lapply(x, f, env=env)
            }
        } else if (is.pairlist(x)) {
            unique(unlist(lapply(x, f, env=env)))
        } else if (is.atomic(x)){
            ## a simple number
            ## 1
            env$theta <- env$theta + 1
            env$df <- rbind(env$df,
                            data.frame(ntheta=env$theta,
                                       neta1=NA,
                                       neta2=NA,
                                       name=NA,
                                       lower= -Inf,
                                       est=as.numeric(eval(x)),
                                       upper=Inf,
                                       fix=FALSE,
                                       err=NA,
                                       label=NA,
                                       condition=NA))
        } else {
            stop("Don't know how to handle type ", typeof(x),
                 call. = FALSE)
        }
    }
    env <- environment(f)
    f(body(fun2), env)
    if (length(df$ntheta) == 0){
        stop(sprintf("Could not find any parameter information."))
    }
    n <- df$name
    n <- n[!is.na(n)]
    if (any(duplicated(n))){
        dups <- unique(n[duplicated(n)])
        stop(sprintf("The following parameter names were duplicated: %s.", paste(dups, collapse=", ")))
    }
    w <- which(is.infinite(df$est));
    if (length(w) > 0){
        stop(sprintf("The following parameters initial estimates are infinite: %s", paste(df$name[w], collapse=", ")))
    }

    w <- which(df$lower == df$est && df$est == df$upper);
    if (length(w) > 0){
        stop(sprintf("The estimate, and upper and lower bounds are the same for the following parameters: %s\nTo fix parameters use %s=fix(%s) instead.",
                     paste(df$name[w], collapse=", "), df$name[w[1]], df$est[w[1]]))
    }
    w <- which(df$lower == df$est || df$est == df$upper);
    if (length(w) > 0){
        tmp <- unique(sort(c(df$est[w[1]], df$lower[w[1]], df$upper[w[1]])));
        tmp <- tmp[!is.infinite(tmp)];
        if (length(tmp) == 1){
            stop(sprintf("The estimate is the same as a boundary for the following parameter: %s\nInstead use %s=%s # est",
                         paste(df$name[w], collapse=", "), df$name[w[1]], tmp))
        } else {
            stop(sprintf("The estimate is the same as a boundary for the following parameters: %s\nInstead use %s=c(%s) # c(lower, est)",
                         paste(df$name[w], collapse=", "), df$name[w[1]], paste(tmp, collapse=", ")))
        }

    }
    w <- which(df$est >= df$upper);
    if (length(w) > 0){
        stop(sprintf("The bounds make no sense for these parameters: %s.\nThey should be ordered as follows: %s=c(%s) # c(lower, est, upper)",
                     paste(df$name[w], collapse=", "), df$name[w[1]], paste(sort(c(df$est[w[1]], df$lower[w[1]], df$upper[w[1]])), collapse=", ")))
}
    w <- which(df$lower >= df$est);
    if (length(w) > 0){
        stop(sprintf("The lower bound is higher than the estimate for these parameters: %s.\nYou can adjust by %s=c(%s, %s) # c(lower, est)",
                     paste(df$name[w], collapse=", "), df$name[w[1]], df$est[w[1]], df$lower[w[1]]))
    }
    class(df) <- c("nlmixrBounds", "data.frame");
    return(df)
}

is.nlmixrBounds <- function(x){
    should <- c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper",  "fix", "err", "label", "condition");
    what <- names(x)
    if (length(should) == length(what)){
        return(all(what == should))
    } else {
        return(FALSE)
    }
}

##' @export
as.data.frame.nlmixrBounds <- function(x, row.names = NULL, optional = FALSE, ...){
    cls <- class(x)
    cls <- cls[cls != "nlmixrBounds"]
    tmp <- x;
    class(tmp) <- cls;
    return(tmp);
}

##' @export
print.nlmixrBounds <- function(x, ...){
    message(paste0(crayon::bold("Fixed Effects")," (", crayon::bold$blue("$theta"), "):"));
    print(x$theta);
    omega <- x$omega;
    if (dim(omega)[1] > 0){
        message(paste0("\n", crayon::bold("Omega")," (", crayon::bold$blue("$omega"), "):"))
        print(omega)
    }
}

##'@export
`$.nlmixrBounds` <- function(obj, arg, exact = TRUE){
    m <- as.data.frame(obj);
    ret <- m[[arg, exact = exact]];
    if (is.null(ret)){
        if (arg == "theta"){
            return(nlmixrBoundsTheta(obj, full=FALSE))
        } else if (arg == "theta.full"){
            return(nlmixrBoundsTheta(obj, full=TRUE))
        } else if (arg == "omega"){
            return(nlmixrBoundsOmega(obj));
        } else if (arg == "random"){
            return(nlmixrBoundsOmega(obj, TRUE));
        } else if (arg == "fixed.form"){
            return(nlmixrBoundsTheta(obj, formula=TRUE))
        } else if (arg == "focei.upper"){
            return(nlmixrBounds.focei.upper.lower(obj, "upper"))
        } else if (arg == "focei.lower"){
            return(nlmixrBounds.focei.upper.lower(obj, "lower"))
        } else if (any(arg == c("theta.names", "focei.names"))){
            return(nlmixrBounds.focei.upper.lower(obj, "name"))
        } else if (arg == "focei.err.type"){
            return(nlmixrBounds.focei.upper.lower(obj, "err"))
        } else if (arg == "eta.names"){
                return(nlmixrBounds.eta.names(obj))
        } else {
            return(NULL)
        }
    } else {
        return(ret)
    }
}

##' Get ETA names
##'
##' @param obj UI object
##' @return ETA names
##' @author Matthew L. Fidler
nlmixrBounds.eta.names <- function(obj){
    df <- as.data.frame(obj);
    df <- df[!is.na(df$neta1), ];
    ## dft.unfixed <- dft[!dft$fix, ];
    return(paste(df[df$neta1 == df$neta2, "name"]))
}


##'@export
str.nlmixrBounds <- function(object, ...){
    str(as.data.frame(object), ...);
    cat(" $ theta     : num ... (theta estimates)\n")
    cat(" $ theta.full: num ... (theta estimates, including error terms)\n")
    cat(" $ omega     : matrix ... (omega matrix)\n")
    cat(" $ random    : matrix class ... (Based on Between Subject Random effects)\n")
    cat(" $ fixed.form: formula  ... (Fixed effect parameters based on theta.)\n")
    cat(" $ focei.upper: Upper bounds for FOCEi\n")
    cat(" $ focei.lower: Lower bounds for FOCEi\n")
    cat(" $ focei.err.type: Residual Error type for FOCEi thetas\n")
    cat(" $ eta.names: Eta names\n")
    cat(" $ focei.names: Theta names for FOCEi\n")
}

##' Get upper/lower/names for THETAs
##'
##' @param obj Bounds object
##' @param type type of object extracted
##' @return lower/upper/name vector
##' @author Matthew L. Fidler
nlmixrBounds.focei.upper.lower <- function(obj, type=c("upper", "lower", "name", "err")){
    type <- match.arg(type);
    df <- as.data.frame(obj);
    dft <- df[!is.na(df$ntheta), ];
    dft.unfixed <- dft[!dft$fix, ];
    ret <- dft.unfixed[[type]]
    if (is(ret, "factor")){
        ret <- paste(ret);
    }
    return(ret)
}

nlmixrBoundsTheta <- function(x, full=TRUE, formula=FALSE){
    if (is.nlmixrBounds(x)){
        x <- as.data.frame(x)
        if (formula) full <- FALSE;
        w <- which(!is.na(x$ntheta));
        tmp <- (x[w, ]);
        nm <- sprintf(ifelse(formula, ".theta.%d", "theta[%d]"), seq_along(w));
        w <- which(!is.na(tmp$name));
        nm[w] <- as.character(tmp$name[w]);
        w <- which(!is.na(tmp$err));
        theta <-  tmp$est
        if (length(w) > 0){
            if (full){
                nm[w] <- sprintf("err[%d]", seq_along(w));
            } else {
                nm <- nm[-w];
                theta <- theta[-w];
            }
        }
        if (formula){
            if (any(duplicated(nm))){
                stop("Duplicated names for thetas; Cannot figure out formula");
            } else {
                return(as.formula(sprintf("%s ~ 1", paste(nm, collapse=" + "))))
            }
        } else {
            names(theta) <- nm;
        }
        return(theta)
    } else {
        return(NULL)
    }
}

nlmixrBoundsOmega <- function(x, nlme=FALSE){
    if (is.nlmixrBounds(x)){
        w <- which(!is.na(x$neta1));
        if (length(w) > 0){
            d <- max(x$neta1);
            df <- x[w, ];
            mx <- max(df$neta1);
            mat <- matrix(rep(0, mx * mx), mx, mx);
            diag <- TRUE;
            for (i in seq_along(df$neta1)){
                neta1 <- df$neta1[i];
                neta2 <- df$neta2[i];
                if (neta1 == neta2){
                    mat[neta1, neta2] <- df$est[i];
                } else {
                    diag <- FALSE
                    mat[neta1, neta2] <- mat[neta2, neta1]<- df$est[i];
                }
            }
            if (is(nlme, "logical") && nlme){
                df.diag <- df[df$neta1 == df$neta2, ];
                n2 <- sprintf(".eta.%d", seq_along(df.diag$name));
                w <- which(!is.na(df.diag$name))
                n2[w] <- as.character(df.diag$name[w])
                frm <- as.formula(paste(paste(n2, collapse=" + "), "~ 1"))
                if (diag){
                    return(nlme::pdDiag(mat, form=frm))
                } else {
                    return(nlme::pdSymm(as.matrix(Matrix::nearPD(mat)$mat), form=frm))
                }
            } else if (is(nlme, "list")){
                df.diag <- df[df$neta1 == df$neta2, ];
                class(df.diag) <- "data.frame"
                n2 <- sprintf(".eta.%d", seq_along(df.diag$name));
                w <- which(!is.na(df.diag$name))
                n2[w] <- sapply(as.character(df.diag$name[w]), function(x){
                    nlme[[x]]
                })
                frm <- as.formula(paste(paste(n2, collapse=" + "), "~ 1"))
                if (diag){
                    return(nlme::pdDiag(mat, form=frm))
                } else {
                    return(nlme::pdSymm(as.matrix(Matrix::nearPD(mat)$mat), form=frm))
                }
            }
            w <- which(df$neta1 == df$neta2);
            dimnames(mat) <- list(paste(df$name[w]), paste(df$name[w]));
            return(mat);
        } else {
            return(matrix(double(), 0, 0))
        }
    } else {
        return(matrix(double(), 0, 0));
    }
}
