##' Extract the Nlmixr bound information from a function.
##'
##' @param fun Function to extract bound information from.
##' @return a dataframe with bound information.
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
nlmixrBounds <- function(fun){
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
    w <- which(regexpr("^ *#.*$", fun2) != -1);
    if (length(w) > 0){
        fun2 <- fun2[-w];
    }
    w <- which(regexpr("#.*", fun2) != -1);
    if (length(w) > 0){
        labels <- gsub(".*# *(.*) *$", "\\1", fun2[w]);
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
            if ((identical(x[[1]], quote(`<-`)) ||
                 identical(x[[1]], quote(`=`))) &&
                is.name(x[[2]])) {
                ## These are theta assignments
                if (length(x[[3]]) > 5){
                    stop(sprintf("%s %s c(%s) syntax is not supported for thetas",
                    as.character(x[[2]]), as.character(x[[1]]), paste(sapply(x[[3]][-1], as.character), collapse=", ")))
                } else if (length(x[[3]]) == 5){
                    if (as.character(x[[3]][[1]]) == "c" &&
                        any(tolower(as.character(x[[3]][[5]])) == c("fix", "fixed"))){
                        ## a = c(1,2,3,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower=as.numeric(x[[3]][[2]]),
                                                   est=as.numeric(x[[3]][[3]]),
                                                   upper=as.numeric(x[[3]][[4]]),
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA))
                    } else {
                        stop(sprintf("%s %s c(%s) syntax is not supported for thetas", as.character(x[[2]]),
                                     as.character(x[[1]]), paste(sapply(x[[3]][-1], as.character), collapse=", ")))
                    }
                } else if (length(x[[3]]) == 4 &&
                    as.character(x[[3]][[1]]) == "c"){
                    if (any(tolower(as.character(x[[3]][[4]])) == c("fix", "fixed"))){
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                         data.frame(ntheta=env$theta,
                                                    neta1=NA,
                                                    neta2=NA,
                                                    name=as.character(x[[2]]),
                                                    lower=as.numeric(x[[3]][[2]]),
                                                    est=as.numeric(x[[3]][[3]]),
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
                                                   lower=as.numeric(x[[3]][[2]]),
                                                   est=as.numeric(x[[3]][[3]]),
                                                   upper=as.numeric(x[[3]][[4]]),
                                                   fix=FALSE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    }
                } else if (length(x[[3]]) == 3 &&
                           as.character(x[[3]][[1]]) == "c"){
                    if  (any(tolower(as.character(x[[3]][[3]])) == c("fix", "fixed"))) {
                        ## a = c(1,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=as.character(x[[2]]),
                                                   lower= -Inf,
                                                   est=as.numeric(x[[3]][[2]]),
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
                                                   lower=as.numeric(x[[3]][[2]]),
                                                   est=as.numeric(x[[3]][[3]]),
                                                   upper=Inf,
                                                   fix=FALSE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    }
                } else if (length(x[[3]]) == 2 &&
                           as.character(x[[3]][[1]]) == "c"){
                    ## a = c(1)
                    env$theta <- env$theta + 1;
                    env$df <- rbind(env$df,
                                    data.frame(ntheta=env$theta,
                                               neta1=NA,
                                               neta2=NA,
                                               name=as.character(x[[2]]),
                                               lower=-Inf,
                                               est=as.numeric(x[[3]][[2]]),
                                               upper=Inf,
                                               fix=FALSE,
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
                                                est=as.numeric(x[[3]][[1]]),
                                                upper=Inf,
                                                fix=FALSE,
                                                err=NA,
                                                label=NA,
                                                condition=NA))
                }
            } else if (identical(x[[1]], quote(`c`))){
                if (length(x) > 5){
                    stop(sprintf("c(%s) syntax is not supported for thetas", paste(sapply(x[-1], as.character), collapse=", ")))
                } else if (length(x) == 5){
                    if (any(tolower(as.character(x[[5]])) == c("fix", "fixed"))){
                        ## c(1,2,3,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=NA,
                                                   lower=as.numeric(x[[2]]),
                                                   est=as.numeric(x[[3]]),
                                                   upper=as.numeric(x[[4]]),
                                                   fix=TRUE,
                                                   err=NA,
                                                   label=NA,
                                                   condition=NA));
                    } else {
                        stop(sprintf("c(%s) syntax is not supported for thetas", paste(sapply(x[-1], as.character), collapse=", ")))
                    }
                } else if (length(x) == 4){
                    ## No assignment...
                    if (any(tolower(as.character(x[[4]])) == c("fix", "fixed"))){
                        ## c(1,2,fixed)
                        env$theta <- env$theta + 1;
                        env$df <- rbind(env$df,
                                        data.frame(ntheta=env$theta,
                                                   neta1=NA,
                                                   neta2=NA,
                                                   name=NA,
                                                   lower=as.numeric(x[[2]]),
                                                   est=as.numeric(x[[3]]),
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
                                                lower=as.numeric(x[[2]]),
                                                est=as.numeric(x[[3]]),
                                                upper=as.numeric(x[[4]]),
                                                fix=FALSE,
                                                err=NA,
                                                label=NA,
                                                condition=NA));
                    }
                } else if (length(x) == 3){
                    if (any(tolower(as.character(x[[3]])) == c("fix", "fixed"))){
                        ## c(1, fixed)
                        env$theta <- env$theta + 1
                        env$df <- rbind(env$df,
                                     data.frame(ntheta=env$theta,
                                                neta1=NA,
                                                neta2=NA,
                                                name=NA,
                                                lower= -Inf,
                                                est=as.numeric(x[[2]]),
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
                                                lower=as.numeric(x[[2]]),
                                                est=as.numeric(x[[3]]),
                                                upper=Inf,
                                                fix=FALSE,
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
                                            est=as.numeric(x[[2]]),
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
                                                est=as.numeric(x[[3]]),
                                                upper=Inf,
                                                fix=FALSE,
                                                err=NA,
                                                label=NA,
                                                condition="ID"));
                    } else {
                        ## et1+et2+et3~c() lower triangular matrix
                        if (as.character(x[[3]][[1]]) == "c"){
                            env$netas <- length(x[[3]]) - 1;
                            num <- sqrt(1+env$netas*8)/2-1/2
                            if (round(num) == num){
                                n <- unlist(strsplit(as.character(x[[2]]), " +[+] +"));
                                n <- n[n != "+"];
                                if(length(n) == num){
                                    r <- x[[3]][-1];
                                    r <- sapply(r, as.numeric);
                                    i <- 0
                                    j <- 1;
                                    for (v in r){
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
                                                                fix=FALSE,
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
                                                       est=as.numeric(x[[3]][[2]]),
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
                                                       est=as.numeric(x[[3]][[2]][[2]]),
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
                                                       est=as.numeric(x[[3]][[3]][[2]]),
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
                                                est=as.numeric(x[[2]]),
                                                upper=Inf,
                                                fix=FALSE,
                                                err=NA,
                                                label=NA,
                                                condition="ID"))
                    } else {
                        ## ~c() lower triangular matrix
                        if (as.character(x[[2]][[1]]) == "c"){
                            env$netas <- length(x[[2]]) - 1;
                            num <- sqrt(1+env$netas*8)/2-1/2
                            if (round(num) == num){
                                r <- x[[2]][-1];
                                r <- sapply(r, as.numeric);
                                i <- 0
                                j <- 1;
                                for (v in r){
                                    i <- i + 1;
                                    env$df <- rbind(env$df,
                                                 data.frame(ntheta=NA,
                                                            neta1=eta1 + j,
                                                            neta2=eta1 + i,
                                                            name=NA,
                                                            lower=-Inf,
                                                            est=v,
                                                            upper=Inf,
                                                            fix=FALSE,
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
                len <- length(env$df$ntheta);
                if (!is.na(env$df$ntheta[len])){
                    tmp <- as.character(env$df$label);
                    tmp[len] <- lab
                    env$df$label <- tmp
                } else {
                    stop("Currently only thetas can be labeled");
                }
            } else if (identical(x[[1]], quote(`condition`))){
                lab <- as.character(x[[2]]);
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
                                    est=as.numeric(x),
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
        stop("Could not find any parameter information.")
    }
    class(df) <- c("nlmixrBounds", "data.frame");
    df
}
