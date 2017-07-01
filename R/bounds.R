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
    labels <- gsub(".*# *(.*) *$", "\\1", fun2[w]);
    labels <- sapply(labels,
                     function(x){
        return(sprintf("label(%s)", paste0(deparse(x))))});
    fun2[w] <- paste0(fun2[w], "\n", labels);
    fun2 <- eval(parse(text=paste0(fun2, collapse = "\n")))
    theta <- 0;
    eta1 <- 0;
    df <- data.frame(theta=numeric(),
                     eta1=numeric(),
                     eta2=numeric(),
                     name=character(),
                     lower=numeric(),
                     est=numeric(),
                     upper=numeric(),
                     fix=logical(),
                     label=character());
    f <- function(x) {
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
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=as.character(x[[2]]),
                                                lower=as.numeric(x[[3]][[2]]),
                                                est=as.numeric(x[[3]][[3]]),
                                                upper=as.numeric(x[[3]][[4]]),
                                                fix=TRUE,
                                                label=NA));
                    } else {
                        stop(sprintf("%s %s c(%s) syntax is not supported for thetas", as.character(x[[2]]),
                                     as.character(x[[1]]), paste(sapply(x[[3]][-1], as.character), collapse=", ")))
                    }
                } else if (length(x[[3]]) == 4 &&
                    as.character(x[[3]][[1]]) == "c"){
                    if (any(tolower(as.character(x[[3]][[4]])) == c("fix", "fixed"))){
                        ## a = c(1,2,fixed)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=as.character(x[[2]]),
                                                lower=as.numeric(x[[3]][[2]]),
                                                est=as.numeric(x[[3]][[3]]),
                                                upper=Inf,
                                                fix=TRUE,
                                                label=NA));
                    } else {
                        ## a = c(1,2,3)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=as.character(x[[2]]),
                                                lower=as.numeric(x[[3]][[2]]),
                                                est=as.numeric(x[[3]][[3]]),
                                                upper=as.numeric(x[[3]][[4]]),
                                                fix=FALSE,
                                                label=NA));
                    }
                } else if (length(x[[3]]) == 3 &&
                           as.character(x[[3]][[1]]) == "c"){
                    if  (any(tolower(as.character(x[[3]][[3]])) == c("fix", "fixed"))) {
                        ## a = c(1,fixed)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=as.character(x[[2]]),
                                                lower= -Inf,
                                                est=as.numeric(x[[3]][[2]]),
                                                upper=Inf,
                                                fix=TRUE,
                                                label=NA));
                    }  else {
                        ## a = c(1,2)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=as.character(x[[2]]),
                                                lower=as.numeric(x[[3]][[2]]),
                                                est=as.numeric(x[[3]][[3]]),
                                                upper=Inf,
                                                fix=FALSE,
                                                label=NA));
                    }
                } else if (length(x[[3]]) == 2 &&
                           as.character(x[[3]][[1]]) == "c"){
                    ## a = c(1)
                    theta <<- theta + 1;
                    df <<- rbind(df,
                                 data.frame(theta=theta,
                                            eta1=NA,
                                            eta2=NA,
                                            name=as.character(x[[2]]),
                                            lower=-Inf,
                                            est=as.numeric(x[[3]][[2]]),
                                            upper=Inf,
                                            fix=FALSE,
                                            label=NA))
                } else if (length(x[[3]]) == 1){
                    ## a = 1
                    theta <<- theta + 1;
                    df <<- rbind(df,
                                 data.frame(theta=theta,
                                            eta1=NA,
                                            eta2=NA,
                                            name=as.character(x[[2]]),
                                            lower=-Inf,
                                            est=as.numeric(x[[3]][[1]]),
                                            upper=Inf,
                                            fix=FALSE,
                                            label=NA))
                }
            } else if (identical(x[[1]], quote(`c`))){
                if (length(x) > 5){
                    stop(sprintf("c(%s) syntax is not supported for thetas", paste(sapply(x[-1], as.character), collapse=", ")))
                } else if (length(x) == 5){
                    if (any(tolower(as.character(x[[5]])) == c("fix", "fixed"))){
                        ## c(1,2,3,fixed)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=NA,
                                                lower=as.numeric(x[[2]]),
                                                est=as.numeric(x[[3]]),
                                                upper=as.numeric(x[[4]]),
                                                fix=TRUE,
                                                label=NA));
                    } else {
                        stop(sprintf("c(%s) syntax is not supported for thetas", paste(sapply(x[-1], as.character), collapse=", ")))
                    }
                } else if (length(x) == 4){
                    ## No assignment...
                    if (any(tolower(as.character(x[[4]])) == c("fix", "fixed"))){
                        ## c(1,2,fixed)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=NA,
                                                lower=as.numeric(x[[2]]),
                                                est=as.numeric(x[[3]]),
                                                upper=Inf,
                                                fix=TRUE,
                                                label=NA));
                    } else {
                        ## c(1,2,3)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=NA,
                                                lower=as.numeric(x[[2]]),
                                                est=as.numeric(x[[3]]),
                                                upper=as.numeric(x[[4]]),
                                                fix=FALSE,
                                                label=NA));
                    }
                } else if (length(x) == 3){
                    if (any(tolower(as.character(x[[3]])) == c("fix", "fixed"))){
                        ## c(1, fixed)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=NA,
                                                lower= -Inf,
                                                est=as.numeric(x[[2]]),
                                                upper=Inf,
                                                fix=TRUE,
                                                label=NA));
                    } else {
                        ## c(1,2)
                        theta <<- theta + 1;
                        df <<- rbind(df,
                                     data.frame(theta=theta,
                                                eta1=NA,
                                                eta2=NA,
                                                name=NA,
                                                lower=as.numeric(x[[2]]),
                                                est=as.numeric(x[[3]]),
                                                upper=Inf,
                                                fix=FALSE,
                                                label=NA));
                    }
                } else if (length(x) == 2){
                    ## c(1)
                    theta <<- theta + 1;
                    df <<- rbind(df,
                                 data.frame(theta=theta,
                                            eta1=NA,
                                            eta2=NA,
                                            name=NA,
                                            lower=-Inf,
                                            est=as.numeric(x[[2]]),
                                            upper=Inf,
                                            fix=FALSE,
                                            label=NA));
                }
            } else if (identical(x[[1]], quote(`~`))){
                if (length(x) == 3){
                    if (length(x[[2]]) == 1){
                        ## et1 ~ 0.2
                        eta1 <<- eta1 + 1;
                        df <<- rbind(df,
                                     data.frame(theta=NA,
                                                eta1=eta1,
                                                eta2=eta1,
                                                name=as.character(x[[2]]),
                                                lower=-Inf,
                                                est=as.numeric(x[[3]]),
                                                upper=Inf,
                                                fix=FALSE,
                                                label=NA));
                    } else {
                        ## et1+et2+et3~c() lower triangular matrix
                        if (as.character(x[[3]][[1]]) == "c"){
                            num <- sqrt(1+(length(x[[3]]) - 1)*8)/2-1/2
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
                                        df <<- rbind(df,
                                                     data.frame(theta=NA,
                                                                eta1=eta1 + j,
                                                                eta2=eta1 + i,
                                                                name=nm,
                                                                lower=-Inf,
                                                                est=v,
                                                                upper=Inf,
                                                                fix=FALSE,
                                                                label=NA))
                                        if (i == j){
                                            j <- j + 1;
                                            i <- 0;
                                        }
                                    }
                                    eta1 <<- eta1 + num;
                                }  else {
                                    stop("The left handed side of the expression must match the number of ETAs in the lower triangular matrix.");
                                }
                            } else {
                                stop(sprintf("%s does not have the right dimensions for a lower triangular matrix.", as.character(x[[2]])))
                            }
                        }
                    }
                } else if (length(x) == 2){
                    if (length(x[[2]]) == 1){
                        ## ~ 0.1
                        eta1 <<- eta1 + 1;
                        df <<- rbind(df,
                                     data.frame(theta=NA,
                                                eta1=eta1,
                                                eta2=eta1,
                                                name=NA,
                                                lower=-Inf,
                                                est=as.numeric(x[[2]]),
                                                upper=Inf,
                                                fix=FALSE,
                                                label=NA))
                    } else {
                        ## ~c() lower triangular matrix
                        if (as.character(x[[2]][[1]]) == "c"){
                            num <- sqrt(1+(length(x[[2]]) - 1)*8)/2-1/2
                            if (round(num) == num){
                                r <- x[[2]][-1];
                                r <- sapply(r, as.numeric);
                                i <- 0
                                j <- 1;
                                for (v in r){
                                    i <- i + 1;
                                    df <<- rbind(df,
                                                 data.frame(theta=NA,
                                                            eta1=eta1 + j,
                                                            eta2=eta1 + i,
                                                            name=NA,
                                                            lower=-Inf,
                                                            est=v,
                                                            upper=Inf,
                                                            fix=FALSE,
                                                            label=NA))
                                    if (i == j){
                                        j <- j + 1;
                                        i <- 0;
                                    }
                                }
                                eta1 <<- eta1 + num;
                            } else {
                                stop(sprintf("%s does not have the right dimensions for a lower triangular matrix.", as.character(x[[2]])))
                            }
                        }
                    }
                }
            } else if (identical(x[[1]], quote(`label`))) {
                lab <- as.character(x[[2]]);
                len <- length(df$theta);
                if (!is.na(df$theta[len])){
                    df$label[len] <<- lab
                } else {
                    stop("Currently only thetas can be labeled");
                }
            } else {
                lapply(x, f)
            }
        } else if (is.pairlist(x)) {
            unique(unlist(lapply(x, f)))
        } else if (is.atomic(x)){
            ## a simple number
            ## 1
            theta <<- theta + 1;
            df <<- rbind(df,
                         data.frame(theta=theta,
                                    eta1=NA,
                                    eta2=NA,
                                    name=NA,
                                    lower= -Inf,
                                    est=as.numeric(x),
                                    upper=Inf,
                                    fix=FALSE,
                                    label=NA))
        } else {
            stop("Don't know how to handle type ", typeof(x),
                 call. = FALSE)
        }
    }
    f(body(fun2))
    if (length(df$theta) == 0){
        stop("Could not find any parameter information.")
    }
    df
}
