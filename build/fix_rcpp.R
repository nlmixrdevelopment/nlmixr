assignInNamespace(".parseLinkingTo",function(linkingTo) {

    if (is.null(linkingTo))
        return (character())

    linkingTo <- strsplit(linkingTo, "\\s*\\,")[[1]]
    result <- gsub("\\s", "", linkingTo)
    ret <- gsub("\\(.*", "", result);
    return(ret[ret != "RcppArmadillo"]);}
   ,"Rcpp");
