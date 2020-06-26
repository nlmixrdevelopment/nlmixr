#' Generic function (should be hidden) for adding a covariate to a given string expression
#'
#' @param funstring the fun.txt file from nlmixr fit object
#' @param covariate the expression for the covariate to be added to the model
#'
#' @author Vipul Mann, Matthew Fidler
#'
addCov <- function(funstring, covariate) {
  gsub("(?<!\\()\\)", paste0(" + ", covariate, ")"), funstring, perl = TRUE)
}

addCov2 <- function(funstring, varName, covariate, theta) {
  f <- function(x, isCov = FALSE) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      if (isCov) {
        if (any(as.character(x) == theta)) {
          return(eval(parse(
            text = paste0("quote(", x, "+", covariate, ")")
          )))
        }
      }
      return(x)
    } else if (is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`{`))) {
        return(paste(unlist(lapply(x[-1], function(x) {
          gsub(" +", "", paste(deparse1(f(x, isCov)), collapse = ""))
        })), collapse = "\n"))
      }
      else if (identical(x[[1]], quote(`~`)) ||
               identical(x[[1]], quote(`=`)) ||
               identical(x[[1]], quote(`<-`))) {
        if (length(x[[2]]) == 1) {
          if (as.character(x[[2]]) == varName) {
            isCov <- TRUE
          }
        }
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      else {
        return(as.call(lapply(x, f, isCov = isCov)))
      }
    }
  }
  
  f(eval(parse(text = paste0(
    "quote({", funstring, "})"
  ))))
}


addCov3 <- function(funstring, varName, covariate, theta) {
  f <- function(x, isCov = FALSE, isLog = FALSE) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      if (isCov && isLog) {
        if (any(as.character(x) == theta)) {
          if (varName %in% list("cl")) {
            if (inherits(covariate, "list")) {
              return(eval(parse(
                text = paste0(
                  "quote(",
                  x,
                  "+",
                  "0.75*",
                  covariate[[2]],
                  "*",
                  "log(",
                  covariate[[1]],
                  ")",
                  ")"
                )
              )))
            }
            else {
              return(eval(parse(
                text = paste0(
                  "quote(",
                  x,
                  "+",
                  "0.75*",
                  "log(",
                  covariate,
                  ")",
                  ")"
                )
              )))
            }
          }
          else {
            if (inherits(covariate, "list")) {
              return(eval(parse(
                text = paste0(
                  "quote(",
                  x,
                  "+",
                  covariate[[2]],
                  "*",
                  "log(",
                  covariate[[1]],
                  ")",
                  ")"
                )
              )))
            }
            
            else {
              return(eval(parse(
                text = paste0("quote(", x, "+", "log(", covariate, ")", ")")
              )))
            }
          }
        }
      }
      # else if(isCov && !isLog){
      #   if (any(as.character(x) == theta)){
      #     return(f2(varName))
      #   }
      # }
      
      return(x)
    } else if (is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`{`))) {
        return(paste(unlist(lapply(x[-1], function(x) {
          gsub(" +", "", paste(deparse1(f(x, isCov)), collapse = ""))
        })), collapse = "\n"))
      }
      else if (identical(x[[1]], quote(`~`)) ||
               identical(x[[1]], quote(`=`)) ||
               identical(x[[1]], quote(`<-`))) {
        if (length(x[[2]]) == 1) {
          if (as.character(x[[2]]) == varName) {
            isCov <- TRUE
            
            # check if exp() form: boolean
            isLog <- grepl("exp", as.character(x[[3]])[[1]])
          }
        }
        return(as.call(lapply(
          x, f, isCov = isCov, isLog = isLog
        )))
      }
      else {
        return(as.call(lapply(
          x, f, isCov = isCov, isLog = isLog
        )))
      }
    }
  }
  
  f2 <- function(varName) {
    funstringLhsRhs <- strsplit(funstring, "(<-|=)")[[1]]
    funstringRhs <- funstringLhsRhs[2]
    funstringLhs <- funstringLhsRhs[1]
    if (inherits(covariate, "list")) {
      if (varName %in% list("cl")) {
        expr <-
          paste0(
            "(",
            funstringRhs,
            ")",
            "*",
            "(",
            covariate[[1]],
            "^",
            "(0.75*",
            covariate[[2]],
            ")",
            ")"
          )
      }
      else {
        expr <-
          paste0("(",
                 funstringRhs,
                 ")",
                 "*",
                 "(",
                 covariate[[1]],
                 "^",
                 covariate[[2]],
                 ")")
      }
    }
    
    else {
      expr <- paste0("(", funstringRhs, ")", "*", "(", covariate, ")")
    }
    
    return(paste0(funstringLhs, "<-", expr))
  }
  
  isLog <- grepl("\\bexp *\\(", funstring, perl = TRUE)
  if (!isLog) {
    f2(varName)
  }
  else {
    f(eval(parse(text = paste0(
      "quote({", funstring, "})"
    ))))
  }
}


addCov4 <- function(funstring, varName, covariate, theta, isLog) {
  f <- function(x, isCov = FALSE) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      if (isCov) {
        if (any(as.character(x) == theta)) {
          return(eval(parse(
            text = paste0("quote(", x, "+", covariate, ")")
          )))
        }
      }
      return(x)
    } else if (is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`{`))) {
        return(paste(unlist(lapply(x[-1], function(x) {
          gsub(" +", "", paste(deparse1(f(x, isCov)), collapse = ""))
        })), collapse = "\n"))
      }
      else if (identical(x[[1]], quote(`~`)) ||
               identical(x[[1]], quote(`=`)) ||
               identical(x[[1]], quote(`<-`))) {
        if (length(x[[2]]) == 1) {
          if (as.character(x[[2]]) == varName) {
            isCov <- TRUE
          }
        }
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      else {
        return(as.call(lapply(x, f, isCov = isCov)))
      }
    }
  }
  
  f(eval(parse(text = paste0(
    "quote({", funstring, "})"
  ))))
  
  f2 <- function(varName) {
    funstringLhsRhs <- strsplit(funstring, "(<-|=)")[[1]]
    funstringRhs <- funstringLhsRhs[2]
    funstringLhs <- funstringLhsRhs[1]
    
    expr <- paste0("(", funstringRhs, ")", "*", "(", covariate, ")")
    
    return(paste0(funstringLhs, "<-", expr))
  }
  
  if (!isLog) {
    f2(varName)
  }
  else{
    f(eval(parse(text = paste0(
      "quote({", funstring, "})"
    ))))
  }
}


#' Add a covariate to the respective variable in a model string
#'
#' @param funstring the fun.txt file from nlmixr fit object
#' @param varName a string of the variable name to which the covariate needs to be added
#' @param covariate the expression for the covariate to be added to the model
#'
#' @author Vipul Mann, Matthew Fidler
#'
addCovVar <- function(fitobject,
                      varName,
                      covariate,
                      norm = c("median", "mean"),
                      norm_type = c("mul", "div", "sub", "add"),
                      categorical = NULL,
                      initialEst = 0,
                      initialEstLB = -Inf,
                      initialEstUB = Inf,
                      assignToFitEnv = TRUE,
                      uidCol) {
  if (!inherits(fitobject, "nlmixrFitCore")) {
    stop("'fit' needs to be a nlmixr fit")
  }
  
  if (inherits(norm, "numeric")) {
    checkmate::assert_numeric(
      norm,
      len = 1,
      lower = .Machine$double.eps,
      null.ok = FALSE,
      any.missing = FALSE
    )
  }
  else {
    norm <- match.arg(norm)
  }
  
  norm_type <- match.arg(norm_type)
  
  if (missing(norm_type)) {
    norm_type <- "sub"
  }
  norm_ops <-
    list(
      "div" = `/`,
      "mul" = `*`,
      "sub" = `-`,
      "add" = `+`
    )
  normOp <- norm_ops[[norm_type]]
  
  funstring <- fitobject$uif$fun.txt
  funstringSplit <-
    unlist(strsplit(funstring, split = "\\\n")) # split the string at \n
  
  # # look for the string that needs to be modified
  idx <- grep(varName, funstringSplit)
  if (length(idx) >= 2) {
    stop("cannot update more than one variable at a time")
  }
  
  # Get the normalization value from data
  data <- getData(fitobject)
  
  if (missing(uidCol)) {
    # search the dataframe for a column name of 'ID'
    colNames <- colnames(data)
    colNamesLower <- tolower(colNames)
    if ("id" %in% colNamesLower) {
      uidCol <- colNames[which("id" %in% colNamesLower)]
    }
    else {
      uidCol <- "ID"
    }
  }
  if (covariate %in% colnames(data)) {
    if (inherits(norm, "numeric")) {
      normValVec <- norm
    }
    else {
      if (norm %in% "mean") {
        # mean of the mean values
        uids <- unlist(unname(data[uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[uidCol] == x,]
          normVal <- mean(unlist(datSlice[covariate]))
        })
        
        normValVec <- mean(unlist(normValVec))
      }
      else {
        # mean of the median values
        uids <- unlist(unname(data[uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[uidCol] == x,]
          normVal <- median(unlist(datSlice[covariate]))
        })
        normValVec <- mean(unlist(normValVec))
      }
    }
    
    # check if log transform required
    
    isLog <- grepl("\\bexp *\\(", funstringSplit[[idx]], perl = TRUE)
    
    res = performNorm(data=data,
                      funstring=funstringSplit[[idx]],
                      covariate=covariate,
                      varName=varName,
                      normValVec=normValVec,
                      normOp=normOp,
                      isLog=isLog,
                      isCat = categorical)
    
    data = res[[1]]
    covNameMod = res[[2]]
    covNames = res[[3]]
    
    funstringSplit[[idx]] <-
      addCov4(funstringSplit[[idx]],
              varName,
              covNameMod,
              names(fitobject$theta),
              isLog)
    
    # if (is.null(categorical)) {
    #   # non-categorical variable
    #   covNames <- list(paste0("cov_", covariate))
    #   
    #   res = performNorm(data=data,
    #                      funstring=funstringSplit[[idx]],
    #                      covariate=covariate,
    #                      varName=varName,
    #                      normValVec=normValVec,
    #                      normOp=normOp,
    #                      isLog=isLog,
    #                      isCat = categorical)
    #   
    #   data = res[[1]]
    #   covNameMod = res[[2]]
    # 
    #   funstringSplit[[idx]] <-
    #     addCov3(funstringSplit[[idx]],
    #             varName,
    #             covNameMod,
    #             names(fitobject$theta))
    # }
    # else {
    #   # categorical variable
    #   datColName <- paste0("categorical_", covariate)
    #   res <- makeDummies(data, covariate = covariate)
    #   data <- res[[1]]
    #   covModExpr <- res[[2]]
    #   covNames <- res[[3]]
    #   
    #   covNameMod <- paste(covModExpr, collapse = "+")
    #   
    #   funstringSplit[[idx]] <-
    #     addCov3(funstringSplit[[idx]],
    #             varName,
    #             covNameMod,
    #             names(fitobject$theta))
    # }
    # cli::cli_alert_success("computed {datColName} performing normalization with {norm}")
    # cli::cli_alert_success(cli::col_blue("added '{cli::style_bold(datColName)}' to data"))
  }
  
  cli::cli_alert_success("added {covNameMod} to {varName}'s equation in the model")
  cli::cli_alert_success("updated value: {funstringSplit[[idx]]}")
  
  updatedMod <- initializeCovars(fitobject,
                                 funstringSplit[[idx]],
                                 covNames,
                                 initialEst,
                                 initialEstLB,
                                 initialEstUB)
  
  list(updatedMod, data) # return updated model and updated data
}


performNorm = function(data,
                       funstring,
                       covariate,
                       varName,
                       normValVec,
                       normOp,
                       isLog = NULL,
                       isCat = NULL) {
  
  if (is.null(isCat)){
    datColNames <- paste0("centered_", covariate)
    data[, datColNames] <-
      normOp(unname(unlist(data[, covariate])), normValVec)
    
    covNameMod1 <- datColNames
    covNameParam1 <- paste0("cov_", covariate)
    covNameMod = paste0(covNameMod1 , '*', covNameParam1)
    covNames= covNameParam1
  }
  
  else{ # categorical variable
    # datColNames <- paste0("categorical_", covariate)
    res <- makeDummies(data, covariate = covariate)
    data <- res[[1]]
    covModExpr <- res[[2]]
    covNames <- res[[3]]
    datColNames <- res[[4]]
    
    covNameMod <- paste(covModExpr, collapse = "+")
    
    return (list(data, covNameMod, covNames))
  }
  
  if (isLog){
    if(varName %in% c('cl')){  # with 0.75 prefactor
      for (datColName in datColNames){  # for loop to handle both non-categorical and categorical vars
        data[,datColName] = 0.75 * log(data[, datColName])
      }
    }
    else{
      for (datColName in datColNames){
        data[,datColName] = log(data[, datColName])
      }
    }
  }
  
  list(data, covNameMod, covNames)
}


updateModel <- function(prevModel, newFuntxt, ini) {
  class(ini) <- c("nlmixrBounds", "data.frame")
  .model <-
    eval(parse(
      text = paste0("function() {", newFuntxt, "}", sep = "\n"),
      keep.source = TRUE
    ))
  
  nlmixr:::.finalizeUiModel(nlmixr:::nlmixrUIModel(.model, ini, NULL),
                            as.list(prevModel$meta)) # add new initial estimates for the updated model?
}

initializeCovars <- function(fitobject,
                             fstring,
                             covNames,
                             initialEst,
                             initialEstLB,
                             initialEstUB) {
  updatedMod <- paste0("model(fitobject,{", fstring, "})")
  updatedMod <- eval(parse(text = updatedMod))
  
  ini2 <- as.data.frame(updatedMod$ini)
  for (covName in covNames) {
    ini2[ini2$name == covName, "est"] <- initialEst
    ini2[ini2$name == covName, "lower"] <- initialEstLB
    ini2[ini2$name == covName, "upper"] <- initialEstUB
  }
  
  class(ini2) <- c("nlmixrBounds", "data.frame")
  updatedMod$ini <- ini2
  
  updatedMod
}

makeDummies <- function(data, covariate) {
  v <- unlist(data[, covariate])
  
  # prefix2 = paste0('cov', ...) # for unknown param in the model
  # paste, collapse
  # expression, initial vals
  
  prefix <- paste0("categorical_", covariate, "_")
  s <- head(sort(unique(v)),-1) # remove the last column
  
  d <- outer(v, s, function(v, s) {
    1L * (v == s)
  })
  colnames(d) <- paste0(prefix, s)
  newdat <- cbind(data, d)
  
  prefix2 <- paste0("cov_", covariate, "_")
  covNames <- paste0(prefix2, s)
  
  covModExpr <- paste0(colnames(d), "*", covNames)
  list(newdat, covModExpr, covNames, colnames(d))
}

addCovMultiple <- function(covInfo, fitobject, parallel = TRUE) {
  # create directory to store 'fit' objects for the covariate search
  
  outputDir <-
    paste0("nlmixrCovariateSearchCache_", as.character(substitute(fitobject)))
  if (!dir.exists(outputDir)) {
    print(outputDir)
    dir.create(outputDir)
  }
  
  # adding covariates sequentially, independent of each other
  if (parallel) {
    mets <- lapply(covInfo, function(x) {
      res <- do.call(addCovVar, c(list(fitobject), x))
      updatedMod <- res[[1]]
      data <- res[[2]]
      
      fit2 <-
        nlmixr(updatedMod, data, est = getFitMethod(fitobject))
      
      fnamefit2 <-
        paste0(
          outputDir,
          "/",
          as.character(substitute(fitobject)),
          "_",
          as.character(x$varName),
          "_",
          as.character(x$covariate),
          ".RData",
          sep = ""
        )
      saveRDS(fit2, fnamefit2)
      cli::cli_h1("Metrics for CovFit: AIC: {fit2$AIC}, BIC: {fit2$BIC}, OBJF: {fit2$OBJF}")
    })
  }
  
  # adding covariates one after the other, appending to the previous model
  else {
    covsAdded <-
      list() # to keep track of covariates added and store in a file
    covsAddedIdx <- 1
    for (x in covInfo) {
      res <- do.call(addCovVar, c(list(fitobject), x))
      updatedMod <- res[[1]]
      data <- res[[2]]
      
      if (length(covsAdded) == 0) {
        covsAdded[[covsAddedIdx]] <- paste0(x$varName, '_', x$covariate)
        covsAddedIdx <- covsAddedIdx + 1
      }
      else {
        covsAdded[[covsAddedIdx]] <-
          paste0(covsAdded[[covsAddedIdx - 1]], "_", x$varName, '_', x$covariate)
        covsAddedIdx <- covsAddedIdx + 1
      }
      
      fit2 <-
        nlmixr(updatedMod, data, est = getFitMethod(fitobject))
      
      print(fit2$fun.txt)
      cli::cli_h1("Metrics for CovFit: AIC: {fit2$AIC}, BIC: {fit2$BIC}, OBJF: {fit2$OBJF}")
      
      fnamefit2 <-
        paste0(outputDir, "/", covsAdded[[covsAddedIdx - 1]],
               ".RData",
               sep = "")
      print(fnamefit2)
      saveRDS(fit2, fnamefit2)
      
      
      fitobject <- fit2
    }
  }
}