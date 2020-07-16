#' Add covariate expression to a function string
#'
#' @param funstring a string giving the expression that needs to be modified
#' @param varName the variable to which the given string corresponds to in the model expression
#' @param covariate the covariate expression that needs to be added (at the appropriate place)
#' @param theta a list of names of the 'theta' parameters in the 'fit' object
#' @param isLog a boolean signifying the presence of log-transformation in the funstring
#' @return returns the modified string with the covariate added to function string
#' @author Vipul Mann, Matthew Fidler
#' @export
#' @noRd
addCovariate <-
  function(funstring,
           varName,
           covariate,
           theta,
           isLog) {
    f <- function(x, isCov = FALSE) {
      if (is.atomic(x)) {
        return(x)
      }
      else if (is.name(x)) {
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
      
      expr <-
        paste0("(", funstringRhs, ")", "*", "(", covariate, ")")
      expr <-
        gsub(" ", "", expr, perl = TRUE) # remove white spaces from the above string
      funstringLhs <-
        gsub(" ", "", funstringLhs, perl = TRUE) # remove white spaces from the above string
      return(paste0(funstringLhs, "<-", expr))
    }
    
    if (!isLog) {
      f2(varName)
    }
    else {
      f(eval(parse(text = paste0(
        "quote({", funstring, "})"
      ))))
    }
  }


#' Remove covariate expression from a function string
#'
#' @param funstring a string giving the expression that needs to be modified
#' @param varName the variable to which the given string corresponds to in the model expression
#' @param covariate the covariate expression that needs to be removed (from the appropriate place)
#' @param theta a list of names of the 'theta' parameters in the 'fit' object
#' @param isLog a boolean signifying the presence of log-transformation in the funstring
#'
#' @return returns the modified string with the covariate removed from the function string
#'
#' @author Vipul Mann, Matthew Fidler
#' @export
#' @noRd
#'
removeCovariate <- function(funstring, varName, covariate, theta) {
  covariate <-
    gsub(" ", "", covariate, perl = TRUE) # remove white spaces from the above string
  covariateSplit <- strsplit(covariate, "\\*|\\+")[[1]]
  
  f <- function(x, isCov = FALSE) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      # x is a name recognized by R eg. in-built functions, etc.
      # if (isCov) {  # the equaitons corresponds to the varName's equation where the covariate is added
      #   # if (any(as.character(x) == theta)) {
      #     if (any(as.character(x) == covariateSplit)){
      #         return(eval(parse(
      #         # text = paste0("quote(", 0,")")
      #           text = ""
      #       )))
      #
      #     }
      #
      #   # }
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
          }
        }
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      else if (isCov && identical(x[[1]], quote(`+`))) {
        # Case 1: + centered_factor*cov_factor
        if (length(x[[3]]) > 1 &&
            identical(quote(`*`), x[[3]][[1]])) {
          if (as.character(x[[3]][[2]]) %in% covariateSplit &&
              as.character(x[[3]][[3]]) %in% covariateSplit) {
            # return(x[[2]])
            if (length(x[[2]]) == 1) {
              return(x[[2]])
            }
            return(as.call(lapply(x[[2]], f, isCov = isCov)))
          }
        }
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      
      else if (isCov && identical(x[[1]], quote(`*`))) {
        # Case 2: *(centered_factor*cov_factor)
        if (length(x[[3]]) > 1 &&
            identical(x[[3]][[1]], quote(`(`))) {
          covExpr <- x[[3]][[2]]
          
          if (length(covExpr) > 1 &&
              identical(quote(`*`), covExpr[[1]])) {
            if (as.character(covExpr[[2]]) %in% covariateSplit &&
                as.character(covExpr[[3]]) %in% covariateSplit) {
              # return(x[[2]])
              return(as.call(lapply(x[[2]], f, isCov = isCov)))
            }
          }
          
          else if (length(covExpr) > 1 &&
                   identical(quote(`+`), covExpr[[1]])) {
            if (length(covExpr[[3]]) > 1 &&
                identical(quote(`*`), covExpr[[3]][[1]])) {
              if (as.character(covExpr[[3]][[2]]) %in% covariateSplit &&
                  as.character(covExpr[[3]][[3]]) %in% covariateSplit) {
                # return(x[[2]])
                return(as.call(lapply(x[[2]], f, isCov = isCov)))
              }
            }
          }
        }
        
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      
      
      else {
        return(as.call(lapply(x, f, isCov = isCov)))
      }
    }
  }
  
  nch <- 0
  ret <- funstring
  
  while (nch != nchar(ret)) {
    nch <- nchar(ret)
    
    ret <- f(eval(parse(text = paste0(
      "quote({", ret, "})"
    ))))
  }
  
  return(ret)
}


#' Adding covariate to a given variable in an nlmixr model expression
#'
#' @param fitobject an nlmixr 'fit' object
#' @param varName a string giving the variable name to which covariate needs to be added
#' @param covariate a string giving the covariate name; must be present in the data used for 'fit'
#' @param norm the kind of normalization to be used while normalizing covariates; must be either 'mean' or 'median'
#' @param norm_type a string defining operator to be used for transforming covariates using 'norm'; must be one among 'mul', 'div', 'sub', 'add'
#' @param categorical a boolean indicating if the 'covariate' is categorical
#' @param isHS a boolean indicating if 'covariate' is of Hockey-stick kind
#' @param initialEst the initial estimate for the covariate parameters to be estimated; default is 0
#' @param initialEstLB a lower bound for the covariate parameters to be estimated; default is -Inf
#' @param initialEstUB an upper bound for the covariate parameters to be estimated; default is Inf
#'
#' @return a list with the updated model expression and data with columns corresponding to normalized covaraite(s) appended
#' @author Vipul Mann, Matthew Fidler
#' @export
#' @noRd
#'
addCovVar <- function(fitobject,
                      varName,
                      covariate,
                      norm = c("median", "mean", 'autoscale'),
                      norm_type = c("mul", "div", "sub", "add", 'autoscale'),
                      categorical = FALSE,
                      isHS = FALSE,
                      initialEst = 0,
                      initialEstLB = -Inf,
                      initialEstUB = Inf) {
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
      "add" = `+`,
      "autoscale" = c(`-`, `/`)
    )
  normOp <- norm_ops[[norm_type]]
  
  funstring <- fitobject$uif$fun.txt
  funstringSplit <-
    unlist(strsplit(funstring, split = "\\\n")) # split the string at \n
  
  # # look for the string that needs to be modified
  # idx <- grep(varName, funstringSplit)
  idx <-
    grep(paste0("(?<!\\w)", varName), funstringSplit, perl = TRUE)  # varName not preceded by any other character
  
  if (length(idx) >= 2) {
    print(funstring)
    stop("cannot update more than one variable at a time")
  }
  
  # Get the normalization value from data
  data <- getData(fitobject)
  
  # search the dataframe for a column name of 'ID'
  colNames <- colnames(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uidCol <- colNames[which("id" %in% colNamesLower)]
  }
  else {
    uidCol <- "ID"
  }
  
  if (covariate %in% colnames(data)) {
    if (inherits(norm, "numeric")) {
      normValVec <- norm
    }
    else {
      if (norm %in% "autoscale") {
        normValVecMean = mean(data[, covariate])
        normValVecSd = sd(data[, covariate])
        
        normValVec <- list(normValVecMean, normValVecSd)
        
        # uids <- unlist(unname(data[uidCol]))
        # normValVecMean <- lapply(uids, function(x) {
        #   datSlice <- data[data[uidCol] == x, ]
        #   normVal <- mean(unlist(datSlice[covariate]))
        # })
        #
        # normValVecSd <- lapply(uids, function(x) {
        #   datSlice <- data[data[uidCol] == x, ]
        #   normVal <- sd(unlist(datSlice[covariate]))
        # })
        
      }
      
      else if (norm %in% "mean") {
        # mean of the mean values
        uids <- unlist(unname(data[uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[uidCol] == x, ]
          normVal <- mean(unlist(datSlice[covariate]))
        })
        
        normValVec <- mean(unlist(normValVec))
      }
      else {
        # mean of the median values
        uids <- unlist(unname(data[uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[uidCol] == x, ]
          normVal <- median(unlist(datSlice[covariate]))
        })
        normValVec <- mean(unlist(normValVec))
      }
    }
    
    # check if log transform required
    isLog <-
      grepl("\\bexp *\\(", funstringSplit[[idx]], perl = TRUE)
    
    res <- performNorm(
      data = data,
      covariate = covariate,
      varName = varName,
      normOp = normOp,
      normValVec = normValVec,
      isLog = isLog,
      isCat = categorical,
      isHS = isHS
    )
    
    data <- res[[1]]
    covNameMod <- res[[2]]
    covNames <- res[[3]]
    
    funstringSplit[[idx]] <-
      addCovariate(funstringSplit[[idx]],
                   varName,
                   covNameMod,
                   names(fitobject$theta),
                   isLog)
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


removeCovVar <- function(fitobject,
                         varName,
                         covariate,
                         categorical = FALSE,
                         isHS = FALSE) {
  if (!inherits(fitobject, "nlmixrFitCore")) {
    stop("'fit' needs to be a nlmixr fit")
  }
  
  funstring <- fitobject$uif$fun.txt
  funstringSplit <-
    unlist(strsplit(funstring, split = "\\\n")) # split the string at \n
  
  # # look for the string that needs to be modified
  idx <-
    grep(paste0("(?<!\\w)", varName), funstringSplit, perl = TRUE)  # varName not preceded by any other character
  
  if (length(idx) >= 2) {
    print(funstring)
    stop("cannot remove more than one covariate at a time")
  }
  
  if (!isHS && !categorical) {
    # not HS, not CAT
    covNameMod <-
      paste0(paste0("centered_", covariate),
             "*",
             paste0("cov_", covariate))
  }
  
  else if (isHS) {
    # HS
    
    prefix <- paste0("centered_", covariate, "_")
    prefix2 <- paste0("cov_", covariate, "_")
    s <- c("lower", "upper")
    covModExpr <-
      paste0(paste0(prefix, s), "*", paste0(prefix2, s))
    covNameMod <- paste(covModExpr, collapse = "+")
    
  }
  
  else if (categorical) {
    # CAT
    prefix <- paste0("categorical_", covariate, "_")
    prefix2 <- paste0("cov_", covariate, "_")
    s <- c("lower", "upper")
    covModExpr <-
      paste0(paste0(prefix, s), "*", paste0(prefix2, s))
    covNameMod <- paste(covModExpr, collapse = "+")
    
  }
  
  else{
    stop("aborting...unknown covariate type", call. = FALSE)
  }
  
  funstringSplit[[idx]] <-
    removeCovariate(funstringSplit[[idx]],
                    varName,
                    covNameMod,
                    names(fitobject$theta))
  
  
  cli::cli_alert_success("removed {covNameMod} from {varName}'s equation in the model")
  cli::cli_alert_success("updated funcition text: {funstringSplit[[idx]]}")
  
  updatedMod = paste0("model(fitobject,{", funstringSplit[[idx]], "})")
  updatedMod <- eval(parse(text = updatedMod))
  
  updatedMod
}

#' Perform normalization of the covariate
#'
#' @param data a dataframe consisting the covariates added
#' @param covariate a string giving the covariate name; must be present in the data used for 'fit'
#' @param varName the variable name to which the covariate is being added
#' @param normOp an operator indicating the kind transformation to be done on the covariate
#' @param normValVec a numeric value to be used for normalization of the covariate
#' @param isLog a boolean indicating the presence of log-transformation in the funstring; default is FALSE
#' @param isCat a boolean indicating if the covariate is categorical; default is FALSE
#' @param isHS a boolean indicating if the covariate is of Hockey-stick kind; default is FALSE
#'
#' @return a list comprising the update dataframe, the expression for covariate, and a list of covariate names
#' @export
#' @author Vipul Mann, Matthew Fidler
#' @noRd
performNorm <- function(data,
                        covariate,
                        varName,
                        normOp,
                        normValVec,
                        isLog = FALSE,
                        isCat = FALSE,
                        isHS = FALSE) {
  if (!(isCat)) {
    # not categorical variable
    if (!(isHS)) {
      # not hockey stick
      datColNames <- paste0("centered_", covariate)
      
      if (length(normOp) > 1) {
        data[, datColNames] <-
          normOp[[1]](unname(unlist(data[, covariate])), normValVec[[1]])
        data[, datColNames] <-
          normOp[[2]](unname(unlist(data[, datColNames])), normValVec[[2]])
      }
      else{
        data[, datColNames] <-
          normOp(unname(unlist(data[, covariate])), normValVec)
      }
      
      covNameMod1 <- datColNames
      covNameParam1 <- paste0("cov_", covariate)
      covNameMod <- paste0(covNameMod1, "*", covNameParam1)
      covNames <- covNameParam1
    }
    
    else {
      res <- makeHockeyStick(data, covariate = covariate)
      
      data <- res[[1]]
      covModExpr <- res[[2]]
      covNames <- res[[3]]
      datColNames <- res[[4]]
      covNameMod <- paste(covModExpr, collapse = "+")
      return(list(data, covNameMod, covNames))
    }
  }
  
  else {
    # categorical variable
    # datColNames <- paste0("categorical_", covariate)
    res <- makeDummies(data, covariate = covariate)
    data <- res[[1]]
    covModExpr <- res[[2]]
    covNames <- res[[3]]
    datColNames <- res[[4]]
    
    covNameMod <- paste(covModExpr, collapse = "+")
    
    return(list(data, covNameMod, covNames))
  }
  
  if (isLog) {
    if (varName %in% c("cl")) {
      # with 0.75 prefactor
      for (datColName in datColNames) {
        # for loop to handle both non-categorical and categorical vars
        data[, datColName] <- 0.75 * log(data[, datColName])
      }
    }
    else {
      for (datColName in datColNames) {
        data[, datColName] <- log(data[, datColName])
      }
    }
  }
  
  list(data, covNameMod, covNames)
}

#' Initializing covariates before estimation
#'
#' @param fitobject an nlmixr 'fit' object
#' @param fstring a string giving the entire expression for the model function string
#' @param covNames  a list of covariate names (parameters) that need to be estimates
#' @param initialEst the initial estimate for the covariate parameters to be estimated; default is 0
#' @param initialEstLB a lower bound for the covariate parameters to be estimated; default is -Inf
#' @param initialEstUB an upper bound for the covariate parameters to be estimated; default is Inf
#'
#' @return updated model object with the modified initial values
#' @export
#' @author Vipul Mann, Matthew Fidler
#' @noRd
#'
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

#' Creating Hockey-stick covariates
#'
#' @param data a dataframe containing the dataset that needs to be used
#' @param covariate the covariate that needs to be converted to hockey-stick; must be present in the data
#'
#' @return a list of updated data with covariates added, an expression that needs to be added to the model expression, the list of covariate names, and the column names corresponding to the hockey-stick covariates
#' @export
#' @author Vipul Mann, Matthew Fidler
#' @noRd
makeHockeyStick <- function(data, covariate) {
  v <- unlist(data[, covariate])
  
  prefix <- paste0("centered_", covariate, "_")
  s <- c("lower", "upper")
  
  # create two columns for below and above the median
  med <- median(v)
  d <- list(v * 1L * (v < med), v * 1L * (v >= med))
  
  names(d) <- paste0(prefix, s)
  newdat <- cbind(data, d)
  
  prefix2 <- paste0("cov_", covariate, "_")
  covNames <- paste0(prefix2, s)
  
  covModExpr <- paste0(names(d), "*", covNames)
  list(newdat, covModExpr, covNames, colnames(d))
}


#' Create categorical covariates
#'
#' @param data a dataframe containing the dataset that needs to be used
#' @param covariate the covariate that needs to be converted to categorical; must be present in the data
#'
#' @return a list of updated data with covariates added, an expression that needs to be added to the model expression, the list of covariate names, and the column names corresponding to the categorical covariates
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
#' @noRd
makeDummies <- function(data, covariate) {
  v <- unlist(data[, covariate])
  
  prefix <- paste0("categorical_", covariate, "_")
  s <- head(sort(unique(v)), -1) # remove the last column
  
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

removeCovMultiple <- function(covInfo, fitobject) {
  covSearchRes = list()  # list to store fitobjects during the search
  
  # removing multiple covariates (independently)
  .env = environment()
  .env$covSearchRes = list()
  lapply(1:length(covInfo), function(idx) {
    x = covInfo[[idx]]
    updatedMod <- do.call(removeCovVar, c(list(fitobject), x))
    data <- getData(fitobject)
    
    fit2 <-
      suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
    
    .env$covSearchRes[[idx]] = list(fit2, paste0(x$covariate , "_", x$varName))
  })
  
  covSearchRes = .env$covSearchRes
  
}


#' Add multiple covariates to a given model, sequentially or all at once
#'
#' @param covInfo a list of lists containing information on the covariates that need to be tested
#' @param fitobject an nlmixr 'fit' object
#' @param indep a boolean indicating if the covariates should be added independently, or sequentially (append to the previous model); default is TRUE
#'
#' @export
#' @author Vipul Mann, Matthew Fidler
#' @noRd
#'
addCovMultiple <- function(covInfo, fitobject, indep = TRUE) {
  covSearchRes = list()  # list to store fitobjects during the search
  
  # create directory to store 'fit' objects for the covariate search
  outputDir <-
    paste0("nlmixrCovariateSearchCache_", as.character(substitute(fitobject)))
  if (!dir.exists(outputDir)) {
    print(outputDir)
    dir.create(outputDir)
  }
  
  # adding covariates independent of each other
  .env = environment()
  .env$covSearchRes = list()
  if (indep) {
    lapply(1:length(covInfo), function(idx) {
      x = covInfo[[idx]]
      res <- do.call(addCovVar, c(list(fitobject), x))
      updatedMod <- res[[1]]
      data <- res[[2]]
      
      fit2 <-
        suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
      
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
      # cli::cli_h1("Metrics for CovFit: AIC: {fit2$AIC}, BIC: {fit2$BIC}, OBJF: {fit2$OBJF}")
      
      .env$covSearchRes[[idx]] = list(fit2, paste0(x$covariate , "_", x$varName))
    })
    
    covSearchRes = .env$covSearchRes
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
        covsAdded[[covsAddedIdx]] <- paste0(x$covariate , "_", x$varName)
        fit2 <-
          suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
        covSearchRes[[covsAddedIdx]] = list(fit2, covsAdded[[covsAddedIdx]])
        covsAddedIdx <- covsAddedIdx + 1
      }
      else {
        covsAdded[[covsAddedIdx]] <-
          paste0(covsAdded[[covsAddedIdx - 1]], "_", x$covariate , "_", x$varName)
        fit2 <-
          suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
        covSearchRes[[covsAddedIdx]] = list(fit2, covsAdded[[covsAddedIdx]])
        covsAddedIdx <- covsAddedIdx + 1
      }
      
      print(fit2$fun.txt)
      # cli::cli_h1("Metrics for CovFit: AIC: {fit2$AIC}, BIC: {fit2$BIC}, OBJF: {fit2$OBJF}")
      
      fnamefit2 <-
        paste0(outputDir, "/", covsAdded[[covsAddedIdx - 1]],
               ".RData",
               sep = "")
      saveRDS(fit2, fnamefit2)
      fitobject <- fit2
    }
  }
  
  covSearchRes
}

#' Stepwise Covariate Model-selection (SCM) method
#'
#' @param fit an nlmixr 'fit' object
#' @param varsVec a list of candidate variables to which the covariates could be added
#' @param covarsVec a list of candidate covariates that need to be tested
#' @param covInformation a list containing information on the variables-covariates pairs
#' @param testAll a boolean indicating if all possible permutations between varsVec and covarsVec need to be tested
#'
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
#' @examples
covarSearchSCM <-
  function(fit,
           varsVec,
           covarsVec,
           covInformation = NULL,
           testAll = TRUE,
           forward = TRUE) {
    if (testAll) {
      possiblePerms <- expand.grid(varsVec, covarsVec)
      possiblePerms <-
        list(as.character(possiblePerms[[1]]),
             as.character(possiblePerms[[2]]))
      names(possiblePerms) <- c("vars", "covars")
    }
    else {
      possiblePerms <-
        list(names(relationsVarCovar), unlist(unname(relationsVarCovar)))
      names(possiblePerms) <- c("vars", "covars")
    }
    
    covInfo <- list()
    for (item in Map(list, possiblePerms$vars, possiblePerms$covars)) {
      listVarName <- paste0(item[[2]], "_", item[[1]])
      if (listVarName %in% names(covInformation)) {
        # search for name of var_covar in the list covInformation
        # covInfo[[length(covInfo) + 1]] <- c(list(varName = item[[1]], covariate = item[[2]]), covInformation[[listVarName]])
        covInfo[[listVarName]] = c(list(varName = item[[1]], covariate = item[[2]]), covInformation[[listVarName]])
      }
      else {
        # covInfo[[length(covInfo) + 1]] <- list(varName = item[[1]], covariate = item[[2]])
        covInfo[[listVarName]] = list(varName = item[[1]], covariate = item[[2]])
      }
    }
    
    if (forward) {
      cli::cli_h1('starting forward search...')
      stepIdx <- 1
      while (length(covInfo) > 0) {
        # forward covariate search
        covSearchRes = addCovMultiple(covInfo, fit, indep = TRUE)
        
        resTable = lapply(covSearchRes, function(res) {
          x = res[[1]]
          nam = res[[2]]
          list(nam, x$objf, x$objf - fit$objf, x$AIC, x$BIC)
        })
        
        resTable = data.frame(do.call(rbind, resTable))
        colnames(resTable) = c('varCovar', 'objf', 'deltObjf', 'AIC', 'BIC')
        
        bestRow = resTable[which.min(resTable$deltObjf),]
        
        if (bestRow$deltObjf < 0) {
          # objf function value improved
          cli::cli_h1('best model at step {stepIdx}: ')
          print(bestRow)
          
          fit <-
            covSearchRes[[which.min(resTable$deltObjf)]][[1]]  # extract fit object corresponding to the best model
          stepIdx <- stepIdx + 1
          covInfo[[as.character(bestRow$varCovar)]] = NULL
          
          cli::cli_h2('removed {bestRow$varCovar}')
        }
        else{
          # objf function value did not improve
          cli::cli_h1('objf value did not improve, exiting the search ...')
          break
        }
      }
      cli::cli_h2(cli::col_red('forward search complete'))
    }
    
    else{
      cli::cli_h1('starting backward search...')

      covSearchRes = addCovMultiple(covInfo, fit, indep = FALSE)
      fit = covSearchRes[[length(covSearchRes)]][[1]]
      
      cli::cli_h2(cli::col_blue('initial function text to remove from:'))
      cli::cli_text(cli::col_red("{fit$fun.txt}"))
      
      # Now remove covars step by step until the objf fun value
      stepIdx <- 1
      while (length(covInfo) > 1) {
        # Remove covars on by one: if objf val increases retain covar; otherwise (objf val decreases), remove the covar
        # At any stage, retain the one that results in highest increase in objf value; exit if removal of none results in increase
        covSearchRes = removeCovMultiple(covInfo, fit)
        
        resTable = lapply(covSearchRes, function(res) {
          x = res[[1]]
          nam = res[[2]]
          list(nam, x$objf, x$objf - fit$objf, x$AIC, x$BIC)
        })
        
        resTable = data.frame(do.call(rbind, resTable))
        colnames(resTable) = c('varCovar', 'objf', 'deltObjf', 'AIC', 'BIC')
        
        bestRow = resTable[which.max(resTable$deltObjf),]
        
        if (bestRow$deltObjf < 0) {
          # objf function value improved
          cli::cli_h1('best model at step {stepIdx}: ')
          print(bestRow)
          
          fit <-
            covSearchRes[[which.max(resTable$deltObjf)]][[1]]  # extract fit object corresponding to the best model
          stepIdx <- stepIdx + 1
          covInfo[[as.character(bestRow$varCovar)]] = NULL
          
          cli::cli_h2('retaining {bestRow$varCovar}')
        }
        else{
          # objf function value did not improve
          cli::cli_h1('objf value did not improve, exiting the search ...')
          break
        }
        
      }

      cli::cli_h2(cli::col_red('backward search complete'))
    }
    resTable
  }
