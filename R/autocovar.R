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
          if (any(as.character(x) == theta) && regexpr("^cov\\_", as.character(x)) == -1) {
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
#'
#' @return returns the modified string with the covariate removed from the function string
#'
#' @author Vipul Mann, Matthew Fidler
#' @export
#'
removeCovariate <- function(funstring, varName, covariate, theta) {
  covariate <-
    gsub(" ", "", covariate, perl = TRUE) # remove white spaces from the above string
  covariateSplit <- strsplit(covariate, "\\*|\\+")[[1]]

  f <- function(x, isCov = FALSE) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
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
#'
addCovVar <- function(fitobject,
                      varName,
                      covariate,
                      norm = c("median", "mean", "autoscale"),
                      norm_type = c("mul", "div", "sub", "add", "autoscale"),
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
    grep(paste0("(?<!\\w)", varName), funstringSplit, perl = TRUE) # varName not preceded by any other character

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

    # infer categorical covariates if not specified
    if (missing(categorical) && is.factor(data[, covariate])) {
      categorical <- TRUE
      cli::cli_alert_info("treating {covariate} as categorical variable ...")
    }

    if (inherits(norm, "numeric")) {
      normValVec <- norm
    }
    else {
      if (norm %in% "autoscale") {
        normValVecMean <- mean(data[, covariate])
        normValVecSd <- sd(data[, covariate])

        if (normValVecSd == 0) {
          # if (!categorical){  # print warning only for non-categorical variables
          #   cli::cli_alert_warning('normalization value for subject ID: {x} is zero; using with 1...')
          # }
          normValVecSd <- 1
        }

        normValVec <- list(normValVecMean, normValVecSd)
      }

      else if (norm %in% "mean") {
        # mean of the mean values
        uids <- unlist(unname(data[, uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[, uidCol] == x, ]
          normVal <- mean(unlist(datSlice[covariate]))
          if (normVal == 0) {
            if (!categorical && !all(unlist(datSlice[covariate]) == 0)) { # print warning only for non-categorical variables
              cli::cli_alert_warning("normalization value for subject ID: {x} is zero; using with 1...")
            }
            normVal <- 1
          }

          normVal
        })

        normValVec <- mean(unlist(normValVec))
      }
      else {
        # mean of the median values
        uids <- unlist(unname(data[, uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[, uidCol] == x, ]
          normVal <- median(unlist(datSlice[covariate]))
          if (normVal == 0) {
            if (!categorical && !all(unlist(datSlice[covariate]) == 0)) { # print warning only for non-categorical variables
              cli::cli_alert_warning("normalization value for subject ID: {x} is zero; using with 1...")
            }
            normVal <- 1
          }
          normVal
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

    # cli::cli_h1('theta: {names(fitobject$theta)}')
    # cli::cli_h1('previous value: {funstringSplit[[idx]]}')
    funstringSplit[[idx]] <-
      addCovariate(
        funstringSplit[[idx]],
        varName,
        covNameMod,
        names(fitobject$theta),
        isLog
      )
  }

  cli::cli_alert_success("added {covNameMod} to {varName}'s equation in the model")
  cli::cli_alert_success("updated value: {funstringSplit[[idx]]}")

  updatedMod <- initializeCovars(
    fitobject,
    funstringSplit[[idx]],
    covNames,
    initialEst,
    initialEstLB,
    initialEstUB
  )

  list(updatedMod, data, covNames) # return updated model and updated data
}


#' Remove covariate from function string
#'
#' Function to remove covariates from a given variable's equation in the function string text
#'
#' @param fitobject an nlmixr 'fit' object
#' @param varName a string giving the variable name to which covariate needs to be added
#' @param covariate a string giving the covariate name; must be present in the data used for 'fit'
#' @param categorical a boolean to represent if the covariate to be added is categorical
#' @param isHS a boolean to represent if the covariate to be added is hockey-stick normalized
#'
#' @return returns a list containing the updated model and the parameter names for the covariates added
#'
#' @author Vipul Mann, Matthew Fidler
#'
#' @export
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
    grep(paste0("(?<!\\w)", varName), funstringSplit, perl = TRUE) # varName not preceded by any other character

  if (length(idx) >= 2) {
    print(funstring)
    stop("cannot remove more than one covariate at a time")
  }

  if (!isHS && !categorical) {
    # not HS, not CAT
    covNameMod <-
      paste0(
        paste0("centered_", covariate),
        "*",
        paste0("cov_", covariate, "_", varName)
      )
    covNames <- paste0("cov_", covariate, "_", varName)
  }

  else if (isHS) {
    # HS
    prefix <- paste0("centered_", covariate, "_")
    prefix2 <- paste0("cov_", covariate, "_")
    s <- c("lower", "upper")
    covModExpr <-
      paste0(paste0(prefix, s), "*", paste0(prefix2, s, "_", varName))
    covNameMod <- paste(covModExpr, collapse = "+")

    covNames <- paste0(prefix2, s, "_", varName)
  }

  else if (categorical) {
    # CAT
    prefix <- paste0("categorical_", covariate, "_")
    prefix2 <- paste0("cov_", covariate, "_")

    v <- unlist(getData(fitobject)[, covariate])
    s <- head(sort(unique(v)), -1) # remove the last column

    covModExpr <-
      paste0(paste0(prefix, s), "*", paste0(prefix2, s, "_", varName))
    covNameMod <- paste(covModExpr, collapse = "+")

    covNames <- paste0(prefix2, s, "_", varName)
  }

  else {
    stop("aborting...unknown covariate type", call. = FALSE)
  }

  funstringSplit[[idx]] <-
    removeCovariate(
      funstringSplit[[idx]],
      varName,
      covNameMod,
      names(fitobject$theta)
    )


  cli::cli_alert_success("removed {covNameMod} from {varName}'s equation in the model")
  cli::cli_alert_success("updated function text: {funstringSplit[[idx]]}")

  updatedMod <- paste0("model(fitobject,{", funstringSplit[[idx]], "})")
  updatedMod <- eval(parse(text = updatedMod))

  list(updatedMod, covNames)
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
        if (normValVec[[2]] == 0) {
          normValVec <- 1
          cli::cli_alert_warning("replacing the normalization value from 0 to 1")
        }

        data[, datColNames] <-
          normOp[[1]](unname(unlist(data[, covariate])), normValVec[[1]])
        data[, datColNames] <-
          normOp[[2]](unname(unlist(data[, datColNames])), normValVec[[2]])
      }
      else {
        if (normValVec == 0) {
          normValVec <- 1
          cli::cli_alert_warning("replacing the normalization value from 0 to 1")
        }

        data[, datColNames] <-
          normOp(unname(unlist(data[, covariate])), normValVec)
      }

      covNameMod1 <- datColNames
      covNameParam1 <- paste0("cov_", covariate, "_", varName)
      covNameMod <- paste0(covNameMod1, "*", covNameParam1)
      covNames <- covNameParam1
    }

    else {
      res <- makeHockeyStick(data, covariate = covariate, varName = varName)

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
    res <- makeDummies(data, covariate = covariate, varName = varName)
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
        if (!all(is.finite(log(data[, datColName])))) {
          stop("non-finite values encountered in log-normalization. aborting...", call. = FALSE)
        }

        # for loop to handle both non-categorical and categorical vars
        data[, datColName] <- 0.75 * log(data[, datColName])
      }
    }
    else {
      for (datColName in datColNames) {
        if (!all(is.finite(log(data[, datColName])))) {
          stop("non-finite values encountered in log-normalization. aborting...", call. = FALSE)
        }

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
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return a list of updated data with covariates added, an expression that needs to be added to the model expression, the list of covariate names, and the column names corresponding to the hockey-stick covariates
#' @export
#' @author Vipul Mann, Matthew Fidler
makeHockeyStick <- function(data, covariate, varName) {
  v <- unlist(data[, covariate])

  prefix <- paste0("centered_", covariate, "_")
  s <- c("lower", "upper")

  # create two columns for below and above the median
  med <- median(v)
  d <- list(v * 1L * (v < med), v * 1L * (v >= med))

  names(d) <- paste0(prefix, s)
  newdat <- cbind(data, d)

  prefix2 <- paste0("cov_", covariate, "_")
  covNames <- paste0(prefix2, s, "_", varName)

  covModExpr <- paste0(names(d), "*", covNames)
  list(newdat, covModExpr, covNames, colnames(d))
}


#' Create categorical covariates
#'
#' @param data a dataframe containing the dataset that needs to be used
#' @param covariate the covariate that needs to be converted to categorical; must be present in the data
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return a list of updated data with covariates added, an expression that needs to be added to the model expression, the list of covariate names, and the column names corresponding to the categorical covariates
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
makeDummies <- function(data, covariate, varName) {
  v <- unlist(data[, covariate])

  prefix <- paste0("categorical_", covariate, "_")
  s <- head(sort(unique(v)), -1) # remove the last column

  d <- outer(v, s, function(v, s) {
    1L * (v == s)
  })
  colnames(d) <- paste0(prefix, s)
  newdat <- cbind(data, d)

  prefix2 <- paste0("cov_", covariate, "_")
  covNames <- paste0(prefix2, s, "_", varName)

  covModExpr <- paste0(colnames(d), "*", covNames)
  list(newdat, covModExpr, covNames, colnames(d))
}


#' Removing multiple covariates
#'
#' @param covInfo a list containing information about each variable-covariate pair
#' @param fitobject an nlmixr 'fit' object
#'
#' @return a list with the updated fit object, the variable-covariate pair string, and the parameter names for the corresponding covaraites removed
#' @export
#'
#' @author Vipul Mann, Matthew Fidler
#'
removeCovMultiple <- function(covInfo, fitobject) {
  covSearchRes <- list() # list to store fitobjects during the search

  # removing multiple covariates (independently)
  .env <- environment()
  .env$covSearchRes <- list()
  lapply(1:length(covInfo), function(idx) {
    x <- covInfo[[idx]]

    res <- do.call(removeCovVar, c(list(fitobject), x))
    updatedMod <- res[[1]]
    data <- getData(fitobject)
    covNames <- res[[2]]

    reassignVars <- rownames(fitobject$parFixedDf)[fitobject$parFixedDf$Estimate != fitobject$parFixedDf[, "Back-transformed"] & fitobject$parFixedDf[, "Back-transformed"] == 0]
    if (length(reassignVars) > 0) {
      ini2 <- as.data.frame(updatedMod$ini)

      for (r in reassignVars) {
        ini2[ini2$name == r, "est"] <- 1.0
        cli::cli_alert_warning("reasssigned initial value for {r} to 1.0")
      }

      class(ini2) <- c("nlmixrBounds", "data.frame")
      updatedMod$ini <- ini2
    }

    # fit2 <-
    #   suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))

    fit2 <- tryCatch(
      {
        fit2 <-
          suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
        fit2 # to return 'fit2'
      },
      error = function(error_message) {
        print("error fitting the model after REMOVING the covariates")
        print(error_message)
        print("returning the previous fitobject")
        return(fitobject) # return NA otherwise (instead of NULL)
      }
    )

    AIC(fit2)

    .env$covSearchRes[[idx]] <- list(fit2, c(x$covariate, x$varName), covNames)
  })

  covSearchRes <- .env$covSearchRes
}


#' Add multiple covariates to a given model, sequentially or all at once
#'
#' @param covInfo a list containing information about each variable-covariate pair
#' @param fitobject an nlmixr 'fit' object
#' @param indep a boolean indicating if the covariates should be added independently, or sequentially (append to the previous model); default is TRUE
#'
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
addCovMultiple <- function(covInfo, fitobject, indep = TRUE) {
  covSearchRes <- list() # list to store fitobjects during the search

  # adding covariates independent of each other
  .env <- environment()
  .env$covSearchRes <- list()
  if (indep) {
    lapply(1:length(covInfo), function(idx) {
      x <- covInfo[[idx]]
      res <- do.call(addCovVar, c(list(fitobject), x))
      updatedMod <- res[[1]]
      data <- res[[2]]
      covNames <- res[[3]]

      reassignVars <- rownames(fitobject$parFixedDf)[fitobject$parFixedDf$Estimate != fitobject$parFixedDf[, "Back-transformed"] & fitobject$parFixedDf[, "Back-transformed"] == 0]
      if (length(reassignVars) > 0) {
        ini2 <- as.data.frame(updatedMod$ini)

        for (r in reassignVars) {
          ini2[ini2$name == r, "est"] <- 1.0
          cli::cli_alert_warning("reasssigned initial value for {r} to 1.0")
        }

        class(ini2) <- c("nlmixrBounds", "data.frame")
        updatedMod$ini <- ini2
      }

      # fit2 <-
      #   suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))

      fit2 <- tryCatch(
        {
          fit2 <-
            suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
          fit2 # to return 'fit2'
        },
        error = function(error_message) {
          print("error fitting the model after ADDING the covariate")
          print(error_message)
          print("returning the previous fitobject")
          return(fitobject) # return NA otherwise (instead of NULL)
        }
      )

      AIC(fit2)

      .env$covSearchRes[[idx]] <- list(fit2, c(x$covariate, x$varName), covNames)
    })

    covSearchRes <- .env$covSearchRes
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
      covNames <- res[[3]]

      if (length(covsAdded) == 0) {
        # covsAdded[[covsAddedIdx]] <- paste0(x$covariate, x$varName)
        covsAdded[[covsAddedIdx]] <- c(x$covariate, x$varName)

        # fit2 <-
        #   suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))

        fit2 <- tryCatch(
          {
            fit2 <-
              suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
            fit2 # to return 'fit2'
          },
          error = function(error_message) {
            print("error fitting the model while SIMULTANEOUSLY ADDING covariates")
            print(error_message)
            print("skipping this covariate")
            return(fitobject) # return NA otherwise (instead of NULL)
          }
        )

        covSearchRes[[covsAddedIdx]] <- list(fit2, covsAdded[[covsAddedIdx]], covNames)
        covsAddedIdx <- covsAddedIdx + 1
      }
      else {
        # covsAdded[[covsAddedIdx]] <-
        #   paste0(covsAdded[[covsAddedIdx - 1]], "_", x$covariate, x$varName)
        covsAdded[[covsAddedIdx]] <-
          c(covsAdded[[covsAddedIdx - 1]], x$covariate, x$varName)

        # fit2 <-
        #   suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))

        fit2 <- tryCatch(
          {
            fit2 <-
              suppressWarnings(nlmixr(updatedMod, data, est = getFitMethod(fitobject)))
            fit2 # to return 'fit2'
          },
          error = function(error_message) {
            print("error fitting the model while SIMULTANEOUSLY ADDING covariates")
            print(error_message)
            print("skipping this covariate")
            return(fitobject) # return NA otherwise (instead of NULL)
          }
        )

        covSearchRes[[covsAddedIdx]] <- list(fit2, covsAdded[[covsAddedIdx]], covNames)
        covsAddedIdx <- covsAddedIdx + 1
      }

      print(fit2$fun.txt)

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
#' @param pVal a named list with names 'fwd' and 'bck' for specifying the p-values for the forward and backward searches, respectively
#' @param covInformation a list containing additionl information on the variables-covariates pairs that should be passed on to addCovMultiple function
#' @param catCovariates a list of covariates that should be treated as categorical
#' @param searchType one of 'scm', 'forward' and 'backward' to specify the covariate search method; default is 'scm'
#' @param restart a boolean that controls if the search should be restarted; default is FALSE
#'
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'     ## You may label each parameter with a comment
#'     tka <- 0.45 # Log Ka
#'     tcl <- log(c(0, 2.7, 100)) # Log Cl
#'     ## This works with interactive models
#'     ## You may also label the preceding line with label("label text")
#'     tv <- 3.45; label("log V")
#'     ## the label("Label name") works with all models
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' fit <- nlmixr(one.cmt, theo_sd,"focei")
#' RxODE::.rxWithWd(tempdir(), {# with temporary directory
#'
#' covarSearchAuto(fit, varsVec = c("ka", "cl"),
#'     covarsVec = c("WT", "SEX"), catCovariates = c("SEX"))
#'
#' })
#'
#' ## Note that this didn't include sex, add it to dataset and restart model
#'
#' d <- theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' fit <- nlmixr(one.cmt, d, "focei")
#'
#' # This would restart if for some reason the search crashed:
#'
#' RxODE::.rxWithWd(tempdir(), {# with temporary directory
#'
#' covarSearchAuto(fit, varsVec = c("ka", "cl"), covarsVec = c("WT", "SEX"),
#'                 catCovariates = c("SEX"), restart = TRUE)
#'
#' covarSearchAuto(fit, varsVec = c("ka", "cl"), covarsVec = c("WT", "SEX"),
#'                 catCovariates = c("SEX"), restart = TRUE,
#'                 searchType = "forward")
#' })
#'
#' }
#'
covarSearchAuto <- # unsuccessful runs info store; check for covInformation before resuming
  function(fit,
           varsVec,
           covarsVec,
           pVal = list(fwd = 0.05, bck = 0.01), # diff default vals for fwd and backward
           covInformation = NULL,
           catCovariates = NULL,
           searchType = c("scm", "forward", "backward"),
           restart = FALSE) {
    if (!is.numeric(AIC(fit))) {
      cli::cli_alert_danger("the 'fit' object needs to have an objective functions value associated with it")
      cli::cli_alert_info("try computing 'AIC(fitobject)' in console to compute and store the corresponding OBJF value")
      stop("aborting...objf value not associated with the current 'fit' object", call. = FALSE)
    }

    .origDat <- getData(fit)
    .new <- intersect(names(.origDat), covarsVec)
    if (length(.new) == 0L) stop("no covariates specified in original dataset")
    if (length(.new) != length(covarsVec)) {
      .old <- covarsVec[(covarsVec %in% .new)]
      warning("some covariates are not in the dataset: '", paste0(.old, "', '"), "'")
      covarsVec <- .new
    }
    if (inherits(catCovariates, "character")) {
      catCovariates <- intersect(catCovariates, covarsVec)
    }

    searchType <- match.arg(searchType)

    if (missing(searchType)) {
      searchType <- "scm"
    }

    if (!all((names(pVal) %in% c("fwd", "bck")))) {
      stop("pVal should be list of two with names  'fwd' and 'bck' ")
    }

    outputDir <-
      paste0("nlmixrCovariateSearchCache_", as.character(substitute(fit)), "_", digest::digest(fit)) # a new directory with this name will be created

    # outputDir <-
    #   paste0("nlmixrCovariateSearchCache_", as.character(substitute(fit)), "_", 'dbdf08b8c1b4cb9cb021c52008d3c343') # a new directory with this name will be created

    if (!dir.exists(outputDir)) {
      dir.create(outputDir)
    }
    else if (dir.exists(outputDir) && restart == TRUE) {
      unlink(outputDir, recursive = TRUE, force = TRUE) # unlink any of the previous directories
      dir.create(outputDir) # create a fresh directory
    }

    possiblePerms <- expand.grid(varsVec, covarsVec)
    possiblePerms <-
      list(
        as.character(possiblePerms[[1]]),
        as.character(possiblePerms[[2]])
      )
    names(possiblePerms) <- c("vars", "covars")


    covInfo <- list() # reversivle listVarName!!
    for (item in Map(list, possiblePerms$vars, possiblePerms$covars)) {
      listVarName <- paste0(item[[2]], item[[1]])

      if (item[[2]] %in% catCovariates) {
        covInformation[[listVarName]]$categorical <- TRUE
      }
      else {
        covInformation[[listVarName]]$categorical <- FALSE
      }

      if (listVarName %in% names(covInformation)) {
        covInfo[[listVarName]] <- c(list(varName = item[[1]], covariate = item[[2]]), covInformation[[listVarName]])
      }
      else {
        covInfo[[listVarName]] <- list(varName = item[[1]], covariate = item[[2]])
      }
    }

    if (searchType %in% "scm") {
      resFwd <- forwardSearch(covInfo, fit, pVal$fwd, outputDir = outputDir, restart = restart)
      resBck <- backwardSearch(covInfo, fitorig = fit, fitupdated = resFwd[[1]], pVal = pVal$bck, reFitCovars = FALSE, outputDir = outputDir, restart = restart)
      summaryTable <- Reduce(rbind, list(resFwd[[2]], resBck[[2]]))

      return(list(summaryTable = summaryTable, resFwd = resFwd, resBck = resBck))
    }

    else if (searchType %in% "forward") {
      resFwd <- forwardSearch(covInfo, fit, pVal = pVal$fwd, outputDir = outputDir, restart = restart)
      summaryTable <- Reduce(rbind, list(resFwd[[2]], NULL))

      return(list(summaryTable = summaryTable, resFwd = resFwd, resBck = NULL))
    }

    else {
      resBck <- backwardSearch(covInfo, fitorig = fit, pVal = pVal$bck, reFitCovars = TRUE, outputDir = outputDir, restart = restart)
      summaryTable <- Reduce(rbind, list(NULL, resBck[[2]]))

      return(list(summaryTable = summaryTable, resFwd = NULL, resBck = resBck))
    }
  }


#' Forward covariate search
#'
#' @param covInfo a list containing information about each variable-covariate pair
#' @param fit  an nlmixr 'fit' object
#' @param pVal p-value that should be used for selecting covariates in the forward search
#' @param outputDir the name of the output directory that stores the covariate search result
#' @param restart a boolean that controls if the search should be restarted; default is FALSE
#'
#' @return returns the updated 'fit' object at the end of the forward search and a table of information for all the covaraites tested
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
forwardSearch <- function(covInfo, fit, pVal = 0.05, outputDir, restart = FALSE) {
  if (missing(outputDir)) {
    stop("please specify output directory to store the results for forward search. aborting ...")
  }

  resTableComplete <- data.frame(matrix(ncol = 14, nrow = 0))
  cli::cli_h1("starting forward search...")
  stepIdx <- 1


  fnameTablePatternForward <-
    paste0("forward_step_", "[0-9]+", "_table_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )
  fnameFitPatternForward <-
    paste0("forward_step_", "[0-9]+", "_fit_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )
  fnameCompleteTablePatternForward <-
    paste0("forward_step_", "[0-9]+", "_completetable_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )

  fileExistsTab <-
    list.files(paste0("./", outputDir), pattern = fnameTablePatternForward)

  fileExistsFit <-
    list.files(paste0("./", outputDir), pattern = fnameFitPatternForward)

  fileExistsCompleteTable <-
    list.files(paste0("./", outputDir), pattern = fnameCompleteTablePatternForward)

  if (length(fileExistsTab) == 0) {
    restart <- TRUE
  }

  if (!restart) {
    resumeTable <- lapply(fileExistsTab, function(x) {
      readRDS(paste0(outputDir, "/", x))
    })

    resumeTable <- data.table::rbindlist(resumeTable)
    fit <- readRDS(paste0(outputDir, "/", fileExistsFit[[length(fileExistsFit)]]))
    resTableComplete <- readRDS(paste0(outputDir, "/", fileExistsCompleteTable[[length(fileExistsCompleteTable)]]))

    # update covInfo and step idx
    testedCovarVars <- paste0(unlist(resumeTable$covar), unlist(resumeTable$var))

    for (x in testedCovarVars) {
      covInfo[[x]] <- NULL
    }

    stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1

    cli::cli_alert_success("loaded forward search data from disk ...")
    cli::cli_alert_success("resuming forward search ...")
  }

  while (length(covInfo) > 0) {
    # forward covariate search
    covSearchRes <- addCovMultiple(covInfo, fit, indep = TRUE)

    resTable <- lapply(covSearchRes, function(res) {
      x <- res[[1]]
      nam_covar <- res[[2]][[1]]
      nam_var <- res[[2]][[2]]
      covNames <- res[[3]]

      # fwd: if deltObjf <0: pchisq=1-pchisq(-deltObjf, dof), else pchisq=1
      # bck: if deltObjf >0: pchisq=1-pchisq(deltObjf, dof), else pchisq=1

      dObjf <- fit$objf - x$objf
      dof <- length(x$uif$ini$est) - length(fit$uif$ini$est)
      if (dObjf < 0) {
        pchisqr <- 1 - pchisq(-dObjf, df = dof)
      }
      else {
        pchisqr <- 1
      }

      l1 <- list(step = stepIdx, covar = nam_covar, var = nam_var, objf = x$objf, deltObjf = dObjf, AIC = x$AIC, BIC = x$BIC, numParams = length(x$uif$ini$est), qchisqr = qchisq(1 - pVal, dof), pchisqr = pchisqr, included = "", searchType = "forward")
      l2 <- list(covNames = covNames, covarEffect = x$parFixedDf[covNames, "Estimate"])

      c(l1, l2)
    })

    resTable <- data.frame(do.call(rbind, resTable))
    colnames(resTable) <- c(names(resTable))
    bestRow <- resTable[which.min(resTable$pchisqr), ]

    colnames(resTableComplete) <- colnames(resTable)

    if (bestRow$pchisqr <= pVal) { # should be based on p-value

      # objf function value improved
      resTable[which.min(resTable$pchisqr), "included"] <- "yes"
      bestRow[, "included"] <- "yes"

      cli::cli_h1("best model at step {stepIdx}: ")
      print(bestRow)

      fit <-
        covSearchRes[[which.min(resTable$pchisqr)]][[1]] # extract fit object corresponding to the best model

      covInfo[[paste0(as.character(bestRow$covar), as.character(bestRow$var))]] <- NULL

      cli::cli_h2("excluding {paste0(as.character(bestRow$covar), as.character(bestRow$var))} from list of covariates ...")

      saveRDS(fit, file = paste0(outputDir, "/", "forward_", "step_", stepIdx, "_", "fit", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      saveRDS(bestRow, file = paste0(outputDir, "/", "forward_", "step_", stepIdx, "_", "table", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      resTableComplete <- rbind(resTableComplete, resTable)
      saveRDS(resTableComplete, file = paste0(outputDir, "/", "forward_", "step_", stepIdx, "_", "completetable", "_", as.character(bestRow$covar), as.character(bestRow$var), ".RData"))

      stepIdx <- stepIdx + 1
    }
    else {
      # objf function value did not improve
      cli::cli_h1("objf value did not improve, exiting the search ...")

      resTableComplete <- rbind(resTableComplete, resTable)
      saveRDS(resTableComplete, file = paste0(outputDir, "/", "forward_", "step_", stepIdx, "_", "completetable", "_", as.character(bestRow$covar), as.character(bestRow$var), ".RData"))

      break
    }
  }

  cli::cli_h2(cli::col_red("forward search complete"))

  list(fit, resTableComplete)
}

#' Backward covariate search
#'
#' @param covInfo a list containing information about each variable-covariate pair
#' @param fitorig the original 'fit' object before forward search
#' @param fitupdated the updatef 'fit' object, if any, after the forward search
#' @param pVal p-value that should be used for selecting covariates in the forward search
#' @param reFitCovars if the covariates should be added before performing backward search - useful for directly performing backward search without forward search; default is FALSE
#' @param outputDir the name of the output directory that stores the covariate search result
#' @param restart a boolean that controls if the search should be restarted; default is FALSE
#'
#' @return returns the updated 'fit' object at the end of the backward search and a table of information for all the covariates tested
#' @export
#'
#' @author Vipul Mann, Matthew Fidler
backwardSearch <- function(covInfo, fitorig, fitupdated, pVal = 0.01, reFitCovars = FALSE, outputDir, restart = FALSE) {
  if (missing(outputDir)) {
    stop("please specify output directory to store the results for backward search. aborting ...")
  }
  cli::cli_h1("starting backward search...")
  resTableComplete <- data.frame(matrix(ncol = 14, nrow = 0))

  stepIdx <- 1

  if (!missing(fitupdated)) {
    if (all(names(fitupdated$ini$theta) %in% names(fitorig$ini$theta))) {
      cli::cli_alert_warning("no covariates added in the forward search, skipping backward search")
      return(list(fitorig, NULL))
    }
    else {
      fit <- fitupdated
    }
  }

  if (reFitCovars) {
    covSearchRes <- addCovMultiple(covInfo, fitorig, indep = FALSE)
    fitupdated <- covSearchRes[[length(covSearchRes)]][[1]] # get the last fit object with all covariates added # DOES NOT ADD $ini
    fit <- fitupdated
  }

  fnameTablePatternBackward <-
    paste0("backward_step_", "[0-9]+", "_table_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )
  fnameFitPatternBackward <-
    paste0("backward_step_", "[0-9]+", "_fit_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )

  fnameCompleteTablePatternBackward <-
    paste0("forward_step_", "[0-9]+", "_completetable_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )


  fileExistsTab <-
    list.files(paste0("./", outputDir), pattern = fnameTablePatternBackward)

  fileExistsFit <-
    list.files(paste0("./", outputDir), pattern = fnameFitPatternBackward)

  fileExistsCompleteTable <-
    list.files(paste0("./", outputDir), pattern = fnameCompleteTablePatternBackward)

  if (length(fileExistsTab) == 0) {
    restart <- TRUE
  }

  if (!restart) {
    resumeTable <- lapply(fileExistsTab, function(x) {
      readRDS(paste0(outputDir, "/", x))
    })

    resumeTable <- data.table::rbindlist(resumeTable)
    fit <- readRDS(paste0(outputDir, "/", fileExistsFit[[length(fileExistsFit)]]))
    resTableComplete <- readRDS(paste0(outputDir, "/", fileExistsCompleteTable[[length(fileExistsCompleteTable)]]))

    # update covInfo and step idx
    testedCovarVars <- paste0(unlist(resumeTable$covar), unlist(resumeTable$var))

    for (x in testedCovarVars) {
      covInfo[[x]] <- NULL
    }

    stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1

    cli::cli_alert_success("loaded backward search data from disk ...")
    cli::cli_alert_success("resuming backward search ...")
  }


  cli::cli_h2(cli::col_blue("initial function text to remove from:"))
  cli::cli_text(cli::col_red("{fit$fun.txt}"))

  # check if covInfo vars in fit; abort if nonoe of the covaraites added in the forward search step



  # Now remove covars step by step until the objf fun value...?
  while (length(covInfo) > 0) {
    # Remove covars on by one: if objf val increases retain covar; otherwise (objf val decreases), remove the covar
    # At any stage, retain the one that results in highest increase in objf value; exit if removal of none results in increase

    covSearchRes <- removeCovMultiple(covInfo, fit)

    resTable <- lapply(covSearchRes, function(res) {
      x <- res[[1]]
      nam_covar <- res[[2]][[1]]
      nam_var <- res[[2]][[2]]
      covNames <- res[[3]]

      # fwd: if deltObjf <0: pchisq=1-pchisq(-deltObjf, dof), else pchisq=1
      # bck: if deltObjf >0: pchisq=1-pchisq(deltObjf, dof), else pchisq=1

      dObjf <- x$objf - fit$objf
      dof <- length(fit$uif$ini$est) - length(x$uif$ini$est)

      if (dObjf > 0) {
        pchisqr <- 1 - pchisq(dObjf, df = dof)
      }
      else {
        pchisqr <- 1
      }

      l1 <- list(step = stepIdx, covar = nam_covar, var = nam_var, objf = x$objf, deltObjf = dObjf, AIC = x$AIC, BIC = x$BIC, numParams = length(x$uif$ini$est), qchisqr = qchisq(1 - pVal, dof), pchisqr = pchisqr, included = "", searchType = "backward")
      l2 <- list(covNames = covNames, covarEffect = fit$parFixedDf[covNames, "Estimate"])

      c(l1, l2)
    })

    resTable <- data.frame(do.call(rbind, resTable))
    colnames(resTable) <- c(names(resTable))

    colnames(resTableComplete) <- colnames(resTable)

    bestRow <- resTable[which.min(resTable$pchisqr), ]


    if (bestRow$pchisqr <= pVal) {
      # objf function value increased after removal of covariate: retain the best covariate at this stage, test for the rest

      resTable[which.min(resTable$pchisqr), "included"] <- "yes"
      bestRow[, "included"] <- "yes"


      cli::cli_h1("best model at step {stepIdx}: ")
      print(bestRow)

      fit <-
        covSearchRes[[which.min(resTable$pchisqr)]][[1]] # extract fit object corresponding to the best model
      covInfo[[paste0(as.character(bestRow$covar), as.character(bestRow$var))]] <- NULL

      saveRDS(fit, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "fit", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      saveRDS(bestRow, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "table", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      cli::cli_h2("retaining {paste0(as.character(bestRow$covar), as.character(bestRow$var))}")
      resTableComplete <- rbind(resTableComplete, resTable)
      saveRDS(resTableComplete, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "completetable", "_", as.character(bestRow$covar), as.character(bestRow$var), ".RData"))

      stepIdx <- stepIdx + 1
    }
    else {
      # objf function value did not improve
      cli::cli_h1("objf value did not improve, exiting the search ...")
      resTableComplete <- rbind(resTableComplete, resTable)
      saveRDS(resTableComplete, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "completetable", "_", as.character(bestRow$covar), as.character(bestRow$var), ".RData"))

      break
    }
  }

  cli::cli_h2(cli::col_red("backward search complete"))

  list(fit, resTableComplete)
}


# fitDapto <- readRDS("daptomycin.Rds")
# covarSearchAuto(fitDapto, c("v1"), c("SEX"), catCovariates = "SEX", restart = TRUE, pVal = list(fwd = 1, bck = 1))

# covarSearchAuto(fitDapto, c('v1', 'v2'), c('WT', 'SEX'), catCovariates = 'SEX', restart = FALSE, pVal=list(fwd=1, bck=1))


# fun.txt: "    cl = exp(tcl+eta.cl)\n    q = exp(tq+eta.q)\n    v1 = exp(tv1+eta.v1)\n    v2=exp(tv2+eta.v2)\n    cp = linCmt()\n    cp ~ add(add.err)"
# data colnames: "id"    "time"  "CL"    "V1"    "Q"     "V2"    "A1"    "A2"    "Cp"    "centr" "peri"  "CRCL"  "WT"    "SEX"   "AGE"   "DV"
