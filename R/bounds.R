#' Extract the nlmixr bound information from a function.
#'
#' @param fun Function to extract bound information from.
#' @return a data.frame with bound information.
#' @author Bill Denney and Matthew L. Fidler
#' @family nlmixrBounds
#' @export
nlmixrBounds <- function(fun) {
  # Prepare the data.frame
  df <- nlmixrBounds_df(fun)
  # Check the data.frame, adjust values as required, and set its class
  as.nlmixrBounds(df)
}

# This is defines the columns and classes for the nlmixrBounds data.frame.
# (Done here to make testing easier and ensure consistency across functions.)
nlmixrBoundsTemplate <-
  data.frame(
    ntheta = NA_real_,
    neta1 = NA_real_,
    neta2 = NA_real_,
    name = NA_character_,
    lower = NA_real_,
    est = NA_real_,
    upper = NA_real_,
    fix = NA,
    err = NA_character_,
    label = NA_character_,
    backTransform = "",
    condition = NA_character_,
    stringsAsFactors = FALSE
  )

# Generate the data.frame form of nlmixrBounds from a function or call
nlmixrBounds_df <- function(fun) {
  df <- nlmixrBoundsTemplate[-1, ]
  funPrepared <- nlmixrBoundsPrepareFun(fun)
  funParsed <- nlmixrBoundsParser(funPrepared)
  currentParse <- 0
  if (!("assign" %in% funParsed[[1]]$operation)) {
    stop("first initialization item must be 'theta', 'omega', or 'sigma'", call. = FALSE)
  }
  for (currentParse in seq_along(funParsed)) {
    if ("assign" %in% funParsed[[currentParse]]$operation) {
      if (funParsed[[currentParse]]$operation[2] == "theta") {
        newRows <-
          nlmixrBoundsParserTheta(x = funParsed[[currentParse]], currentData = nlmixrBoundsTemplate)
        newRows$ntheta <-
          if (nrow(df) == 0) {
            1
          } else if (all(is.na(df$ntheta))) {
            1
          } else {
            max(df$ntheta, na.rm = TRUE) + 1
          }
      } else if (funParsed[[currentParse]]$operation[2] == "omega") {
        newRows <-
          nlmixrBoundsParserOmega(x = funParsed[[currentParse]], currentData = nlmixrBoundsTemplate)
        maxPreviousEta <-
          if (nrow(df) == 0) {
            0
          } else if (all(is.na(df$neta1))) {
            0
          } else {
            max(df$neta1, na.rm = TRUE)
          }
        newRows$neta1 <- maxPreviousEta + newRows$neta1
        newRows$neta2 <- maxPreviousEta + newRows$neta2
      } else {
        stop(paste0("report a bug.  unknown assignment for '", funParsed[[currentParse]]$operation[2], "'"), call. = FALSE) # nocov
      }
      currentRows <- nrow(df) + seq_len(nrow(newRows))
      df <- rbind(df, newRows)
    } else if (funParsed[[currentParse]]$operation %in% "attribute") {
      df[currentRows, ] <-
        nlmixrBoundsParserAttribute(
          x = funParsed[[currentParse]],
          currentData = df[currentRows, ]
        )
    } else {
      stop(paste0("report a bug.  unknown nlmixrBounds operation: '", funParsed[[currentParse]]$operation, "'"), call. = FALSE) # nocov
    }
  }
  df
}

nlmixrBoundsPrepareFun <- function(fun) {
  ret <- fun
  if (!is.null(attr(fun, "srcref"))) {
    # Check for comments, and if comments exist, try to convert them to
    # 'label(%s)'.
    hasComments <-
      any(getParseData(
        parse(text = as.character(attr(fun, "srcref")), keep.source = TRUE),
      )$token == "COMMENT")
    if (hasComments) {
      cli::cli_alert_info("parameter labels from comments will be replaced by 'label()'")
      ret <- nlmixrBoundsPrepareFunComments(as.character(attr(fun, "srcref")))
    }
  }
  ret
}

#' Prepare an bounds function with comments by extracting the comments and
#' converting them to label()
#'
#' @param fun_char The function as a vector of character strings
#' @return A function body with comments converted to \code{label()} and pipes
#'   converted to \code{condition()} calls.
#' @noRd
nlmixrBoundsPrepareFunComments <- function(fun_char) {
  # drop comment-only lines
  w <- which(regexpr("^ *#+.*", fun_char) == 1)
  if (length(w) > 0) {
    fun_char <- fun_char[-w]
  }
  # convert comments to 'label()' values
  w <- which(regexpr("^ *[^\n\"]+ *#+.*", fun_char) != -1)
  if (length(w) > 0) {
    labels <- gsub(x = fun_char[w], pattern = "^ *[^\n\"]+ *#+ *(.*) *$", replacement = "\\1")
    labels <- sapply(
      labels,
      function(x) {
        # Ensure that quotes and other special characters are correctly escaped.
        return(sprintf("label(%s)", paste0(deparse(x))))
      }
    )
    # Remove the comment and then insert the labels as new lines.
    fun_char[w] <- gsub(x = fun_char[w], pattern = "#.*$", replacement = "")
    fun_char <- c(fun_char, labels)[order(c(seq_along(fun_char), w))]
  }
  # Perform final parsing of the modified function
  fun_parsed <-
    try(
      eval(parse(text = paste0(fun_char, collapse = "\n"))),
      silent = TRUE
    )
  if (inherits(fun_parsed, "try-error")) {
    stop("error parsing bounds: possible (unsupported) comment/condition inside bounds", call. = FALSE)
  }
  # The current environment is not valid for the function as parsed here.
  environment(fun_parsed) <- emptyenv()
  fun_parsed
}

nlmixrBoundsSuggest <- function(varname, lower, est, upper, fixed) {
  varnameC <- na.omit(varname)
  maskDupVarname <- duplicated(varnameC)
  if (any(maskDupVarname)) {
    stop(
      paste0(
        "duplicated parameter names: '",
        paste(unique(varnameC[maskDupVarname]), collapse = "', '"),
        "'"
      ),
      call. = FALSE
    )
  }
  maskNAEst <- is.na(est)
  maskNALower <- is.na(lower)
  maskNAUpper <- is.na(upper)
  if (any(maskNAEst | maskNALower | maskNAUpper)) {
    stop(
      "NA values\n",
      if (any(maskNAEst)) {
        paste0("  estimates: '", paste(varname[maskNAEst], collapse = "', '"), "'\n")
      },
      if (any(maskNALower)) {
        paste0("  lower bounds: '", paste(varname[maskNALower], collapse = "', '"), "'\n")
      },
      if (any(maskNAUpper)) {
        paste0("  upper bounds: '", paste(varname[maskNAUpper], collapse = "', '"), "'\n")
      },
      call. = FALSE
    )
  }
  maskGood <-
    !is.infinite(est) &
      lower < est &
      est < upper
  maskInfEst <-
    is.infinite(est)
  maskSuggestFixed <-
    !maskInfEst &
      (lower == upper |
        lower == est |
        upper == est)
  maskSuggestReorder <-
    !(maskInfEst | maskSuggestFixed) &
      (lower > est |
        est > upper |
        lower > upper)
  maskUnknown <-
    !(maskGood |
      maskInfEst |
      maskSuggestFixed |
      maskSuggestReorder)
  # Generate the messages
  messageInfEst <- NULL
  messageFixed <- NULL
  messageReorder <- NULL
  messageUnknown <- NULL
  if (any(maskInfEst)) {
    messageInfEst <-
      paste(
        "  infinite estimates:",
        paste(varname[maskInfEst], collapse = ", ")
      )
  }
  if (any(maskSuggestFixed)) {
    messageFixed <-
      paste(
        "  consider fixing these:\n",
        paste(
          sprintf("    %s = fixed(%g)", varname[maskSuggestFixed], est[maskSuggestFixed]),
          collapse = "\n"
        )
      )
  }
  if (any(maskSuggestReorder)) {
    hasLower <- is.finite(lower[maskSuggestReorder])
    hasUpper <- is.finite(upper[maskSuggestReorder])
    bestOrder <-
      # matrix is necessary if there is a single item that needs to be reordered
      # to prevent it from returning as a vector.
      matrix(
        apply(
          cbind(
            lower[maskSuggestReorder],
            est[maskSuggestReorder],
            upper[maskSuggestReorder]
          ),
          MARGIN = 1,
          FUN = sort
        ),
        nrow = sum(maskSuggestReorder)
      )
    messageReorderLower <-
      ifelse(
        hasLower | hasUpper,
        paste0(format(bestOrder[, 1]), ", "),
        ""
      )
    messageReorderUpper <-
      ifelse(
        hasUpper,
        paste0(", ", format(bestOrder[, 3])),
        ""
      )
    messageReorderConcatStart <-
      ifelse(
        fixed[maskSuggestReorder],
        "fixed(",
        ifelse(
          hasLower | hasUpper,
          "c(",
          ""
        )
      )
    messageReorderConcatEnd <-
      ifelse(
        messageReorderConcatStart == "",
        "",
        ")"
      )
    messageReorder <-
      paste(
        "  reorder bounds:\n",
        paste0(
          "    ",
          varname[maskSuggestReorder],
          " = ",
          messageReorderConcatStart,
          messageReorderLower,
          format(bestOrder[, 2]),
          messageReorderUpper,
          messageReorderConcatEnd,
          collapse = "\n"
        )
      )
  }
  if (any(maskUnknown)) {
    messageUnknown <-
      paste0(paste(
        "  unknown issue with parameters: '",
        paste(varname[maskUnknown], collapse = "', '")
      ), "'")
  }
  if (any(!is.null(messageFixed), !is.null(messageReorder), !is.null(messageUnknown))) {
    stop(
      "initial conditions error:\n",
      messageInfEst, messageFixed, messageReorder, messageUnknown,
      call. = FALSE
    )
  }
  NULL
}

#' Verify the accuracy of a nlmixrBounds object and update initial conditions,
#' as required.
#'
#' @param df The data.frame to check and convert to an nlmixrBounds object.
#' @param addMissingCols Should missing columns be added to the object?  (Should
#'   typically be FALSE except for testing.)
#' @return An nlmixrBounds object with data confirmed to be consistent.
#'
#' @noRd
as.nlmixrBounds <- function(df, addMissingCols = FALSE) {
  # Ensure that the format is data.frame (instead of data.table, tibble, etc.)
  df <- as.data.frame(df)
  if (nrow(df) == 0) {
    stop("no parameter information", call. = FALSE)
  }
  extraColumns <- setdiff(names(df), names(nlmixrBoundsTemplate))
  if (length(extraColumns)) {
    stop(paste0("extra columns found: '", paste(extraColumns, collapse = "', '"), "'"), call. = FALSE)
  }
  missingColumns <- setdiff(names(nlmixrBoundsTemplate), names(df))
  if (length(missingColumns)) {
    if (!addMissingCols) {
      stop(paste0("columns missing: '", paste(missingColumns, collapse = "', '"), "'"), call. = FALSE)
    } else {
      # Add in the missing columns, if requested.  This is mostly for ensuring
      # that testing works even when new columns are added.
      for (nm in missingColumns) {
        df[[nm]] <- nlmixrBoundsTemplate[[nm]]
      }
    }
  }
  # Ensure that the columns are in the expected order (mainly for simplification
  # of testing)
  df <- df[, names(nlmixrBoundsTemplate)]
  nlmixrBoundsSuggest(
    varname = df$name, lower = df$lower, est = df$est, upper = df$upper, fixed = df$fix
  )
  .w <- which(!is.finite(df$est))
  if (length(.w) > 0) stop("infinite/NA initial parameters: '",
                           paste(df$name[.w], collapse="', '"),
                           "'", call.=FALSE)
  w <- which(df$lower == 0)
  if (length(w) > 0) df$lower[w] <- sqrt(.Machine$double.eps)
  w <- which(df$upper == 0)
  if (length(w) > 0) df$upper[w] <- -sqrt(.Machine$double.eps)
  class(df) <- c("nlmixrBounds", "data.frame")
  df
}

# nlmixrBoundsParser #####

#' Functions to assist with setting initial conditions and boundaries
#'
#' These functions are not intended to be called by a user.  They are intended
#' to be internal to nlmixr
#'
#' @param x the object to attempt extraction from
#' @return A list with how the object will be used
#' @family nlmixrBounds
#' @export
nlmixrBoundsParser <- function(x) {
  UseMethod("nlmixrBoundsParser")
}
# When nested assignments occur (like '({({a = 1})})'), unnest assignments so that
# the result is a flat list.
nlmixrBoundsParserUnnest <- function(x) {
  if ("operation" %in% names(x)) {
    list(x)
  } else if (length(x) == 1) {
    nlmixrBoundsParserUnnest(x[[1]])
  } else {
    ret <- list()
    for (idx in seq_along(x)) {
      ret <- append(ret, nlmixrBoundsParserUnnest(x[idx]))
    }
    ret
  }
}
#' @export
nlmixrBoundsParser.default <- function(x) {
  stop(
    "cannot parse initial condition: '", deparse(x), "', class: ", class(x),
    call. = FALSE
  )
}
# @describeIn nlmixrBoundsParser For functions, apply to the function body
#' @export
nlmixrBoundsParser.function <- function(x) {
  nlmixrBoundsParser(body(x))
}
# @describeIn nlmixrBoundsParser For function bodies and similar.
#' @export
`nlmixrBoundsParser.{` <- function(x) {
  # Recurse; there is nothing more to do
  nlmixrBoundsParserUnnest(
    lapply(x[-1], nlmixrBoundsParser)
  )
}
#' @describeIn nlmixrBoundsParser For function bodies and similar.
#' @export
`nlmixrBoundsParser.(` <- function(x) {
  `nlmixrBoundsParser.{`(x)
}
# @describeIn nlmixrBoundsParser Assignments to thetas with names
#' @export
`nlmixrBoundsParser.<-` <- function(x) {
  list(
    operation = c("assign", "theta"),
    varname = as.character(x[[2]]),
    value = x[[3]]
  )
}
# @describeIn nlmixrBoundsParser Assignments to thetas with names
#' @export
`nlmixrBoundsParser.=` <- function(x) {
  `nlmixrBoundsParser.<-`(x)
}
# @describeIn nlmixrBoundsParser Assignments to thetas without names,
#   assignment to omegas with or without names, or setting of attributes.
#' @export
nlmixrBoundsParser.call <- function(x) {
  if (as.character(x[[1]]) == "c" |
    grepl(x = as.character(x[[1]]), pattern = "^fix(ed)?$", ignore.case = TRUE)) {
    # unnamed assignment can happen either with 'c()', with 'fix()', or with
    # 'fixed()' (where fix and fixed are case insensitive).
    list(
      operation = c("assign", "theta"),
      varname = NA_character_,
      value = x
    )
  } else if (as.character(x[[1]]) %in% "~") {
    if (length(x) == 2) {
      # one-sided formula
      list(
        operation = c("assign", "omega"),
        varname = NA_character_,
        value = x
      )
    } else if (length(x) == 3) {
      list(
        operation = c("assign", "omega"),
        varname = x[[2]],
        value = x[c(1, 3)]
      )
    } else {
      # This should never be possible because formula only parse with length=2
      # or 3.
      stop("report a bug.  invalid 'omega' assignment: ", deparse(x), call. = FALSE) # nocov
    }
  } else if (as.character(x[[1]]) %in% c("label", "backTransform")) {
    list(
      operation = "attribute",
      varname = as.character(x[[1]]),
      value = x
    )
  } else {
    # This may be a valid R expression that could be evaluated
    stop("invalid call in initial conditions: ", deparse(x), call. = FALSE)
  }
}
# @describeIn nlmixrBoundsParser Assignments of numbers to thetas without names
#' @export
nlmixrBoundsParser.numeric <- function(x) {
  list(
    operation = c("assign", "theta"),
    varname = NA_character_,
    value = x
  )
}
# @describeIn nlmixrBoundsParser Assignments of numbers to thetas without names
#' @export
nlmixrBoundsParser.integer <- function(x) {
  nlmixrBoundsParser.numeric(x)
}

# Convert a parsed theta assignment to data.frame form.
#' Forms of fixed that are allowed are:
#'
#' * \code{fixed(1)}
#' * \code{fixed(1, 2)}
#' * \code{fixed(1, 2, 3)}
#' * \code{c(1, fixed)}
#' * \code{c(1, 2, fixed)}
#' * \code{c(1, 2, 3, fixed)}
#' * \code{c(1, fixed(2))}
#' * \code{c(1, fixed(2), 3)}
#'
#' Where 'fixed' can be 'FIX', 'FIXED', 'fix', or 'fixed'.
nlmixrBoundsParserTheta <- function(x, currentData) {
  currentData$name <- x$varname
  valueFix <- nlmixrBoundsValueFixed(x$value)
  # Set the lower bound, estimate, and upper bound for theta
  value <- valueFix$value
  if (length(value) == 1) {
    currentData$lower <- -Inf
    currentData$est <- value
    currentData$upper <- Inf
  } else if (length(value) == 2) {
    currentData$lower <- value[1]
    currentData$est <- value[2]
    currentData$upper <- Inf
  } else if (length(value) == 3) {
    currentData$lower <- value[1]
    currentData$est <- value[2]
    currentData$upper <- value[3]
  } else {
    stop("Syntax is not supported for thetas: ", deparse(x$value), call. = FALSE)
  }
  if (all(valueFix$fixed)) {
    currentData$fix <- TRUE
  } else if (!any(valueFix$fixed)) {
    currentData$fix <- FALSE
  } else {
    currentData$fix <- valueFix$fixed[2]
    if (length(valueFix$fixed) == 2 && valueFix$fixed[1]) {
      stop(
        "cannot fix theta lower bound without fixed estimate: ",
        deparse(x),
        call. = FALSE
      )
    } else if (length(valueFix$fixed) == 3 &&
      (valueFix$fixed[1] | valueFix$fixed[3])) {
      stop(
        "cannot fix lower or upper bounds without fixing estimate: ",
        deparse(x),
        call. = FALSE
      )
    }
  }
  currentData
}

# Convert a parsed omega assignment to data.frame form.
nlmixrBoundsParserOmega <- function(x, currentData) {
  if (x$value[[2]][[1]] == as.name("|")) {
    # The formula is conditional (i.e. ~a|b)
    conditionValue <- all.vars(x$value[[2]][[3]])
    if (length(conditionValue) != 1) {
      stop(
        "invalid conditional expression, cannot parse: ", deparse(x$value),
        call. = FALSE
      )
    }
    valueCor <- nlmixrBoundsValueCor(x$value[[2]][[2]])
  } else {
    # The condition is not specified, set to "ID"
    conditionValue <- "ID"
    valueCor <- nlmixrBoundsValueCor(x$value[[2]])
  }
  if (length(valueCor$value) == 1 & all(valueCor$cor)) {
    warning("'cor(...)' with a single value is ignored: ", deparse(x$value))
  } else if (any(valueCor$cor) & !all(valueCor$cor)) {
    stop("'cor(...)' must enclose all or none of the block: ", deparse(x$value))
  }
  # Make it a matrix so that decorrelation can occur, if necessary.
  tmpValue <- valueCor$value
  valueMatrix <- matrix(numeric(), nrow = 0, ncol = 0)
  currentDiagIdx <- 0
  while (length(tmpValue)) {
    currentDiagIdx <- currentDiagIdx + 1
    if (length(tmpValue) < currentDiagIdx) {
      stop(
        "incorrect lower triangular matrix dimensions: ", deparse(x$value),
        call. = FALSE
      )
    }
    valueMatrix <-
      rbind(
        cbind(
          valueMatrix,
          tmpValue[seq_len(currentDiagIdx - 1)]
        ),
        tmpValue[seq_len(currentDiagIdx)]
      )
    tmpValue <- tmpValue[-seq_len(currentDiagIdx)]
  }
  if (all(valueCor$cor)) {
    # Confirm that all or none are fixed
    if (!(all(valueCor$fixed) | all(!valueCor$fixed))) {
      stop(
        "either all or none of the elements may be fixed with cor(...): ",
        deparse(x$value)
      )
    }
    # Decorrelate the matrix
    sdValues <- diag(valueMatrix)
    diag(valueMatrix) <- 1
    # Convert to variance-covariance matrix
    valueMatrix <- sweep(sweep(valueMatrix, 1, sdValues, "*"), 2, sdValues, "*")
  }
  currentDiagIdx <- 0
  est <- numeric()
  fix <- logical()
  neta1 <- numeric()
  neta2 <- numeric()
  for (currentDiagIdx in seq_len(nrow(valueMatrix))) {
    nextValues <- seq_len(currentDiagIdx)
    est <- c(est, valueMatrix[currentDiagIdx, nextValues])
    fix <- c(fix, valueCor$fixed[nextValues])
    neta1 <- c(neta1, rep(currentDiagIdx, currentDiagIdx))
    neta2 <- c(neta2, nextValues)
    # Drop the values of 'fixed' that had just been used
    valueCor$fixed <- valueCor$fixed[-nextValues]
  }
  if (class(x$varname) != "character") {
    # It has a name
    nameVec <- all.vars(x$varname)
    if (length(nameVec) != currentDiagIdx) {
      stop(
        "omega assignment left handed side must match lower triangular matrix size",
        call. = FALSE
      )
    }
    name1 <- nameVec[neta1]
    name2 <- nameVec[neta2]
  } else {
    name1 <- name2 <- NA_character_
  }
  # rep() ensures that the default values are set for each row
  finalData <- currentData[rep(1, length(neta1)), ]
  rownames(finalData) <- NULL
  finalData$lower <- -Inf
  finalData$est <- est
  finalData$upper <- Inf
  finalData$fix <- fix
  finalData$neta1 <- neta1
  finalData$neta2 <- neta2
  finalData$name <-
    ifelse(
      is.na(name1) | name1 == name2,
      name1,
      sprintf("(%s,%s)", name1, name2)
    )
  finalData$condition <- conditionValue
  finalData
}

# Add a parsed attribute to the existing theta or omega data.frame.
nlmixrBoundsParserAttribute <- function(x, currentData) {
  if (x$varname == "label") {
    if (!length(x$value) == 2) {
      stop("only apply a single label: ", deparse(x$value), call. = FALSE)
    } else if (!is.character(x$value[[2]])) {
      # We could try to coerce it to a character string, but that will
      # be more likely to yield a bug at some point.
      stop("label must be a character string: ", deparse(x$value), call. = FALSE)
    } else if (!all(is.na(currentData$label))) {
      warning("only last label used: ", deparse(x$value))
    }
    currentData$label <- x$value[[2]]
  } else if (x$varname == "backTransform") {
    if (!all(currentData$backTransform %in% "")) {
      warning("only last backTransform used: ", deparse(x$value))
    }
    if (length(x$value) == 1) {
      # It is `backTransform()` without an argument, set to ""
      currentData$backTransform <- ""
    } else if (length(x$value) == 2) {
      # `backTransform()` with a single argument, the expected way to set a backTransform
      if (is.name(x$value[[2]])) {
        # Make sure that a function name is a call instead of a bare name
        x$value[[2]] <- as.call(list(x$value[[2]]))
      } else if (is.character(x$value[[2]])) {
        # Convert a character string to a call
        x$value[[2]] <- as.call(list(as.name(x$value[[2]])))
      } else if (!is.call(x$value[[2]])) {
        stop(
          "'backTransform()' argument must be a call, function name, or character string giving a function name: ",
          deparse(x$value)
        )
      }
      currentData$backTransform <- deparse(x$value[[2]])
    } else {
      # `backTransform()` with a multiple arguments, cannot handle
      stop("'backTransform()' must have zero or one arguments: ", deparse(x$value))
    }
  } else {
    # We could ignore this, but it is likely to indicate some form of accidental
    # misspecification.
    stop("cannot handle attribute: ", x$varname, call. = FALSE)
  }
  currentData
}

#' Determine the values and what is fixed
#'
#' @param x A call that may have fixed values
#' @return A list with elements of:
#' * value: the numeric value of evaluating the expression
#' * all_fixed: Are all values from the expression fixed ?
#' * fixed: Which value(s) from \code{x} are fixed?
#' @seealso \code{\link{nlmixrBoundsReplaceFixed}},
#'   \code{\link{nlmixrBoundsValueCor}}
#' @author Bill Denney
#' @noRd
nlmixrBoundsValueFixed <- function(x) {
  valueFixed <- nlmixrBoundsReplaceFixed(x, replacementName = NULL)
  # determine the numeric value after removing 'fixed' names and using 'fixed()'
  # like 'c()'
  value <- try(eval(valueFixed$call, list(fixed = c)), silent=TRUE)
  if (inherits(value, "try-error")) {
    stop(
      "error parsing initial condition '", deparse(x), "': ", attr(value, "condition")$message,
      call. = FALSE
    )
  } else if (!is.numeric(value)) {
    stop(
      "non-numeric values in initial condition: ", deparse(x),
      call. = FALSE
    )
  } else if (any(is.nan(value))) {
    stop("NaN values in initial condition: ", deparse(x), call. = FALSE)
  }
  isFixed <- valueFixed$fixed
  if (length(isFixed) != 1) {
    stop( # nocov
      "report as a bug.  'length(isFixed)' > 1 in nlmixrBoundsValueFixed: ", # nocov
      length(isFixed), ", '", deparse(x), "'", # nocov
      call. = FALSE # nocov
    ) # nocov
  }
  if (!isFixed) {
    # Determine which specific values are fixed by assigning them to NaN
    fixedVec <-
      is.nan(
        eval(
          valueFixed$call,
          envir =
            list(
              fixed = function(...) {
                rep(NaN, length(sapply(list(...), FUN = eval)))
              }
            )
        )
      )
  } else {
    fixedVec <- rep(FALSE, length(value))
  }
  list(
    value = value,
    fixed = fixedVec | any(isFixed)
  )
}

#' Find all \code{fixed} names and calls in a call or other language object and
#' detect the fixed status and replace them for later evaluation.
#'
#' @details Note that fixed calls, like \code{fixed(1)} do not make the return
#'   list element of \code{fixed=TRUE}.  When the call form of fixed (e.g.
#'   \code{fixed(1)}) is used, then fixed may only apply to a subset of
#'   elements, and that is accounted for by \code{eval(ret$call)} outside of
#'   this function.
#'
#'   \code{replacementName} can either be numeric (often NaN) which is usable in
#'   a call, or it can be NULL which will remove \code{fixed} from the call
#'   arguments, or it can be something that will be converted to a name.  If an
#'   actual character string replacement is desired, make that into a name,
#'   first.
#'
#' @param x The object to search
#' @param replacementFun The function name (coerced using \code{as.name()}, if
#'   necessary) to use to replace \code{fixed()}
#' @param replacementName \code{NULL}, a number (including \code{NaN}) or the
#'   name (coerced using \code{as.name()}) to replace bare uses of \code{fixed},
#'   such as \code{c(1, fixed)}.
#' @return A list with names of \code{call} which is the modified call version
#'   of \code{x} and \code{fixed} indicating if a \code{replacementName} was
#'   used within.
#' @seealso \code{\link{nlmixrBoundsValueFixed}}
#' @author Bill Denney
#' @noRd
nlmixrBoundsReplaceFixed <- function(x, replacementFun = "fixed", replacementName = NULL) {
  fixedNames <- sapply(c("fix", "FIX", "fixed", "FIXED"), as.name)
  if (!is.name(replacementFun)) {
    if (length(replacementFun) != 1) {
      stop("'replacementFun' must be scalar", call. = FALSE)
    }
    replacementFun <- as.name(replacementFun)
  }
  if (!is.numeric(replacementName) & !is.null(replacementName) & !is.name(replacementName)) {
    if (length(replacementName) != 1) {
      stop("'replacementName' must be scalar or NULL", call. = FALSE)
    }
    replacementName <- as.name(replacementName)
  }
  ret <- x
  fixed <- FALSE
  if (is.call(x)) {
    if (identical(x[[1]], as.name("{")) |
      identical(x[[1]], as.name("("))) {
      # Note that while these types of calls are allowed, the function is only
      # intended to be called on a single initial condition assignment at a time
      # and not to be called on the entire initial condition assignment block.
      retPrep <- lapply(X = x[-1], nlmixrBoundsReplaceFixed)
      x[-1] <- retPrep$call
      ret <- x
      # Potential fragile code: Are there any cases where some but not all of
      # the entries may be fixed?  (Those cases should all be caught by
      # 'fixed()' used as a function which happens elsewhere. When it is used as
      # a name, it should likely always be all fixed.)
      fixed <- any(retPrep$fixed)
    } else if (any(sapply(fixedNames, identical, x[[1]]))) {
      # Replace fixed used as a function, like 'fixed(1)', to the
      # 'replacementFun'
      x[[1]] <- replacementFun
      ret <- x
      # 'fixed' is not set to TRUE here as the setting to TRUE will be based on
      # the evaluation of the function.  (Evaluation happens outside of this
      # function.)
    } else {
      # Find fixed used as a name, like 'c(1, fixed)'.  'fixed' is only valid at
      # the end when used as a name, 'c(1, fixed)' is valid while 'c(fixed, 1)'
      # is invalid.
      fixedPrep <-
        lapply(
          x,
          FUN = nlmixrBoundsReplaceFixed,
          replacementFun = replacementFun,
          replacementName = replacementName
        )
      fixed <- sapply(fixedPrep, "[[", "fixed")
      if (any(fixed)) {
        if (length(fixed) == 1) {
          # This should not be possible because this should either be a call
          # (like 'fixed()') and caught above (x[[1]] %in% fixedNames) or a name
          # (which would not be here because of the outer if block here as
          # is.call(x)).
          stop(
            "report a bug.  Invalid detection of scalar 'fixed' within a call: ", deparse(x),
            call. = FALSE
          )
        }
        if (any(fixed[-length(fixed)])) {
          stop(
            "'fixed' may only be the last item in a list: ", deparse(x),
            call. = FALSE
          )
        }
        # When 'fixed' is at the end of a vector, it applies to the entire
        # vector, and it should be dropped (or modified by replacementName) for
        # later evaluation.
        fixed <- TRUE
      } else {
        fixed <- FALSE
      }
      for (idx in rev(seq_along(fixedPrep))) {
        # Reversing the assignment is needed to ensure that if a $call is NULL,
        # it doesn't mess up the subsequent assignments.
        ret[[idx]] <- fixedPrep[[idx]]$call
      }
    }
  } else if (is.name(x)) {
    fixed <- any(sapply(fixedNames, FUN = identical, x))
    if (fixed) {
      ret <- replacementName
    }
  }
  # No 'else' is required.  Other classes including name, numeric, character,
  # and logical that are likely valid within a call but not fixed.
  list(call = ret, fixed = fixed)
}

#' Detect \code{cor()} in omega blocks
#'
#' @param x A call to check for correlation
#' @return A list with elements for 'value', 'fixed', and 'cor'
#' @seealso \code{\link{nlmixrBoundsValueFixed}}
#' @noRd
nlmixrBoundsValueCor <- function(x) {
  x_nofixed_call <-
    replaceCallName(
      x = x,
      replacementFun = "c",
      sourceNames = c("fix", "fixed", "FIX", "FIXED")
    )
  x_nofixed <-
    replaceNameName(
      x_nofixed_call,
      replacementName = NULL,
      sourceNames = c("fix", "fixed", "FIX", "FIXED")
    )
  # Handle both correlation and
  ret <-
    nlmixrBoundsValueFixed(
      replaceCallName(x, replacementFun = "c", sourceNames = "cor")
    )
  isCor <-
    is.nan(eval(
      x_nofixed,
      envir =
        list(
          cor = function(...) {
            rep(NaN, length(sapply(list(...), FUN = eval)))
          }
        )
    ))
  ret$cor <- isCor
  ret
}

#' Find all calls (i.e. function calls) and replace them
#'
#' This does not apply to names that are not calls.
#'
#' @param x A call to replace calls within
#' @param replacementFun The name (or character string) to use as a replacement
#' @param sourceNames The scalar or vector of names (or character strings) to
#'   replace.
#' @return \code{x} with calls to \code{sourceNames} replaced with
#'   \code{replacementFun}
#' @seealso \code{\link{replaceNameName}}
#' @examples
#' replaceCallName(x = a ~ b(), replacementFun = "c", sourceNames = "b")
#' # names that are not calls are not replaced
#' replaceCallName(x = a ~ b(b), replacementFun = "c", sourceNames = "b")
#' @noRd
replaceCallName <- function(x, replacementFun, sourceNames) {
  if (length(replacementFun) != 1) {
    stop("'replacementFun' must be a scalar")
  } else if (!is.name(replacementFun)) {
    replacementFun <- as.name(replacementFun)
  }
  if (!all(sapply(sourceNames, is.name))) {
    sourceNames <- sapply(X = sourceNames, FUN = as.name)
  }
  if (is.call(x)) {
    ret <- x
    if (any(sapply(X = sourceNames, FUN = identical, y = x[[1]]))) {
      # If one of the source names is the name of the call, replace it.
      ret[[1]] <- replacementFun
    }
    for (idx in rev(seq_len(length(ret) - 1) + 1)) {
      # Recurse through all the parts of the call, if there are more than one.
      ret[[idx]] <-
        replaceCallName(
          ret[[idx]],
          replacementFun = replacementFun,
          sourceNames = sourceNames
        )
    }
  } else {
    # Almost any other class can be part of a call, so just return anything
    # other than a call as-is.
    ret <- x
  }
  ret
}

#' Replace a name that is not used as a function call with a new name.
#'
#' This does not apply to names that are the function name of calls.
#'
#' @inheritParams replaceCallName
#' @param replacementName The name (or character string) to use as a
#'   replacement; if \code{NULL}, the name is removed (which can make an invalid
#'   call).
#' @return \code{x} with calls to \code{sourceNames} replaced with
#'   \code{replacementFun}
#' @seealso \code{\link{replaceCallName}}
#' @noRd
replaceNameName <- function(x, replacementName, sourceNames) {
  if (is.null(replacementName)) {
    # do nothing, NULL is allowed and it removes a name
  } else if (length(replacementName) != 1) {
    stop("'replacementName' must be a scalar")
  } else if (!is.name(replacementName)) {
    replacementName <- as.name(replacementName)
  }
  if (!all(sapply(sourceNames, is.name))) {
    sourceNames <- sapply(X = sourceNames, FUN = as.name)
  }
  if (is.call(x)) {
    ret <- x
    # Recurse through all the parts of the call after the first (i.e. do not
    # replace function calls), if there are more than one.
    for (idx in rev(seq_len(length(ret) - 1) + 1)) {
      ret[[idx]] <-
        replaceNameName(
          ret[[idx]],
          replacementName = replacementName,
          sourceNames = sourceNames
        )
    }
  } else if (is.name(x)) {
    if (any(sapply(X = sourceNames, FUN = identical, y = x))) {
      ret <- replacementName
    } else {
      ret <- x
    }
  } else {
    # Almost any other class can be part of a call, so just return anything
    # other than a call or a name as-is.
    ret <- x
  }
  ret
}

# nlmixrBounds helpers ####

is.nlmixrBounds <- function(x) {
  should <- names(nlmixrBoundsTemplate)
  what <- names(x)
  if (length(should) == length(what)) {
    return(all(what == should))
  } else {
    return(FALSE)
  }
}

##' @export
as.data.frame.nlmixrBounds <- function(x, row.names = NULL, optional = FALSE, ...) {
  cls <- class(x)
  cls <- cls[cls != "nlmixrBounds"]
  tmp <- x
  class(tmp) <- cls
  return(tmp)
}

##' @export
print.nlmixrBounds <- function(x, ...) {
  cat(paste0(crayon::bold("Fixed Effects"), " (", crayon::bold$blue("$theta"), "):"), "\n")
  print(x$theta)
  omega <- x$omega
  if (dim(omega)[1] > 0) {
    cat(paste0("\n", crayon::bold("Omega"), " (", crayon::bold$blue("$omega"), "):"), "\n")
    print(omega)
  }
}

##' @export
`$.nlmixrBounds` <- function(obj, arg, exact = TRUE) {
  m <- as.data.frame(obj, stringsAsFactors = FALSE)
  ret <- m[[arg, exact = exact]]
  if (is.null(ret)) {
    if (arg == "theta") {
      return(nlmixrBoundsTheta(obj, full = FALSE))
    } else if (arg == "theta.full") {
      return(nlmixrBoundsTheta(obj, full = TRUE))
    } else if (arg == "omega") {
      return(nlmixrBoundsOmega(obj))
    } else if (arg == "random") {
      return(nlmixrBoundsOmega(obj, TRUE))
    } else if (arg == "fixed.form") {
      return(nlmixrBoundsTheta(obj, formula = TRUE))
    } else if (arg == "focei.upper") {
      return(nlmixrBounds.focei.upper.lower(obj, "upper"))
    } else if (arg == "focei.lower") {
      return(nlmixrBounds.focei.upper.lower(obj, "lower"))
    } else if (any(arg == c("theta.names", "focei.names"))) {
      return(nlmixrBounds.focei.upper.lower(obj, "name"))
    } else if (arg == "focei.err.type") {
      return(nlmixrBounds.focei.upper.lower(obj, "err"))
    } else if (arg == "eta.names") {
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
nlmixrBounds.eta.names <- function(obj) {
  df <- as.data.frame(obj, stringsAsFactors = FALSE)
  df <- df[!is.na(df$neta1), ]
  ## dft.unfixed <- dft[!dft$fix, ];
  return(paste(df[df$neta1 == df$neta2, "name"]))
}


##' @export
str.nlmixrBounds <- function(object, ...) {
  str(as.data.frame(object), ..., stringsAsFactors = FALSE)
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
nlmixrBounds.focei.upper.lower <- function(obj, type = c("upper", "lower", "name", "err")) {
  type <- match.arg(type)
  df <- as.data.frame(obj, stringsAsFactors = FALSE)
  dft <- df[!is.na(df$ntheta), ]
  ret <- dft[[type]]
  if (is(ret, "factor")) {
    ret <- paste(ret)
  }
  return(ret)
}

nlmixrBoundsTheta <- function(x, full = TRUE, formula = FALSE) {
  if (is.nlmixrBounds(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (formula) full <- FALSE
    w <- which(!is.na(x$ntheta))
    tmp <- (x[w, ])
    nm <- sprintf(ifelse(formula, ".theta.%d", "theta[%d]"), seq_along(w))
    w <- which(!is.na(tmp$name))
    nm[w] <- as.character(tmp$name[w])
    w <- which(!is.na(tmp$err))
    theta <- tmp$est
    if (length(w) > 0) {
      if (full) {
        nm[w] <- sprintf("err[%d]", seq_along(w))
      } else {
        nm <- nm[-w]
        theta <- theta[-w]
      }
    }
    if (formula) {
      if (any(duplicated(nm))) {
        stop("duplicated theta names", call. = FALSE)
      } else {
        return(as.formula(sprintf("%s ~ 1", paste(nm, collapse = " + "))))
      }
    } else {
      names(theta) <- nm
    }
    return(theta)
  } else {
    return(NULL)
  }
}

nlmixrBoundsOmega <- function(x, nlme = FALSE) {
  if (is.nlmixrBounds(x)) {
    w <- which(!is.na(x$neta1))
    if (length(w) > 0) {
      d <- max(x$neta1)
      df <- x[w, ]
      mx <- max(df$neta1)
      mat <- matrix(rep(0, mx * mx), mx, mx)
      diag <- TRUE
      for (i in seq_along(df$neta1)) {
        neta1 <- df$neta1[i]
        neta2 <- df$neta2[i]
        if (neta1 == neta2) {
          mat[neta1, neta2] <- df$est[i]
        } else {
          diag <- FALSE
          mat[neta1, neta2] <- mat[neta2, neta1] <- df$est[i]
        }
      }
      if (is(nlme, "logical") && nlme) {
        df.diag <- df[df$neta1 == df$neta2, ]
        n2 <- sprintf(".eta.%d", seq_along(df.diag$name))
        w <- which(!is.na(df.diag$name))
        n2[w] <- as.character(df.diag$name[w])
        frm <- as.formula(paste(paste(n2, collapse = " + "), "~ 1"))
        if (diag) {
          return(nlme::pdDiag(mat, form = frm))
        } else {
          return(nlme::pdSymm(as.matrix(Matrix::nearPD(mat)$mat), form = frm))
        }
      } else if (is(nlme, "list")) {
        df.diag <- df[df$neta1 == df$neta2, ]
        class(df.diag) <- "data.frame"
        n2 <- sprintf(".eta.%d", seq_along(df.diag$name))
        w <- which(!is.na(df.diag$name))
        n2[w] <- sapply(as.character(df.diag$name[w]), function(x) {
          nlme[[x]]
        })
        frm <- as.formula(paste(paste(n2, collapse = " + "), "~ 1"))
        if (diag) {
          return(nlme::pdDiag(mat, form = frm))
        } else {
          return(nlme::pdSymm(as.matrix(Matrix::nearPD(mat)$mat), form = frm))
        }
      }
      w <- which(df$neta1 == df$neta2)
      dimnames(mat) <- list(paste(df$name[w]), paste(df$name[w]))
      return(mat)
    } else {
      return(matrix(double(), 0, 0))
    }
  } else {
    return(matrix(double(), 0, 0))
  }
}
