#' Format confidence bounds for a variable into bracketed notation using string formatting
#'
#' @param var a list of values for the varaible
#' @param confLower the lower bounds for each of the values
#' @param confUpper the upper bounds for each of the values
#' @param sigdig the number of significant digits
#'
#' @noRd
addConfboundsToVar <-
  function(var, confLower, confUpper, sigdig = 3) {
    res <- lapply(seq_along(var), function(idx) {
      paste0(
        signif(var[idx], sigdig),
        " (",
        signif(confLower[idx], sigdig),
        ", ",
        signif(confUpper[idx], sigdig),
        ")"
      )
    })

    unlist(res)
  }

#' Bootstrap nlmixr fit
#'
#' Bootstrap input dataset and rerun the model to get confidence bounds and aggregated parameters
#'
#' @param fit the nlmixr fit object
#' @param nboot an integer giving the number of bootstrapped models to be fit; default value is 100
#' @param nSampIndiv an integer specifying the number of samples in each bootstrapped sample; default is the number of unique subjects in the original dataset
#' @param pvalues a vector of pvalues indicating the probability of each subject to get selected; default value is NULL implying that probability of each subject is the same
#' @param resume a boolean that indicates if a previous session has to be resumed; default value is TRUE
#'
#'
#' @author Vipul Mann, Matthew Fidler
#' @export
#'
#' @examples
#' bootstrapFit(fit)
#' bootstrapFit(fit, nboot = 5, resume = FALSE) # overwrites any of the existing data or model files
#' bootstrapFit(fit, nboot = 7) # resumes fitting using the stored data and model files
bootstrapFit <- function(fit,
                         nboot = 100,
                         nSampIndiv,
                         pvalues = NULL,
                         resume = TRUE) {
  modelsList <-
    modelBootstrap(fit, nboot, nSampIndiv, pvalues, resume) # multiple models
  bootSummary <-
    getBootstrapSummary(modelsList) # aggregate values/summary

  # modify the fit object
  nrws <- nrow(bootSummary$parFixedDf$mean)
  sigdig <- fit$control$sigdig

  newParFixedDf <- fit$parFixedDf
  newParFixed <- fit$parFixed

  # Add Estimate_boot
  est <- unname(bootSummary$parFixedDf$mean[1:nrws, 1])
  cLower <- unname(bootSummary$parFixedDf$confLower[1:nrws, 1])
  cUpper <- unname(bootSummary$parFixedDf$confUpper[1:nrws, 1])
  estEst <- est

  estimateBoot <- addConfboundsToVar(est, cLower, cUpper, sigdig)

  # Add SE_boot
  seBoot <- unname(bootSummary$parFixedDf$stdDev[1:nrws, 1])

  # Add Back-transformed
  est <- unname(bootSummary$parFixedDf$mean[1:nrws, 4])
  cLowerBT <- unname(bootSummary$parFixedDf$confLower[1:nrws, 4])
  cUpperBT <- unname(bootSummary$parFixedDf$confUpper[1:nrws, 4])
  backTransformed <-
    addConfboundsToVar(est, cLowerBT, cUpperBT, sigdig)
  estBT <- est


  newParFixedDf["Bootstrap Estimate"] <- estEst
  newParFixedDf["Bootstrap SE"] <- seBoot
  newParFixedDf["Bootstrap %RSE"] <- seBoot / estEst * 100
  newParFixedDf["Bootstrap CI Lower"] <- cLowerBT
  newParFixedDf["Bootstrap CI Upper"] <- cUpperBT
  newParFixedDf["Bootstrap Back-transformed"] <- estBT

  newParFixed["Bootstrap Estimate"] <- estimateBoot
  newParFixed["Bootstrap SE"] <- signif(seBoot, sigdig)
  newParFixed["Bootstrap %RSE"] <-
    signif(seBoot / estEst * 100, sigdig)
  newParFixed["Bootstrap Back-transformed(95%CI)"] <-
    backTransformed

  assign("parFixedDf", newParFixedDf, envir = fit$env)
  assign("parFixed", newParFixed, envir = fit$env)

  assign("omegaSummary", bootSummary$omega, envir = fit$env)
}


#' Perform bootstrap-sampling from a given dataframe
#'
#' @param data the original dataframe object to sample from for bootstrapping
#' @param nsamp an integer specifying the number of samples in each bootstrapped sample; default is the number of unique subjects in the original dataset
#' @param uid_colname a string representing the unique ID of each subject in the data; default values is 'ID'
#' @param pvalues a vector of pvalues indicating the probability of each subject to get selected; default value is NULL implying that probability of each subject is the same
#'
#' @return returns a bootstrap sampled dataframe object
#' @author Vipul Mann, Matthew Fidler
#'
#' @examples
#' sampling(data)
#' sampling(data, 10)
#' @noRd
sampling <- function(data,
                     nsamp,
                     uid_colname,
                     pvalues = NULL) {
  checkmate::assert_data_frame(data)
  if (missing(nsamp)) {
    nsamp <- length(unique(data[, uid_colname]))
  }
  else {
    checkmate::assert_integerish(nsamp,
      len = 1,
      any.missing = FALSE,
      lower = 2
    )
  }

  checkmate::assert_integerish(nsamp,
    lower = 2,
    len = 1,
    any.missing = FALSE
  )

  if (missing(uid_colname)) {
    # search the dataframe for a column name of 'ID'
    colNames <- colnames(data)
    colNamesLower <- tolower(colNames)
    if ("id" %in% colNames) {
      uid_colname <- colNames[which("id" %in% colNamesLower)]
    }
    else {
      uid_colname <- "ID"
    }
  }
  else {
    checkmate::assert_character(uid_colname)
  }

  uids <- unique(data[, uid_colname])
  uids_samp <- sample(uids,
    size = nsamp,
    replace = TRUE,
    prob = pvalues
  )

  sampled_df <-
    data.frame(samp_dat)[0, ] # initialize an empty dataframe with the same col names

  # populate dataframe based on sampled uids
  # new_id = 1
  .env <- environment()
  .env$new_id <- 1

  do.call(rbind, lapply(uids_samp, function(u) {
    data_slice <- data[data[, uid_colname] == u, ]
    start <- NROW(sampled_df) + 1
    end <- start + NROW(data_slice) - 1

    data_slice[uid_colname] <-
      .env$new_id # assign a new ID to the sliced dataframe
    .env$new_id <- .env$new_id + 1
    # assign("new_id", .env$new_id+1, .env)
    data_slice
  }))
}


#' Fitting multiple bootstrapped models without aggregaion; called by the function bootstrapFit()
#'
#' @param fit the nlmixr fit object
#' @param nboot an integer giving the number of bootstrapped models to be fit; default value is 100
#' @param nSampIndiv an integer specifying the number of samples in each bootstrapped sample; default is the number of unique subjects in the original dataset
#' @param pvalues a vector of pvalues indicating the probability of each subject to get selected; default value is NULL implying that probability of each subject is the same
#' @param resume a boolean that indicates if a previous session has to be resumed; default value is TRUE
#'
#' @return a list of lists containing the different attributed of the fit object for each of the bootstrapped models
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' modelBootstrap(fit)
#' modelBootstrap(fit, 5)
#' modelBootstrap(fit, 5, 20)
#' @noRd
modelBootstrap <- function(fit,
                           nboot = 100,
                           nSampIndiv,
                           pvalues = NULL,
                           resume = FALSE) {
  fitName <- as.character(substitute(fit))

  if (!inherits(fit, "nlmixrFitCore")) {
    stop("'fit' needs to be a nlmixr fit", call. = FALSE)
  }

  checkmate::assert_integerish(nboot,
    len = 1,
    any.missing = FALSE,
    lower = 1
  )

  data <- getData(fit)

  if (missing(nSampIndiv)) {
    nSampIndiv <- length(unique(data[, uidCol]))
  }
  else {
    checkmate::assert_integerish(
      nSampIndiv,
      len = 1,
      any.missing = FALSE,
      lower = 2
    )
  }

  # infer the ID column from data
  colNames <- colnames(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNames) {
    uid_colname <- colNames[which("id" %in% colNamesLower)]
  }
  else {
    stop("cannot find the 'ID' column! aborting ...", call. = FALSE)
  }


  # if (missing(uidCol)) {
  #   # search the dataframe for a column name of 'ID'
  #   colNames <- colnames(data)
  #   colNamesLower <- tolower(colNames)
  #   if ("id" %in% colNames) {
  #     uid_colname <- colNames[which("id" %in% colNamesLower)]
  #   }
  #   else {
  #     uid_colname <- "ID"
  #   }
  # }
  # else {
  #   checkmate::assert_character(uid_colname)
  # }


  uif <- fit$uif
  fitMeth <- getFitMethod(fit)

  bootData <- vector(mode = "list", length = nboot)

  output_dir <-
    paste0("nlmixrBootstrapCache_", fitName) # a new directory with this name will be created
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  fnameBootData <- paste0(output_dir, "/",
    as.character(substitute(boot_data)),
    ".RData",
    sep = ""
  )

  if (resume) {
    if (file.exists(fnameBootData)) {
      cli::cli_alert_success("resuming bootstrap data sampling using data at {fnameBootData}")
      bootData <- readRDS(fnameBootData)
      # assign('.Random.seed', attr(bootData, 'randomSeed'), envir = .GlobalEnv)

      startCtr <- sum(!sapply(bootData, is.null)) + 1

      if (startCtr > nboot) {
        cli::cli_alert_danger(
          cli::col_red(
            "the data file already has {startCtr-1} datasets when max datasets is {nboot}"
          )
        )
        stop(
          "aborting ... the number of models in resume file exceeds the number of models to be estimated",
          call. = FALSE
        )
      }
    }
    else {
      cli::cli_alert_danger(cli::col_red(
        "need the file at {paste0(getwd(), '/', fnameBootData)} to resume"
      ))
      stop("aborting...resume file missing", call. = FALSE)
    }
  }

  else {
    startCtr <- 1
  }

  for (mod_idx in startCtr:nboot) {
    bootData[[mod_idx]] <- sampling(data,
      nsamp = nSampIndiv,
      uid_colname = uidCol,
      pvalues = pvalues
    )


    # save bootData in curr directory: read the file using readRDS()
    attr(bootData, "randomSeed") <- .Random.seed
    saveRDS(bootData, file = fnameBootData)
  }

  bootData <- readRDS(fnameBootData)

  # check if number of samples in stored file is the same as required number of samples
  if (length(bootData) == nboot) {
    cli::cli_alert_success(
      cli::col_silver(
        "sampling complete! saved data is at {paste0(getwd(), '/', fnameBootData)}"
      )
    )
  }
  else {
    cli::cli_alert_danger(
      cli::col_red(
        "could not save all data. resume bootstrapping using saved data at {paste0(getwd(), '/', fnameBootData)}"
      )
    )
    stop(
      "aborting... could not save all the data; resume using data saved at {paste0(getwd(), '/', fnameBootData)}",
      call. = FALSE
    )
  }


  # Fitting models to bootData now
  .env <- environment()
  fnameModelsEnsemble <- paste0(output_dir, "/",
    as.character(substitute(modelsEnsemble)),
    ".RData",
    sep = ""
  )

  if (resume) {
    if (file.exists(fnameModelsEnsemble) &&
      (file.exists(fnameBootData))) {
      cli::cli_alert_success(
        "resuming bootstrap model fitting using data at {fnameModelsEnsemble} and {fnameBootData}"
      )
      bootData <- readRDS(fnameBootData)
      modelsEnsembleLoaded <- readRDS(fnameModelsEnsemble)

      .env$mod_idx <- length(modelsEnsembleLoaded) + 1

      if (.env$mod_idx > nboot) {
        cli::cli_alert_danger(
          cli::col_red(
            "the model file already has {.env$mod_idx-1} models when max models is {nboot}"
          )
        )
        stop(
          "aborting...the number of models in resume file exceeds the number of models to be estimated",
          call. = FALSE
        )
      }
    }
    else {
      cli::cli_alert_danger(
        cli::col_red(
          "need both the files: {paste0(getwd(), '/', fnameModelsEnsemble)} and {paste0(getwd(), '/', fnameBootData)} to resume"
        )
      )
      stop(
        "aborting...data and model files missing: {paste0(getwd(), '/', fnameModelsEnsemble)} and {paste0(getwd(), '/', fnameBootData)}",
        call. = FALSE
      )
    }
  }

  else {
    .env$mod_idx <- 1
  }


  # get control settings for the 'fit' object and save computation effort by not computing the tables
  .ctl <- fit$origControl
  .ctl$print <- 0
  .ctl$covMethod <- ""
  .ctl$calcTables <- FALSE

  modelsEnsemble <-
    lapply(bootData[.env$mod_idx:nboot], function(boot_data) {
      cli::cli_h1("Running nlmixr for model index: {.env$mod_idx}")
      assign("mod_idx", .env$mod_idx + 1, .env)

      fit <- suppressWarnings(nlmixr(
        uif,
        boot_data,
        est = fitMeth,
        control = .ctl
      ))

      .env$multipleFits <- list(
        objf = fit$OBJF,
        aic = fit$AIC,
        omega = fit$omega,
        parFixedDf = fit$parFixedDf,
        method = fit$method,
        message = fit$message,
        warnings = fit$warnings
      )
    })

  if (resume) {
    modelsEnsemble <- c(modelsEnsembleLoaded, modelsEnsemble)
  }

  saveRDS(modelsEnsemble, file = fnameModelsEnsemble)

  modelsEnsemble <- readRDS(fnameModelsEnsemble)

  if (length(modelsEnsemble) == nboot) {
    cli::cli_alert_success(
      cli::col_silver(
        "fitting complete! saved models at {paste0(getwd(), '/', fnameModelsEnsemble)}"
      )
    )
  }
  else {
    cli::cli_alert_danger(
      cli::col_red(
        "all models not saved. resume bootstrapping using saved models at {paste0(getwd(),'/', fnameModelsEnsemble)}"
      )
    )
    stop(
      "aborting...could not save all the models; resume bootstrapping useing models saved at {paste0(getwd(),'/', fnameModelsEnsemble)}",
      call. = FALSE
    )
  }

  modelsEnsemble
}

#' Get the nlmixr method used for fitting the model
#'
#' @param fit the nlmixr fit object
#'
#' @return returns a string representing the method used by nlmixr for fitting the given model
#'
#' @author Vipul Mann, Matthew Fidler
#'
#' @examples
#' getFitMethod(fit)
#' @noRd
getFitMethod <- function(fit) {
  methodsList <-
    c(
      "nlmixrFOCEi",
      "nlmixrNlmeUI",
      "nlmixrSaem",
      "nlmixrFOCE",
      "nlmixrFOi",
      "nlmixrFO",
      "nlmixrPosthoc"
    )
  methodsListMap <-
    c("focei", "nlme", "saem", "foce", "foi", "fo", "posthoc")

  if (!(inherits(fit, "nlmixrFitCore"))) {
    stop("'fit' needs to be a nlmixr fit", call. = FALSE)
  }

  res <- lapply(methodsList, function(met) {
    inherits(fit, met)
  })

  methodsListMap[which(res == TRUE)]
}

#' Extract all the relevant variables from a set of bootstrapped models
#'
#' @param fitlist a list of lists containing information on the multiple bootstrapped models; similar to the output of modelsBootstrap() function
#' @param id a character representing the variable of interest: OBJF, AIC, omega, parFixedDf, method, message, warnings
#'
#' @return returns a vector or list across of the variable of interest from all the fits/bootstrapped models
#'
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' extractVars(fitlist, 1) # returns a vector of OBJF values
#' extractVars(fitlist, 4) # returns a list of dataframes containing parFixedDf values
#' @noRd
extractVars <- function(fitlist, id = "objf") {
  if (id == "method") {
    # no lapply for 'method'
    unlist(unname(fitlist[[1]][id]))
  }
  else {
    # if id not equal to 'method'
    res <- lapply(fitlist, function(x) {
      x[[id]]
    })


    if (!(id == "omega" ||
      id == "parFixedDf")) {
      # check if all message strings are empty
      if (id == "message") {
        prev <- TRUE
        for (i in length(res)) {
          status <- (res[[i]] == "") && prev
          prev <- status
        }
        if (status == TRUE) {
          c("")
        }
        else {
          # if non-empty 'message'
          unlist(res)
        }
      }

      else {
        # if id does not equal 'message'
        unlist(res)
      }
    }
    else {
      # if id equals 'omega' or 'parFixedDf
      res
    }
  }
}

#' Summarize the bootstrapped fits/models
#'
#' @param fitList a list of lists containing information on the multiple bootstrapped models; similar to the output of modelsBootstrap() function
#' @return returns aggregated quantities (mean, median, standard deviation, and variance) as a list for all the quantities
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' getBootstrapSummary(fitlist)
#' @noRd
getBootstrapSummary <- function(fitList, ci = 0.95) {
  if (!(ci < 1 && ci > 0)) {
    stop("'ci' needs to be between 0 and 1", call. = FALSE)
  }
  quantLevels <-
    c(0.5, (1 - ci) / 2, 1 - (1 - ci) / 2) # median, (1-ci)/2, 1-(1-ci)/2

  varIds <-
    names(fitList[[1]]) # number of different variables present in fitlist
  summaryList <- lapply(varIds, function(id) {
    if (!(id %in% c("omega", "parFixedDf", "method", "message", "warnings"))) {
      varVec <- extractVars(fitList, id)
      mn <- mean(varVec)
      median <- median(varVec)
      sd <- sd(varVec)

      c(
        mean = mn,
        median = median,
        stdDev = sd
      )
    }
    else if (id == "omega") {
      # omega estimates
      varVec <- simplify2array(extractVars(fitList, id))
      mn <- apply(varVec, 1:2, mean)
      sd <- apply(varVec, 1:2, sd)

      quants <- apply(varVec, 1:2, function(x) {
        unname(quantile(x, quantLevels))
      })

      median <- quants[1, , ]
      confLower <- quants[2, , ]
      confUpper <- quants[3, , ]

      lst <- list(
        mean = mn,
        median = median,
        stdDev = sd,
        confLower = confLower,
        confUpper = confUpper
      )
    }

    else if (id == "parFixedDf") {
      # parameter estimates (dataframe)
      varVec <- extractVars(fitList, id)
      mn <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, mean, na.rm = TRUE)
      sd <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, sd, na.rm = TRUE)

      quants <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, function(x) {
          unname(quantile(x, quantLevels, na.rm = TRUE))
        })

      median <- quants[1, , ]
      confLower <- quants[2, , ]
      confUpper <- quants[3, , ]

      lst <- list(
        mean = mn,
        median = median,
        stdDev = sd,
        confLower = confLower,
        confUpper = confUpper
      )
    }

    else {
      # if id equals method, message, or warning
      extractVars(fitList, id)
    }
  })

  names(summaryList) <- varIds
  summaryList

  class(summaryList) <- "nlmixrBoostrapSummary"
}

#' Print a well-formatted summary for the bootstrap models
#'
#' @param x the summary object returned by the getBootstrapSummary() function
#'
#' @author Vipul Mann, Matthew Fidler
#'
#' @noRd
print.nlmixrBootstrapSummary <- function(x, fitObj) {
  if (!inherits(fitObj, "nlmixrFitCore")) {
    stop("'fit' needs to be a nlmixr fit", call. = FALSE)
  }
  sigdig <- fitObj$control$sigdig

  objf <- x$objf
  aic <- x$aic
  method <- x$method
  message <- x$message
  warnings <- x$warnings

  omega <- x$omega
  parFixedDf <- x$parFixedDf

  cli::cli_h1(
    cli::col_red(
      "Summary of the bootstrap models using method: {cli::col_yellow(method)}"
    )
  )
  cli::cli_ol()
  cli::cli_li(cli::col_blue(
    cli::style_bold("Objective function"),
    cli::col_yellow(" (summary$objf)")
  ))
  print(objf)

  cli::cli_li(cli::col_blue(cli::style_bold("AIC"), cli::col_yellow(" (summary$aic)")))
  print(aic)

  cli::cli_li(cli::col_magenta(
    cli::style_bold(
      "Omega matrices: mean, median, standard deviation, and confidence bousnds"
    ),
    cli::col_yellow(" (summary$omega)")
  ))

  lapply(seq_along(omega), function(x) {
    cli::cli_text(cli::col_green(paste0("$", names(omega)[x])))
    print(signif(omega[[x]], sigdig))
  })

  cli::cli_li(cli::col_magenta(
    cli::style_bold(
      "Estimated parameters: mean, median, standard deviation, and confidence bounds"
    ),
    cli::col_yellow(" (summary$parFixedDf)")
  ))
  lapply(seq_along(parFixedDf), function(x) {
    cli::cli_text(cli::col_yellow(paste0("$", names(parFixedDf)[x])))
    print(signif(parFixedDf[[x]], sigdig))
  })

  cli::cli_li(cli::cli_text(
    cli::bg_green(cli::style_bold("Method")),
    cli::col_yellow(" (summary$method)")
  ))
  print(method)

  cli::cli_li(cli::cli_text(
    cli::bg_yellow(cli::style_bold("Messages")),
    cli::col_yellow(" (summary$message)")
  ))
  print(message)

  cli::cli_li(cli::cli_text(
    cli::bg_red(cli::style_bold(cli::col_white("Warnings"))),
    cli::col_yellow(" (summary$warnings)")
  ))
  print(warnings)

  cli::cli_h1("end")
}


#' Assign a set of variables to the nlmixr fit environment
#'
#' @param namedVars a named list of variables that need to be assigned to the given environment
#' @param fitobject the nlmixr fit object that contains its environment information
#' @noRd
#'
assignToEnv <- function(namedVars, fitobject) {
  if (!inherits(fitobject, "nlmixrFitCore")) {
    stop("'fit' needs to be a nlmixr fit", call. = FALSE)
  }

  if (is.null(names(namedVars))) {
    stop("'namedVars needs to be a named list", call. = FALSE)
  }

  if (length(namedVars) != length(names(namedVars))) {
    stop("'namedVars does not have all the elements named", call. = FALSE)
  }
  env <- fitobject$env
  lapply(names(namedVars), function(x) {
    assign(x, namedVars[[x]], envir = env)
  })
}
