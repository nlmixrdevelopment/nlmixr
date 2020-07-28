#' Format confidence bounds for a variable into bracketed notation using string formatting
#'
#' @param var a list of values for the varaible
#' @param confLower the lower bounds for each of the values
#' @param confUpper the upper bounds for each of the values
#' @param sigdig the number of significant digits
#'
#' @author Vipul Mann
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
#' @param restart a boolean that indicates if a previous session has to be restarted; default value is FALSE
#'
#'
#' @author Vipul Mann, Matthew Fidler
#' @export
#'
#' @examples
#' \dontrun{
#' one.cmt <- function() {
#'   ini({
#'     ## You may label each parameter with a comment
#'     tka <- 0.45 # Log Ka
#'     tcl <- 1 # Log Cl
#'     ## This works with interactive models
#'     ## You may also label the preceding line with label("label text")
#'     tv <- 3.45
#'     label("log V")
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
#' fit <- nlmixr(one.cmt, theo_sd, "focei")
#'
#' bootstrapFit(fit)
#' bootstrapFit(fit, nboot = 5, restart = TRUE) # overwrites any of the existing data or model files
#' bootstrapFit(fit, nboot = 7) # resumes fitting using the stored data and model files
# }
bootstrapFit <- function(fit,
                         nboot = 500,
                         nSampIndiv,
                         stratVar,
                         stdErrType = c("perc", "se"),
                         ci = 0.95,
                         pvalues = NULL,
                         restart = FALSE,
                         plotHist = TRUE) {
  stdErrType <- match.arg(stdErrType)
  if (missing(stdErrType)) {
    stdErrType <- "perc"
  }

  if (!(ci < 1 && ci > 0)) {
    stop("'ci' needs to be between 0 and 1", call. = FALSE)
  }

  if (missing(stratVar)) {
    performStrat <- FALSE
  }
  else {
    if (!(stratVar %in% colnames(getData(fit)))) {
      cli::cli_alert_danger("{stratVar} not in data")
      stop("aborting ...stratifying variable not in data", call. = FALSE)
    }
    performStrat <- TRUE
  }

  fitName <- as.character(substitute(fit))

  if (is.null(fit$bootstrapMd5)) {
    bootstrapMd5 <- digest::digest(fit)
    assign("bootstrapMd5", bootstrapMd5, envir = fit$env)
  }

  if (performStrat) {
    resBootstrap <-
      modelBootstrap(
        fit,
        nboot = nboot,
        nSampIndiv = nSampIndiv,
        stratVar = stratVar,
        pvalues = pvalues,
        restart = restart,
        fitName = fitName
      ) # multiple models
    
    modelsList = resBootstrap[[1]]
    fitList = resBootstrap[[2]]
  }
  else {
    resBootstrap <-
      modelBootstrap(
        fit,
        nboot = nboot,
        nSampIndiv = nSampIndiv,
        pvalues = pvalues,
        restart = restart,
        fitName = fitName
      ) # multiple models
    modelsList = resBootstrap[[1]]
    fitList = resBootstrap[[2]]
  }


  bootSummary <-
    getBootstrapSummary(modelsList, ci, stdErrType) # aggregate values/summary

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
  
  # compute bias
  bootstrapBias <-
    bootSummary$objf[[1]] - fit$objf # 1 corresponds to the mean value, 2 corresponds to the median

  # compute covariance matrix
  covMatrix <- cov(getData(fit), getData(fit))
  corMatrix <- cor(getData(fit), getData(fit))

  # assign("deltOBJF", deltOBJF, envir = fit$env)
  assign("bootstrapBias", bootstrapBias, envir = fit$env)
  assign("covMatrix", covMatrix, envir = fit$env)
  assign("corMatrix", corMatrix, envir = fit$env)
  assign("parFixedDf", newParFixedDf, envir = fit$env)
  assign("parFixed", newParFixed, envir = fit$env)
  assign("omegaSummary", bootSummary$omega, envir = fit$env)

  # plot histogram
  if (plotHist) {
    
    # compute delta objf values for each of the models
    origData = getData(fit)
    
    if (is.null(fit$bootstrapMd5)) {
      bootstrapMd5 <- digest::digest(fit)
      assign("bootstrapMd5", bootstrapMd5, envir = fit$env)
    }
    
    # already exists
    output_dir <- paste0("nlmixrBootstrapCache_", as.character(substitute(fit)), "_", fit$bootstrapMd5)
    
    deltOBJFloaded = NULL
    deltOBJF = NULL
    if (!restart){
      deltOBJFloaded = readRDS(paste0("./", output_dir,"/",'deltOBJF',".RData"))
      deltOBJF = c(deltOBJFloaded, deltOBJF)
    }
    else{
      deltOBJF <- lapply(fitList, function(x) {
        xPosthoc = suppressWarnings(nlmixr(x, data=origData, est='posthoc'))
        xPosthoc$objf - fit$objf
      })
      deltOBJF = c(deltOBJFloaded, deltOBJF)
      saveRDS(deltOBJF, file = paste0("./", output_dir,"/",'deltOBJF',".RData"))
    }

    .deltaO <- sort(abs(unlist(deltOBJF)))

    .deltaN <- length(.deltaO)

    .df <- length(fit$ini$est)

    .chisq <- rbind(data.frame(deltaofv=qchisq(seq(0,0.99,0.01),df=.df),
                               quantiles=seq(0,0.99,0.01),
                               Distribution=1L,
                               stringsAsFactors = FALSE),
                    data.frame(deltaofv=.deltaO,
                               quantiles=seq(.deltaN) / .deltaN,
                               Distribution=2L,
                               stringsAsFactors = FALSE))

    .fdelta <- approxfun(seq(.deltaN) / .deltaN, .deltaO)

    .df2 <- round(mean(.deltaO, na.rm=TRUE))

    .dfD <- data.frame(label=paste(c("df\u2248", "df="), c(.df2, .df)),
                       Distribution=c(2L, 1L),
                       quantiles=0.7,
                       deltaofv=c(.fdelta(0.7), qchisq(0.7, df=.df))
                       )

    .dfD$Distribution <- factor(.dfD$Distribution, c(1L, 2L),
                                c("Reference distribution", "\u0394 objective function"))

    .chisq$Distribution <- factor(.chisq$Distribution, c(1L, 2L),
                                  c("Reference distribution", "\u0394 objective function"))

    assign(".chisq", .chisq, globalenv())

    .plot <- ggplot2::ggplot(.chisq, aes(quantiles, deltaofv, color=Distribution)) +
      ggplot2::geom_line() + ggplot2::ylab("\u0394 objective function") +
      ggplot2::geom_text(data=.dfD, aes(label=label), hjust=0) +
      ggplot2::xlab("Distribution quantiles") +
      ggplot2::scale_color_manual(values=c("red", "blue")) +
      RxODE::rxTheme() +
      ggplot2::theme(legend.position="bottom",legend.box="horizontal")

    if (requireNamespace("ggtext", quietly = TRUE)) {
      .plot <- .plot +
        ggplot2::theme(plot.title = ggtext::element_markdown(),
                       legend.position="none") +
        ggplot2::labs(
          title = paste0(
            'Bootstrap <span style="color:blue; opacity: 0.2;">\u0394 objective function (', .deltaN,
            ' models, df\u2248', .df2, ')</span> vs <span style="color:red; opacity: 0.2;">reference \u03C7\u00B2(df=',
            length(fit$ini$est), ")</style>"
          ),
          caption = "\u0394 objective function curve should be on or below the reference distribution curve"
        )
    } else {
      .plot <- ggplot2::labs(
        title = paste0("Distribution of \u0394 objective function values for ", .deltaN, ' df=', .df2, " models"),
        caption = "\u0394 objective function curve should be on or below the reference distribution curve"
      )
    }



    ## df <-
    ##   data.frame(
    ##     vals = c(.deltaO, rchisq(2000, df = length(fit$ini$est))),
    ##     Distribution=factor(c(rep(1, length(.deltaO)), rep(2, 2000)),
    ##                         levels=c(1, 2), labels=c("delta objective function", "reference distribution"))
    ##   )

    ## .plot <- ggplot2::ggplot(df, ggplot2::aes(vals, fill=Distribution)) +
    ##   ggplot2::geom_density( color = NA, alpha = 0.2) +
    ##   ggplot2::xlab("\u0394 Objective function") +
    ##   ggplot2::ylab("Density") +
    ##   ggplot2::scale_fill_manual(values=c("red", "blue"))
    ## .plot <- .plot + RxODE::rxTheme(grid=FALSE)


    ## .plot
    ## +
    ##   ggplot2::scale_color_manual(name='', values = c('delta objective function' = 'blue',
    ##                                                   'reference distribution' = 'red')) +
    ##   RxODE::rxTheme()
    .plot
  }
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
                     pvalues = NULL,
                     performStrat = FALSE,
                     stratVar) {
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

  if (performStrat && missing(stratVar)) {
    print("stratVar is required for stratifying")
    stop("aborting... stratVar not specified", call. = FALSE)
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


  if (performStrat) {
    stratLevels <-
      as.character(unique(data[, stratVar])) # char to access freq. values

    dataSubsets <- lapply(stratLevels, function(x) {
      data[data[, stratVar] == x, ]
    })

    names(dataSubsets) <- stratLevels

    tab <- table(theo_sd[stratVar])
    nTab <- sum(tab)

    sampledDataSubsets <- lapply(names(dataSubsets), function(x) {
      dat <- dataSubsets[[x]]

      uids <- unique(dat[, uid_colname])
      uids_samp <- sample(
        list(uids),
        size = ceiling(nsamp * unname(tab[x]) / nTab),
        replace = TRUE,
        prob = pvalues
      )

      sampled_df <-
        data.frame(dat)[0, ] # initialize an empty dataframe with the same col names

      # populate dataframe based on sampled uids
      # new_id = 1
      .env <- environment()
      .env$new_id <- 1
      do.call(rbind, lapply(uids_samp, function(u) {
        data_slice <- dat[dat[, uid_colname] == u, ]
        start <- NROW(sampled_df) + 1
        end <- start + NROW(data_slice) - 1

        data_slice[uid_colname] <-
          .env$new_id # assign a new ID to the sliced dataframe
        .env$new_id <- .env$new_id + 1
        data_slice
      }))
    })

    do.call("rbind", sampledDataSubsets)
  }

  else {
    uids <- unique(data[, uid_colname])
    uids_samp <- sample(uids,
      size = nsamp,
      replace = TRUE,
      prob = pvalues
    )

    sampled_df <-
      data.frame(data)[0, ] # initialize an empty dataframe with the same col names

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
      data_slice
    }))
  }
}


#' Fitting multiple bootstrapped models without aggregaion; called by the function bootstrapFit()
#'
#' @param fit the nlmixr fit object
#' @param nboot an integer giving the number of bootstrapped models to be fit; default value is 100
#' @param nSampIndiv an integer specifying the number of samples in each bootstrapped sample; default is the number of unique subjects in the original dataset
#' @param pvalues a vector of pvalues indicating the probability of each subject to get selected; default value is NULL implying that probability of each subject is the same
#' @param restart a boolean that indicates if a previous session has to be restarted; default value is FALSE
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
                           stratVar,
                           pvalues = NULL,
                           restart = FALSE,
                           fitName = "fit") {
  if (!inherits(fit, "nlmixrFitCore")) {
    stop("'fit' needs to be a nlmixr fit", call. = FALSE)
  }

  if (missing(stratVar)) {
    performStrat <- FALSE
    stratVar <- NULL
  }
  else {
    performStrat <- TRUE
  }

  data <- getData(fit)

  .w <- tolower(names(data)) == "id"
  uidCol <- names(data)[.w]

  checkmate::assert_integerish(nboot,
    len = 1,
    any.missing = FALSE,
    lower = 1
  )

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
  colNames <- names(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uid_colname <- colNames[which("id" %in% colNamesLower)]
  }
  else {
    stop("cannot find the 'ID' column! aborting ...", call. = FALSE)
  }

  uif <- fit$uif
  fitMeth <- getFitMethod(fit)

  bootData <- vector(mode = "list", length = nboot)

  if (is.null(fit$bootstrapMd5)) {
    bootstrapMd5 <- digest::digest(fit)
    assign("bootstrapMd5", bootstrapMd5, envir = fit$env)
  }

  output_dir <-
    paste0("nlmixrBootstrapCache_", fitName, "_", fit$bootstrapMd5) # a new directory with this name will be created
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  else if (dir.exists(output_dir) && restart == TRUE) {
    unlink(output_dir, recursive = TRUE, force = TRUE) # unlink any of the previous directories
    dir.create(output_dir) # create a fresh directory
  }

  fnameBootDataPattern <-
    paste0(as.character(substitute(boot_data)),
      "_", "[0-9]+", ".RData",
      sep = ""
    )
  fileExists <-
    list.files(paste0("./", output_dir), pattern = fnameBootDataPattern)

  if (length(fileExists) == 0) {
    restart <- TRUE
  }

  if (!restart) {
    # read saved bootData from boot_data files on disk
    if (length(fileExists) > 0) {
      cli::cli_alert_success("resuming bootstrap data sampling using data at {paste0('./', output_dir)}")
      
      bootData <- lapply(fileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })

      startCtr <- length(bootData) + 1
    }
    else {
      cli::cli_alert_danger(
        cli::col_red(
          "need the data files at {.file {paste0(getwd(), '/', output_dir)}} to resume"
        )
      )
      stop("aborting...resume file missing", call. = FALSE)
    }
  }

  else {
    startCtr <- 1
  }

  # Generate additional samples (if nboot>startCtr)
  if (nboot >= startCtr) {
    for (mod_idx in startCtr:nboot) {
      bootData[[mod_idx]] <- sampling(
        data,
        nsamp = nSampIndiv,
        uid_colname = uidCol,
        pvalues = pvalues,
        performStrat = performStrat,
        stratVar = stratVar
      )

      # save bootData in curr directory: read the file using readRDS()
      attr(bootData, "randomSeed") <- .Random.seed
      saveRDS(bootData[[mod_idx]],
        file = paste0(
          "./",
          output_dir,
          "/",
          as.character(substitute(boot_data)),
          "_",
          mod_idx,
          ".RData"
        )
      )
    }
  }

  # check if number of samples in stored file is the same as the required number of samples
  fileExists <-
    list.files(paste0("./", output_dir), pattern = fnameBootDataPattern)
  bootData <- lapply(fileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })

  currBootData <- length(bootData)

  # Fitting models to bootData now
  .env <- environment()
  fnameModelsEnsemblePattern <-
    paste0(as.character(substitute(modelsEnsemble)), "_", "[0-9]+",
      ".RData",
      sep = ""
    )
  modFileExists <-
    list.files(paste0("./", output_dir), pattern = fnameModelsEnsemblePattern)
  
  fnameFitEnsemblePattern <-
    paste0(as.character(substitute(fitEnsemble)), "_", "[0-9]+",
           ".RData",
           sep = ""
    )
  fitFileExists <- list.files(paste0("./", output_dir), pattern = fnameFitEnsemblePattern)

  if (!restart) {
    if (length(modFileExists) > 0 &&
      (length(fileExists) > 0)) {
      
      # read bootData and modelsEnsemble files from disk
      cli::cli_alert_success(
        "resuming bootstrap model fitting using data and models stored at {paste0(getwd(), '/', output_dir)}"
      )

      bootData <- lapply(fileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })
      modelsEnsembleLoaded <- lapply(modFileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })
      
      fitEnsembleLoaded <- lapply(fitFileExists, function(x){
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })

      .env$mod_idx <- length(modelsEnsembleLoaded) + 1

      currNumModels <- .env$mod_idx - 1

      if (currNumModels > nboot) {
        cli::cli_alert_danger(
          cli::col_red(
            "the model file already has {.env$mod_idx-1} models when max models is {nboot}; using only the first {nboot} model(s)"
          )
        )
        return(list(modelsEnsembleLoaded[1:nboot], fitEnsembleLoaded[1:nboot]))
        
        # return(modelsEnsembleLoaded[1:nboot])
      }

      else if (currNumModels == nboot) {
        cli::col_red(
          "the model file already has {.env$mod_idx-1} models when max models is {nboot}; loading from {nboot} models already saved on disk"
        )
        return(list(modelsEnsembleLoaded, fitEnsembleLoaded))
        
        # return(modelsEnsembleLoaded)
      }

      else if (currNumModels < nboot) {
        cli::col_red("estimating the additional models ... ")
      }
    }

    else {
      cli::cli_alert_danger(
        cli::col_red(
          "need both the data and the model files at: {paste0(getwd(), '/', output_dir)} to resume"
        )
      )
      stop(
        "aborting...data and model files missing at: {paste0(getwd(), '/', output_dir)}",
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

      fit <- tryCatch(
        {
          fit <- suppressWarnings(nlmixr(uif,
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

          fit # to return 'fit'
        },
        error = function(error_message) {
          print("error fitting the model")
          print(error_message)
          print("storing the models as NA ...")
          return(NA) # return NA otherwise (instead of NULL)
        }
      )

      saveRDS(
        .env$multipleFits,
        file = paste0(
          "./",
          output_dir,
          "/",
          as.character(substitute(modelsEnsemble)),
          "_",
          .env$mod_idx,
          ".RData"
        )
      )
      
      saveRDS(
        fit,
        file = paste0(
          "./",
          output_dir,
          "/",
          as.character(substitute(fitEnsemble)),
          "_",
          .env$mod_idx,
          ".RData"
        )
      )
      
      assign("mod_idx", .env$mod_idx + 1, .env)
    })
  
  fitEnsemble <- NULL

  if (!restart) {
    modelsEnsemble <- c(modelsEnsembleLoaded, modelsEnsemble)
    fitEnsemble <- c(fitEnsembleLoaded, fitEnsemble)
  }

  modFileExists <-
    list.files(paste0("./", output_dir), pattern = fnameModelsEnsemblePattern)
  
  modelsEnsemble <- lapply(modFileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })
  
  fitFileExists <- list.files(paste0("./", output_dir), pattern = fnameFitEnsemblePattern)
  fitEnsemble <- lapply(fitFileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })
  
  

  list(modelsEnsemble, fitEnsemble)
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
      "nlmixrFOCEi" = "focei",
      "nlmixrNlmeUI" = "nlme",
      "nlmixrSaem" = "saem",
      "nlmixrFOCE" = "foce",
      "nlmixrFOi" = "foi",
      "nlmixrFO" = "fo",
      "nlmixrPosthoc" = "posthoc"
    )

  if (!(inherits(fit, "nlmixrFitCore"))) {
    stop("'fit' needs to be a nlmixr fit", call. = FALSE)
  }

  res <- sapply(names(methodsList), function(met) {
    inherits(fit, met)
  })
  .w <- which(res == TRUE)
  if (length(.w) != 1) {
    stop("cannot determine the method the nlmixr fit used, please submit a bug report",
      call. = FALSE
    )
  }
  setNames(methodsList[.w], NULL)
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
getBootstrapSummary <-
  function(fitList,
           ci = 0.95,
           stdErrType = "perc") {
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

        if (stdErrType != "perc") {
          confLower <- mn - qnorm(quantLevels[[2]]) * sd
          confUpper <- mn + qnorm(quantLevels[[3]]) * sd
        }

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

        if (stdErrType != "perc") {
          confLower <- mn - qnorm(quantLevels[[2]]) * sd
          confUpper <- mn + qnorm(quantLevels[[3]]) * sd
        }

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

    class(summaryList) <- "nlmixrBoostrapSummary"
    summaryList
  }

#' @export
print.nlmixrBootstrapSummary <- function(x, ..., sigdig = NULL) {
  if (is.null(sigdig)) {
    if (any(names(x) == "sigdig")) {
      sigdig <- x$sigdig
    } else {
      sigdig <- 3
    }
  }

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
  invsibile(x)
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
