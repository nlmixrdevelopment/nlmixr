.nmMuTable <- function(x) {
  .mu <- x$nmodel$mu.ref
  if (length(.mu) == 0) {
    return(NULL)
  }
  .ret <- do.call(rbind, lapply(names(.mu), function(x) {
    .data.frame(theta = .mu[[x]], eta = x)
  }))
  if (length(x$nmodel$cov.ref) > 0) {
    .covRef <- x$cov.ref
    .covRef <- do.call(`c`, lapply(names(.covRef), function(x) {
      .ret <- .covRef[[x]]
      setNames(.ret, paste0(x, "*", names(.ret)))
    }))
    .env <- new.env(parent = emptyenv())
    .env$ref <- c()
    lapply(names(.covRef), function(x) {
      .n <- .covRef[x]
      if (any(.n == names(.env$ref))) {
        .ref <- .env$ref
        .ref[.n] <- paste0(.env$ref[.n], ", ", x)
      } else {
        .ref <- .env$ref
        .ref[.n] <- x
      }
      assign("ref", .ref, envir = .env)
    })
    if (length(.env$ref) > 0) {
      .ret <- merge(.ret, .data.frame(theta = names(.env$ref), covariates = as.character(.env$ref)),
        by = "theta", all.x = TRUE
      )
      .ret$covariates <- paste(.ret$covariates)
      .ret$covariates[.ret$covariates == "NA"] <- ""
    }
  }
  if (requireNamespace("huxtable", quietly = TRUE)) {
    .ret <- huxtable::hux(.ret) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_position("center") %>%
      huxtable::set_all_borders(TRUE)
  }
  return(.ret)
}

.nmEstMethod <- function(x) {
  .cls1 <- class(x)[1]
  if (.cls1 == "nlmixrNlmeUI") {
    return("nlme")
  }
  if (.cls1 == "nlmixrPosthoc") {
    return("posthoc")
  }
  if (.cls1 == "nlmixrSaem") {
    return("saem")
  }
  if (.cls1 == "nlmixrFOCEi") {
    return("FOCEi")
  }
  if (.cls1 == "nlmixrFOCE") {
    return("FOCE")
  }
  if (.cls1 == "nlmixrFOi") {
    return("FOi")
  }
  if (.cls1 == "nlmixrFO") {
    return("FO")
  }
  if (.cls1 == "nlmixrPosthoc") {
    return("posthoc")
  }
}

.yamlHeader <- function(x, bound = "bound") {
  .cor <- x$omega
  diag(.cor) <- 0
  .lst <- list(
    "Full nlmixr Version:" = .getFullNlmixrVersion(),
    "Full RxODE Version:" = .getFullNlmixrVersion("RxODE"),
    "R Model Function Name ($modelName):" = x$modelName,
    "R Fit Object:" = bound,
    "R Data Name ($dataName):" = x$dataName,
    "Estimation Method:" = .nmEstMethod(x),
    "All variables mu-referenced:" = dim(x$omega)[1] == length(x$mu.ref),
    "Covariance Method ($covMethod):" = x$covMethod,
    "Modeled Correlations:" = !all(.cor == 0),
    "CWRES:" = ifelse(inherits(x, "nlmixrFitData"), any(names(x) == "CWRES"), FALSE),
    "NPDE:" = ifelse(inherits(x, "nlmixrFitData"), any(names(x) == "NPDE"), FALSE),
    "Goodness Of Fit Metrics:" = as.list(x$objDf),
    "Timing Metrics" = as.list(x$time)
  )
  return(yaml::as.yaml(.lst))
}



.nmHuxHeader <- function(x, bound) {
  .cor <- x$omega
  diag(.cor) <- 0
  .tab <- tibble::tribble(
    ~Description, ~Value,
    "Full nlmixr Version:", .getFullNlmixrVersion(),
    "Full RxODE Version:", .getFullNlmixrVersion("RxODE"),
    "R Model Function Name ($modelName):", x$modelName,
    "R Fit Object:", bound,
    "R Data Name ($dataName):", x$dataName,
    "Estimation Method:", .nmEstMethod(x),
    "All variables mu-referenced:", paste(dim(x$omega)[1] == length(x$mu.ref)),
    "Covariance Method ($covMethod):", x$covMethod,
    "Modeled Correlations:", paste0(!all(.cor == 0)),
    "CWRES:", paste0(ifelse(inherits(x, "nlmixrFitData"), any(names(x) == "CWRES"), FALSE)),
    "NPDE:", paste0(ifelse(inherits(x, "nlmixrFitData"), any(names(x) == "NPDE"), FALSE)),
  )
  if (requireNamespace("huxtable", quietly = TRUE)) {
    huxtable::as_hux(.tab) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = huxtable::everywhere, col = 1, value = TRUE) %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_all_borders(TRUE) ->
    .tab
  }
  return(.tab)
}


.nmHuxTime <- function(x) {
  .tab <- x$time
  if (requireNamespace("huxtable", quietly = TRUE)) {
    huxtable::as_hux(.tab) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_all_borders(TRUE) -> .tab
  }
  return(.tab)
}


.nmHuxObjf <- function(x) {
  .tmp <- x$objDf
  .tmp <- .data.frame(type = row.names(.tmp), .tmp, check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  if (requireNamespace("huxtable", quietly = TRUE)) {
    huxtable::as_hux(.tmp) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_all_borders(TRUE) ->
    .tmp
  }
  return(.tmp)
}

.getFullNlmixrVersion <- function(pkg = "nlmixr") {
  .x <- tryCatch(utils::packageDescription(pkg, lib.loc = .libPaths()),
    error = function(e) NA, warning = function(e) NA
  )
  .pv <- utils::packageVersion(pkg)
  if (identical(.x, NA) || length(.x$RemoteSha) == 0) {
    return(paste(.pv))
  } else {
    return(sprintf("%s %s", .pv, substr(paste(.x$RemoteSha), 1, 10)))
  }
}

.nmGetPrintLines <- function(x) {
  .tempfile <- tempfile()
  .zz <- file(.tempfile, open = "wt")
  .unsink <- FALSE
  sink(.zz, type = "output")
  sink(.zz, type = "message")
  on.exit({
    if (!.unsink) {
      sink(type = "output")
      sink(type = "message")
      close(.zz)
    }
    unlink(.tempfile)
  })
  if (is.null(x)) {
    print(sessionInfo())
  } else {
    print(x)
  }
  sink(type = "output")
  sink(type = "message")
  close(.zz)
  .unsink <- TRUE
  suppressWarnings(readLines(.tempfile))
}

.nmIterHist <- function(x) {
  if (inherits(x, "nlmixrSaem")) {
    .tmp <- x$parHist
    .nit <- x$mcmc$niter
    return(list(
      "SA iterations" = .tmp[seq(1, .nit[1]), ],
      "EM iterations" = .tmp[-seq(1, .nit[1]), ]
    ))
  } else if (inherits(x, "nlmixrNlmeUI") || inherits(x, "posthoc")) {
    return(NULL)
  } else if (exists("parHistData", x$env)) {
    .tmp <- x$parHistData$type
    class(.tmp) <- NULL
    .w <- which(.tmp <= 4)
    if (length(.w) > 0) {
      .grad <- x$parHistData[.w, ]
      .grad <- .grad[, names(.grad) != "objf"]
      .est <- x$parHistData[-.w, ]
    } else {
      .grad <- NULL
      .est <- x$parHistData
    }
    .scale <- .est[.est$type == "Scaled", ]
    .scale <- .scale[, names(.scale) != "type"]
    .unscale <- .est[.est$type == "Unscaled", ]
    .unscale <- .unscale[, names(.unscale) != "type"]
    .back <- .est[.est$type == "Back-Transformed", ]
    .back <- .back[, names(.back) != "type"]
    .ret <- list("Scaled" = .scale, "Unscaled" = .unscale, "Back-Transformed" = .back)
    if (!is.null(.grad)) {
      .ret <- c(.ret, list("Gradient" = .grad))
    }
    return(.ret)
  }
}

.nmIterHistDoc <- function(x, doc, headerStyle, normalStyle, preformatStyle, width) {
  .x <- .nmIterHist(x)
  .doc <- doc
  .nrows <- getOption("datatable.print.nrows")
  options(datatable.print.nrows = 11)
  on.exit(options(datatable.print.nrows = .nrows))
  if (length(.x) > 0) {
    .doc <- .doc %>% officer::body_add_par("Iteration History", style = headerStyle)
  }
  for (.i in seq_along(.x)) {
    .doc <- .doc %>%
      officer::body_add_par(names(.x)[.i], style = normalStyle)
    .doc <- .nmDocxPreformat(
      data.table::data.table(.x[[.i]]), doc,
      preformatStyle, width
    )
  }
  return(.doc)
}

.nmDocxPreformat <- function(x, doc, preformatStyle, width) {
  .cw <- Sys.getenv("RSTUDIO_CONSOLE_WIDTH")
  .width <- getOption("width")
  options(width = width)
  Sys.setenv("RSTUDIO_CONSOLE_WIDTH" = width)
  on.exit({
    options(width = .width)
    if (.cw == "") {
      Sys.unsetenv("RSTUDIO_CONSOLE_WIDTH")
    } else {
      Sys.setenv("RSTUDIO_CONSOLE_WIDTH" = .cw)
    }
  })
  .doc <- doc
  .lines <- .nmGetPrintLines(x)
  for (.i in seq_along(.lines)) {
    .doc <- officer::body_add_par(.doc, .lines[.i], style = preformatStyle)
  }
  .doc
}

##' Create a run summary word document
##'
##' @param x nlmixr fit object.
##'
##' @param docxOut Output file for run information document.  If not
##'     specified it is the name of R object where the fit is located
##'     with the \code{-YEAR-MONTH-DAY.docx} appended.  If it is
##'     \code{NULL} the document is not saved, but the \code{officer}
##'     object is returned.
##'
##' @param docxTemplate This is the document template.  If not
##'     specified it defaults to
##'     \code{option("nlmixr.docx.template")}.  If
##'     \code{option("nlmixr.docx.template")} is not specified it uses
##'     the included nlmixr document template.  When
##'     \code{docxTemplate=NULL} it uses the \code{officer} blank
##'     document.
##'
##' @param plot Boolean indicating if the default goodness of fit
##'     plots are added to the document.  By default \code{TRUE}
##'
##' @param titleStyle This is the word style name for the nlmixr
##'     title; Usually this is \code{nlmixr version (R
##'     object)}. Defaults to \code{option("nlmixr.docx.title")} or
##'     \code{Title}
##'
##' @param subtitleStyle This is the word style for the subtitle which
##'     is \code{nlmixr model name and date}. Defaults to
##'     \code{option("nlmixr.docx.subtitle")} or \code{Subtitle}
##'
##' @param normalStyle This is the word style for normal text. Defaults to
##'     \code{option("nlmixr.docx.normal")} or \code{Normal}
##'
##' @param headerStyle This is the word style for heading text. Defaults to
##'     \code{option("nlmixr.docx.heading1")} or \code{Heading 1}
##'
##' @param centeredStyle This is the word style for centered text
##'     which is used for the figures. Defaults to
##'     \code{option("nlmixr.docx.centered")} or \code{centered}
##'
##' @param preformattedStyle This is the preformatted text style for R
##'     output lines.  Defaults to
##'     \code{option("nlmixr.docx.preformatted")} or \code{HTML Preformatted}
##'
##' @param width Is an integer representing the number of characters
##'     your preformatted style supports.  By default this is
##'     \code{option("nlmixr.docx.width")} or \code{69}
##'
##' @param save Should the docx be saved in a zip file with the R rds
##'     data object for the fit?  By default this is \code{FALSE} with
##'     \code{nmDocx} and \code{TRUE} with \code{nmSave}
##'
##' @param ... when using `nmSave` these arguments are passed to `nmDocx`
##'
##' @return An officer docx object
##'
##' @author Matthew Fidler
##'
##' @examples
##' \dontrun{
##'  library(nlmixr)
##' pheno <- function() {
##'     # Pheno with covariance
##'   ini({
##'     tcl <- log(0.008) # typical value of clearance
##'     tv <-  log(0.6)   # typical value of volume
##'     ## var(eta.cl)
##'     eta.cl + eta.v ~ c(1,
##'                        0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
##'                       # interindividual variability on clearance and volume
##'     add.err <- 0.1    # residual variability
##'   })
##'   model({
##'     cl <- exp(tcl + eta.cl) # individual value of clearance
##'     v <- exp(tv + eta.v)    # individual value of volume
##'     ke <- cl / v            # elimination rate constant
##'     d/dt(A1) = - ke * A1    # model differential equation
##'     cp = A1 / v             # concentration in plasma
##'     cp ~ add(add.err)       # define error model
##'   })
##' }
##'
##' fit.s <- nlmixr(pheno, pheno_sd, "saem")
##'
##' ## Save output information into a word document
##' nmDocx(fit.s)
##'
##' }
##' @export
nmDocx <- function(x,
                   docxOut = NULL,
                   docxTemplate = NULL,
                   plot = TRUE,
                   titleStyle = getOption("nlmixr.docx.title", "Title"),
                   subtitleStyle = getOption("nlmixr.docx.subtitle", "Subtitle"),
                   normalStyle = getOption("nlmixr.docx.normal", "Normal"),
                   headerStyle = getOption("nlmixr.docx.heading1", "Heading 1"),
                   centeredStyle = getOption("nlmixr.docx.centered", "centered"),
                   preformattedStyle = getOption("nlmixr.docx.preformatted", "HTML Preformatted"),
                   width = getOption("nlmixr.docx.width", 69),
                   save = FALSE) {
  RxODE::rxReq("officer")
  RxODE::rxReq("flextable")
  RxODE::rxReq("data.table")
  RxODE::rxReq("huxtable")
  if (!inherits(x, "nlmixrFitCore")) {
    stop("This only applies to nlmixr fit objects")
  }
  .parent <- globalenv()
  .bound <- ""
  .bound2 <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], x)) {
      return(.cur)
    }
    return(NULL)
  }))
  if (length(.bound2) > 0) {
    .bound <- .bound2[order(sapply(.bound2, nchar))]
    .bound <- .bound[1]
  }
  if (missing(docxTemplate)) {
    docxTemplate <- getOption(
      "nlmixr.docx.template",
      file.path(
        system.file(package = "nlmixr"),
        "nlmixr-template.docx"
      )
    )
  }
  if (missing(docxOut)) {
    docxOut <- file.path(getwd(), paste0(
      ifelse(any(.bound == c("", ".")), "nlmixr", .bound),
      "-", .nmEstMethod(x), "-", x$modelName,
      format(Sys.time(), "-%Y-%m-%d.docx")
    ))
  }
  .doc <- officer::read_docx(docxTemplate)
  .style <- officer::styles_info(.doc)
  if (!any(titleStyle == .style$style_name)) stop("Need to specify a valid titleStyle")
  if (!any(subtitleStyle == .style$style_name)) stop("Need to specify a valid subtitleStyle")
  if (!any(headerStyle == .style$style_name)) stop("Need to specify a valid headerStyle")
  if (!any(preformattedStyle == .style$style_name)) stop("Need to specify a valid preformattedStyle")
  if (.bound != "") {
    .hreg <- eval(parse(text = sprintf('huxtable::huxreg("%s"=x)', .bound)))
  } else {
    .hreg <- huxtable::huxreg(x)
  }
  .doc <- .doc %>%
    officer::body_add_par(sprintf(
      "nlmixr %s (%s)", utils::packageVersion("nlmixr"),
      .bound
    ), style = titleStyle) %>%
    officer::body_add_par(sprintf(
      "function: %s; %s", x$modelName,
      format(Sys.time(), "%a %b %d %X %Y")
    ),
    style = subtitleStyle
    ) %>%
    officer::body_add_par("", style = normalStyle) %>%
    officer::body_add_toc()
  .doc <- .doc %>%
    flextable::body_add_flextable(flextable::autofit(.asFlx(.nmHuxHeader(x, .bound))))
  .doc <- .doc %>%
    officer::body_add_par("Timing ($time, in seconds)", style = headerStyle) %>%
    flextable::body_add_flextable(flextable::autofit(.asFlx(.nmHuxTime(x))))
  .doc <- .doc %>%
    officer::body_add_par(sprintf("Parameter Information huxreg(%s)", .bound), style = headerStyle) %>%
    flextable::body_add_flextable(flextable::autofit(.asFlx(.hreg)))
  doc <- .doc %>%
    officer::body_add_par("UI Model information ($uif):", style = headerStyle)
  .doc <- .nmDocxPreformat(x$uif, .doc, preformattedStyle, width)

  if (x$message != "") {
    .doc <- .doc %>%
      officer::body_add_par("Minimization Information ($mesage):", style = headerStyle) %>%
      officer::body_add_par(paste0("$message: ", x$message), style = preformattedStyle)
  }
  if (length(x$warnings) > 0) {
    .doc <- .doc %>%
      officer::body_add_par("Warning(s) ($warnings):", style = headerStyle)
    for (.i in seq_along(x$warnings)) {
      .doc <- officer::body_add_par(.doc, x$warnings[.i], style = preformattedStyle)
    }
  }
  .mu <- .nmMuTable(x)
  print
  if (length(.mu) != 0) {
    .doc <- .doc %>%
      officer::body_add_par("Parsed Mu-referencing:", style = headerStyle) %>%
      flextable::body_add_flextable(flextable::autofit(.asFlx(.mu)))
  }
  .doc <- .nmIterHistDoc(x, .doc, headerStyle, normalStyle, preformattedStyle, width)
  if (any(.nmEstMethod(x) == c("FOCEi", "FOCE", "FO", "FOi", "posthoc"))) {
    .doc <- .doc %>%
      officer::body_add_par("Scaling and gradient information ($scaleInfo):", style = headerStyle)
    .si <- x$scaleInfo
    .t1 <- c("est", "scaleC")
    .t2 <- c("est", "Initial Gradient", "Forward aEps", "Forward rEps", "Central aEps", "Central rEps")
    .t3 <- c("est", "Covariance Gradient", "Covariance aEps", "Covariance rEps")
    .t1 <- .si[, .t1]
    .t2 <- .si[, .t2]
    .t3 <- .si[, .t3]
    .t1 <- huxtable::hux(.t1) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_position("center") %>%
      huxtable::set_all_borders(TRUE)
    .t2 <- huxtable::hux(.t2) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_position("center") %>%
      huxtable::set_all_borders(TRUE)
    .t3 <- huxtable::hux(.t3) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_position("center") %>%
      huxtable::set_all_borders(TRUE)
    .doc <- .doc %>%
      officer::body_add_par("Scaling information (est/scaleC):", style = normalStyle) %>%
      flextable::body_add_flextable(flextable::autofit(.asFlx(.t1))) %>%
      officer::body_add_par("Gradient Information, used in estimation:", style = normalStyle) %>%
      flextable::body_add_flextable(flextable::autofit(.asFlx(.t2))) %>%
      officer::body_add_par("Gradient Information, used in covariance:", style = normalStyle) %>%
      flextable::body_add_flextable(flextable::autofit(.asFlx(.t3)))
  }
  .doc <- .doc %>%
    officer::body_add_par("Original arguments to control ($origControl):", style = headerStyle)
  .doc <- .nmDocxPreformat(x$origControl, .doc, preformattedStyle, width)
  if (plot && inherits(x, "nlmixrFitData")) {
    .doc <- .doc %>%
      officer::body_add_par(sprintf("Basic Goodness of fit, ie plot(%s)", .bound), style = headerStyle)
    .lst <- plot(x)
    for (.i in seq_along(.lst)) {
      .doc <- officer::body_add_gg(.doc, value = .lst[[.i]], style = centeredStyle)
    }
  }
  .doc <- .doc %>%
    officer::body_add_par("sessionInfo():", style = headerStyle)
  .doc <- .nmDocxPreformat(NULL, .doc, preformattedStyle, width)
  if (!is.null(docxOut)) {
    print(.doc, target = docxOut)
    message(sprintf("Printed document information to: %s", docxOut))
    if (save) {
      .rdsOut <- sub("[.]docx", ".rds", docxOut)
      saveRDS(x, file = .rdsOut)
      message(sprintf("Saved R object to: %s", .rdsOut))
      .zipOut <- sub("[.]docx", ".zip", docxOut)
      if (file.exists(.zipOut)) {
        warning(sprintf("Overwriting %s", .zipOut))
        unlink(.zipOut)
      }
      utils::zip(.zipOut, c(docxOut, .rdsOut), flags = "-9Xj")
      if (file.exists(.zipOut)) {
        message(sprintf("Zipped to: %s", .zipOut))
        unlink(.rdsOut)
        unlink(docxOut)
        message("Removed rds/docx")
      }
    }
    return(invisible(.doc))
  } else {
    return(.doc)
  }
}

##' Create a large output based on a nlmixr fit
##'
##' @param x nlmixr fit object
##'
##' @param lst Listing file.  If not specified, it is determined by
##'     the day and the model/R-object name.  If it is specified as
##'     \code{NULL} the listing output is displayed on the screen.
##'
##' @return invisibly returns fit
##'
##' @author Matthew Fidler
##' @export
nmLst <- function(x,
                  lst = NULL) {
  RxODE::rxReq("data.table")
  RxODE::rxReq("yaml")
  RxODE::rxReq("huxtable")
  if (!inherits(x, "nlmixrFitCore")) {
    stop("This only applies to nlmixr fit objects")
  }
  .parent <- globalenv()
  .bound <- ""
  .bound2 <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], x)) {
      return(.cur)
    }
    return(NULL)
  }))
  if (length(.bound2) > 0) {
    .bound <- .bound2[order(sapply(.bound2, nchar))]
    .bound <- .bound[1]
  }
  if (.bound != "") {
    .hreg <- eval(parse(text = sprintf('huxtable::huxreg("%s"=x)', .bound)))
  } else {
    .hreg <- huxtable::huxreg(x)
  }
  .rule <- function(msg) {
    message("")
    cli::rule(msg)
    message("")
  }
  if (missing(lst)) {
    lst <- file.path(getwd(), paste0(
      ifelse(any(.bound == c("", ".")), "nlmixr", .bound),
      "-", .nmEstMethod(x), "-", x$modelName,
      format(Sys.time(), "-%Y-%m-%d.lst")
    ))
  }
  if (!is.null(lst)) {
    .zz <- file(lst, open = "wt")
    .unsink <- FALSE
    sink(.zz, type = "output")
    sink(.zz, type = "message")
    on.exit({
      if (!.unsink) {
        sink(type = "output")
        sink(type = "message")
        close(.zz)
      }
    })
  }
  message("---")
  message(.yamlHeader(x, .bound))
  message("---")
  message("")
  .rule(sprintf("Parameter Information as_hux(%s)", .bound))
  .hreg %>% huxtable::print_screen(colnames = FALSE)
  .rule(sprintf("Parsed model information %s$uif", .bound))
  print(x$uif)
  if (x$message != "") {
    .rule("Minimization Information ($mesage):")
    message(paste("$message: ", x$message, collapse = " "))
  }
  if (length(x$warnings) > 0) {
    .rule("Warning(s) ($warnings):")
    for (.i in seq_along(x$warnings)) {
      message(x$warnings[.i])
    }
  }
  .mu <- .nmMuTable(x)
  if (length(.mu) != 0) {
    .rule("Parsed Mu-referencing:")
    .mu %>% huxtable::print_screen(colnames = FALSE)
  }
  .it <- .nmIterHist(x)
  if (length(.it) > 0) {
    .rule("Iteration History")
    for (.i in seq_along(.it)) {
      message(paste0(names(.it)[.i], ":"))
      print(data.table::data.table(.it[[.i]]), nrows = 11)
    }
  }
  if (any(.nmEstMethod(x) == c("FOCEi", "FOCE", "FO", "FOi", "posthoc"))) {
    .rule("Scaling and gradient information ($scaleInfo):")
    print(x$scaleInfo)
  }
  if (!is.null(lst)) {
    .rule("Original arguments to control ($origControl):")
    print(x$origControl)
    .rule("sessionInfo():")
    print(sessionInfo())
    sink(type = "output")
    sink(type = "message")
    close(.zz)
    .unsink <- TRUE
    message(sprintf("Wrote listing to %s", lst))
  }
  invisible(x)
}

##' @rdname nmDocx
##' @export
nmSave <- function(x, ..., save = TRUE) {
  nmDocx(x, ..., save = save)
}

## Ugly hacks to make word conversion work without loading flextable, officer, etc
.huxNS <- NULL
.asFlx <- function(x, ...) {
  RxODE::rxReq("huxtable")
  if (is.null(.huxNS)) {
    ## ugly hack
    .huxNS <- loadNamespace("huxtable")
    assignInMyNamespace(".huxNS", .huxNS)
  }
  .env <- new.env(parent = .huxNS)
  .env$x <- x
  with(.env, as_flextable.huxtable(x))
}
