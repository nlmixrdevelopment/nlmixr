##' Change a nlmixr fit object to a huxtable
##'
##' This is a thin layer that differs from
##' \code{\link[hustable]{huxreg}} to make it easier to produce
##' reasonable huxtables with nlmixr fit objects.
##'
##' \itemize{
##'
##' \item Drops \code{NA} values from tables, so (\code{NA}) is not
##' shown for estimates without a standard error.
##'
##' \item Filters out any blank rows.
##'
##' \item Removes R.squared from the statistics and replaces with the
##' objective function value.
##'
##' \item Adjust broom separators for an easier to read format.
##'
##' }
##'
##' @inheritParams huxtable::huxreg
##'@export
asHux.nlmixrFitCore  <- function(...,
                                 error_format    = "({std.error})",
                                 error_style     = c("stderr", "ci", "statistic", "pvalue"),
                                 error_pos       = c("below", "same", "right"),
                                 number_format   = "%.3f",
                                 align           = ".",
                                 pad_decimal     = ".",
                                 ci_level        = NULL,
                                 tidy_args       = NULL,
                                 stars           = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                                 bold_signif     = NULL,
                                 borders         = 0.4,
                                 outer_borders   = 0.8,
                                 note            = if (is.null(stars)) NULL else "{stars}.",
                                 statistics      = c("N" = "nobs", "Objective Function" = "OBJF", "logLik", "AIC"),
                                 coefs           = NULL,
                                 omit_coefs      = NULL,
                                 na_omit         = c("all", "any", "none")){
    RxODE::rxReq("huxtable")
    .opt  <- options()
    on.exit(options(.opt));
    options(broom.mixed.sep1=": ",
            broom.mixed.sep2=", ")
    .lst <- as.list(match.call(expand.dots=TRUE))[-1];
    ## Extract all of the error format columns
    .glueReg  <- "[{][^}]*[}]"
    .isComplex  <- FALSE
    .errCols <- sapply(unlist(stringr::str_extract_all(error_format, .glueReg)),
                       function(x){
        .err  <- substr(x, 2, nchar(x) - 1);
        if (is.name(eval(parse(text=sprintf("quote(%s)", .err))))){
            return(.err)
        } else {
            .isComplex <<- TRUE
            return(NULL)
        }
    })
    if (!.isComplex){
        .lst$error_format <- paste0("{ifelse(!sapply(seq_along(",.errCols[1],"), function(x){any(is.na(c(",
                                    paste(paste0(.errCols,"[x]"),collapse=","),
                                    ")))}), paste0(\"",
                                    gsub(rex::rex("{",capture(or(.errCols)),"}"), "\",\\1,\"",error_format),
                                    "\"),\"\")}");
    }
    .lst$statistics  <- statistics
    .ret <- do.call(huxtable::huxreg,.lst)
    .w  <- which(sapply(seq_along(.ret$names),function(i){!all(as.character(.ret[i,])=="")}))
    return(.ret[.w])
}

##'@importFrom huxtable as_hux
##'@export
as_hux  <- huxtable::as_hux

##'@importFrom huxtable as_huxtable
##'@export
as_huxtable  <- huxtable::as_huxtable

# This shouldn't be necessary, but it is :(
##'@S3method as_huxtable nlmixrFitCore
##'@export as_huxtable.nlmixrFitCore
as_huxtable.nlmixrFitCore  <- function(x,...){
    if (missing(x)){
        .args <- list(...)
    } else {
        .args  <- list(x,...);
    }
    do.call(asHux.nlmixrFitCore, .args)
}


.nmHuxHeader  <- function(x, bound){
    huxtable::tribble_hux(
                  ~ Description, ~ Value,
                  "Full nlmixr Version:", .getFullNlmixrVersion(),
                  "R Model Function Name ($modelName):",  x$modelName,
                  "R Fit Object", bound,
                  "R Data Name ($dataName):", x$dataName
              ) %>%
        huxtable::add_colnames() %>%
            huxtable::set_bold(row = huxtable::everywhere, col = 1, value = TRUE) %>%
            huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
            huxtable::set_all_borders(TRUE)
}


.nmHuxTime  <- function(x){
    huxtable::as_hux(x$time) %>%
        huxtable::add_colnames() %>%
        huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
        huxtable::set_all_borders(TRUE)
}


.nmHuxObjf  <-function(x){
    .tmp  <- x$objDf;
    .tmp  <- data.frame(type=row.names(.tmp),.tmp, check.rows=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
    huxtable::as_hux(.tmp) %>%
        huxtable::add_colnames() %>%
        huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
        huxtable::set_all_borders(TRUE)
}

.getFullNlmixrVersion  <- function(){
    .x <- tryCatch(utils::packageDescription("nlmixr", lib.loc = .libPaths()),
                   error = function(e) NA, warning = function(e) NA)
    .pv  <- packageVersion("nlmixr");
    if (identical(.x, NA)) {
        return(paste(.pv))
    } else {
        return(sprintf("%s %s", .pv, substr(paste(.x$RemoteSha),1,10)));
    }
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
nmDocx  <- function(x,
                    docxOut=NULL,
                    docxTemplate=NULL,
                    titleStyle=getOption("nlmixr.docx.title", "Title"),
                    subtitleStyle=getOption("nlmixr.docx.subtitle","Subtitle"),
                    normalStyle=getOption("nlmixr.docx.normal", "Normal"),
                    headerStyle=getOption("nlmixr.docx.heading1", "Heading 1"),
                    centeredStyle=getOption("nlmixr.docx.centered","centered")
                    ){
    RxODE::rxReq("officer");
    RxODE::rxReq("flextable");
    if (!inherits(x, "nlmixrFitCore")){
        stop("This only applies to nlmixr fit objects");
    }
    .parent <- globalenv();
    .bound  <- "";
    .bound2 <- do.call("c", lapply(ls(.parent), function(.cur){
                                if (identical(.parent[[.cur]], x)){
                                    return(.cur)
                                }
                                return(NULL);
                            }))
    if (length(.bound2) > 0){
        .bound <- .bound2[order(sapply(.bound2, nchar))];
        .bound <- .bound[1]
    }
    if (missing(docxTemplate)){
        docxTemplate  <- getOption("nlmixr.docx.template",
                                   file.path(system.file(package="nlmixr"),
                                             "nlmixr-template.docx"))
    }
    if (missing(docxOut)){
        docxOut  <- file.path(getwd(),paste0(ifelse(any(.bound==c("",".")),"nlmixr",.bound),
                                             format(Sys.time(), "-%Y-%m-%d.docx")))
    }
    .doc <- officer::read_docx(docxTemplate);
    .style  <- officer::styles_info(.doc)
    if (!any(titleStyle==.style$style_name)) stop("Need to specify a valid titleStyle");
    if (!any(subtitleStyle==.style$style_name)) stop("Need to specify a valid subtitleStyle");
    if (!any(headerStyle==.style$style_name)) stop("Need to specify a valid headerStyle");
    if (.bound !=""){
        .hreg  <- eval(parse(text=sprintf('nlmixr::as_hux("%s"=x)', .bound)));
    } else {
        .hreg  <- nlmixr::as_hux(x);
    }
    .doc <- .doc %>%
        officer::body_add_par(sprintf("nlmixr %s (%s)", packageVersion("nlmixr"),
                                      .bound), style=titleStyle) %>%
        officer::body_add_par(sprintf("function: %s; %s", x$modelName,
                                      format(Sys.time(), "%a %b %d %X %Y")),
                              style=subtitleStyle) %>%
        officer::body_add_par("",style=normalStyle) %>%
        flextable::body_add_flextable(flextable::autofit(huxtable::as_flextable(.nmHuxHeader(x, .bound)))) %>%
        officer::body_add_par("Run Time ($time)", style=headerStyle) %>%
        flextable::body_add_flextable(flextable::autofit(huxtable::as_flextable(.nmHuxTime(x)))) %>%
        officer::body_add_par(sprintf("Parameter Information as_hux(%s)", .bound), style=headerStyle) %>%
        flextable::body_add_flextable(flextable::autofit(huxtable::as_flextable(.hreg))) %>%
        officer::body_add_par(sprintf("Basic Goodness of fit, ie plot(%s)",.bound), style=headerStyle)
    if (inherits(x, "nlmixrFitData")){
        .lst <- plot(x);
        for (.i in seq_along(.lst)){
            .doc <- officer::body_add_gg(.doc, value = .lst[[.i]], style = "centered")
        }
    }
    if (!is.null(docxOut)){
        print(.doc, target=docxOut);
        message(sprintf("Printed document information to: %s", docxOut))
        return(invisible(.doc));
    } else {
        return(.doc)
    }
}
