#' Add CWRES
#'
#' This returns a new fit object with CWRES attached
#'
#' @param fit nlmixr fit without WRES/CWRES
#' @param updateObject Boolean indicating if the original fit object
#'     should be updated. By default this is true.
#' @param envir Environment that should be checked for object to
#'     update.  By default this is the global environment.
#' @return fit with CWRES
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
#' f <- try(nlmixr(one.cmt, theo_sd, "saem"))
#'
#' print(f)
#'
#' # even though you may have forgotten to add the cwres, you can add it to the data.frame:
#'
#' if (!inherits(f, "try-error")) {
#'   f <- try(addCwres(f))
#'   print(f)
#' }
#'
#' # Note this also adds the FOCEi objective function
#' }
#' @author Matthew L. Fidler
#' @export
addCwres <- function(fit, updateObject = TRUE, envir = parent.frame(1)) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .pt <- proc.time()
  .oTime <- fit$env$time
  .objName <- substitute(fit)
  if (any(names(fit) == "CWRES")) {
    ## warning("Already contains CWRES")
    return(fit)
  }
  .uif <- fit$uif
  .saem <- fit$saem
  .od <- fit$origData
  if (!is.null(.saem)) {
    .calcResid <- inherits(fit, "nlmixrFitData")
    if (.calcResid) {
      assign("saem", NULL, fit$env)
      on.exit({
        assign("saem", .saem, fit$env)
      })
      .newFit <- as.focei.saemFit(.saem, .uif,
        data = .nmGetData(fit), calcResid = TRUE, obf = NA,
        calcCov = fit$cov, covMethod = fit$covMethod,
        calcCovTime = as.vector(fit$time[["covariance"]])
      )
      assign("saem", .saem, fit$env)
      assign("adjObj", fit$env$adjObj, .newFit$env)
      .df <- .newFit[, c("WRES", "CRES", "CWRES", "CPRED")]
      ## Add CWRES timing to fit.
      .new <- cbind(fit, .df)
    } else {
      .newFit <- nlmixr(fit, getData(fit), "focei",
        control = foceiControl(
          maxOuterIterations = 0L, maxInnerIterations = 0L,
          etaMat = as.matrix(fit$eta[, -1]), calcResid = FALSE,
          covMethod = ""
        )
      )
    }
    .ob1 <- .newFit$objDf
    .ob2 <- fit$objDf
    if (any(names(.ob2) == "Condition Number")) {
      .cn <- unique(.ob2[["Condition Number"]])
      .cn <- .cn[!is.na(.cn)]
      if (length(.cn) == 1) {
        .ob1[, "Condition Number"] <- .cn
      } else {
        .ob1[, "Condition Number"] <- NA
      }
    }
    .ob1 <- rbind(.ob1, .ob2)
    .ob1 <- .ob1[order(row.names(.ob1)), ]
    .ob1 <- .ob1[!is.na(.ob1$OBJF), ]
    if (.calcResid) {
      assign("objDf", .ob1, envir = .newFit$env)
    } else {
      assign("objDf", .ob1, envir = fit$env)
      setOfv(fit, "FOCEi")
      .env <- fit$env
      .env$time <- .data.frame(.oTime, foceiLik = (proc.time() - .pt)["elapsed"], check.names = FALSE)
      return(fit)
    }
  }
  .nlme <- fit$nlme
  if (!is.null(.nlme)) {
    .newFit <- as.focei(.nlme, .uif, data = .nmGetData(fit), calcResid = TRUE)
    assign("adjObj", fit$env$adjObj, .newFit$env)
    .df <- .newFit[, c("WRES", "CRES", "CWRES", "CPRED")]
    .new <- cbind(fit, .df)
  }
  class(.new) <- class(.newFit)
  if (!is.null(.saem)) {
    setOfv(.new, "FOCEi")
  }
  if (updateObject) {
    .parent <- envir
    .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
      if (.cur == .objName && identical(.parent[[.cur]]$env, fit$env)) {
        return(.cur)
      }
      return(NULL)
    }))
    if (length(.bound) == 1) {
      if (exists(.bound, envir = .parent)) {
        assign(.bound, .new, envir = .parent)
      }
    }
  }
  .env <- .new$env
  .env$time <- .data.frame(.oTime, cwres = (proc.time() - .pt)["elapsed"], check.names = FALSE)
  assign("origData", .od, .env)
  return(.new)
}

.setOfvFo <- function(fit, type = c("focei", "foce", "fo")) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .pt <- proc.time()
  .type <- match.arg(type)
  .oTime <- fit$env$time
  .etaMat <- as.matrix(fit$eta[, -1]) ## with fo should be NULL
  .fo <- FALSE
  if (.type == "focei") {
    .interaction <- TRUE
    .rn <- "FOCEi"
  } else if (.type == "foce") {
    .interaction <- FALSE
    .rn <- "FOCE"
  } else {
    .interaction <- FALSE
    .fo <- TRUE
    .etaMat <- NULL
    .rn <- "FO"
  }
  .inObjDf <- fit$objDf
  if (any(rownames(.inObjDf) == .rn)) {
    return(fit)
  }
  .newFit <- nlmixr(fit, getData(fit), "focei",
    control = foceiControl(
      maxOuterIterations = 0L, maxInnerIterations = 0L,
      interaction = .interaction, fo = .fo,
      etaMat = .etaMat, calcResid = FALSE,
      covMethod = ""
    )
  )
  .env <- fit$env
  .ob1 <- .newFit$objDf
  if (any(names(.inObjDf) == "Condition Number")) {
    .cn <- unique(.inObjDf[["Condition Number"]])
    .cn <- .cn[!is.na(.cn)]
    if (length(.cn) == 1) {
      .ob1[, "Condition Number"] <- .cn
    } else {
      .ob1[, "Condition Number"] <- NA
    }
  }
  .ob1 <- rbind(.ob1, .inObjDf)
  .ob1 <- .ob1[order(row.names(.ob1)), ]
  .ob1 <- .ob1[!is.na(.ob1$OBJF), ]
  assign("objDf", .ob1, fit$env)
  .env$time <- .data.frame(foceiLik = (proc.time() - .pt)["elapsed"], .oTime, check.names = FALSE)
  names(.env$time)[1] <- paste0(.type, "Lik")
  setOfv(fit, .type)
  return(invisible(fit))
}

##' Set/get Objective function type for a nlmixr object
##'
##' @param x nlmixr fit object
##' @param type Type of objective function to use for AIC, BIC, and
##'     $objective
##' @return Nothing
##' @author Matthew L. Fidler
##' @export
setOfv <- function(x, type) {
  if (inherits(x, "nlmixrFitCore") || inherits(x, "nlmixrFitCoreSilent")) {
    .objDf <- x$objDf
    .w <- which(tolower(row.names(.objDf)) == tolower(type))
    if (length(.w) == 1) {
      .env <- x$env
      .objf <- .objDf[.w, "OBJF"]
      .lik <- .objDf[.w, "Log-likelihood"]
      attr(.lik, "df") <- attr(get("logLik", .env), "df")
      attr(.lik, "nobs") <- attr(get("logLik", .env), "nobs")
      class(.lik) <- "logLik"
      .bic <- .objDf[.w, "BIC"]
      .aic <- .objDf[.w, "AIC"]
      assign("OBJF", .objf, .env)
      assign("objf", .objf, .env)
      assign("objective", .objf, .env)
      assign("logLik", .lik, .env)
      assign("AIC", .aic, .env)
      assign("BIC", .bic, .env)
      if (!is.null(x$saem)) {
        .setSaemExtra(.env, type)
      }
      .env$ofvType <- type
      invisible(x)
    } else {
      if (any(tolower(type) == c("focei", "foce", "fo"))) {
        return(.setOfvFo(x, tolower(type)))
      } else if (!is.null(x$saem)) {
        .ret <- x$saem
        .reg <- rex::rex(start, "laplace", capture(.regNum), end)
        .regG <- rex::rex(start, "gauss", capture(.regNum), "_", capture(.regNum), end)
        if (regexpr(.reg, type, perl = TRUE) != -1) {
          .nnode <- 1
          .nsd <- as.numeric(sub(.reg, "\\1", type, perl = TRUE))
        } else if (regexpr(.regG, type, perl = TRUE) != -1) {
          .nnode <- as.numeric(sub(.regG, "\\1", type, perl = TRUE))
          .nsd <- as.numeric(sub(.regG, "\\2", type, perl = TRUE))
        } else {
          stop("cannot switch objective function to '", type, "' type", call. = FALSE)
        }
        .likTime <- proc.time()
        .saemObf <- calc.2LL(x$saem, nnodes.gq = .nnode, nsd.gq = .nsd)
        .likTime <- proc.time() - .likTime
        .likTime <- .likTime["elapsed"]
        .env <- x$env
        .time <- .env$time
        if (any(names(.time) == "logLik")) {
          .time$logLik <- .time$logLik + .likTime
        } else {
          .time <- data.frame(.time, logLik = .likTime, check.names = FALSE)
        }
        .env$time <- .time
        .llik <- -.saemObf / 2
        .nobs <- .env$nobs
        attr(.llik, "df") <- attr(get("logLik", .env), "df")
        .objf <- ifelse(.env$adjObj, .saemObf - .nobs * log(2 * pi), .saemObf)
        .tmp <- data.frame(
          OBJF = .objf, AIC = .saemObf + 2 * attr(get("logLik", .env), "df"),
          BIC = .saemObf + log(.env$nobs) * attr(get("logLik", .env), "df"),
          "Log-likelihood" = as.numeric(.llik), check.names = FALSE
        )
        if (any(names(.env$objDf) == "Condition Number")) {
          .cn <- unique(.env$objDf[["Condition Number"]])
          .cn <- .cn[!is.na(.cn)]
          if (length(.cn) == 1) {
            .tmp <- data.frame(.tmp, "Condition Number" = .cn, check.names = FALSE)
          } else {
            .tmp <- data.frame(.tmp, "Condition Number" = NA, check.names = FALSE)
          }
        }
        .rn <- row.names(.env$objDf)
        .env$objDf <- rbind(
          .env$objDf,
          .tmp
        )
        row.names(.env$objDf) <- c(.rn, type)
        .env$objDf <- .env$objDf[order(row.names(.env$objDf)), ]
        .env$objDf <- .env$objDf[!is.na(.env$objDf$OBJF), ]
        return(setOfv(x, type))
      }
      stop("cannot switch objective function to '", type, "' type",
        call. = FALSE
      )
    }
  } else {
    stop("wrong type of object",
      call. = FALSE
    )
  }
}

##' @rdname setOfv
##' @export
getOfvType <- function(x) {
  return(x$ofvType)
}
