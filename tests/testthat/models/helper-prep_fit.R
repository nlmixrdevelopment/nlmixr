# This header variable assignment simplifies testing and comparison across
# platforms
if (!exists("fit")) {
  fit <- list()
}
if (!exists("expected_values")){
    expected_values <- list()
}
if (!exists("verbose_minimization")) verbose_minimization <- FALSE
default_control <-
  nlme::nlmeControl(
    returnObject=TRUE,
    msMaxiter = 100,
    maxIter = 100,
    pnlsMaxIter = 100,
    msVerbose = verbose_minimization
    )


defaultControl <- function(x){
    if (x == "nlme"){
        return(nlme::nlmeControl(
                   returnObject=TRUE,
                   msMaxiter = 100,
                   maxIter = 100,
                   pnlsMaxIter = 100,
                   msVerbose = verbose_minimization))
    } else if (x == "saem") {
        return(saemControl())
    } else {
        return(foceiControl())
    }
}

defaultTable <- tableControl(cwres=TRUE); ## CWRES needed for FOCEi likelihood

generate_expected_values <- function(x=FALSE) {
  ret <-
    list(
      lik=round(c(logLik(fit[[runno]]), AIC(fit[[runno]]), BIC(fit[[runno]])), 2),
      param=unname(signif(fixef(fit[[runno]]), 5)),
      stdev_param=unname(signif(as.numeric(VarCorr(fit[[runno]])[1:length(fixef(fit[[runno]])), "StdDev"]), 5)),
      sigma=signif(fit[[runno]]$sigma, 5)
    )
  if (x){
      sink(paste0("values-", runno, "-", .Platform$OS.type, ".R"))
      on.exit(sink());
  }
  cat(
    "expected_values[[runno]] <-\n",
    "  list(\n",
    "    lik=c(", paste(ret$lik, collapse=", "), "),\n",
    "    param=c(", paste(ret$param, collapse=", "), "),\n",
    "    stdev_param=c(", paste(ret$stdev_param, collapse=", "), "),\n",
    "    sigma=c(", paste(ret$sigma, collapse=", "), ")\n",
    "  )\n",
    sep=""
  )
  if (!x){
      clip_con <- file(description="clipboard", open="w")
      cat(
          "expected_values[[runno]] <-\n",
          "  list(\n",
          "    lik=c(", paste(ret$lik, collapse=", "), "),\n",
          "    param=c(", paste(ret$param, collapse=", "), "),\n",
          "    stdev_param=c(", paste(ret$stdev_param, collapse=", "), "),\n",
          "    sigma=c(", paste(ret$sigma, collapse=", "), ")\n",
          "  )\n",
          file=clip_con,
          sep=""
      )
      close(con=clip_con)
  }
  invisible(ret)
}


genIfNeeded <- function(){
    .ret <- paste0("values-",runno, "-", .Platform$OS.type, ".R");
    if (!file.exists(.ret)){
        generate_expected_values(TRUE);
    }
    return(.ret)
}
