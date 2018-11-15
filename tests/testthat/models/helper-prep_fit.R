# This header variable assignment simplifies testing and comparison across
# platforms
if (!exists("fit")) {
  fit <- list()
}
if (!exists("verbose_minimization")) verbose_minimization <- FALSE
default_control <-
  nlmeControl(
    returnObject=TRUE,
    msMaxiter = 100,
    maxIter = 100,
    pnlsMaxIter = 100,
    msVerbose = verbose_minimization
  )
  # nlmeControl(
  #   returnObject=TRUE,
  #   msMaxIter=1L,
  #   maxIter=1L,
  #   pnlsMaxIter=1L
  # )

generate_expected_values <- function(x) {
  ret <-
    list(
      lik=round(c(logLik(fit[[runno]]), AIC(fit[[runno]]), BIC(fit[[runno]])), 2),
      param=unname(signif(fixef(fit[[runno]]), 5)),
      stdev_param=unname(signif(as.numeric(VarCorr(fit[[runno]])[1:length(fixef(fit[[runno]])), "StdDev"]), 5)),
      sigma=signif(fit[[runno]]$sigma, 5)
    )
  cat(
    "expected_values <-\n",
    "  list(\n",
    "    lik=c(", paste(ret$lik, collapse=", "), "),\n",
    "    param=c(", paste(ret$param, collapse=", "), "),\n",
    "    stdev_param=c(", paste(ret$stdev_param, collapse=", "), "),\n",
    "    sigma=c(", paste(ret$sigma, collapse=", "), ")\n",
    "  )\n",
    sep=""
  )
  invisible(ret)
}
