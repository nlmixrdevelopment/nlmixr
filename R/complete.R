##' @importFrom utils .DollarNames
##' @export
.DollarNames.nlmixrBounds <- function(x, pattern) {
  grep(pattern,
    c(
      names(x),
      "theta", "theta.full", "omega", "random", "fixed.form",
      "focei.upper", "focei.lower", "theta.names", "focei.names",
      "focei.err.type", "eta.names"
    ),
    value = TRUE
  )
}
##' @export
.DollarNames.nlmixrUI <- function(x, pattern) {
  ## ui
  ## uiCompletions
  .cmp <- c(
    "ini", "nmodel", "model", "nlme.fun.mu", "dynmodel.fun", "nlme.fun",
    "nlme.fun.mu.cov", "nlme.specs", "nlme.specs.mu", "nlme.specs.mu.cov",
    "nlme.var", "rxode.pred", "theta.pars", "focei.inits", "focei.fixed",
    "focei.mu.ref", "saem.fixed", "saem.theta.name", "saem.eta.trans",
    "saem.model.omega", "saem.res.mod", "saem.ares", "saem.bres", "saem.log.eta",
    "saem.fit", "saem.model", "saem.init.theta", "saem.init.omega", "saem.init",
    "saem.omega.name", "saem.res.name", "model.desc", "meta", "saem.distribution",
    "lincmt.dvdx", "notfixed_bpop", "poped.notfixed_bpop",
    "poped.ff_fun", "poped.d", "poped.sigma", "logThetasList", "random.mu",
    "bpop", "multipleEndpoint", "muRefTable",
    ## bounds
    .DollarNames.nlmixrBounds(x$ini, ""),
    names(x$nmodel),
    ls(x$meta)
  )
  grep(pattern, .cmp, value = TRUE)
}
##' @export
.DollarNames.nlmixrFitCore <- function(x, pattern) {
  .env <- x$env
  .cmp <- c(
    names(x),
    "posthoc", "notes",
    "logLik", "value", "obf", "ofv", "objf", "OBJF", "objective", "AIC", "BIC",
    "value", "obf", "ofv", "objf", "sigma", "coefficients", "parHist", "par.hist",
    "parHistStacked", "par.hist.stacked", "omegaR", "omega.R", "par.fixed", "eta",
    "ranef", "theta", "fixef", "varFix", "thetaMat", "cov",
    "env", ls(.env)
  )
  if (exists("saem", .env)) {
    .cmp <- c(.cmp, "seed", "saem.cfg")
  }
  if (exists("uif", .env)) {
    .cmp <- c(
      .cmp, "model.name", "modelName", "dataName", "data.name",
      .DollarNames.nlmixrUI(.env$uif, "")
    )
  }
  .cmp <- c(.cmp, "env")
  grep(pattern, .cmp, value = TRUE)
}
