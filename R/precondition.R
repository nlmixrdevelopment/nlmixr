preconditionFit <- function(fit, estType = c("full", "posthoc", "none")) {
  pre <- preCondInv(fit$R)
  P <- symengine::Matrix(pre)
  d0 <- dimnames(fit$R)[[1]]
  d <- paste0("nlmixrPre_", dimnames(fit$R)[[1]])
  d2 <- symengine::Matrix(d)
  modExtra <- paste(paste0(d0, "=", sapply(as.vector(P %*% d2), as.character)), collapse = "\n")
  preInv <- solve(pre)

  .ini <- as.data.frame(fit$uif$ini)
  newEst <- setNames(as.vector(preInv %*% matrix(fit$theta[d0])), d0)
  for (v in d0) {
    .w <- which(.ini$name == v)
    .ini$lower[.w] <- -Inf
    .ini$upper[.w] <- Inf
    .ini$est[.w] <- newEst[v]
    .ini$name[.w] <- paste0("nlmixrPre_", v)
  }
  .newModel <- eval(parse(text = paste0("function(){", modExtra, "\n", fit$uif$fun.txt, "}"), keep.source = TRUE))
  class(ini) <- c("nlmixrBounds", "data.frame")
  newModel <- .finalizeUiModel(
    nlmixrUIModel(.newModel, .ini, NULL),
    as.list(fit$uif$meta)
  )

  .ctl <- fit$origControl
  estType <- match.arg(estType)
  if (estType == "none") {
    .ctl$maxInnerIterations <- 0
    .ctl$maxOuterIterations <- 0
    .ctl$boundTol <- 0
    .ctl$etaMat <- as.matrix(fit$eta[, -1])
    .ctl$calcTables <- FALSE
  } else if (estType == "posthoc") {
    .ctl$maxOuterIterations <- 0
    .ctl$boundTol <- 0
    .ctl$calcTables <- FALSE
  } else if (estType == "full") {
    .ctl$boundTol <- 0
    .ctl$calcTables <- FALSE
  }
  ## FIXME use vipul's get estimate routine
  newFit <- nlmixr(newModel, getData(fit), est = "focei", control = .ctl)
  cov <- pre %*% newFit$cov %*% t(pre)
  dimnames(cov) <- dimnames(pre)
  ## FIXME use back-transformation implemented by bill to update se to have full model definition
  return(cov)
}
