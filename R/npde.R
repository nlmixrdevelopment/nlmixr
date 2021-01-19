
.npdeSim <- function(object, nsim = 300, ties = TRUE, seed = 1009, updateObject = TRUE,
                     cholSEtol = (.Machine$double.eps)^(1 / 3), ..., addDosing=FALSE, subsetNonmem=TRUE) {
  set.seed(seed)
  .si <- object$simInfo
  .rx <- .si$rx
  .rx <- gsub(rex::rex(capture("ipred"), or("=", "~"), capture(except_any_of("\n;")), any_of("\n;")), "ipred~\\2;\n", .rx)
  .rx <- gsub(rex::rex("d/dt(", capture(except_any_of("\n;)")), ")", or("=", "~")), "d/dt(\\1)~", .rx)
  .rx <- gsub(
    rex::rex("sim", or("=", "~"), "rxTBSi(", capture(except_any_of(",)")), ",", anything, any_of("\n;")),
    "sim=\\1", .rx
  )
  .si$rx <- .rx
  .si$object <- object
  .si$returnType <- "data.frame.TBS"
  .si$nsim <- nsim
  .si <- c(.si, list(...))
  .si$modelName <- "NPDE"
  .pt <- proc.time()
  .si$dfObs <- 0
  .si$dfSub <- 0
  .si$thetaMat <- NA
  .si$addDosing <- addDosing
  .si$subsetNonmem <- subsetNonmem
  do.call("nlmixrSim", .si)
}
