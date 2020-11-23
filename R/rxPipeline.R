##' @export
et.nlmixrFitData <- function(x, ..., envir = parent.frame()) {
  .si <- x$simInfo
  RxODE::.clearPipe(
    rx = RxODE::RxODE(.si$rx),
    ## events=nlmixrData(getData(x)),
    params = .si$params,
    thetaMat = .si$thetaMat,
    dfObs = .si$dfObs,
    omega = .si$omega,
    dfSub = .si$dfSub,
    sigma = .si$sigma
  )
  do.call(RxODE::et, c(list(...), list(envir = envir)), envir = envir)
}

##' @export
rxParams.nlmixrFitData <- function(obj, ...) {
  .si <- obj$simInfo
  RxODE::.clearPipe(
    rx = RxODE::RxODE(.si$rx),
    events = nlmixrData(getData(obj)),
    thetaMat = .si$thetaMat,
    dfObs = .si$dfObs,
    omega = .si$omega,
    dfSub = .si$dfSub,
    sigma = .si$sigma
  )
  do.call(RxODE::rxParams, list(...))
}
