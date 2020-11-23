##' Get the optimal forward difference interval by Gill83 method
##'
##' @inheritParams base::do.call
##' @param which Which parameters to calculate the forward difference
##'     and optimal forward difference interval
##' @inheritParams foceiControl
##' @author Matthew Fidler
##'
##' @return
##'
##' A data frame with the following columns:
##' \itemize{
##' \item{info}{Gradient evaluation/forward difference information}
##' \item{hf}{Forward difference final estimate}
##' \item{df}{Derivative estimate}
##' \item{df2}{2nd Derivative Estimate}
##' \item{err}{Error of the final estimate derivative}
##' \item{aEps}{Absolute difference for forward numerical differences}
##' \item{rEps}{Relative Difference for backward numerical differences}
##' \item{aEpsC}{Absolute difference for central numerical differences}
##' \item{rEpsC}{Relative difference for central numerical differences}
##' }
##'
##' The \code{info} returns one of the following:
##' \itemize{
##' \item{Not Assessed}{Gradient wasn't assessed}
##' \item{Good}{Success in Estimating optimal forward difference interval}
##' \item{High Grad Error}{Large error; Derivative estimate error \code{fTol} or more of the derivative}
##' \item{Constant Grad}{Function constant or nearly constant for this parameter}
##' \item{Odd/Linear Grad}{Function odd or nearly linear, df = K, df2 ~ 0}
##' \item{Grad changes quickly}{df2 increases rapidly as h decreases}
##' }
##' @examples
##'
##' ## These are taken from the numDeriv::grad examples to show how
##' ## simple gradients are assessed with nlmixrGill83
##'
##' nlmixrGill83(sin, pi)
##'
##' nlmixrGill83(sin, (0:10)*2*pi/10)
##'
##' func0 <- function(x){ sum(sin(x))  }
##' nlmixrGill83(func0 , (0:10)*2*pi/10)
##'
##' func1 <- function(x){ sin(10*x) - exp(-x) }
##' curve(func1,from=0,to=5)
##'
##' x <- 2.04
##' numd1 <- nlmixrGill83(func1, x)
##' exact <- 10*cos(10*x) + exp(-x)
##' c(numd1$df, exact, (numd1$df - exact)/exact)
##'
##' x <- c(1:10)
##' numd1 <- nlmixrGill83(func1, x)
##' exact <- 10*cos(10*x) + exp(-x)
##' cbind(numd1=numd1$df, exact, err=(numd1$df - exact)/exact)
##'
##' sc2.f <- function(x){
##'   n <- length(x)
##'    sum((1:n) * (exp(x) - x)) / n
##' }
##'
##' sc2.g <- function(x){
##'   n <- length(x)
##'   (1:n) * (exp(x) - 1) / n
##' }
##'
##' x0 <- rnorm(100)
##' exact <- sc2.g(x0)
##'
##' g <- nlmixrGill83(sc2.f, x0)
##'
##' max(abs(exact - g$df)/(1 + abs(exact)))
##'
##' @export
nlmixrGill83 <- function(what, args, envir = parent.frame(),
                         which, gillRtol = sqrt(.Machine$double.eps), gillK = 10L, gillStep = 2, gillFtol = 0) {
  if (missing(which)) {
    which <- rep(TRUE, length(args))
  }
  return(nlmixrGill83_(what, args, envir, which,
    gillRtol = sqrt(.Machine$double.eps), gillK = 10L, gillStep = 2, gillFtol = 0
  ))
}

.nlmixrGradInfo <- new.env(parent = emptyenv())
##' Create a gradient function based on gill numerical differences
##'
##' @param thetaNames Names for the theta parameters
##' @inheritParams nlmixrGill83
##' @inheritParams foceiControl
##' @param theta for the internal functions theta is the parameter
##'     values
##' @param md5 the md5 identifier for the internal gradient function
##'     information.
##' @examples
##'
##' func0 <- function(x){ sum(sin(x))  }
##'
##' ## This will printout every interation or when print=X
##' gf <- nlmixrGradFun(func0)
##'
##' ## x
##' x <- (0:10)*2*pi/10;
##' gf$eval(x)
##' gf$grad(x)
##'
##' ## x2
##' x2 <- x+0.1
##' gf$eval(x2)
##' gf$grad(x2)
##'
##' ## Gives the parameter history as a data frame
##' gf$hist()
##'
##' @export
nlmixrGradFun <- function(what, envir = parent.frame(), which, thetaNames,
                          gillRtol = sqrt(.Machine$double.eps), gillK = 10L, gillStep = 2, gillFtol = 0,
                          useColor = crayon::has_color(),
                          printNcol = floor((getOption("width") - 23) / 12),
                          print = 1) {
  .md5 <- digest::digest(list(what, gillRtol, gillK, gillStep, gillFtol))
  .nlmixrGradInfo[["printNcol"]] <- printNcol
  .nlmixrGradInfo[["useColor"]] <- useColor
  .nlmixrGradInfo[["isRstudio"]] <- (Sys.getenv("RSTUDIO") == "1")
  .nlmixrGradInfo[["print"]] <- print
  if (!missing(which)) {
    .nlmixrGradInfo[[paste0(.md5, ".w")]] <- which
  }
  if (!missing(thetaNames)) {
    .nlmixrGradInfo[["thetaNames"]] <- thetaNames
  }
  .nlmixrGradInfo[[paste0(.md5, ".n")]] <- 0L
  .nlmixrGradInfo[[paste0(.md5, ".f")]] <- what
  .nlmixrGradInfo[[paste0(.md5, ".e")]] <- envir
  .nlmixrGradInfo[[paste0(.md5, ".rtol")]] <- gillRtol
  .nlmixrGradInfo[[paste0(.md5, ".k")]] <- gillK
  .nlmixrGradInfo[[paste0(.md5, ".s")]] <- gillStep
  .nlmixrGradInfo[[paste0(.md5, ".ftol")]] <- gillFtol
  .eval <- eval(parse(text = paste0("function(theta){
        nlmixrEval_(theta, \"", .md5, "\");
    }")))
  .grad <- eval(parse(text = paste0("function(theta){
        nlmixrGrad_(theta, \"", .md5, "\");
    }")))
  .hist <- eval(parse(text = paste0("function(){
        nlmixrParHist_(md5=\"", .md5, "\");
    }")))
  .unscaled <- eval(parse(text = paste0("function(theta){
        nlmixrUnscaled_(theta,md5=\"", .md5, "\");
    }")))
  return(list(eval = .eval, grad = .grad, hist = .hist, unscaled = .unscaled))
}

##' Calculate Hessian
##'
##' Unlike `stats::optimHess` which assumes the gradient is accurate,
##' nlmixrHess does not make as strong an assumption that the gradient
##' is accurate but takes more function evaluations to calculate the
##' Hessian.  In addition, this procedures optimizes the forward
##' difference interval by \code{\link{nlmixrGill83}}
##'
##' If you have an analytical gradient function, you should use
##' `stats::optimHess`
##'
##' @inheritParams stats::optimHess
##' @param ... Extra arguments sent to \code{\link{nlmixrGill83}}
##' @inheritParams base::do.call
##' @author Matthew Fidler
##' @export
##' @seealso \code{\link{nlmixrGill83}}, \code{\link{optimHess}}
##' @references
##'
##' \url{https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm}
##'
##' @examples
##'  func0 <- function(x){ sum(sin(x))  }
##'  x <- (0:10)*2*pi/10
##'  nlmixrHess(x, func0)
##'
##' fr <- function(x) {   ## Rosenbrock Banana function
##'     x1 <- x[1]
##'     x2 <- x[2]
##'     100 * (x2 - x1 * x1)^2 + (1 - x1)^2
##' }
##' grr <- function(x) { ## Gradient of 'fr'
##'     x1 <- x[1]
##'     x2 <- x[2]
##'     c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
##'        200 *      (x2 - x1 * x1))
##' }
##'
##' h1 <- optimHess(c(1.2,1.2), fr, grr)
##'
##' h2 <- optimHess(c(1.2,1.2), fr)
##'
##' ## in this case h3 is closer to h1 where the gradient is known
##'
##' h3 <- nlmixrHess(c(1.2,1.2), fr)
nlmixrHess <- function(par, fn, ..., envir = parent.frame()) {
  .gill <- nlmixrGill83(fn, par, envir = envir, ...)
  return(nlmixrHess_(par, fn, envir, .gill))
}
