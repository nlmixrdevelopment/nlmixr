.amtTrans <- function(data) {
  if (!any(names(data) == "amt")) {
    stop("need 'amt' aesthetic")
  }
  if (!any(names(data) == "x")) {
    if (any(names(data) == "time")) {
      data$x <-data$time
    } else {
      stop("need 'x' aesthetic")
    }
  }
  #.range <- range(data$y,na.rm=TRUE)
  .dat <- data
  .dat <- data[,c("x","amt")]
  .dat <- data[!is.na(.dat$amt),]
  ret <- data.frame(x=.dat$x,xend=.dat$x,y=-Inf,yend=Inf,amt=.dat$amt)
  return(ret)
}

GeomAmt <- ggplot2::ggproto("GeomAmt", ggplot2::GeomSegment,
                            required_aes = c("x", "y", "xend", "yend"),
                            default_aes = ggplot2::aes(colour="black", linetype="dotted",size=0.5, alpha = 1, fill="black"))

StatAmt <- ggplot2::ggproto("StatAmt", ggplot2::Stat,
                            compute_group = function(data, scales) {
                              .amtTrans(data)
                            },
                            required_aes=c("x","amt"))


##' Dosing/Amt geom/stat
##'
##' This is a dosing geom that shows the vertical lines where a dose occurs
##'
##' Requires the following aesthetics:
##'
##' \itemize{
##' \item x representing the x values, usually time
##' \item amt representing the dosing values;  They are missing or zero when no dose is given
##' }
##'
##' @export
##' @inheritParams ggplot2::stat_identity
stat_amt <- function(mapping = NULL, data = NULL,
                     position = "identity", show.legend = NA,
                     inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatAmt, data = data, mapping = mapping, geom = GeomAmt,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = TRUE,  ...)
  )
}

##' @rdname stat_amt
##' @export
geom_amt <- stat_amt
