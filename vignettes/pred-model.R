library(nlmixr)
one.cmt <- function() {
     ini({
         tka <- 0.47   # Log Ka
         tcl <- -3.22 # Log Cl
         tke <- -2.45    # Log Ke
         eta.ka ~ 1
         eta.cl ~ 2
         add.err <- c(0, 0.1,1)
     })
     model({
         ka <- exp(tka + eta.ka)
         cl <- exp(tcl + eta.cl)
         ke <- exp(tke)
         cp = dose * ke * ka * (exp(-ke * time)-exp(-ka * time))/cl/(ka-ke)
         cp ~ add(add.err)
     })
 }


fit <- nlmixr(one.cmt, theo_sd, est="nlme")
