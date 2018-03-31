multi2 = function (mu, vmat, n) 
{
    eta <- matrix(rnorm(length(mu) * n), ncol = n, nrow = length(mu))
    Q <- chol(vmat, pivot = TRUE)
    pivot <- attr(Q, "pivot")
    oo <- order(pivot)
    para <- t(Q[, oo]) %*% eta
    sweep(para, 1, mu, "+")
}

rmvnorm = function(n, mu, vmat) multi2(mu, vmat, n) 

vpc.saemFit = function(fit, dat, nsim = 100, by=NULL, ...) {
  saem.cfg = attr(fit, "saem.cfg")
  dopred <- attr(fit, "dopred")
  red.mod = sum((fit$sig2 != 0) * 1:2)

  if (!is.null(by)) {
    if (by %in% names(dat)) {
      dat$grp = eval(parse(text=paste0("dat$",by)))
    } else {
      msg = paste0(by, " not found in data")
      stop(msg)
    }
  }
  else dat$grp = T

  xd = subset(dat, EVID==0)
  nsub = length(unique(xd$ID))
  ntim = dim(xd)[1]
  ord=rep(1:ntim, nsim)
  sim=rep(1:nsim, each=ntim)

  mpost_rand = t(matrix(fit$Plambda, length(fit$Plambda), nsub))
  s = lapply(1:nsim, function(k) {
    mpost_rand1 = t(rmvnorm(nsub, fit$Plambda[saem.cfg$i1+1], fit$Gamma2_phi1))
    mpost_rand[, saem.cfg$i1+1] = mpost_rand1
    p = dopred(mpost_rand, saem.cfg$evt, saem.cfg$opt)
    if      (red.mod==1) res = rnorm(ntim,0,sqrt(fit$sig2[1]))
    else if (red.mod==2) res = p*fit$sig2[2]*rnorm(ntim,0,1)
    else if (red.mod==3) res = (fit$sig2[1]+p*fit$sig2[2])*rnorm(ntim,0,1)
    p+res
  })
  xs = do.call("cbind",s)
  df = cbind(xd[ord, c("ID", "TIME", "grp")], DV=as.vector(xs), SIM=sim)

  if (!is.null(by)) {
      p = vpc::vpc(sim = df, obs = dat, strat=c("grp"), facet="wrap", ...)
  } else {
      p = vpc::vpc(sim = df, obs = dat, ...)
  }
  print(p)

  invisible(df)
}
