## FIXME: g by endpint

##' Log-likelihood using Gaussian Quadrature
##'
##' Estimate the log-likelihood using Gaussian Quadrature (multidimensional
##' grid)
##'
##'
##' @param fit saemFit fit
##'
##' @inheritParams saemControl
##'
##' @return log-likelihood calculated by Gaussian Quadrature
##'
##' @references  Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear
##' mixed effects models. Computational Statistics and Data Analysis 49, 4
##' (2005), 1020-1038.
##'
##' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
##' 20th meeting of the Population Approach Group in Europe, Athens, Greece
##' (2011), Abstr 2173.
##'
##' @export
calc.2LL <- function(fit, nnodes.gq = 8, nsd.gq = 4) {
  ## nnodes.gq=8, nsd.gq=4
  saem.cfg <- attr(fit, "saem.cfg")
  .evtM <- saem.cfg$evtM
  dopred <- attr(fit, "dopred")
  resMat <- fit$resMat
  transMat <- fit$transMat
  ares <- resMat[, 1]
  bres <- resMat[, 2]
  lambda <- transMat[, 1]
  yj <- transMat[, 2]
  low <- transMat[, 3]
  hi <- transMat[, 4]
  i1 <- saem.cfg$i1 + 1
  nphi1 <- saem.cfg$nphi1
  nphi0 <- saem.cfg$nphi0
  nphi <- nphi0 + nphi1
  N <- saem.cfg$N
  ntotal <- saem.cfg$ntotal
  ix_endpnt <- saem.cfg$ix_endpnt[1:ntotal] + 1
  ares <- ares[ix_endpnt]
  bres <- bres[ix_endpnt]
  lambda <- lambda[ix_endpnt]
  yj <- yj[ix_endpnt]
  low <- low[ix_endpnt]
  hi <- hi[ix_endpnt]
  i1 <- saem.cfg$i1 + 1

  phi <- fit$mpost_phi
  IOmega.phi1 <- solve(fit$Gamma2_phi1)
  Omega <- fit$Gamma2_phi1
  yobs <- yobss <- saem.cfg$y
  evt <- saem.cfg$evt
  id <- evt[evt[, "EVID"] == 0, "ID"]
  nb_measures <- table(id)
  mlen <- max(nb_measures)
  io <- t(sapply(nb_measures, function(x) rep(1:0, c(x, mlen - x))))
  ind.io <- grep(1, t(io))
  DYF <- matrix(0, mlen, N)

  phiM <- matrix(scan(saem.cfg$phiMFile, quiet = TRUE), byrow = TRUE, ncol = nphi)
  dim(phiM) <- c(N, saem.cfg$nmc, saem.cfg$niter, nphi)
  ## print(head(phiM))
  ## print(dim(phiM))
  cond.mean.phi <- apply(phiM, c(1, 4), mean)
  var.all <- lapply(1:N, function(k) {
    x <- phiM[k, , , ]
    dim(x) <- c(saem.cfg$nmc * saem.cfg$niter, nphi)
    var(x)
  })
  condsd.eta <- t(sapply(var.all, function(x) sqrt(diag(x))))

  y <- gqg.mlx(nphi1, nnodes.gq)
  # xform [0,1] => [-1, 1]
  x <- (y$nodes - 0.5) * 2
  w <- (y$weights) * (2^nphi1)
  nx <- dim(x)[1]
  ## nx.1 = max(as.integer(nx/10), 1)
  xmin <- cond.mean.phi[, i1] - nsd.gq * condsd.eta[, i1]
  xmax <- cond.mean.phi[, i1] + nsd.gq * condsd.eta[, i1]
  a <- (xmin + xmax) / 2
  dim(a) <- c(N, nphi1)
  b <- (xmax - xmin) / 2
  dim(b) <- c(N, nphi1)

  Q <- 0
  if (nnodes.gq == 1) {
    message(sprintf("Calculating Laplace -2LL (nsd=%s)", nsd.gq))
  } else {
    message(sprintf("Calculating -2LL by Gaussian quadrature (nnodes=%s,nsd=%s)", nnodes.gq, nsd.gq))
  }
  RxODE::rxProgress(nx)
  ysave <- yobs
  yobs <- .Call(`_nlmixr_powerD`, yobs, lambda, as.integer(yj), as.double(low), as.double(hi))
  on.exit(RxODE::rxProgressAbort("Error calculating likelihood"))
  for (j in 1:nx) {
    phi[, i1] <- a + b * matrix(rep(x[j, ], N), ncol = nphi1, byrow = TRUE)
    f <- fsave <- as.vector(dopred(phi, saem.cfg$evt, saem.cfg$opt))
    f <- .Call(`_nlmixr_powerD`, f, lambda, as.integer(yj), as.double(low), as.double(hi))
    g <- ares + bres * abs(fsave)
    g[g < 1.0e-200] <- 1.0e-200
    DYF[ind.io] <- -0.5 * ((yobs - f) / g)^2 - log(g)
    ly <- colSums(DYF)
    dphi1 <- phi[, i1] - fit$mprior_phi[, i1]
    lphi1 <- -0.5 * rowSums((dphi1 %*% IOmega.phi1) * dphi1)
    ltot <- ly + lphi1
    ltot[is.na(ltot)] <- -Inf
    Q <- Q + w[j] * exp(ltot)
    RxODE::rxTick()
  }
  RxODE::rxProgressStop()
  # - 2 * saem.cfg$extraLL
  ll2 <- 2 * sum(log(Q) + rowSums(log(b))) - N * log(det(Omega)) - (N * nphi1 + ntotal) * log(2 * pi) -
    2 * .Call(`_nlmixr_powerL`, ysave, lambda, as.integer(yj), as.double(low), as.double(hi))
  -ll2
}

#' Plot an SAEM model fit
#'
#' Plot an SAEM model fit
#'
#' @param x a saemFit object
#' @param ... others
#' @return a list
#' @export
plot.saemFit <- function(x, ...) {
  CMT <- RES <- NULL
  fit <- x
  saem.cfg <- attr(fit, "saem.cfg")
  ## .env$model$assignPtr()
  .evtM <- saem.cfg$evtM
  dat <- .as.data.frame(saem.cfg$evt)
  dat <- cbind(dat[dat$EVID == 0, ], DV = saem.cfg$y)
  df <- rbind(cbind(dat, grp = 1), cbind(dat, grp = 2), cbind(dat, grp = 3))
  dopred <- attr(fit, "dopred")
  yp <- dopred(fit$mprior_phi, saem.cfg$evt, saem.cfg$opt)
  yi <- dopred(fit$mpost_phi, saem.cfg$evt, saem.cfg$opt)
  df$DV[df$grp == 2] <- yp
  df$DV[df$grp == 3] <- yi
  df0 <- df

  m <- fit$par_hist
  df <- .data.frame(val = as.vector(m), par = rep(1:ncol(m),
    each = nrow(m)
  ), iter = rep(1:nrow(m), ncol(m)))
  p1 <- ggplot2::ggplot(df, aes(iter, val)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~par, scales = "free_y")
  print(p1)

  for (cmt in sort(unique(df0$CMT))) {
    df <- subset(df0, CMT == cmt)
    p6 <- ggplot(subset(df, grp == 1), aes(TIME, DV)) +
      geom_point() +
      facet_wrap(~ID) +
      geom_line(aes(TIME, DV), subset(df, grp == 2), col = "blue") +
      geom_line(aes(TIME, DV), subset(df, grp == 3), col = "red")

    df <- cbind(subset(df, grp == 1), PRED = subset(df, grp == 2)[, "DV"])
    df$RES <- df$DV - df$PRED
    p2 <- ggplot(df, aes(PRED, DV)) +
      geom_point() +
      geom_abline(
        intercept = 0,
        slope = 1, col = "red"
      )
    p3 <- ggplot(df, aes(PRED, RES)) +
      geom_point() +
      geom_abline(
        intercept = 0,
        slope = 0, col = "red"
      )

    df <- subset(df0, CMT == cmt)
    df <- cbind(subset(df, grp == 1), IPRED = subset(df, grp == 3)[, "DV"])
    df$IRES <- df$DV - df$IPRED
    p4 <- ggplot(df, aes(IPRED, DV)) +
      geom_point() +
      geom_abline(
        intercept = 0,
        slope = 1, col = "red"
      )
    p5 <- ggplot(df, aes(IPRED, IRES)) +
      geom_point() +
      geom_abline(
        intercept = 0,
        slope = 0, col = "red"
      )

    print(p2)
    print(p4)
    print(p3)
    print(p5)
    print(p6)
  }
  invisible(NULL)
}

# this fn is from saemix
# import it once saem is back on CRAN
gqg.mlx <- function(dim, nnodes.gq) {
  # GQG.MLX Nodes and weights for numerical integration on grids
  # (multidimensional Gaussian Quadrature)
  #    dim  : dimension of the integration problem
  #    nnodes.gq   : number of points on any 1-D grid
  #
  #    x    = matrix of nodes with dim columns
  #    w    = row vector of corresponding weights
  #
  if (nnodes.gq > 25) {
    cat("The number of nodes for Gaussian Quadrature should be less than 25.\n")
    return(list(nodes = NULL, weights = c()))
  }
  if (nnodes.gq == 1) {
    n <- c(5.0000000000000000e-001)
    w <- c(1.0000000000000000e+000)
  }
  if (nnodes.gq == 2) {
    n <- c(7.8867513459481287e-001)
    w <- c(5.0000000000000000e-001)
  }
  if (nnodes.gq == 3) {
    n <- c(5.0000000000000000e-001, 8.8729833462074170e-001)
    w <- c(4.4444444444444570e-001, 2.7777777777777712e-001)
  }
  if (nnodes.gq == 4) {
    n <- c(6.6999052179242813e-001, 9.3056815579702623e-001)
    w <- c(3.2607257743127516e-001, 1.7392742256872484e-001)
  }
  if (nnodes.gq == 5) {
    n <- c(5.0000000000000000e-001, 7.6923465505284150e-001, 9.5308992296933193e-001)
    w <- c(2.8444444444444655e-001, 2.3931433524968501e-001, 1.1846344252809174e-001)
  }
  if (nnodes.gq == 6) {
    n <- c(6.1930959304159849e-001, 8.3060469323313235e-001, 9.6623475710157603e-001)
    w <- c(2.3395696728634746e-001, 1.8038078652407072e-001, 8.5662246189581834e-002)
  }
  if (nnodes.gq == 7) {
    n <- c(5.0000000000000000e-001, 7.0292257568869854e-001, 8.7076559279969723e-001, 9.7455395617137919e-001)
    w <- c(2.0897959183673620e-001, 1.9091502525256090e-001, 1.3985269574463935e-001, 6.4742483084431701e-002)
  }
  if (nnodes.gq == 8) {
    n <- c(5.9171732124782495e-001, 7.6276620495816450e-001, 8.9833323870681348e-001, 9.8014492824876809e-001)
    w <- c(1.8134189168918213e-001, 1.5685332293894469e-001, 1.1119051722668793e-001, 5.0614268145185180e-002)
  }
  if (nnodes.gq == 9) {
    n <- c(5.0000000000000000e-001, 6.6212671170190451e-001, 8.0668571635029518e-001, 9.1801555366331788e-001, 9.8408011975381304e-001)
    w <- c(1.6511967750063075e-001, 1.5617353852000226e-001, 1.3030534820146844e-001, 9.0324080347429253e-002, 4.0637194180784583e-002)
  }
  if (nnodes.gq == 10) {
    n <- c(5.7443716949081558e-001, 7.1669769706462361e-001, 8.3970478414951222e-001, 9.3253168334449232e-001, 9.8695326425858587e-001)
    w <- c(1.4776211235737713e-001, 1.3463335965499873e-001, 1.0954318125799158e-001, 7.4725674575290599e-002, 3.3335672154342001e-002)
  }
  if (nnodes.gq == 11) {
    n <- c(5.0000000000000000e-001, 6.3477157797617245e-001, 7.5954806460340585e-001, 8.6507600278702468e-001, 9.4353129988404771e-001, 9.8911432907302843e-001)
    w <- c(1.3646254338895086e-001, 1.3140227225512388e-001, 1.1659688229599563e-001, 9.3145105463867520e-002, 6.2790184732452625e-002, 2.7834283558084916e-002)
  }
  if (nnodes.gq == 12) {
    n <- c(5.6261670425573451e-001, 6.8391574949909006e-001, 7.9365897714330869e-001, 8.8495133709715235e-001, 9.5205862818523745e-001, 9.9078031712335957e-001)
    w <- c(1.2457352290670189e-001, 1.1674626826917781e-001, 1.0158371336153328e-001, 8.0039164271673444e-002, 5.3469662997659276e-002, 2.3587668193254314e-002)
  }
  if (nnodes.gq == 13) {
    n <- c(5.0000000000000000e-001, 6.1522915797756739e-001, 7.2424637551822335e-001, 8.2117466972017006e-001, 9.0078904536665494e-001, 9.5879919961148907e-001, 9.9209152735929407e-001)
    w <- c(1.1627577661543741e-001, 1.1314159013144903e-001, 1.0390802376844462e-001, 8.9072990380973202e-002, 6.9436755109893875e-002, 4.6060749918864378e-002, 2.0242002382656228e-002)
  }
  if (nnodes.gq == 14) {
    n <- c(5.5402747435367183e-001, 6.5955618446394482e-001, 7.5762431817907705e-001, 8.4364645240584268e-001, 9.1360065753488251e-001, 9.6421744183178681e-001, 9.9314190434840621e-001)
    w <- c(1.0763192673157916e-001, 1.0259923186064811e-001, 9.2769198738969161e-002, 7.8601583579096995e-002, 6.0759285343951711e-002, 4.0079043579880291e-002, 1.7559730165874574e-002)
  }
  if (nnodes.gq == 15) {
    n <- c(5.0000000000000000e-001, 6.0059704699871730e-001, 6.9707567353878175e-001, 7.8548608630426942e-001, 8.6220886568008503e-001, 9.2410329170521366e-001, 9.6863669620035298e-001, 9.9399625901024269e-001)
    w <- c(1.0128912096278091e-001, 9.9215742663556039e-002, 9.3080500007781286e-002, 8.3134602908497196e-002, 6.9785338963077315e-002, 5.3579610233586157e-002, 3.5183023744054159e-002, 1.5376620998057434e-002)
  }
  if (nnodes.gq == 16) {
    n <- c(5.4750625491881877e-001, 6.4080177538962946e-001, 7.2900838882861363e-001, 8.0893812220132189e-001, 8.7770220417750155e-001, 9.3281560119391593e-001, 9.7228751153661630e-001, 9.9470046749582497e-001)
    w <- c(9.4725305227534431e-002, 9.1301707522462000e-002, 8.4578259697501462e-002, 7.4797994408288562e-002, 6.2314485627767105e-002, 4.7579255841246545e-002, 3.1126761969323954e-002, 1.3576229705875955e-002)
  }
  if (nnodes.gq == 17) {
    n <- c(5.0000000000000000e-001, 5.8924209074792389e-001, 6.7561588172693821e-001, 7.5634526854323847e-001, 8.2883557960834531e-001, 8.9075700194840068e-001, 9.4011957686349290e-001, 9.7533776088438384e-001, 9.9528773765720868e-001)
    w <- c(8.9723235178103419e-002, 8.8281352683496447e-002, 8.4002051078225143e-002, 7.7022880538405308e-002, 6.7568184234262890e-002, 5.5941923596702053e-002, 4.2518074158589644e-002, 2.7729764686993612e-002, 1.2074151434273140e-002)
  }
  if (nnodes.gq == 18) {
    n <- c(5.4238750652086765e-001, 6.2594311284575277e-001, 7.0587558073142131e-001, 7.7988541553697377e-001, 8.4584352153017661e-001, 9.0185247948626157e-001, 9.4630123324877791e-001, 9.7791197478569880e-001, 9.9578258421046550e-001)
    w <- c(8.4571191481571939e-002, 8.2138241872916504e-002, 7.7342337563132801e-002, 7.0321457335325452e-002, 6.1277603355739306e-002, 5.0471022053143716e-002, 3.8212865127444665e-002, 2.4857274447484968e-002, 1.0808006763240719e-002)
  }
  if (nnodes.gq == 19) {
    n <- c(5.0000000000000000e-001, 5.8017932282011264e-001, 6.5828204998181494e-001, 7.3228537068798050e-001, 8.0027265233084055e-001, 8.6048308866761469e-001, 9.1135732826857141e-001, 9.5157795180740901e-001, 9.8010407606741501e-001, 9.9620342192179212e-001)
    w <- c(8.0527224924391946e-002, 7.9484421696977337e-002, 7.6383021032929960e-002, 7.1303351086803413e-002, 6.4376981269668232e-002, 5.5783322773667113e-002, 4.5745010811225124e-002, 3.4522271368820669e-002, 2.2407113382849821e-002, 9.7308941148624341e-003)
  }
  if (nnodes.gq == 20) {
    n <- c(5.3826326056674867e-001, 6.1389292557082253e-001, 6.8685304435770977e-001, 7.5543350097541362e-001, 8.1802684036325757e-001, 8.7316595323007540e-001, 9.1955848591110945e-001, 9.5611721412566297e-001, 9.8198596363895696e-001, 9.9656429959254744e-001)
    w <- c(7.6376693565363113e-002, 7.4586493236301996e-002, 7.1048054659191187e-002, 6.5844319224588346e-002, 5.9097265980759248e-002, 5.0965059908620318e-002, 4.1638370788352433e-002, 3.1336024167054569e-002, 2.0300714900193556e-002, 8.8070035695753026e-003)
  }
  if (nnodes.gq == 21) {
    n <- c(5.0000000000000000e-001, 5.7278092708044759e-001, 6.4401065840120053e-001, 7.1217106010371944e-001, 7.7580941794360991e-001, 8.3356940209870611e-001, 8.8421998173783889e-001, 9.2668168229165859e-001, 9.6004966707520034e-001, 9.8361341928315316e-001, 9.9687608531019478e-001)
    w <- c(7.3040566824845346e-002, 7.2262201994985134e-002, 6.9943697395536658e-002, 6.6134469316668845e-002, 6.0915708026864350e-002, 5.4398649583574356e-002, 4.6722211728016994e-002, 3.8050056814189707e-002, 2.8567212713428641e-002, 1.8476894885426285e-002, 8.0086141288864491e-003)
  }
  if (nnodes.gq == 22) {
    n <- c(5.3486963665986109e-001, 6.0393021334411068e-001, 6.7096791044604209e-001, 7.3467791899337853e-001, 7.9382020175345580e-001, 8.4724363159334137e-001, 8.9390840298960406e-001, 9.3290628886015003e-001, 9.6347838609358694e-001, 9.8503024891771429e-001, 9.9714729274119962e-001)
    w <- c(6.9625936427816129e-002, 6.8270749173007697e-002, 6.5586752393531317e-002, 6.1626188405256251e-002, 5.6466148040269712e-002, 5.0207072221440600e-002, 4.2970803108533975e-002, 3.4898234212260300e-002, 2.6146667576341692e-002, 1.6887450792407110e-002, 7.3139976491353280e-003)
  }
  if (nnodes.gq == 23) {
    n <- c(5.0000000000000000e-001, 5.6662841214923310e-001, 6.3206784048517251e-001, 6.9515051901514546e-001, 7.5475073892300371e-001, 8.0980493788182306e-001, 8.5933068156597514e-001, 9.0244420080942001e-001, 9.3837617913522076e-001, 9.6648554341300807e-001, 9.8627123560905761e-001, 9.9738466749877608e-001)
    w <- c(6.6827286093053176e-002, 6.6231019702348404e-002, 6.4452861094041150e-002, 6.1524542153364815e-002, 5.7498320111205814e-002, 5.2446045732270824e-002, 4.6457883030017563e-002, 3.9640705888359551e-002, 3.2116210704262994e-002, 2.4018835865542369e-002, 1.5494002928489686e-002, 6.7059297435702412e-003)
  }
  if (nnodes.gq == 24) {
    n <- c(5.3202844643130276e-001, 5.9555943373680820e-001, 6.5752133984808170e-001, 7.1689675381302254e-001, 7.7271073569441984e-001, 8.2404682596848777e-001, 8.7006209578927718e-001, 9.1000099298695147e-001, 9.4320776350220048e-001, 9.6913727600136634e-001, 9.8736427798565474e-001, 9.9759360999851066e-001)
    w <- c(6.3969097673376246e-002, 6.2918728173414318e-002, 6.0835236463901793e-002, 5.7752834026862883e-002, 5.3722135057982914e-002, 4.8809326052057039e-002, 4.3095080765976693e-002, 3.6673240705540205e-002, 2.9649292457718385e-002, 2.2138719408709880e-002, 1.4265694314466934e-002, 6.1706148999928351e-003)
  }
  if (nnodes.gq == 25) {
    n <- c(5.0000000000000000e-001, 5.6143234630535521e-001, 6.2193344186049426e-001, 6.8058615290469393e-001, 7.3650136572285752e-001, 7.8883146512061142e-001, 8.3678318423673415e-001, 8.7962963151867890e-001, 9.1672131438041693e-001, 9.4749599893913761e-001, 9.7148728561448716e-001, 9.8833196072975871e-001, 9.9777848489524912e-001)
    w <- c(6.1588026863357799e-002, 6.1121221495155122e-002, 5.9727881767892461e-002, 5.7429129572855862e-002, 5.4259812237131867e-002, 5.0267974533525363e-002, 4.5514130991481903e-002, 4.0070350167500532e-002, 3.4019166906178545e-002, 2.7452347987917691e-002, 2.0469578350653148e-002, 1.3177493307516108e-002, 5.6968992505125535e-003)
  }
  n1 <- 1 - n
  if (nnodes.gq %% 2 == 0) {
    x <- c(rev(n1), n)
    w <- c(rev(w), w)
  } else {
    x <- c(rev(n1[-1]), n)
    w <- c(rev(w[-1]), w)
  }
  mw <- nodes <- matrix(0, nrow = nnodes.gq**dim, ncol = dim)
  for (j in 1:dim) {
    nodes[, j] <- rep(rep(x, each = nnodes.gq**(dim - j)), nnodes.gq**(j - 1))
    mw[, j] <- rep(rep(w, each = nnodes.gq**(dim - j)), nnodes.gq**(j - 1))
  }
  weights <- apply(mw, 1, prod)
  return(list(nodes = nodes, weights = weights))
}

cutoff <- function(x, cut = .Machine$double.xmin) {
  x[x < cut] <- cut
  x
}

##' Covariance matrix by Fisher Information Matrix via linearization
##'
##' Get the covariance matrix of fixed effect estimates via calculating Fisher Information Matrix by linearization
##'
##' @param fit0 saemFit fit
##'
##' @return standard error of fixed effects
##'
##' @references Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
##' 20th meeting of the Population Approach Group in Europe, Athens, Greece
##' (2011), Abstr 2173.
##'
calc.COV <- function(fit0) {
  message("Calculating covariance matrix")
  if (is(fit0, "saemFit")) {
    fit <- fit0
  } else {
    fit <- as.saem(fit0)
  }
  .env <- attr(fit, "env")
  saem.cfg <- attr(fit, "saem.cfg")
  .evtM <- saem.cfg$evtM
  .rx <- .env$model
  .pars <- .rx$params
  .pars <- setNames(rep(1.1, length(.pars)), .pars)
  dopred <- attr(fit, "dopred")
  resMat <- fit$resMat
  ares <- resMat[, 1]
  bres <- resMat[, 2]
  transMat <- fit$transMat
  ares <- resMat[, 1]
  bres <- resMat[, 2]
  lambda <- transMat[, 1]
  yj <- transMat[, 2]
  low <- transMat[, 3]
  hi <- transMat[, 4]
  i1 <- saem.cfg$i1 + 1
  nphi1 <- saem.cfg$nphi1
  nphi0 <- saem.cfg$nphi0
  nphi <- nphi0 + nphi1
  N <- saem.cfg$N
  RxODE::rxProgress(N + nphi)
  on.exit(RxODE::rxProgressAbort("Error calculating covariance via linearization"))

  ntotal <- saem.cfg$ntotal
  ix_endpnt <- saem.cfg$ix_endpnt[1:ntotal] + 1
  ares <- ares[ix_endpnt]
  bres <- bres[ix_endpnt]
  lambda <- lambda[ix_endpnt]
  yj <- yj[ix_endpnt]
  low <- low[ix_endpnt]
  hi <- hi[ix_endpnt]
  if (is.null(names(saem.cfg$inits$theta))) {
    names(saem.cfg$inits$theta) <- rep("", length(saem.cfg$inits$theta))
  }
  s <- matrix(names(saem.cfg$inits$theta), ncol = nphi, byrow = TRUE)
  fix.ix <- c(s) == "FIXED"
  s <- matrix(saem.cfg$inits$theta, ncol = nphi, byrow = TRUE)
  cov.ix <- !is.na(c(s))
  cov.est.ix <- !fix.ix & cov.ix

  phi <- hat.phi <- fit$mpost_phi
  omega <- fit$Gamma2_phi1
  evt <- saem.cfg$evt
  id <- evt[evt[, "EVID"] == 0, "ID"] + 1

  dphi <- cutoff(abs(phi) * 1e-4, 1e-10)
  f1 <- sapply(1:nphi, function(j) {
    phi[, j] <- hat.phi[, j] + dphi[, j]
    ret <- .Call(
      `_nlmixr_powerD`, as.vector(dopred(phi, saem.cfg$evt, saem.cfg$opt)),
      lambda, as.integer(yj), as.double(low), as.double(hi)
    )
    RxODE::rxTick()
    return(ret)
  })
  f0 <- f0s <- as.vector(dopred(hat.phi, saem.cfg$evt, saem.cfg$opt))
  f0 <- .Call(`_nlmixr_powerD`, f0, lambda, as.integer(yj), as.double(low), as.double(hi))
  DF <- (f1 - f0) / dphi[id, ]
  g <- ares + bres * abs(f0s)

  # spectral decom for invVi, idea from saemix
  Xi <- lapply(1:N, function(i) {
    ix <- id == i
    nobs <- sum(ix)
    DFi <- DF[ix, ]
    dim(DFi) <- c(nobs, nphi)
    DFi.i1 <- DFi[, i1]
    dim(DFi.i1) <- c(nobs, nphi1)
    m <- diag(nobs)
    diag(m) <- g[ix]^2
    Vi <- DFi.i1 %*% omega %*% t(DFi.i1) + m
    VD <- try(eigen(Vi, symmetric = TRUE))
    if (inherits(VD, "try-error")) stop("Spectral decom failure when computing FIM")
    D <- Re(VD$values)
    V <- Re(VD$vectors)
    m <- diag(length(D))
    diag(m) <- 1 / sqrt(D)
    invVi.5 <- m %*% t(V) # backsolve(chol(Vi), diag(11)); chol() is worse
    Ai <- kronecker(diag(nphi), saem.cfg$Mcovariables[i, ])
    DFAi <- DFi %*% t(Ai[cov.est.ix, ]) # CHECK!
    ret <- invVi.5 %*% DFAi
    RxODE::rxTick()
    return(ret)
  })
  npar <- sum(cov.est.ix)
  X <- do.call("rbind", Xi)
  Ri <- backsolve(qr.R(qr(X)), diag(npar))
  ret <- crossprod(t(Ri))
  RxODE::rxProgressStop()
  return(ret)
}
