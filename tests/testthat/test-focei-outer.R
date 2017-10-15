rxPermissive({
    context("Outer Problem Gradient Tests")

    mat.indices <- function(nETA){
        idx = do.call("rbind",
                      lapply(1:nETA, function(k) cbind(k:nETA, k)))
        H = matrix(1:(nETA^2), nETA, nETA)
        Hlo.idx = row(H)>=col(H)
        lo.idx = H[row(H)>col(H)]
        hi.idx = t(H)[row(H)>col(H)]

        list(idx=idx,                       # (r, c) of lo-half
             Hlo.idx=Hlo.idx,       # index of lo-half
             lo.idx=lo.idx,         # index of strict lo-half
             hi.idx=hi.idx)         # index of strict hi-half
    }


    ev <- eventTable() %>%
        add.sampling(c(95.99, 119.99, 143.99, 144.25, 144.5, 144.75,
                       145, 145.5, 146, 146.5, 147, 148, 150, 152, 156, 160, 164, 167.99,
                       191.99, 215.99, 216.25, 216.5, 216.75, 217, 217.5, 218, 218.5, 219,
                       220, 222, 224, 228, 232, 236, 240, 252, 264, 276, 288)) %>%
        add.dosing(dose=60000, start.time=72, nbr.doses=7, dosing.interval=24)

    dv <- c(263.6, 164.7, 287.3, 1248.7, 1211.5, 1017.7, 1690.1, 1029.8,
            890.7, 598.4, 1009.3, 1159.8, 742.2, 724.6, 728.2, 509.7, 243.1,
            259.9, 242.2, 281.4, 1500.1, 1281.4, 1200.2, 1378.8, 1373.2,
            582.9, 960.2, 720.3, 852.6, 950.3, 654.7, 402.5, 456, 346.5,
            268.2, 134.2, 42.6, 25.9, 14.6)


    m1 <- RxODE({
        C2 = centr/V;
        d/dt(centr) = - CL*C2
    })


    pred = function() C2


    mypar1 = function ()
    {
        CL = exp(THETA[1] + ETA[1])
        V = exp(THETA[2] + ETA[2])
    }


    m2ag <- rxSymPySetupPred(m1, pred, mypar1, function(){err ~ prop(0.1)}, grad=TRUE)

    omega <- matrix(c(0.1, 0, 0, 0.1), nrow=2)

    symo <- rxSymInvCreate(omega);

    symenv <- rxSymInv(symo, c(sqrt(0.1), sqrt(0.1)))

    THETA <- c(1.6, 4.5, sqrt(0.1));
    ETA <- c(-0.147736086922763, -0.294637022436797)


    tmp1 <- m2ag$outer %>% solve(ev, theta=THETA, eta=ETA)

    tmp2 <- m2ag %>% rxFoceiTheta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv)

    tmp2.nm <- m2ag %>% rxFoceiTheta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, nonmem=TRUE)

    tmp5.g <- m2ag %>%rxFoceiInner(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, invisible=1, scale.to=1)

    tmp5.g2 <- m2ag %>%rxFoceiInner(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, invisible=1,
                                    inits.vec=rep(0.5, 5), scale.to=1)

    test_that("rxFoceiTheta makes sense", {

        expect_equal(tmp1$rx_pred_, tmp2$f); ## F
        err <- matrix(tmp1$rx_pred_ - dv, ncol=1)
        expect_equal(err, tmp2$err) ## Err
        R <- matrix(tmp1$rx_r_, ncol=1)
        expect_equal(R, tmp2$R) ## R (Varinace)
        m <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_", "_sens_rx_pred__ETA_2_")])
        dimnames(m) <- list(NULL, NULL)
        m2 <- tmp2$dErr
        dimnames(m2) <- list(NULL, NULL)
        expect_equal(m, m2) ## Check dErr
        m <- as.matrix(tmp1[, c("_sens_rx_r__ETA_1_", "_sens_rx_r__ETA_2_")]);
        dimnames(m) <- list(NULL, NULL);
        m2 <- tmp2$dR
        dimnames(m2) <- list(NULL, NULL)
        expect_equal(m, m2) ## Check dR
        c <- list(matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1),matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1))
        expect_equal(c, tmp2$c) ## c
        B <- matrix(2 / tmp1[["rx_r_"]],ncol=1)
        expect_equal(B, tmp2$B)

        a <- list(matrix(tmp1[["_sens_rx_pred__ETA_1_"]],ncol=1) - err/R*matrix(tmp1[["_sens_rx_r__ETA_1_"]]),
                  matrix(tmp1[["_sens_rx_pred__ETA_2_"]],ncol=1)- err/R*matrix(tmp1[["_sens_rx_r__ETA_2_"]]));

        expect_equal(a, tmp2$a);

        a <- list(matrix(tmp2.nm$dErr[, 1], ncol=1),
                  matrix(tmp2.nm$dErr[, 2], ncol=1))
        expect_equal(tmp2.nm$a, a)

        ## Does not include the matrix mult part (done with RcppArmadillo)
        lp <- matrix(c(NA, NA), ncol=1)
        c <- matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1)
        expect_equal(c, tmp2$c[[1]])
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_" )])
        lp[1, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)
        c <- matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1)
        expect_equal(c, tmp2$c[[2]])
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_2_" )])
        lp[2, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)

        expect_equal(lp, tmp2$lp)

        llik <- -0.5 * sum(err ^ 2 / R + log(R));

        expect_equal(llik, tmp2$llik);

        ##
        expect_equal(tmp2$dErr.dTheta,
                     matrix(c(tmp1[["_sens_rx_pred__THETA_1_"]],tmp1[["_sens_rx_pred__THETA_2_"]],tmp1[["_sens_rx_pred__THETA_3_"]],
                              rep(0, dim(tmp2$dErr.dTheta)[1] * 2)),ncol=5))

        expect_equal(tmp2$dR.dTheta,
                     matrix(c(tmp1[["_sens_rx_r__THETA_1_"]],tmp1[["_sens_rx_r__THETA_2_"]],tmp1[["_sens_rx_r__THETA_3_"]],
                              rep(0, dim(tmp2$dErr.dTheta)[1] * 2)),ncol=5))

        expect_equal(tmp2$dR2,
                     list(matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__ETA_1_"]],tmp1[["_sens_rx_r__BY_ETA_1__ETA_2_"]]),ncol=2),
                          matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__ETA_2_"]],tmp1[["_sens_rx_r__BY_ETA_2__ETA_2_"]]),ncol=2)
                          ))

        expect_equal(tmp2$dErr2,
                     list(matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__ETA_1_"]],tmp1[["_sens_rx_pred__BY_ETA_1__ETA_2_"]]),ncol=2),
                          matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__ETA_2_"]],tmp1[["_sens_rx_pred__BY_ETA_2__ETA_2_"]]),ncol=2)));

        ## Now Test dErr.dEta.dTheta
        expect_equal(tmp2$dErr.dEta.dTheta,
                     list(matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__THETA_1_"]],
                                   tmp1[["_sens_rx_pred__BY_ETA_2__THETA_1_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__THETA_2_"]],
                                   tmp1[["_sens_rx_pred__BY_ETA_2__THETA_2_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__THETA_3_"]],
                                   tmp1[["_sens_rx_pred__BY_ETA_2__THETA_3_"]]), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2)));

        ## Now test dR.dEta.dTheta
        expect_equal(tmp2$dR.dEta.dTheta,
                     list(matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__THETA_1_"]],
                                   tmp1[["_sens_rx_r__BY_ETA_2__THETA_1_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__THETA_2_"]],
                                   tmp1[["_sens_rx_r__BY_ETA_2__THETA_2_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__THETA_3_"]],
                                   tmp1[["_sens_rx_r__BY_ETA_2__THETA_3_"]]), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2)))

        ## Now test H2.
        h2f <- function(k, l){
            ## Equation 13
            dErr.k <- tmp2$dErr[, k];
            dErr.l <- tmp2$dErr[, l];
            dR.k <- tmp2$dR[, k];
            dR.l <- tmp2$dR[, l];
            dErr.k.l <- tmp2$dErr2[[k]][, l];
            dR.k.l <- tmp2$dR2[[k]][, l];
            err <- tmp2$err;
            R <- tmp2$R;
            return(-0.5 * sum(2 * dErr.l * dErr.k / R - 2 * err * dR.l * dErr.k / (R * R) +
                              2 * err * dErr.k.l / R - err * err * dR.k.l / (R * R) +
                              2 * err * err * dR.k * dR.l / (R * R * R) -
                              2 * err * dR.k * dErr.l / (R * R) -
                              dR.k * dR.l / (R * R) + dR.k.l / R) - symenv$omegaInv[k, l]);
        }

        h2 <- matrix(c(h2f(1, 1), h2f(2, 1), h2f(1, 2), h2f(2, 2)), ncol=2)

        expect_equal(h2, tmp2$H2);

        ## Now test l.dEta.dTheta

        lEH <- function(k, m){
            ## Equation 47
            dErr.m <- tmp2$dErr.dTheta[, m];
            dErr.k <- tmp2$dErr[, k];
            dErr.k.m <- tmp2$dErr.dEta.dTheta[[m]][, k];
            ##
            dR.m <- tmp2$dR.dTheta[, m];
            dR.k <- tmp2$dR[, k];
            dR.k.m <- tmp2$dR.dEta.dTheta[[m]][, k];
            ##
            err <- tmp2$err;
            R <- tmp2$R;
            nomega <- length(tmp2$dOmega)
            ntheta <- tmp2$ntheta;
            if (m > ntheta){
                ome <- tmp2$omega.47[k, m - ntheta];
            } else {
                ome <- 0;
            }
            return(-0.5 * sum(2 * dErr.m * dErr.k / R -
                              2 * err * dR.m * dErr.k / (R * R) +
                              2 * err * dErr.k.m / R -
                              err * err * dR.k.m / (R * R) +
                              2 * err * err * dR.k * dR.m / (R * R * R) -
                              2 * err * dR.k * dErr.m / (R * R) +
                              dR.m * dR.k / (R * R) +
                              dR.k.m / R) - ome);
        }

        df <- expand.grid(k=c(1, 2), theta=1:5);

        tmp2a <- matrix(as.vector(apply(df, 1, function(x) {return(lEH(x[1], x[2]))})), nrow=2)

        expect_equal(tmp2a, tmp2$l.dEta.dTheta)

        ## Now test  #46
        expect_equal(tmp2$dEta.dTheta, -solve(tmp2$H2) %*% tmp2$l.dEta.dTheta);

        ## Now test dErr.dTheta. (Equation #33)
        f <- function(m){
            tmp2$dErr.dTheta[, m] + tmp2$dErr %*% tmp2$dEta.dTheta[, m]
        }

        expect_equal(tmp2$dErr.dTheta., cbind(f(1), f(2), f(3), f(4), f(5)))

        f <- function(m){
            tmp2$dR.dTheta[, m] + tmp2$dR %*% tmp2$dEta.dTheta[, m]
        }

        expect_equal(tmp2$dR.dTheta., cbind(f(1), f(2), f(3), f(4), f(5)));

        ## Now #37
        f <- function(k, m){
            dErr.k.m <- tmp2$dErr.dEta.dTheta[[m]][, k];
            err1 <- tmp2$dErr2[[k]];
            dNdH <- tmp2$dEta.dTheta[, m];
            return(dErr.k.m + err1 %*% dNdH);
        }

        expect_equal(tmp2$dErr.dEta.dTheta., list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                                                  matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)))

        ## Now #37 for dR
        f <- function(k, m){
            dR.k.m <- tmp2$dR.dEta.dTheta[[m]][, k];
            err1 <- tmp2$dR2[[k]];
            dNdH <- tmp2$dEta.dTheta[, m];
            return(dR.k.m + err1 %*% dNdH);
        }

        expect_equal(tmp2$dR.dEta.dTheta., list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                                                matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)))

        ## Now dc/dTheta
        f <- function(k, m){
            dR.m <- tmp2$dR.dTheta.[, m];
            dR.k <- tmp2$dR[, k];
            dR.k.m <- tmp2$dR.dEta.dTheta.[[k]][, m];
            R <- tmp2$R;
            return(-dR.m * dR.k / (R * R) + dR.k.m / R);
        }

        expect_equal(tmp2$dc.dTheta, list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                                          matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)))

        ## Now dB/dTheta
        f <- function(m){
            dR.m <- tmp2$dR.dTheta.[, m];
            R <- tmp2$R;
            return(-2 * dR.m / (R * R));
        }

        expect_equal(tmp2$dB.dTheta, matrix(c(f(1), f(2), f(3), f(4), f(5)), ncol=5))

        ## Now da/dTheta
        f <- function(k, m){
            dErr.m.k <- tmp2$dErr.dEta.dTheta.[[k]][, m];
            dErr.m <- tmp2$dErr.dTheta.[, m];
            dR.k <- tmp2$dR[, k];
            dR.m <- tmp2$dR.dTheta.[, m];
            dR.m.k <- tmp2$dR.dEta.dTheta.[[k]][, m];
            err <- tmp2$err;
            R <- tmp2$R;
            return(dErr.m.k - dErr.m * dR.k / R +
                   err * dR.m * dR.k / (R * R) -
                   err * dR.m.k / R);
        }

        expect_equal(tmp2$da.dTheta,
                     list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                          matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)));

        f <- function(k, m){
            dErr.m.k <- tmp2$dErr.dEta.dTheta.[[k]][, m];
            return(dErr.m.k);
        }

        expect_equal(tmp2.nm$da.dTheta,
                     list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                          matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)));

        ## Test dH/dTheta

        f <- function(m, k, l){
            da.l.m <- tmp2$da.dTheta[[l]][, m];
            B <- tmp2$B
            a.k <- tmp2$a[[k]];
            a.l <- tmp2$a[[l]];
            B.m <- tmp2$dB.dTheta[, m];
            da.k.m <- tmp2$da.dTheta[[k]][, m];
            dc.l.m <- tmp2$dc.dTheta[[l]][, m];
            dc.k.m <- tmp2$dc.dTheta[[k]][, m];
            c.k <- tmp2$c[[k]];
            c.l <- tmp2$c[[l]];
            if (m > 3){
                ome <- symenv$dOmegaInv[[m - 3]][k, l];
            } else {
                ome <- 0;
            }
            return(-0.5 * sum(da.l.m * B * a.k +
                              a.l * B.m * a.k +
                              a.l * B * da.k.m -
                              dc.l.m * c.k -
                              c.l * dc.k.m) - ome);
        }

        df <- expand.grid(m=1:5, k=1:2, l=1:2);
        df <- df[order(df$m, df$k, df$l), ];

        v <- apply(df, 1, function(x) {return(f(x[1], x[2], x[3]))})

        expect_equal(tmp2$dH.dTheta,
                     list(matrix(v[1:4], 2),
                          matrix(v[5:8], 2),
                          matrix(v[9:12], 2),
                          matrix(v[13:16], 2),
                          matrix(v[17:20], 2)))

        f <- function(m, k, l){
            da.l.m <- tmp2.nm$da.dTheta[[l]][, m];
            B <- tmp2.nm$B
            a.k <- tmp2.nm$a[[k]];
            a.l <- tmp2.nm$a[[l]];
            B.m <- tmp2.nm$dB.dTheta[, m];
            da.k.m <- tmp2.nm$da.dTheta[[k]][, m];
            dc.l.m <- tmp2.nm$dc.dTheta[[l]][, m];
            dc.k.m <- tmp2.nm$dc.dTheta[[k]][, m];
            c.k <- tmp2.nm$c[[k]];
            c.l <- tmp2.nm$c[[l]];
            if (m > 3){
                ome <- symenv$dOmegaInv[[m - 3]][k, l];
            } else {
                ome <- 0;
            }
            return(-0.5 * sum(da.l.m * B * a.k +
                              a.l * B.m * a.k +
                              a.l * B * da.k.m +
                              dc.l.m * c.k +
                              c.l * dc.k.m) - ome);
        }

        v <- apply(df, 1, function(x) {return(f(x[1], x[2], x[3]))})

        expect_equal(tmp2.nm$dH.dTheta,
                     list(matrix(v[1:4], 2),
                          matrix(v[5:8], 2),
                          matrix(v[9:12], 2),
                          matrix(v[13:16], 2),
                          matrix(v[17:20], 2)));

        ## Test Inverses
        expect_equal(solve(tmp2$H), tmp2$Hinv);
        expect_equal(solve(tmp2.nm$H), tmp2.nm$Hinv);

        f <- function(m){
            dErr.m <- tmp2$dErr.dTheta[, m];
            dR.m <- tmp2$dR.dTheta[, m];
            err <- tmp2$err;
            R <- tmp2$R;
            eta <- tmp2$eta.mat
            if (m > 3){
                ome <- tmp2$omega.28[m - 3]
            } else {
                ome <- 0;
            }
            Hinv <- tmp2$Hinv;
            dH <- tmp2$dH.dTheta[[m]];
            return(-0.5 * sum(2 * err * dErr.m / R -
                              err * err * dR.m / (R * R) +
                              dR.m / R)  + ome - 0.5 * sum(diag(Hinv %*% dH)));
        }

        expect_equal(tmp2$l.dTheta, c(f(1), f(2), f(3), f(4), f(5)))

        f <- function(m){
            dErr.m <- tmp2.nm$dErr.dTheta[, m];
            dR.m <- tmp2.nm$dR.dTheta[, m];
            err <- tmp2.nm$err;
            R <- tmp2.nm$R;
            eta <- tmp2.nm$eta.mat
            if (m > 3){
                ome <- tmp2.nm$omega.28[m - 3]
            } else {
                ome <- 0;
            }
            Hinv <- tmp2.nm$Hinv;
            dH <- tmp2.nm$dH.dTheta[[m]];
            return(-0.5 * sum(2 * err * dErr.m / R -
                              err * err * dR.m / (R * R) +
                              dR.m / R)  + ome - 0.5 * sum(diag(Hinv %*% dH)));
        }

        expect_equal(tmp2.nm$l.dTheta, c(f(1), f(2), f(3), f(4), f(5)))

        ## Now test omega.28
        f <- function(m){
            eta <- tmp2$eta.mat;
            omegaInv <- symenv$omegaInv;
            dOmega <- symenv$dOmega[[m]];
            return(0.5 * t(eta) %*% omegaInv %*% dOmega %*% omegaInv %*% eta - 0.5 * sum(diag(omegaInv %*% dOmega)))
        }

        expect_equal(tmp2$omega.28, c(f(1), f(2)))

        ## Now test omega.47

        f <- function(m, k){
            eta <- tmp2$eta.mat;
            eta1 <- rep(0, length(as.vector(eta)));
            eta1[k] <- 1;
            eta1 <- matrix(eta1, ncol=1);
            omegaInv <- symenv$omegaInv;
            dOmega <- symenv$dOmega[[m]];
            return(t(eta) %*% omegaInv %*% dOmega %*% omegaInv %*% eta1);
        }

        df <- expand.grid(m=1:2, k=1:2);
        df <- df[order(df$m), ];

        v <- matrix(apply(df, 1, function(x) {return(f(x[1], x[2]))}), nrow=2);

        expect_equal(v, tmp2$omega.47);

        ## Test scaling
        expect_equal(attr(tmp5.g,"grad") * 2, attr(tmp5.g2,"grad"));

        expect_equal(attr(tmp5.g,"dEta.dTheta"), attr(tmp5.g2,"dEta.dTheta"));

    })
}, on.validate="NLMIXR_VALIDATION", silent=TRUE)
