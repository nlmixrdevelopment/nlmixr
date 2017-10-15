rxPermissive({
    context("Focei Inner");

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

    context("Inner Problem Tests")

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

    m2a <- rxSymPySetupPred(m1, pred, mypar1, function(){err ~ prop(0.1)})

    omega <- matrix(c(0.1, 0, 0, 0.1), nrow=2)

    symo <- rxSymInvCreate(omega);

    symenv <- rxSymInv(symo, c(sqrt(0.1), sqrt(0.1)))

    THETA <- c(1.6, 4.5, sqrt(0.1));
    ETA <- c(-0.147736086922763, -0.294637022436797)

    tmp1 <- m2a$inner %>% rxSolve(ev, theta=THETA, eta=ETA)

    tmp2 <- m2a %>% rxFoceiEta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv)

    tmp2.nm <- m2a %>% rxFoceiEta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, nonmem=TRUE)

    tmp3 <- m2a %>% rxFoceiLik(ev, theta=THETA, eta=ETA,dv=dv, inv.env=tmp2)

    tmp4 <- m2a %>% rxFoceiLp(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv)

    tmp5 <- m2a %>%rxFoceiInner(ev, theta=THETA, eta=c(10, 10),dv=dv, inv.env=symenv, invisible=1)

    tmp5 <- m2a %>%rxFoceiInner(ev, theta=THETA, eta=c(10, 10),dv=dv, inv.env=symenv, invisible=1, pred.minus.dv=FALSE)

    tmp5a <- m2a %>%rxFoceiInner(ev, theta=THETA, eta=c(0, 0),dv=dv, inv.env=symenv, invisible=1, nonmem=TRUE)

    test_that("rxFoceiEta makes sense", {

        expect_equal(tmp1$rx_pred_, tmp2$f); ## F
        err <- matrix(tmp1$rx_pred_ - dv , ncol=1)
        expect_equal(err, tmp2$err) ## Err
        R <- matrix(tmp1$rx_r_, ncol=1) ## check
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
        expect_equal(c, tmp2$c) ## check
        B <- matrix(2 / tmp1[["rx_r_"]],ncol=1)
        expect_equal(B, tmp2$B) ## check

        ## This is different...  According to the Almquist paper 2015
        a <- list(matrix(tmp1[["_sens_rx_pred__ETA_1_"]],ncol=1) - err/R*matrix(tmp1[["_sens_rx_r__ETA_1_"]]),
                  matrix(tmp1[["_sens_rx_pred__ETA_2_"]],ncol=1)- err/R*matrix(tmp1[["_sens_rx_r__ETA_2_"]]));
        expect_equal(a, tmp2$a);

        a <- list(matrix(tmp2.nm$dErr[, 1], ncol=1),
                  matrix(tmp2.nm$dErr[, 2], ncol=1))
        expect_equal(tmp2.nm$a, a)

        ## Does not include the matrix multilpaction part (done with RcppArmadillo)
        lp <- matrix(c(NA, NA), ncol=1)
        c <- matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1)
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_" )])
        lp[1, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)
        c <- matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1)
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_2_" )])
        lp[2, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)
        expect_equal(lp, tmp2$lp)

        llik <- -0.5 * sum(err ^ 2 / R + log(R));

        expect_equal(llik, tmp2$llik);

        eta <- ETA
        omegaInv <- symenv$omegaInv
        llik <- -0.5 * sum(err ^ 2 / R + log(R)) - 0.5 * t(matrix(eta,ncol=1)) %*% omegaInv %*% matrix(eta,ncol=1);
        llik <- -llik;

        expect_equal(as.vector(llik), tmp3);

        lp2 <- lp - omegaInv %*% matrix(eta, ncol=1)

        expect_equal(-lp2, tmp4);

        ## Now Test the Hessian calculation
        Hidx <- mat.indices(2);
        H <- matrix(1:(2^2), 2, 2)
        k = Hidx$idx[,1];
        l = Hidx$idx[,2];
        a = cbind(tmp2$a[[1]], tmp2$a[[2]])
        B = as.vector(tmp2$B);
        c = cbind(tmp2$c[[1]], tmp2$c[[2]])
        NONMEM <- 0;
        H[Hidx$Hlo.idx] = apply(a[, k]*B*a[, l] + c(-1, 1)[NONMEM+1] *c[,k]*c[,l], 2, sum)
        H[Hidx$hi.idx] = H[Hidx$lo.idx]
        H = -.5*H - omegaInv;

        rxHessian(tmp2);
        dimnames(tmp2$H) <- list(NULL, NULL);
        dimnames(H) <- list(NULL, NULL);
        expect_equal(tmp2$H, H);

        H.neg.5 = tryCatch({
            chol(-H)
        }, error = function(e) {
            cat("Warning: Hessian not positive definite\n")
            print(-H)
            .m <- -H
            .md = matrix(0, nETA, nETA)
            diag(.md) = abs(diag(.m)) * 1.1 + .001
            chol(.md)
        })
        ##H.neg.5 = chol(-H)
        log.det.H.neg.5 = sum(log(diag(H.neg.5)))

        RxODE_focei_eta_lik(tmp2$eta, tmp2)
        log.det.OMGAinv.5 <- symenv$log.det.OMGAinv.5
        llik = -tmp2$llik2 + log.det.OMGAinv.5             #note no -1/2 with OMGAinv.5
        llik.lapl =llik - log.det.H.neg.5               #note no 1/2
        llik.lapl = as.vector(llik.lapl);
        attr(llik.lapl, "fitted") <- tmp2$f;
        attr(llik.lapl, "posthoc") <- tmp2$eta;
        attr(llik.lapl, "corrected") <- 0L

        lik2 <- RxODE_focei_finalize_llik(tmp2);

        attr(lik2, "Vi") <- NULL; # FIXME test
        attr(lik2, "Vfo") <- NULL; # FIXME test
        attr(lik2, "dErr_dEta") <- NULL # FIXME test
        attr(lik2, "omega.28") <- NULL

        expect_equal(lik2, llik.lapl);

        ## Now test NONMEM approximation.
        NONMEM <- 1;
        Hidx <- mat.indices(2);
        H <- matrix(1:(2^2), 2, 2)
        a = cbind(tmp2.nm$a[[1]], tmp2.nm$a[[2]])
        B = as.vector(tmp2.nm$B);
        c = cbind(tmp2.nm$c[[1]], tmp2.nm$c[[2]])
        H[Hidx$Hlo.idx] = apply(a[, k]*B*a[, l] + c(-1, 1)[NONMEM+1] *c[,k]*c[,l], 2, sum)
        H[Hidx$hi.idx] = H[Hidx$lo.idx]
        H = -.5*H - omegaInv;

        rxHessian(tmp2.nm);
        dimnames(tmp2.nm$H) <- list(NULL, NULL);
        dimnames(H) <- list(NULL, NULL);
        expect_equal(tmp2.nm$H, H);

        H.neg.5 = tryCatch({
            chol(-H)
        }, error = function(e) {
            cat("Warning: Hessian not positive definite\n")
            print(-H)
            .m <- -H
            .md = matrix(0, nETA, nETA)
            diag(.md) = abs(diag(.m)) * 1.1 + .001
            chol(.md)
        })
        ##H.neg.5 = chol(-H)
        log.det.H.neg.5 = sum(log(diag(H.neg.5)))

        RxODE_focei_eta_lik(tmp2.nm$eta, tmp2.nm)

        llik = -tmp2.nm$llik2 + log.det.OMGAinv.5             #note no -1/2 with OMGAinv.5
        llik.lapl =llik - log.det.H.neg.5               #note no 1/2
        llik.lapl = as.vector(llik.lapl);
        attr(llik.lapl, "fitted") <- tmp2.nm$f;
        attr(llik.lapl, "posthoc") <- tmp2.nm$eta;
        attr(llik.lapl, "corrected") <- 0L

        lik2 <- RxODE_focei_finalize_llik(tmp2.nm);

        attr(lik2, "Vi") <- NULL; # FIXME test
        attr(lik2, "Vfo") <- NULL; # FIXME test
        attr(lik2, "dErr_dEta") <- NULL # FIXME test
        attr(lik2, "omega.28") <- NULL # FIXME test

        expect_equal(lik2, llik.lapl);

        expect_equal(as.vector(lik2), -209.467627968204);

        expect_equal(as.vector(tmp5a), -209.467627968204);

    })


}, on.validate="NLMIXR_VALIDATION", silent=TRUE)
