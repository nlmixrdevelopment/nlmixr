nlmixrTest({
  test_that("issue 501", {
    nlmixr_threecmt_mm_no_add_wtcl_pdtg_kout_delay <- function() {
      ini({
        tf_sc <- log(999)
        tf_infilt <- log(999)
        tka_sc <- log(999)
        tka_infilt <- log(999)
        tcl_low <- log(999)
        tcl_high <- log(999)
        tcl_c50 <- log(3000)
        e_wt_cl <- fixed(999)
        tv <- log(999)
        tq1 <- log(999)
        tvp1 <- log(10)
        tq2 <- log(999)
        tvp2 <- log(20)
        eta_cl~999
        eta_v~999
        prop_err <- 999

        tg_bl <- log(999)
        eta_tg_bl~999
        tg_kel <- log(999)
        tg_ec50 <- log(5000)
        tg_emax_kel <- log(2)
        ktr_tg <- log(999)
        prop_err_tg <- 999
      })
      model({
        # PK setup
        f_sc <- exp(tf_sc)
        f_infilt <- exp(tf_infilt)
        ka_sc <- exp(tka_sc)
        ka_infilt <- exp(tka_infilt)
        cl_low <- exp(tcl_low + eta_cl)*(WEIGHT_BL/85)^e_wt_cl
        cl_high <- exp(tcl_high + eta_cl)*(WEIGHT_BL/85)^e_wt_cl
        cl_c50 <- exp(tcl_c50)
        v <- exp(tv + eta_v)
        q1 <- exp(tq1)
        vp1 <- exp(tvp1)
        q2 <- exp(tq2)
        vp2 <- exp(tvp2)

        # PK micro-parameters
        ke_low <- cl_low/v
        ke_high <- cl_high/v
        kc_p1 <- q1/v
        kp1_c <- q1/vp1
        kc_p2 <- q2/v
        kp2_c <- q2/vp2

        # differential equations
        cp <- CENTRAL/v*1e3 # 1e3 is for unit conversion
        ke <- ke_low + (ke_high - ke_low)*cp/(cp + cl_c50)
        d/dt(IVINFILT) =             - ka_infilt * IVINFILT
        d/dt(SC)       = -ka_sc * SC
        d/dt(CENTRAL)  =  ka_sc * SC + ka_infilt * IVINFILT - ke*CENTRAL - kc_p1*CENTRAL + kp1_c*P1 - kc_p2*CENTRAL + kp2_c*P2
        d/dt(P1) =                                                         kc_p1*CENTRAL - kp1_c*P1
        d/dt(P2) =                                                                                    kc_p2*CENTRAL - kp2_c*P2

        f(SC) <- f_sc
        f(IVINFILT) <- f_infilt

        # TG Emax model
        tgbl <- exp(tg_bl + eta_tg_bl)
        kin_tg <- tgbl*exp(tg_kel)
        TG(0) <- tgbl
        ktr_TG <- exp(ktr_tg)
        d/dt(TG_TR) = ktr_tg*cp - ktr_tg*TG_TR
        kout_tg <- exp(tg_kel) + exp(tg_emax_kel)*TG_TR/(TG_TR + exp(tg_ec50))
        d/dt(TG) = kin_tg - kout_tg*TG

        # Residual error models
        cp ~ prop(prop_err)
        TG ~ prop(prop_err_tg)
      })
    }

    expect_error(nlmixr(nlmixr_threecmt_mm_no_add_wtcl_pdtg_kout_delay), NA)
  })
},
test = "lvl2")
