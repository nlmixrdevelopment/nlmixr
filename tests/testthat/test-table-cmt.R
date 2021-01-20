nlmixrTest({
  test_that("proper table outputs", {

    df <- data.frame(
      ID = c(123, 123, 123, 124, 124, 124),
      MDV = c(0,1,0, 0,1,0),
      CMT = c(2,1,1, 2,1,1),
      DV = c(9,NA,15, 10,NA,14),
      AMT = c(NA,15,NA, NA,15,NA),
      RATE = c(NA,0,NA, NA,0,NA),
      ADDL = c(NA,1,NA, NA,1,NA),
      II = c(NA,NA,NA, NA,NA,NA),
      TIME = c(0,0,1, 0,0,1),
      PCA = c(NA,NA,32, NA,NA,32),
      WT = c(NA,NA,2, NA,NA,2),
      CRPZERO = 5
    )



    f <- function(){
      ini({
        # thetas
        theta1  <- c(0.387)     # Add err PK
        theta2  <- c(log(50.9))      # Vd
        theta3  <- c(log(3.42))      # Cl
        theta4  <- c(31.2)      # TM50
        theta5  <- c(3.68)      # Hill
        theta6  <- c(1.46)      # K growth
        theta7  <- c(0.187)     # K death
        theta8  <- c(1.52)      # Emax bact
        theta9  <- c(0.304)     # EC50 bact
        theta10 <- c(4.99)     # Gamma bact
        theta11 <- c(0.162)    # Add err PD
        theta12 <- c(0.274)    # Prop err PD
        theta13 <- c(log(0.276))    # base
        theta14 <- c(log(0.0431))   # Kout
        theta15 <- c(1.22)     # EC50/10^3
        theta16 <- c(0.134)    # Emax
        # ETAs
        eta1 ~ c(1.33)                     # Base
        eta2 + eta3 ~ c(0.194,
                        -0.0159, 0.306) # Vd, Cl
        eta4 + eta5 ~ c(0.521,
                       -0.435, 0.83)   # Emax, Kout
      })
      model({
        # maturation params
        TM50 <- theta4
        HILL <- theta5
        # central compartmental params
        V <- (WT/70)*exp(theta2 + eta2)
        Cl <- ((WT/70)**0.75) * (PCA**HILL) / ((PCA**HILL)+(TM50**HILL))*exp(theta3 + eta3)
        # bacterial params
        KGROWTH <- theta6
        KDEATH <- theta7
        EMAXBACT <- theta8
        EC50BACT <- theta9
        GAMMABACT <- theta10
        BACINIT <- 10**6
        # crp level params
        BASE <- exp(theta13 + eta1)
        KOUT <- exp(theta14 + eta5)
        KIN <- BASE*KOUT
        EC50 <- theta15*1000
        EMAX <- theta16*(1 + eta4)
        # initial conditions
        bact(0) <- BACINIT
        crp(0) <- CRPZERO
        # algebraic expressions
        K = Cl/V
        cp = centr/V
        # differential equations
        d/dt(centr) = -K*centr
        d/dt(crp)   = KIN + EMAX*bact/(EC50 + bact) - KOUT*crp
        d/dt(bact)  = KGROWTH*bact - KDEATH*bact - EMAXBACT*(centr/V)**GAMMABACT/(EC50BACT**GAMMABACT+(centr/V)**GAMMABACT)*bact
        # error model
        cp ~ add(theta1) | centr
        crp ~ prop(theta12) + add(theta11) | crp
      })
    }


    fit.s <- nlmixr(
      object = f,
      data = df,
      est='focei',
      control = foceiControl(
        covMethod="r,s",
        interaction = TRUE,
        maxOuterIterations = 0,
        iter.max=0, calcTables=FALSE)
    )

    tab1 <- addTable(fit.s, table=tableControl(cwres=FALSE, npde=FALSE))
    expect_true(all(c("CMT", "CRPZERO","WT", "PCA") %in% names(tab1)))
    expect_true(all(!is.na(tab1$CMT)))
    expect_true(inherits(tab1$CMT, "factor"))

    tab2 <- addTable(fit.s, table=tableControl(cwres=TRUE, npde=FALSE))
    expect_true(all(c("CMT", "CRPZERO","WT", "PCA") %in% names(tab2)))
    expect_true(all(!is.na(tab2$CMT)))
    expect_true(inherits(tab2$CMT, "factor"))

    tab3 <- addTable(fit.s, table=tableControl(cwres=FALSE, npde=TRUE))
    expect_true(all(c("CMT", "CRPZERO","WT", "PCA") %in% names(tab3)))
    expect_true(all(!is.na(tab3$CMT)))
    expect_true(inherits(tab3$CMT, "factor"))

    tab4 <- addTable(fit.s, table=tableControl(cwres=TRUE, npde=TRUE))
    expect_true(all(c("CMT", "CRPZERO","WT", "PCA") %in% names(tab4)))
    expect_true(all(!is.na(tab4$CMT)))
    expect_true(inherits(tab4$CMT, "factor"))

  })

},
test = "table"
)
