rxPermissive({
    ## From https://raw.githubusercontent.com/bbolker/broom.mixed/master/tests/testthat/helper-checkers.R

    ##' test the basics of tidy/augment/glance output: is a data frame, no row names
    check_tidiness <- function(o) {
        testthat::expect_is(o, "tbl_df")
        testthat::expect_equal(rownames(o), as.character(seq_len(nrow(o))))
    }


    #' check the output of a tidy function
    check_tidy <- function(o, exp.row = NULL, exp.col = NULL, exp.names = NULL) {
        check_tidiness(o)

        if (!is.null(exp.row)) {
            testthat::expect_equal(nrow(o), exp.row)
        }
        if (!is.null(exp.col)) {
            testthat::expect_equal(ncol(o), exp.col)
        }
        if (!is.null(exp.names)) {
            testthat::expect_true(all(exp.names %in% colnames(o)))
        }
    }

    library(testthat)

    options(nlmixr.save=TRUE,
            nlmixr.save.dir=system.file(package="nlmixr"));

    one.compartment <- function() {
        ini({
            tka <- 0.45 # Log Ka
            tcl <- 1 # Log Cl
            tv <- 3.45    # Log V
            eta.ka ~ 0.6
            eta.cl ~ 0.3
            eta.v ~ 0.1
            add.err <- 0.7
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d/dt(depot) = -ka * depot
            d/dt(center) = ka * depot - cl / v * center
            cp = center / v
            cp ~ add(add.err)
        })
    }


    context("broom tidy nlmixr SAEM")

    fitS <- nlmixr(one.compartment, theo_sd, est="saem")
    test_that("tidy works on nlmixr fit SAEM fits", {
        td <- tidy(fitS)
        check_tidy(td,7,5,c("effect", "group", "term", "estimate", "std.error"));
        expect_equal(
            td$term,
            c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
              "add.err")
        )
        td <- tidy(fitS, conf.level=0.9)
        check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error",
                               "conf.low", "conf.high"))
        expect_equal(
            td$term,
            c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
              "add.err")
        )
        expect_equal(td$estimate,c(1.57029661020455, 2.76582308297812, 31.4669109412755, 0.645283223343858,
                                   0.264933148852423, 0.138974316163357, 0.691576073984559))
        expect_equal(td$std.error, c(0.308138001334383, 0.231308197834056, 1.47478964883279, NA,
                                     NA, NA, NA))
        expect_equal(td$conf.low, c(1.13711717110817, 2.41036388512486, 29.1322447510814, NA, NA,
                                    NA, NA))

        td <- tidy(fitS, conf.level=0.9,exponentiate=FALSE)
        check_tidy(td)
        expect_equal(td$estimate, c(0.451264525213545, 1.01733826987306, 3.44893654741748, 0.645283223343858,
                                    0.264933148852423, 0.138974316163357, 0.691576073984559))
        expect_equal(td$std.error, c(0.196229170547751, 0.0836308725809727, 0.0468679512769808,
                                     NA, NA, NA, NA),tolerance=1e-5)
        expect_equal(td$conf.low, c(0.128496262324398, 0.879777725783136, 3.37184562777175, NA,
                                    NA, NA, NA))

        for (ef in c("ran_vals", "random")){
            td  <- tidy(fitS, effects=ef)
            td1  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            td  <- tidy(fitS, effects=ef, exponentiate=FALSE)
            td2  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            td  <- tidy(fitS, effects=ef, exponentiate=TRUE)
            td3  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            expect_equal(td1, td2)
            expect_equal(td2, td3)
        }

        td  <- tidy(fitS, effects="ran_coef")
        td1  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(td1, c(1.74332254976563, 1.95893430769561, 2.26961708166631, 1.1902950696036,
                            1.50007267530839, 1.07145178365103, 0.711788830300308, 1.31089126598457,
                            6.49753490073617, 0.750856047789294, 3.46062641589219, 0.926766706324023,
                            1.69997844437216, 3.19499780816707, 2.85330242232265, 2.72475980105706,
                            2.36872877586276, 3.99217464726057, 3.24641810409366, 3.27189857789171,
                            2.88537953056852, 1.87135452302142, 3.66470901400725, 2.43719470204722,
                            29.0464831288899, 32.1223417486136, 33.3443913346158, 31.2402834906231,
                            27.3050749891848, 38.5330188945709, 32.9677569899538, 34.6297085084578,
                            31.8509966680452, 26.5593807055713, 36.4042269585165, 25.8134435928474))

        td  <- tidy(fitS, effects="ran_coef", exponentiate=FALSE)
        td2  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td  <- tidy(fitS, effects="ran_coef", exponentiate=TRUE)
        td3  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(log(td1),td2)
        expect_equal(td2,log(td3))

        td  <- tidy(fitS, effects="ran_pars")
        td1  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
        expect_equal(td1, c(0.645283223343858, 0.264933148852423, 0.138974316163357, 0.691576073984559))

        td  <- tidy(fitS, effects="ran_pars", exponentiate=FALSE)
        td2  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

        td  <- tidy(fitS, effects="ran_pars", exponentiate=TRUE)
        td3  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

        expect_equal(td1, td2)
        expect_equal(td2, td3)
    })

    for (f in c("focei", "foce")){
        context(sprintf("broom tidy nlmixr %s", f))
        fitF <- nlmixr(one.compartment, theo_sd, est=f)
        test_that(sprintf("tidy works on nlmixr fit %s fits",f), {
            td <- tidy(fitF)
            check_tidy(td,7,5,c("effect", "group", "term", "estimate", "std.error"));
            expect_equal(
                td$term,
                c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
                  "add.err")
            )
            td <- tidy(fitF, conf.level=0.9)
            check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error",
                                   "conf.low", "conf.high"))
            expect_equal(
                td$term,
                c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
                  "add.err")
            )
            expect_equal(td$estimate, c(1.59887298634974, 2.751052394453, 31.7991981878092, 0.633040351157386,
                                        0.262526810818508, 0.13901259012342, 0.694382258169284))
            expect_equal(td$std.error,c(0.316777567222603, 0.636930199174004, 2.02365165436791, NA,
                                        NA, NA, NA))
            expect_equal(td$conf.low, c(1.15420464906127, 1.87979569886401, 28.6388773164848, NA, NA,
                                        NA, NA))

            td <- tidy(fitF, conf.level=0.9,exponentiate=FALSE)
            check_tidy(td)
            expect_equal(td$estimate, c(0.469298997519686, 1.01198352736376, 3.45944107524841, 0.633040351157386,
                                        0.262526810818508, 0.13901259012342, 0.694382258169284))
            expect_equal(td$std.error, c(0.198125535878753, 0.231522380474562, 0.0636384490708233, NA,
                                         NA, NA, NA),tolerance=1e-5)
            expect_equal(td$conf.low, c(0.143411491237816, 0.631163100119741, 3.3547651414807, NA,
                                        NA, NA, NA))
            for (ef in c("ran_vals", "random")){
                td  <- tidy(fitF, effects=ef)
                td1  <- td$estimate
                check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
                ##
                td  <- tidy(fitF, effects=ef, exponentiate=FALSE)
                td2  <- td$estimate
                check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
                ##
                td  <- tidy(fitF, effects=ef, exponentiate=TRUE)
                td3  <- td$estimate
                check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
                ##
                expect_equal(td1, td2)
                expect_equal(td2, td3)
            }
            ##
            td  <- tidy(fitF, effects="ran_coef")
            td1  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
            ##
            expect_equal(td1, c(1.73333605124159, 1.93923492222759, 2.30743914124927, 1.20006316113323,
                                1.52351562065568, 1.08314630146302, 0.735683438046524, 1.34571208216029,
                                6.2155586111074, 0.764076492994444, 3.38332449848861, 0.945602258131827,
                                1.71270724282176, 3.17778074697626, 2.82946200590297, 2.70153237014438,
                                2.36803819562017, 3.98738977025344, 3.1830250393921, 3.24250654748063,
                                2.87891653916599, 1.87860041579119, 3.66735073724982, 2.42767276009263,
                                29.0286578425194, 31.9701679855541, 33.5413148355593, 31.3937950732769,
                                27.5493232024115, 38.6357092611212, 33.6218848424434, 34.9365714093507,
                                31.7987634177591, 26.7756665146766, 36.467713097871, 25.9856375895782))
            ##
            td  <- tidy(fitF, effects="ran_coef", exponentiate=FALSE)
            td2  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
            ##
            td  <- tidy(fitF, effects="ran_coef", exponentiate=TRUE)
            td3  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
            ##
            expect_equal(log(td1),td2)
            expect_equal(td2,log(td3))
            ##
            td  <- tidy(fitF, effects="ran_pars")
            td1  <- td$estimate
            check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
            expect_equal(td1, c(0.633040351157386, 0.262526810818508, 0.13901259012342, 0.694382258169284))
            ##
            td  <- tidy(fitF, effects="ran_pars", exponentiate=FALSE)
            td2  <- td$estimate
            check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
            ##
            td  <- tidy(fitF, effects="ran_pars", exponentiate=TRUE)
            td3  <- td$estimate
            check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
            ##
            expect_equal(td1, td2)
            expect_equal(td2, td3)
        })
    }

    for (f in c("foi", "fo")){
        context(sprintf("broom tidy nlmixr %s", f))
        fitF <- nlmixr(one.compartment, theo_sd, est=f)
        test_that(sprintf("tidy works on nlmixr fit %s fits",f), {
            td <- tidy(fitF)
            check_tidy(td,7,5,c("effect", "group", "term", "estimate", "std.error"));
            expect_equal(
                td$term,
                c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
                  "add.err")
            )
            td <- tidy(fitF, conf.level=0.9)
            check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error",
                                   "conf.low", "conf.high"))
            expect_equal(
                td$term,
                c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
                  "add.err")
            )
            expect_equal(td$estimate, c(2.89710832634736, 2.77084291838898, 33.267220888755, 0.913975198223592,
                                        0.301151147027046, 0.138493033269316, 0.614809287631279))
            expect_equal(td$std.error, c(0.570538001778867, 0.252416964062574, 1.42574896058229, NA, NA, NA, NA))
            expect_equal(td$conf.low, c(2.09548734182092, 2.38526317259775, 31.0028237043988, NA, NA, NA, NA))

            td <- tidy(fitF, conf.level=0.9,exponentiate=FALSE)
            check_tidy(td)
            expect_equal(td$estimate, c(1.063713110683, 1.01915157657608, 3.50457255449315, 0.913975198223592,
                                        0.301151147027046, 0.138493033269316, 0.614809287631279))
            expect_equal(td$std.error, c(0.196933610176115, 0.0910975365609443, 0.0428574711831196,
                                         NA, NA, NA, NA),tolerance=1e-5)
            expect_equal(td$conf.low, c(0.739786147716175, 0.869309463157471, 3.43407828757563, NA,
                                        NA, NA, NA))
            for (ef in c("ran_vals", "random")){
                td  <- tidy(fitF, effects=ef)
                td1  <- td$estimate
                check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
                ##
                td  <- tidy(fitF, effects=ef, exponentiate=FALSE)
                td2  <- td$estimate
                check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
                ##
                td  <- tidy(fitF, effects=ef, exponentiate=TRUE)
                td3  <- td$estimate
                check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
                ##
                expect_equal(td1, td2)
                expect_equal(td2, td3)
            }
            ##
            td  <- tidy(fitF, effects="ran_coef")
            td1  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
            ##
            expect_equal(td1, c(1.78607811545834, 1.99879379710882, 2.4457152378203, 1.21852345600953,
                                1.54300845032712, 1.115112941314, 0.736117352263772, 1.38629389306969,
                                7.75235800311953, 0.765254586219227, 3.69885613949955, 0.941784323453709,
                                1.65146288892371, 3.17449221616754, 2.79453246697996, 2.68064354736301,
                                2.34544110189374, 4.03684598011593, 3.19498289650118, 3.23564374066479,
                                2.83543925845555, 1.84492455327811, 3.66817114592915, 2.42085213933379,
                                29.4003244785328, 32.3320477116677, 34.1400211763302, 31.6826596146577,
                                27.698268411408, 39.4145711619581, 33.8798088650667, 35.4916389115595,
                                32.4289892337905, 26.8861472959111, 37.2953913429352, 25.953067911439))
            ##
            td  <- tidy(fitF, effects="ran_coef", exponentiate=FALSE)
            td2  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
            ##
            td  <- tidy(fitF, effects="ran_coef", exponentiate=TRUE)
            td3  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))
            ##
            expect_equal(log(td1),td2)
            expect_equal(td2,log(td3))
            ##
            td  <- tidy(fitF, effects="ran_pars")
            td1  <- td$estimate
            check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
            expect_equal(td1, c(0.913975198223592, 0.301151147027046, 0.138493033269316, 0.614809287631279))
            ##
            td  <- tidy(fitF, effects="ran_pars", exponentiate=FALSE)
            td2  <- td$estimate
            check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
            ##
            td  <- tidy(fitF, effects="ran_pars", exponentiate=TRUE)
            td3  <- td$estimate
            check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
            ##
            expect_equal(td1, td2)
            expect_equal(td2, td3)
        })
    }


    context("broom tidy nlmixr nlme")
    fitN <- nlmixr(one.compartment, theo_sd, est="nlme", control=nlmeControl(pnlsTol=0.6))
    test_that("tidy works on nlmixr fit nlme fits", {
        td <- tidy(fitN)
        check_tidy(td,7,5,c("effect", "group", "term", "estimate", "std.error"));
        expect_equal(
            td$term,
            c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
              "add.err")
        )
        td <- tidy(fitN, conf.level=0.9)
        check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error",
                               "conf.low", "conf.high"))
        expect_equal(
            td$term,
            c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
              "add.err")
        )
        expect_equal(td$estimate,c(1.57024196423777, 2.73781705816092, 32.0487403004864, 0.585672058408116,
                                   0.261936118166243, 0.151850683238502, 0.679682049379696))
        expect_equal(td$std.error, c(0.280687718615492, 0.231106270788711, 1.59868259994279, NA,
                                     NA, NA, NA))
        expect_equal(td$conf.low, c(1.17023530086566, 2.38289117235058, 29.5241297546193, NA, NA,
                                    NA, NA))

        td <- tidy(fitN, conf.level=0.9,exponentiate=FALSE)
        check_tidy(td)
        expect_equal(td$estimate,c(0.451229724834257, 1.00716090876872, 3.46725787839726, 0.585672058408116,
                                   0.261936118166243, 0.151850683238502, 0.679682049379696))
        expect_equal(td$std.error, c(0.178754437219326, 0.0844126053272356, 0.0498828529593886,
                                     NA, NA, NA, NA),tolerance=1e-5)
        expect_equal(td$conf.low, c(0.157204840440379, 0.86831452873579, 3.38520788678432, NA,
                                    NA, NA, NA))

        for (ef in c("ran_vals", "random")){
            td  <- tidy(fitN, effects=ef)
            td1  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            td  <- tidy(fitN, effects=ef, exponentiate=FALSE)
            td2  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            td  <- tidy(fitN, effects=ef, exponentiate=TRUE)
            td3  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            expect_equal(td1, td2)
            expect_equal(td2, td3)
        }

        td  <- tidy(fitN, effects="ran_coef")
        td1  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(td1, c(1.72226403760418, 1.41329990751741, 2.30632275730473, 1.15465139318611,
                            1.51696020368894, 1.1181717241747, 0.885333058618201, 1.35698639312647,
                            6.13921581288724, 0.868994981850219, 3.35175338509676, 0.95132044333853,
                            1.75845332854771, 3.23512284853291, 2.82383851218207, 2.72153935333835,
                            2.41600047651454, 3.89382547226046, 3.02456192246528, 3.22443496749178,
                            2.78349382482623, 1.86124997462384, 3.6180351182618, 2.35895500533327,
                            28.7885643987419, 30.8926359025941, 33.6008055989365, 31.2885290005653,
                            27.5692272091284, 40.5644092052199, 36.0460989861211, 35.128058850313,
                            33.3510634258588, 27.6431872364855, 36.6729317194168, 26.2298710258424
                            ))

        td  <- tidy(fitN, effects="ran_coef", exponentiate=FALSE)
        td2  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td  <- tidy(fitN, effects="ran_coef", exponentiate=TRUE)
        td3  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(log(td1),td2)
        expect_equal(td2,log(td3))

        td  <- tidy(fitN, effects="ran_pars")
        td1  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
        expect_equal(td1, c(0.585672058408116, 0.261936118166243, 0.151850683238502, 0.679682049379696))

        td  <- tidy(fitN, effects="ran_pars", exponentiate=FALSE)
        td2  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

        td  <- tidy(fitN, effects="ran_pars", exponentiate=TRUE)
        td3  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

        expect_equal(td1, td2)
        expect_equal(td2, td3)
    })


    context("broom tidy nlmixr posthoc")
    fitP <- nlmixr(one.compartment, theo_sd, est="posthoc")

    test_that("tidy works on posthoc fit fits", {
        td <- tidy(fitP)
        check_tidy(td,7,5,c("effect", "group", "term", "estimate", "std.error"));
        expect_equal(
            td$term,
            c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
              "add.err")
        )
        td <- tidy(fitP, conf.level=0.9)
        check_tidy(td, 7, 7, c("effect", "group", "term", "estimate", "std.error",
                               "conf.low", "conf.high"))
        expect_equal(
            td$term,
            c("tka", "tcl", "tv", "sd__eta.ka", "sd__eta.cl", "sd__eta.v",
              "add.err")
        )
        expect_equal(td$estimate, c(1.56831218549017, 2.71828182845905, 31.5003923087479, 0.774596669241483,
                                    0.547722557505166, 0.316227766016838, 0.7))
        expect_equal(td$std.error, c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                                     NA_real_))
        expect_equal(td$conf.low, c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                                    NA_real_))

        td <- tidy(fitP, conf.level=0.9,exponentiate=FALSE)
        check_tidy(td)
        expect_equal(td$estimate,c(0.45, 1, 3.45, 0.774596669241483, 0.547722557505166, 0.316227766016838,
                                   0.7))
        expect_equal(td$std.error, c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                                     NA_real_))
        expect_equal(td$conf.low, c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                                    NA_real_))

        for (ef in c("ran_vals", "random")){
            td  <- tidy(fitP, effects=ef)
            td1  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            td  <- tidy(fitP, effects=ef, exponentiate=FALSE)
            td2  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            td  <- tidy(fitP, effects=ef, exponentiate=TRUE)
            td3  <- td$estimate
            check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

            expect_equal(td1, td2)
            expect_equal(td2, td3)
        }

        td  <- tidy(fitP, effects="ran_coef")
        td1  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(td1, c(1.75611262206464, 1.92886734474446, 2.36872361414073, 1.18776870402381,
                            1.48421770349779, 1.1603466406079, 0.728548560764397, 1.37496750947703,
                            6.63504345280614, 0.728113370661147, 3.56192283422217, 0.87667797314966,
                            1.62338691901516, 3.22874000491057, 2.81088025117325, 2.7068918775575,
                            2.37701456492688, 4.04231470515951, 3.23451222726721, 3.26106788182658,
                            2.88069603771061, 1.86803407932834, 3.73237924876953, 2.49388025078765,
                            29.2305913871308, 31.8379070004274, 33.9077584685152, 31.2454399174589,
                            27.0633487903421, 40.6859508443866, 33.6487106457298, 35.5118803877647,
                            31.940988239621, 26.0728299589223, 37.2139995031868, 24.7515885649664))

        td  <- tidy(fitP, effects="ran_coef", exponentiate=FALSE)
        td2  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        td  <- tidy(fitP, effects="ran_coef", exponentiate=TRUE)
        td3  <- td$estimate
        check_tidy(td, 36, 5, c("effect", "group", "level", "term", "estimate"))

        expect_equal(log(td1),td2)
        expect_equal(td2,log(td3))

        td  <- tidy(fitP, effects="ran_pars")
        td1  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))
        expect_equal(td1, c(0.774596669241483, 0.547722557505166, 0.316227766016838, 0.7))

        td  <- tidy(fitP, effects="ran_pars", exponentiate=FALSE)
        td2  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

        td  <- tidy(fitP, effects="ran_pars", exponentiate=TRUE)
        td3  <- td$estimate
        check_tidy(td, 4, 4, c("effect", "group", "term", "estimate"))

        expect_equal(td1, td2)
        expect_equal(td2, td3)
    })
}, cran=TRUE)
