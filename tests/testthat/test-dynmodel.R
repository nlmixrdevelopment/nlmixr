# Error Model Test Script -------------------------------------------------
# NOTES:
# This script is used to test all the error models in the dynmodel()
# Nelder-Mead is used for optimization for all tests
# Each test is a written within a function, just call the function to test
# Each test function can input custom tolerances, and error model parameters
# Tolerances are defined by relative fraction difference between reference values and estiamted values
# #########################################################################

# Library -----------------------------------------------------------------
# devtools::use_testthat() # deprecated
## usethis::use_testthat()
library(RxODE)
library(devtools)
library(latticeExtra)
## devtools::use_testthat()
# devtools::load_all("C:/Users/Mason/OneDrive/Novartis_2019_internship/dynmodel_github/nlmixr/")
## devtools::install("C:/Users/Mason/OneDrive/Novartis_2019_internship/dynmodel_github/nlmixr", dependencies = FALSE, quick=TRUE)
library(nlmixr)

# Additive and Proportional Error Model Test ------------------------------
additiveProportionalErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                                          sigma.add = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL,
                                          method = NULL, covMethod = "nlmixrHess", nlmixrOutput = NULL) {


  # library -----------------------------------------------------------------
  library(gridExtra)
  library(lattice)
  library(testthat)
  graphics.off()

  # RxODE 1CM Code  ---------------------------------------------------------
  # Structural Model ---
  ode <- "
      k10 = CL1/V1;
      d/dt(centr) = -k10*centr;
      C1=centr/V1;
      PRED = C1
      "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # default values ----------------------------------------------------------
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.add)) sigma.add <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(method)) method <- "bobyqa"
  if (is.null(nlmixrOutput)) nlmixrOutput <- F

  # error model -------------------------------------------------------------
  # Define error terms ---
  eps.add <- rnorm(length(output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.add

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  et$dv <- c(0, data$cp)
  et$cp <- c(0, data$cp)


  # dynmodel ----------------------------------------------------------------
  control <- dynmodelControl(method = method, nlmixrOutput = nlmixrOutput, covMethod = covMethod) # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, data = et, inits = inits, control = control))

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10))
  )
  # yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working Tests
# Check proc time outside of the function proc.time()
additiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa")
additiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "PORT")
additiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "Nelder-Mead")
additiveProportionalErrorTest(print.results = T, print.plot = T, sigma.prop = 0, seed = 1, method = "bobyqa", covMethod = "optimHess")

additiveProportionalErrorTest(print.results = T, print.plot = T, sigma.add = 0, seed = 1)
additiveProportionalErrorTest(print.results = T, print.plot = T, sigma.add = 0, sigma.prop = 0, seed = 1)

# #########################################################################

# Power Error Model Test --------------------------------------------------
powerErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, pow.tol = NULL, pow2.tol = NULL,
                           sigma.add = NULL, sigma.pow = NULL, pow2 = NULL,
                           print.results = NULL, print.plot = NULL, seed = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
      k10 = CL1/V1;
      d/dt(centr) = -k10*centr;
      C1=centr/V1;
      PRED = C1
      "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25

  if (is.null(add.tol)) add.tol <- 0.75
  if (is.null(pow.tol)) pow.tol <- 0.75
  if (is.null(pow2.tol)) pow2.tol <- 0.75

  if (is.null(sigma.add)) sigma.add <- 0
  if (is.null(sigma.pow)) sigma.pow <- 0.01
  if (is.null(pow2)) pow2 <- 0.75

  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE

  if (!is.null(seed)) {
    set.seed(seed)
  }

  graphics.off()

  # Define error terms ---
  eps.add <- if (sigma.pow == 0) {
    0
  } else {
    rnorm(length(output$C1), 0, sigma.add)
  }

  eps.pow <- rnorm(length(output$C1), 0, sigma.pow)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.pow > 0) c(sigma, pow = sigma.pow) else c(sigma)
  sigma <- if (pow2 > 0) c(sigma, pow2 = pow2) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + eps.pow * .F^(pow2) + eps.add

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }
  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 2, V1 = 75)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "pow") & any(names(sigma) %in% "pow2")) {
    error.model <- cp ~ C1 + add(0.01) + pow(0.01) + pow2(0.75)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, pow = pow.tol, pow2 = pow2.tol)
  } else if (any(names(sigma) %in% "pow") & any(names(sigma) %in% "pow2")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, pow = pow.tol, pow2 = pow2.tol)
  } else if (any(names(sigma) %in% "pow2")) {
    error.model <- cp ~ C1 + pow2(0.75)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, pow2 = pow2.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = "Nelder-Mead") # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT
  # control = dynmodelControl(method = "Nelder-Mead", normType = "constant", scaleType = "norm")

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10))
  )
  # yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })


  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Power Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# #########################################################################

# logn Error Model Test ---------------------------------------------------
lognErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                          sigma.add = NULL, sigma.prop = NULL,
                          print.results = NULL, print.plot = NULL, seed = NULL,
                          pow.tol = NULL, pow2.tol = NULL, sigma.pow = NULL, pow2 = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
      k10 = CL1/V1;
      d/dt(centr) = -k10*centr;
      C1=centr/V1;
      PRED = C1
      "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.5
  if (is.null(prop.tol)) prop.tol <- 0.5

  if (is.null(sigma.add)) sigma.add <- 0.025
  if (is.null(sigma.prop)) sigma.prop <- 0.025

  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE

  if (is.null(pow.tol)) pow.tol <- 0.3
  if (is.null(pow2.tol)) pow2.tol <- 0.3
  if (is.null(sigma.pow)) sigma.pow <- 0
  if (is.null(pow2)) pow2 <- 1

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # load library for BoxCox transform and inverse
  ## library(forecast)

  # Error model ---
  .sim.output <- output
  .output <- data.frame(time = .sim.output$time)

  # Define error terms ---

  # resample make above -1/lambda
  eps.add <- rnorm(length(.sim.output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(.sim.output$C1), 0, sigma.prop)
  eps.pow <- rnorm(length(.sim.output$C1), 0, sigma.pow)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)
  sigma <- if (sigma.pow > 0) c(sigma, pow = sigma.pow, pow2 = pow2) else c(sigma)

  # transform observed data
  h.y <- boxCox(.sim.output$C1, 0)

  # add error on transformed data
  if (sigma.pow > 0) {
    h.y.error <- h.y + eps.pow * .sim.output$C1^(pow2) + eps.add
  } else {
    h.y.error <- h.y + eps.prop * .sim.output$C1 + eps.add
  }

  # inverse the transform
  .output$PRED <- iBoxCox(h.y.error, 0)

  # resample instead of deleting (optional)
  if (any(.output$PRED < 0)) stop("negative observations not permitted with log normal testing")
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }
  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  .sim.output <- if (!is.null(r)) .sim.output[-c(r), ] else .sim.output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01) + logn(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop") & !any(names(sigma) %in% "pow")) {
    error.model <- cp ~ C1 + add(0.01) + logn(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01) + logn(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "pow") & !any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + logn(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, pow = pow.tol, pow2 = pow2.tol)
  } else if (any(names(sigma) %in% "pow") & any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + add(0.01) + logn(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, pow = pow.tol, pow2 = pow2.tol, add = sigma.add)
  }

  et <- eventTable()
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = "bobyqa") # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = .sim.output$time, Y = .sim.output$C1, Type = rep("Errorless  Sim", length(.sim.output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g")
  ) # , scales = list(y = list(log = 10)), yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = .sim.output$time, RES = .sim.output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = .sim.output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)


  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)


  test_that("log Normal, Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(rel.diff[i] < ref.tol[i], paste(
        names(rel.diff[i]),
        "Estimation out of tolerance range: +/-", ref.tol[i],
        "(", rel.diff[i], ")"
      ))
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working tests
# additive error test
lognErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, seed = 1)

# additive + proportional error test
lognErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, seed = 1)

# proportional error test
lognErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, seed = 1) # fails with seed=4

# power model
lognErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, sigma.pow = 0.025, pow2 = 2, seed = 1)

# #########################################################################

# dlnorm Error Model Test ---------------------------------------------------
dlnormErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                            sigma.add = NULL, sigma.prop = NULL,
                            print.results = NULL, print.plot = NULL, seed = NULL,
                            pow.tol = NULL, pow2.tol = NULL, sigma.pow = NULL, pow2 = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
      k10 = CL1/V1;
      d/dt(centr) = -k10*centr;
      C1=centr/V1;
      PRED = C1
      "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.5
  if (is.null(prop.tol)) prop.tol <- 0.5

  if (is.null(sigma.add)) sigma.add <- 0.025
  if (is.null(sigma.prop)) sigma.prop <- 0.025

  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE

  if (is.null(pow.tol)) pow.tol <- 0.3
  if (is.null(pow2.tol)) pow2.tol <- 0.3
  if (is.null(sigma.pow)) sigma.pow <- 0
  if (is.null(pow2)) pow2 <- 1

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # load library for BoxCox transform and inverse

  # Error model ---
  .sim.output <- output
  .output <- data.frame(time = .sim.output$time)

  # Define error terms ---

  # resample make above -1/lambda
  eps.add <- rnorm(length(.sim.output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(.sim.output$C1), 0, sigma.prop)
  eps.pow <- rnorm(length(.sim.output$C1), 0, sigma.pow)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)
  sigma <- if (sigma.pow > 0) c(sigma, pow = sigma.pow, pow2 = pow2) else c(sigma)

  # transform observed data
  h.y <- boxCox(.sim.output$C1, 0)

  # add error on transformed data
  if (sigma.pow > 0) {
    h.y.error <- h.y + eps.pow * .sim.output$C1^(pow2) + eps.add
  } else {
    h.y.error <- h.y + eps.prop * .sim.output$C1 + eps.add
  }

  # inverse the transform
  .output$PRED <- iBoxCox(h.y.error, 0)

  # resample instead of deleting (optional)
  if (any(.output$PRED < 0)) stop("negative observations not permitted with log normal testing")
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }
  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  .sim.output <- if (!is.null(r)) .sim.output[-c(r), ] else .sim.output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01) + dlnorm(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop") & !any(names(sigma) %in% "pow")) {
    error.model <- cp ~ C1 + add(0.01) + dlnorm(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01) + dlnorm(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "pow") & !any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + dlnorm(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, pow = pow.tol, pow2 = pow2.tol)
  } else if (any(names(sigma) %in% "pow") & any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + add(0.01) + dlnorm(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, pow = pow.tol, pow2 = pow2.tol, add = sigma.add)
  }

  et <- eventTable()
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = "bobyqa") # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = .sim.output$time, Y = .sim.output$C1, Type = rep("Errorless  Sim", length(.sim.output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g")
  ) # , scales = list(y = list(log = 10)), yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = .sim.output$time, RES = .sim.output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = .sim.output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)


  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)


  test_that("lol Normal (dlnorm), Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(rel.diff[i] < ref.tol[i], paste(
        names(rel.diff[i]),
        "Estimation out of tolerance range: +/-", ref.tol[i],
        "(", rel.diff[i], ")"
      ))
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working tests
# additive error test
dlnormErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, seed = 1)

# additive + proportional error test
dlnormErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, seed = 1)

# proportional error test
dlnormErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, seed = 1) # fails with seed=4

# power model
dlnormErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, sigma.pow = 0.025, pow2 = 2, seed = 1)

# #########################################################################

# norm Error Model Test ---------------------------------------------------
normErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, norm.tol = NULL, prop.tol = NULL,
                          sigma.norm = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
    k10 = CL1/V1;
    d/dt(centr) = -k10*centr;
    C1=centr/V1;
    PRED = C1
    "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  library(gridExtra)

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(norm.tol)) norm.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.norm)) sigma.norm <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }

  graphics.off()

  # Define error terms ---
  eps.norm <- rnorm(length(output$C1), 0, sigma.norm)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.norm > 0) c(sigma, norm = sigma.norm) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.norm

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "norm") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + norm(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, norm = norm.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "norm") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + norm(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, norm = norm.tol)
  } else if (!any(names(sigma) %in% "norm") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + norm(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = "Nelder-Mead") # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Testing ---

  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  if (any("prop" %in% names(sigma)) & any("norm" %in% names(sigma))) {
    ref.fit.pars <- c(ref.fit.pars[1], ref.fit.pars[2], ref.fit.pars[4], ref.fit.pars[3])
  }
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10)),
    yscale.components = yscale.components.log10ticks
  )

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Norm Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      if (any("prop" %in% names(sigma)) & any("norm" %in% names(sigma))) {
        fit$res <- cbind(actual = c(parameters, c(sigma[2], sigma[1])), fit$res)
      } else {
        fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      }
      print(fit$res)
    }
  })
}
# Working Tests
normErrorTest(print.results = T, print.plot = T, seed = 1)
normErrorTest(print.results = T, print.plot = T, sigma.prop = 0, seed = 1)
normErrorTest(print.results = T, print.plot = T, sigma.norm = 0, seed = 1)

# #########################################################################

# dnorm Error Model Test --------------------------------------------------
dnormErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, dnorm.tol = NULL, prop.tol = NULL,
                           sigma.dnorm = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
    k10 = CL1/V1;
    d/dt(centr) = -k10*centr;
    C1=centr/V1;
    PRED = C1
    "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  library(gridExtra)

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(dnorm.tol)) dnorm.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.dnorm)) sigma.dnorm <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }

  graphics.off()

  # Define error terms ---
  eps.dnorm <- rnorm(length(output$C1), 0, sigma.dnorm)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.dnorm > 0) c(sigma, dnorm = sigma.dnorm) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.dnorm

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "dnorm") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + dnorm(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, dnorm = dnorm.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "dnorm") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + dnorm(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, dnorm = dnorm.tol)
  } else if (!any(names(sigma) %in% "dnorm") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + dnorm(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = "Nelder-Mead") # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Testing ---

  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  if (any("prop" %in% names(sigma)) & any("dnorm" %in% names(sigma))) {
    ref.fit.pars <- c(ref.fit.pars[1], ref.fit.pars[2], ref.fit.pars[4], ref.fit.pars[3])
  }
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10)),
    yscale.components = yscale.components.log10ticks
  )

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Norm Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      if (any("prop" %in% names(sigma)) & any("dnorm" %in% names(sigma))) {
        fit$res <- cbind(actual = c(parameters, c(sigma[2], sigma[1])), fit$res)
      } else {
        fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      }
      print(fit$res)
    }
  })
}
# Working Tests
dnormErrorTest(print.results = T, print.plot = T, seed = 1)
dnormErrorTest(print.results = T, print.plot = T, sigma.prop = 0, seed = 1)
dnormErrorTest(print.results = T, print.plot = T, sigma.dnorm = 0, seed = 1)
# #########################################################################

# Optimzation Test --------------------------------------------------
optimizationAdditiveProportionalErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                                                      sigma.add = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL, method = NULL, nlmixrOutput = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
    k10 = CL1/V1;
    d/dt(centr) = -k10*centr;
    C1=centr/V1;
    PRED = C1
    "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  library(gridExtra)
  library(lattice)

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.add)) sigma.add <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(method)) method <- "bobyqa"
  if (is.null(nlmixrOutput)) nlmixrOutput <- F

  graphics.off()

  # Define error terms ---
  eps.add <- rnorm(length(output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.add

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  et$dv <- c(0, data$cp)
  et$cp <- c(0, data$cp)

  control <- dynmodelControl(method = method, nlmixrOutput = nlmixrOutput) # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, data = et, inits = inits, control = control))

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10))
  )
  # yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
      print(fit$value)
    }
  })
}

# Working Tests
optimizationAdditiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "Nelder-Mead")
optimizationAdditiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa")
optimizationAdditiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "lbfgsb3c")
optimizationAdditiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "PORT")

# #########################################################################

# boxCox Error Model Test -------------------------------------------------
boxCoxErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL, lambda.tol = NULL,
                            sigma.add = NULL, sigma.prop = NULL, lambda = NULL,
                            print.results = NULL, print.plot = NULL, seed = NULL,
                            pow.tol = NULL, pow2.tol = NULL, sigma.pow = NULL, pow2 = NULL, method = NULL) {
  library(gridExtra)

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
    k10 = CL1/V1;
    d/dt(centr) = -k10*centr;
    C1=centr/V1;
    PRED = C1
    "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.5
  if (is.null(prop.tol)) prop.tol <- 0.5
  if (is.null(lambda.tol)) lambda.tol <- 0.5

  if (is.null(sigma.add)) sigma.add <- 0 # 0.025
  if (is.null(sigma.prop)) sigma.prop <- 0 # 0.025
  if (is.null(lambda)) lambda <- 0

  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE

  if (is.null(pow.tol)) pow.tol <- 0.25
  if (is.null(pow2.tol)) pow2.tol <- 0.25
  if (is.null(sigma.pow)) sigma.pow <- 0 # 0.01
  if (is.null(pow2)) pow2 <- 1

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(method)) method <- "bobyqa"

  # load library for BoxCox transform and inverse
  library(forecast)

  # Error model ---
  .sim.output <- output
  .output <- data.frame(time = .sim.output$time)

  # Define error terms ---

  # resample make above -1/lambda
  eps.add <- rnorm(length(.sim.output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(.sim.output$C1), 0, sigma.prop)
  eps.pow <- rnorm(length(.sim.output$C1), 0, sigma.pow)



  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)
  sigma <- if (sigma.pow > 0) c(sigma, pow = sigma.pow, pow2 = pow2) else c(sigma)
  sigma <- if (lambda != 0 & !is.null(sigma)) c(sigma, boxCox = lambda) else c(sigma)

  # transform observed data
  h.y <- boxCox(.sim.output$C1, lambda)

  # add error on transformed data
  if (sigma.pow > 0 & sigma.prop == 0) {
    h.y.error <- h.y + eps.pow * .sim.output$C1^(pow2) + eps.add
  } else {
    h.y.error <- h.y + eps.prop * .sim.output$C1 + eps.add
  }

  # inverse the transform
  .output$PRED <- iBoxCox(h.y.error, lambda)

  # resample instead of deleting (optional)
  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }
  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  .sim.output <- if (!is.null(r)) .sim.output[-c(r), ] else .sim.output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  if (!any(names(sigma) %in% "boxCox") & any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01) + boxCox(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol, lambda = lambda.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop") & !any(names(sigma) %in% "pow")) {
    error.model <- cp ~ C1 + add(0.01) + boxCox(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, lambda = lambda.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01) + boxCox(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol, lambda = lambda.tol)
  } else if (any(names(sigma) %in% "pow") & !any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + boxCox(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol)
  } else if (any(names(sigma) %in% "pow") & any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + add(0.01) + boxCox(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol, add = sigma.add)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = method) # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = .sim.output$time, Y = .sim.output$C1, Type = rep("Errorless  Sim", length(.sim.output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g")
  ) # , scales = list(y = list(log = 10)), yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = .sim.output$time, RES = .sim.output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = .sim.output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  test_that("boxCox, Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(rel.diff[i] < ref.tol[i], paste(
        names(rel.diff[i]),
        "Estimation out of tolerance range: +/-", ref.tol[i],
        "(", rel.diff[i], ")"
      ))
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working tests
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, sigma.pow = .125, pow2 = 0.75, lambda = 1.2, seed = 1)

boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, sigma.pow = .125, pow2 = 0.75, lambda = 1.2, seed = 5) # 4 parameters may be too many to estimate for the error model

boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = 1, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = 1, seed = 2)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = 1, seed = 5)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = 1, seed = 1)

boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = 1.2, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = 1.2, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = 1.2, seed = 2)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = 1.2, seed = 1)

boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = 0.75, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = 0.75, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = 0.75, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = 0.75, seed = 2) # fails: 1, 3, 4, 6, 7, 9, 10

boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = -1.2, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = -1.2, seed = 2)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = -1.2, seed = 7) # Hessian
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = -1.2, seed = 1) # Hessian

boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = -0.75, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = -0.75, seed = 1)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = -0.75, seed = 2)
boxCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = -0.75, seed = 5) # fails: 1, 3, 4, 6, 7, 9, 10

# #########################################################################

# tbs Error Model Test ----------------------------------------------------
tbsErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL, lambda.tol = NULL,
                         sigma.add = NULL, sigma.prop = NULL, lambda = NULL,
                         print.results = NULL, print.plot = NULL, seed = NULL,
                         pow.tol = NULL, pow2.tol = NULL, sigma.pow = NULL, pow2 = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
    k10 = CL1/V1;
    d/dt(centr) = -k10*centr;
    C1=centr/V1;
    PRED = C1
    "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.5
  if (is.null(prop.tol)) prop.tol <- 0.5
  if (is.null(lambda.tol)) lambda.tol <- 0.5

  if (is.null(sigma.add)) sigma.add <- 0 # 0.025
  if (is.null(sigma.prop)) sigma.prop <- 0 # 0.025
  # if (is.null(lambda)) lambda = 0

  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE

  if (is.null(pow.tol)) pow.tol <- 0.25
  if (is.null(pow2.tol)) pow2.tol <- 0.25
  if (is.null(sigma.pow)) sigma.pow <- 0 # 0.01
  if (is.null(pow2)) pow2 <- 1

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # load library for BoxCox transform and inverse
  library(forecast)

  # Error model ---
  .sim.output <- output
  .output <- data.frame(time = .sim.output$time)

  # Define error terms ---

  # resample make above -1/lambda
  eps.add <- rnorm(length(.sim.output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(.sim.output$C1), 0, sigma.prop)
  eps.pow <- rnorm(length(.sim.output$C1), 0, sigma.pow)



  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)
  sigma <- if (sigma.pow > 0) c(sigma, pow = sigma.pow, pow2 = pow2) else c(sigma)
  sigma <- if (lambda != 0 & !is.null(sigma)) c(sigma, tbs = lambda) else c(sigma)

  # transform observed data
  h.y <- boxCox(.sim.output$C1, lambda)

  # add error on transformed data
  if (sigma.pow > 0 & sigma.prop == 0) {
    h.y.error <- h.y + eps.pow * .sim.output$C1^(pow2) + eps.add
  } else {
    h.y.error <- h.y + eps.prop * .sim.output$C1 + eps.add
  }

  # inverse the transform
  .output$PRED <- iBoxCox(h.y.error, lambda)

  # resample instead of deleting (optional)
  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }
  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  .sim.output <- if (!is.null(r)) .sim.output[-c(r), ] else .sim.output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  if (!any(names(sigma) %in% "tbs") & any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01) + tbs(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol, lambda = lambda.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop") & !any(names(sigma) %in% "pow")) {
    error.model <- cp ~ C1 + add(0.01) + tbs(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, lambda = lambda.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01) + tbs(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol, lambda = lambda.tol)
  } else if (any(names(sigma) %in% "pow") & !any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + tbs(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol)
  } else if (any(names(sigma) %in% "pow") & any(names(sigma) %in% "add")) {
    error.model <- cp ~ C1 + pow(0.01) + pow2(0.75) + add(0.01) + tbs(2)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol, add = sigma.add)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = "bobyqa") # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = .sim.output$time, Y = .sim.output$C1, Type = rep("Errorless  Sim", length(.sim.output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g")
  ) # , scales = list(y = list(log = 10)), yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = .sim.output$time, RES = .sim.output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = .sim.output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  test_that("boxCox, Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(rel.diff[i] < ref.tol[i], paste(
        names(rel.diff[i]),
        "Estimation out of tolerance range: +/-", ref.tol[i],
        "(", rel.diff[i], ")"
      ))
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working tests
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, sigma.pow = .125, pow2 = 0.75, lambda = 1.2, seed = 1)

tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, sigma.pow = .125, pow2 = 0.75, lambda = 1.2, seed = 5) # 4 parameters may be too many to estimate for the error model

tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = 1, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = 1, seed = 2)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = 1, seed = 5)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = 1, seed = 1)

tbsCoxErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = 1.2, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = 1.2, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = 1.2, seed = 2)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = 1.2, seed = 1)

tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = 0.75, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = 0.75, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = 0.75, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = 0.75, seed = 2) # fails: 1, 3, 4, 6, 7, 9, 10

tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = -1.2, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = -1.2, seed = 2)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = -1.2, seed = 7) # Hessian
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = -1.2, seed = 1) # Hessian

tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, lambda = -0.75, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0, lambda = -0.75, seed = 1)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0.025, lambda = -0.75, seed = 2)
tbsErrorTest(print.plot = T, print.results = T, sigma.add = 0.025, sigma.prop = 0.025, lambda = -0.75, seed = 5) # fails: 1, 3, 4, 6, 7, 9, 10

# #########################################################################

# Boundarries Test --------------------------------------------------
boudaryTestErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                                 sigma.add = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL, method = NULL,
                                 lower = NULL, upper = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
    k10 = CL1/V1;
    d/dt(centr) = -k10*centr;
    C1=centr/V1;
    PRED = C1
    "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  library(gridExtra)

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.add)) sigma.add <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(lower)) .lower <- -Inf else .lower <- lower
  if (is.null(upper)) .upper <- Inf else .upper <- upper

  graphics.off()

  # Define error terms ---
  eps.add <- rnorm(length(output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.add

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = method, lower = .lower, upper = .upper) # bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10)),
    yscale.components = yscale.components.log10ticks
  )

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Boundary and Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working Tests
boudaryTestErrorTest(print.results = T, print.plot = T, sigma.prop = 0.01, seed = 1, method = "bobyqa", lower = c(-5, -5, -5, -5))
# #########################################################################

# nlmixrDynmodelConvert Additive + Proportional ---------------------------
nlmixrDynmodelConvertErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                                           sigma.add = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  f <- function() {
    ini({
      iCL1 <- 3 # Cl (L/hr)
      iV1 <- 75 # Vc (L)
      prop.err <- c(0, 0.2, 1)
    })
    model({
      CL1 <- iCL1
      V1 <- iV1
      kel <- CL1 / V1
      d / dt(centr) <- -kel * centr
      C1 <- centr / V1
      C1 ~ prop(prop.err)
    })
  }

  nmf <- nlmixr(f)

  dynNlmixr <- nlmixrDynmodelConvert(nmf)

  model <- dynNlmixr$system

  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )

  parameters <- dynNlmixr$inits

  output <- as.data.frame(rxSolve(model, et, parameters))
  # #########################################################################

  library(gridExtra)

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.add)) sigma.add <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }

  graphics.off()

  # Define error terms ---
  eps.add <- rnorm(length(output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.add

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = "Nelder-Mead") # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$nlmixr_pred, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10)),
    yscale.components = yscale.components.log10ticks
  )

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working Tests
nlmixrDynmodelConvertErrorTest(print.results = T, print.plot = T, seed = 1)
nlmixrDynmodelConvertErrorTest(print.results = T, print.plot = T, sigma.prop = 0, seed = 1)
nlmixrDynmodelConvertErrorTest(print.results = T, print.plot = T, sigma.add = 0, seed = 1)
nlmixrDynmodelConvertErrorTest(print.results = T, print.plot = T, sigma.add = 0, sigma.prop = 0, seed = 1)

# #########################################################################

# Scaling Test ------------------------------
ScalingTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                        sigma.add = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL,
                        method = NULL, nlmixrOutput = NULL, control = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
    k10 = CL1/V1;
    d/dt(centr) = -k10*centr;
    C1=centr/V1;
    PRED = C1
    "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))

  # #########################################################################

  library(gridExtra)
  library(lattice)

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.add)) sigma.add <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(method)) method <- "bobyqa"
  if (is.null(nlmixrOutput)) nlmixrOutput <- F
  if (is.null(control)) dynmodelControl()

  graphics.off()

  # Define error terms ---
  eps.add <- rnorm(length(output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.add

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  et$dv <- c(0, data$cp)
  et$cp <- c(0, data$cp)

  control$method <- method
  control$nlmixrOutput <- nlmixrOutput
  assign("data", et, envir = .GlobalEnv)
  (fit <- dynmodel(system = model, model = error.model, data = et, inits = inits, control = control))

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10))
  )
  # yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

#
# normType=c("constant","rescale2", "mean", "rescale", "std", "len"),
# scaleType=c("norm","nlmixr", "mult", "multAdd"),
# scaleCmax=1e5,
# scaleCmin=1e-5,
# scaleC=NULL,
# scaleC0=1e5,

# Notes
# normType - Feature scaling all parameters
# scaleType - The scaling scheme for nlmixr
# scaleTo - Scale the initial parameter estimate to this value

control <- dynmodelControl(scaleType = "nlmixr", normType = "constant")
ScalingTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa", control = control)

control <- dynmodelControl(scaleType = "nlmixr", normType = "rescale2")
ScalingTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa", control = control)

control <- dynmodelControl(scaleType = "nlmixr", normType = "mean")
ScalingTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa", control = control)

control <- dynmodelControl(scaleType = "nlmixr", normType = "rescale")
ScalingTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa", control = control)

control <- dynmodelControl(scaleType = "nlmixr", normType = "std")
ScalingTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa", control = control)

control <- dynmodelControl(scaleType = "nlmixr", normType = "len")
ScalingTest(print.results = T, print.plot = T, seed = 1, method = "bobyqa", control = control)




## Not Working ##

# power model not working with yeoJohnson. May need to fix all parameters except error terms?

# #########################################################################

# yeoJohnson Error Model Test ---------------------------------------------
yeoJohnsonErrorTest <- function(kin.tol = NULL, kout.tol = NULL, Imax.tol = NULL, IC50.tol = NULL,
                                add.tol = NULL, prop.tol = NULL, lambda.tol = NULL,
                                sigma.add = NULL, sigma.prop = NULL, lambda = NULL,
                                print.results = NULL, print.plot = NULL, seed = NULL,
                                pow.tol = NULL, pow2.tol = NULL, sigma.pow = NULL, pow2 = NULL, method = NULL) {
  library(gridExtra)

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
        k10 = CL1/V1;
        d/dt(centr) = -k10*centr;
        C1=centr/V1;
        PRED = C1

        eff(0) = kin/kout
        d/dt(eff) = kin*(1-(C1*Imax)/(C1+IC50)) - kout*eff
        "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75, kin = 0.1, kout = 0.2, Imax = 5, IC50 = 2)
  output <- as.data.frame(rxSolve(model, et, parameters))
  # plot(output$time, output$eff)

  # #########################################################################

  # Set default tolerances and error model parameters
  if (is.null(kin.tol)) kin.tol <- 0.6
  if (is.null(kout.tol)) kout.tol <- 0.6
  if (is.null(Imax.tol)) Imax.tol <- 0.6
  if (is.null(IC50.tol)) IC50.tol <- 0.6

  if (is.null(add.tol)) add.tol <- 0.5
  if (is.null(prop.tol)) prop.tol <- 0.5
  if (is.null(lambda.tol)) lambda.tol <- 0.5

  if (is.null(sigma.add)) sigma.add <- 0 # 0.025
  if (is.null(sigma.prop)) sigma.prop <- 0 # 0.025
  if (is.null(lambda)) lambda <- 0

  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE

  if (is.null(pow.tol)) pow.tol <- 0.25
  if (is.null(pow2.tol)) pow2.tol <- 0.25
  if (is.null(sigma.pow)) sigma.pow <- 0 # 0.01
  if (is.null(pow2)) pow2 <- 1

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(method)) method <- "bobyqa"

  # Error model ---
  .sim.output <- output
  .output <- data.frame(time = .sim.output$time)

  # Define error terms ---

  # resample make above -1/lambda
  eps.add <- rnorm(length(.sim.output$eff), 0, sigma.add)
  eps.prop <- rnorm(length(.sim.output$eff), 0, sigma.prop)
  eps.pow <- rnorm(length(.sim.output$eff), 0, sigma.pow)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)
  sigma <- if (sigma.pow > 0) c(sigma, pow = sigma.pow, pow2 = pow2) else c(sigma)
  sigma <- if (lambda != 0 & !is.null(sigma)) c(sigma, yeoJohnson = lambda) else c(sigma)

  # transform observed data
  h.y <- yeoJohnson(.sim.output$eff, lambda)

  # add error on transformed data
  if (sigma.pow > 0 & sigma.prop == 0) {
    h.y.error <- h.y + eps.pow * .sim.output$eff^(pow2) + eps.add
  } else {
    h.y.error <- h.y + eps.prop * .sim.output$eff + eps.add
  }

  # inverse the transform
  .output$PRED <- iYeoJohnson(h.y.error, lambda)

  # Parameter Estimation ---
  inits <- c(kin = 0.1, kout = 0.05, Imax = 1, IC50 = 1)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  if (!any(names(sigma) %in% "yeoJohnson") & any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ eff + add(0.01) + prop(0.01)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      add = add.tol, prop = prop.tol
    )
  } else if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ eff + add(0.01) + prop(0.01) + yeoJohnson(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      prop = prop.tol, lambda = lambda.tol
    )
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop") & !any(names(sigma) %in% "pow")) {
    error.model <- cp ~ eff + add(0.01) + yeoJohnson(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      lambda = lambda.tol
    )
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ eff + prop(0.01) + yeoJohnson(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      prop = prop.tol, lambda = lambda.tol
    )
  } else if (any(names(sigma) %in% "pow") & !any(names(sigma) %in% "add")) {
    error.model <- cp ~ eff + pow(0.01) + pow2(0.75) + yeoJohnson(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol
    )
  } else if (any(names(sigma) %in% "pow") & any(names(sigma) %in% "add")) {
    error.model <- cp ~ eff + pow(0.01) + pow2(0.75) + add(0.01) + yeoJohnson(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol, add = sigma.add
    )
  } else {
    error.model <- cp ~ eff
    ref.tol <- c(kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = method, fixPars = c(CL1 = 3, V1 = 75)) # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = .sim.output$time, Y = .sim.output$eff, Type = rep("Errorless  Sim", length(.sim.output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))

  sim.parameters <- c(3, 75, fit$res[c(1, 2, 3, 4)])
  names(sim.parameters) <- c("CL1", "V1", "kin", "kout", "Imax", "IC50")
  fit.data <- as.data.frame(rxSolve(model, et, sim.parameters))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$eff, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g")
  ) # , scales = list(y = list(log = 10)), yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = .sim.output$time, RES = .sim.output$eff - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = .sim.output$eff)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters[-c(1, 2)], sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  test_that("yeoJohnson Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(rel.diff[i] < ref.tol[i], paste(
        names(rel.diff[i]),
        "Estimation out of tolerance range: +/-", ref.tol[i],
        "(", rel.diff[i], ")"
      ))
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters[-c(1, 2)], sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working Tests
yeoJohnsonErrorTest(print.results = T, print.plot = T, sigma.add = 0, sigma.prop = 0, seed = 1)
yeoJohnsonErrorTest(print.results = T, print.plot = T, sigma.add = 0.025, sigma.prop = 0, lambda = 1, seed = 1)
yeoJohnsonErrorTest(print.results = T, print.plot = T, sigma.add = 0.025, sigma.prop = 0, lambda = 0.75, seed = 1)
yeoJohnsonErrorTest(print.results = T, print.plot = T, sigma.add = 0.45, sigma.prop = 0.5, lambda = 1.4, seed = 7)

# Non-working Tests
yeoJohnsonErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, sigma.pow = 1, pow2 = 0.75, lambda = 1, seed = 5) # 4 parameters may be too many to estimate for the error model


# #########################################################################

# tbsYJ Error Model Test --------------------------------------------------
tbsYjErrorTest <- function(kin.tol = NULL, kout.tol = NULL, Imax.tol = NULL, IC50.tol = NULL,
                           add.tol = NULL, prop.tol = NULL, lambda.tol = NULL,
                           sigma.add = NULL, sigma.prop = NULL, lambda = NULL,
                           print.results = NULL, print.plot = NULL, seed = NULL,
                           pow.tol = NULL, pow2.tol = NULL, sigma.pow = NULL, pow2 = NULL, method = NULL) {
  library(gridExtra)

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
      k10 = CL1/V1;
      d/dt(centr) = -k10*centr;
      C1=centr/V1;
      PRED = C1

      eff(0) = kin/kout
      d/dt(eff) = kin*(1-(C1*Imax)/(C1+IC50)) - kout*eff
      "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75, kin = 0.1, kout = 0.2, Imax = 5, IC50 = 2)
  output <- as.data.frame(rxSolve(model, et, parameters))
  # plot(output$time, output$eff)

  # #########################################################################

  # Set default tolerances and error model parameters
  if (is.null(kin.tol)) kin.tol <- 0.6
  if (is.null(kout.tol)) kout.tol <- 0.6
  if (is.null(Imax.tol)) Imax.tol <- 0.6
  if (is.null(IC50.tol)) IC50.tol <- 0.6

  if (is.null(add.tol)) add.tol <- 0.5
  if (is.null(prop.tol)) prop.tol <- 0.5
  if (is.null(lambda.tol)) lambda.tol <- 0.5

  if (is.null(sigma.add)) sigma.add <- 0 # 0.025
  if (is.null(sigma.prop)) sigma.prop <- 0 # 0.025
  if (is.null(lambda)) lambda <- 0

  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE

  if (is.null(pow.tol)) pow.tol <- 0.25
  if (is.null(pow2.tol)) pow2.tol <- 0.25
  if (is.null(sigma.pow)) sigma.pow <- 0 # 0.01
  if (is.null(pow2)) pow2 <- 1

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(method)) method <- "bobyqa"

  # Error model ---
  .sim.output <- output
  .output <- data.frame(time = .sim.output$time)

  # Define error terms ---

  # resample make above -1/lambda
  eps.add <- rnorm(length(.sim.output$eff), 0, sigma.add)
  eps.prop <- rnorm(length(.sim.output$eff), 0, sigma.prop)
  eps.pow <- rnorm(length(.sim.output$eff), 0, sigma.pow)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)
  sigma <- if (sigma.pow > 0) c(sigma, pow = sigma.pow, pow2 = pow2) else c(sigma)
  sigma <- if (lambda != 0 & !is.null(sigma)) c(sigma, tbsYj = lambda) else c(sigma)

  print(sigma)

  # transform observed data
  h.y <- yeoJohnson(.sim.output$eff, lambda)

  # add error on transformed data
  if (sigma.pow > 0 & sigma.prop == 0) {
    h.y.error <- h.y + eps.pow * .sim.output$eff^(pow2) + eps.add
  } else {
    h.y.error <- h.y + eps.prop * .sim.output$eff + eps.add
  }

  # inverse the transform
  .output$PRED <- iYeoJohnson(h.y.error, lambda)

  # Parameter Estimation ---
  inits <- c(kin = 0.1, kout = 0.05, Imax = 1, IC50 = 1)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  if (!any(names(sigma) %in% "tbsYj") & any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ eff + add(0.01) + prop(0.01)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      add = add.tol, prop = prop.tol
    )
  } else if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ eff + add(0.01) + prop(0.01) + tbsYj(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      prop = prop.tol, lambda = lambda.tol
    )
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop") & !any(names(sigma) %in% "pow")) {
    error.model <- cp ~ eff + add(0.01) + tbsYj(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      lambda = lambda.tol
    )
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ eff + prop(0.01) + tbsYj(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      prop = prop.tol, lambda = lambda.tol
    )
  } else if (any(names(sigma) %in% "pow") & !any(names(sigma) %in% "add")) {
    error.model <- cp ~ eff + pow(0.01) + pow2(0.75) + tbsYj(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol
    )
  } else if (any(names(sigma) %in% "pow") & any(names(sigma) %in% "add")) {
    error.model <- cp ~ eff + pow(0.01) + pow2(0.75) + add(0.01) + tbsYj(2)
    ref.tol <- c(
      kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol,
      lambda = lambda.tol, pow = pow.tol, pow2 = pow2.tol, add = sigma.add
    )
  } else {
    error.model <- cp ~ eff
    ref.tol <- c(kin = kin.tol, kout = kout.tol, Imax = Imax.tol, IC50 = IC50.tol)
  }

  print(error.model)

  et <- eventTable()
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = data$time
  )

  control <- dynmodelControl(method = method, fixPars = c(CL1 = 3, V1 = 75)) # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, evTable = et, inits = inits, data = data, control = control))

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = .sim.output$time, Y = .sim.output$eff, Type = rep("Errorless  Sim", length(.sim.output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))

  sim.parameters <- c(3, 75, fit$res[c(1, 2, 3, 4)])
  names(sim.parameters) <- c("CL1", "V1", "kin", "kout", "Imax", "IC50")
  fit.data <- as.data.frame(rxSolve(model, et, sim.parameters))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$eff, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g")
  ) # , scales = list(y = list(log = 10)), yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = .sim.output$time, RES = .sim.output$eff - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = .sim.output$eff)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters[-c(1, 2)], sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  test_that("tbdYJ Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(rel.diff[i] < ref.tol[i], paste(
        names(rel.diff[i]),
        "Estimation out of tolerance range: +/-", ref.tol[i],
        "(", rel.diff[i], ")"
      ))
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters[-c(1, 2)], sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working Tests
tbsYjErrorTest(print.results = T, print.plot = T, sigma.add = 0, sigma.prop = 0, seed = 1)
tbsYjErrorTest(print.results = T, print.plot = T, sigma.add = 0.025, sigma.prop = 0, lambda = 1, seed = 1)
yeoJohnsonErrorTest(print.results = T, print.plot = T, sigma.add = 0.025, sigma.prop = 0, lambda = 0.75, seed = 1)
yeoJohnsonErrorTest(print.results = T, print.plot = T, sigma.add = 0.45, sigma.prop = 0.5, lambda = 1.4, seed = 7)

# Non-working Tests
yeoJohnsonErrorTest(print.plot = T, print.results = T, sigma.add = 0, sigma.prop = 0, sigma.pow = 1, pow2 = 0.75, lambda = 1, seed = 5) # 4 parameters may be too many to estimate for the error model


# #########################################################################


# nlmixr input ------------------------------------------------------------
# Additive and Proportional Error Model Test ------------------------------
additiveProportionalErrorTest <- function(CL1.tol = NULL, V1.tol = NULL, add.tol = NULL, prop.tol = NULL,
                                          sigma.add = NULL, sigma.prop = NULL, print.results = NULL, print.plot = NULL, seed = NULL,
                                          method = NULL, covMethod = "nlmixrHess", nlmixrOutput = NULL) {

  # RxODE 1CM Code --------------------------------------------------------------
  # Structural Model ---
  ode <- "
      k10 = CL1/V1;
      d/dt(centr) = -k10*centr;
      C1=centr/V1;
      PRED = C1
      "
  model <- RxODE(model = ode)
  et <- eventTable()
  dose <- 250
  et$add.dosing(
    dose = dose
  )
  et$add.sampling(
    time = 1:96
  )
  parameters <- c(CL1 = 3, V1 = 75)
  output <- as.data.frame(rxSolve(model, et, parameters))


  # nlmixr 1CM Code ---------------------------------------------------------
  one.compartment.IV.model.solve <- function() {
    ini({ # Where initial conditions/variables are specified
      # '<-' or '=' defines population parameters
      # Simple numeric expressions are supported
      Cl <- 1.6 # log Cl (L/hr)
      Vc <- 4.5 # log V (L)
      # Bounds may be specified by c(lower, est, upper), like NONMEM:
      # Residuals errors are assumed to be population parameters
      prop.err <- c(0, 0.3, 1)
      # Between subject variability estimates are specified by '~'
      # Semicolons are optional
      # eta.Vc ~ 0.1   #IIV V
      # eta.Cl ~ 0.1   #IIV Cl
    })
    model({ # Where the model is specified
      # The model uses the ini-defined variable names
      # Vc <- exp(lVc + eta.Vc)
      # Cl <- exp(lCl + eta.Cl)
      linCmt() ~ prop(prop.err)
    })
  }





  # #########################################################################

  library(gridExtra)
  library(lattice)
  library(testthat)

  # Set default tolerances and error model parameters
  if (is.null(CL1.tol)) CL1.tol <- 0.25
  if (is.null(V1.tol)) V1.tol <- 0.25
  if (is.null(add.tol)) add.tol <- 0.75
  if (is.null(prop.tol)) prop.tol <- 0.75
  if (is.null(sigma.add)) sigma.add <- 0.02
  if (is.null(sigma.prop)) sigma.prop <- 0.1
  if (is.null(print.results)) print.results <- FALSE
  if (is.null(print.plot)) print.plot <- FALSE
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(method)) method <- "bobyqa"
  if (is.null(nlmixrOutput)) nlmixrOutput <- F


  graphics.off()

  # Define error terms ---
  eps.add <- rnorm(length(output$C1), 0, sigma.add)
  eps.prop <- rnorm(length(output$C1), 0, sigma.prop)

  # combine sigma terms (for testing below)
  sigma <- c()
  sigma <- if (sigma.add > 0) c(sigma, add = sigma.add) else c(sigma)
  sigma <- if (sigma.prop > 0) c(sigma, prop = sigma.prop) else c(sigma)

  # Error model ---
  .output <- data.frame(time = output$time)
  .F <- output$C1
  .output$PRED <- .F + .F * eps.prop + eps.add

  # resample instead of deleting (optional)

  any(.output$PRED < 0)
  r <- c()
  for (i in 1:length(.output$PRED)) {
    if (.output$PRED[i] < 0) {
      r <- c(r, i)
    }
  }

  .output <- if (!is.null(r)) .output[-c(r), ] else .output
  output <- if (!is.null(r)) output[-c(r), ] else output

  # Parameter Estimation ---
  inits <- c(CL1 = 1, V1 = 10)
  data <- data.frame(time = .output$time, cp = .output$PRED)

  # Used to edit error model for combination of parameters
  if (any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01) + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol, prop = prop.tol)
  } else if (any(names(sigma) %in% "add") & !any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + add(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, add = add.tol)
  } else if (!any(names(sigma) %in% "add") & any(names(sigma) %in% "prop")) {
    error.model <- cp ~ C1 + prop(0.01)
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol, prop = prop.tol)
  } else {
    error.model <- cp ~ C1
    ref.tol <- c(CL1 = CL1.tol, V1 = V1.tol)
  }

  et <- eventTable()
  et$add.dosing(
    dose = 250
  )
  et$add.sampling(
    time = data$time
  )

  et$dv <- c(0, data$cp)
  et$cp <- c(0, data$cp)

  control <- dynmodelControl(method = method, nlmixrOutput = nlmixrOutput, covMethod = covMethod) # , lower=c(0,0,0), upper = c(5,100,1)) #bobyqa, Nelder-Mead, lbfgsb3c, PORT

  (fit <- dynmodel(system = model, model = error.model, data = et, inits = inits, control = control))

  # Testing ---
  test.fit.pars <- fit$res[, 1]
  ref.fit.pars <- c(parameters, sigma)
  rel.diff <- abs((test.fit.pars - ref.fit.pars) / ref.fit.pars)

  # Simulation and Estimation Plot ---
  errorless.data <- data.frame(Time = output$time, Y = output$C1, Type = rep("Errorless  Sim", length(output$time)))
  error.data <- data.frame(Time = .output$time, Y = .output$PRED, Type = rep("Error Sim", length(.output$time)))
  fit.data <- as.data.frame(rxSolve(model, et, fit$res[c(1, 2)]))
  fit.data <- data.frame(Time = fit.data$time, Y = fit.data$PRED, Type = rep("Error Fit", length(fit.data$time)))
  plot1.data <- rbind(errorless.data, error.data, fit.data)

  p1 <- xyplot(Y ~ Time,
    data = plot1.data, groups = factor(Type, labels = c("Errorless", "Error", "Fit")), main = "Simulated and Estimation",
    pch = 20, auto.key = list(columns = 3), type = c("p", "g"), scales = list(y = list(log = 10))
  )
  # yscale.components = yscale.components.log10ticks)

  # Residual Plot ---
  residual.data <- data.frame(Time = output$time, RES = output$C1 - .output$PRED)
  p2 <- xyplot(RES ~ Time, data = residual.data, main = "Residuals: Errorless - Error", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 0)
  })

  # Histogram Plot of Error vs. Non-error
  p3 <- histogram(~RES, data = residual.data, main = "Residuals: Errorless - Error", type = "density", breaks = dim(residual.data)[1])

  # Error vs. Non-error
  obs.err.plot <- data.frame(Error = .output$PRED, noError = output$C1)
  p4 <- xyplot(Error ~ noError, data = obs.err.plot, main = "Error vs. Errorless", xlab = "Errorless", panel = function(x, y) {
    panel.xyplot(x, y)
    panel.abline(a = 0, b = 1)
  })

  options(warn = -1) # mute warnings here
  if (print.plot) {
    library(latticeExtra)
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  options(warn = 0)

  test_that("Additive and Proportional Error Model Test", {
    for (i in 1:length(ref.tol)) {
      expect(
        rel.diff[i] < ref.tol[i],
        paste(names(rel.diff[i]), "Estimation out of tolerance range: +/-", ref.tol[i], "(", rel.diff[i], ")")
      )
    }
    if (print.results) {
      fit$res <- cbind(actual = c(parameters, sigma), fit$res)
      print(fit$res)
    }
  })
}

# Working Tests
# Check proc time outside of the function proc.time()
additiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "PORT")
additiveProportionalErrorTest(print.results = T, print.plot = T, seed = 1, method = "Nelder-Mead")
additiveProportionalErrorTest(print.results = T, print.plot = T, sigma.prop = 0, seed = 1, method = "bobyqa", covMethod = "optimHess")

additiveProportionalErrorTest(print.results = T, print.plot = T, sigma.add = 0, seed = 1)
additiveProportionalErrorTest(print.results = T, print.plot = T, sigma.add = 0, sigma.prop = 0, seed = 1)

# #########################################################################
