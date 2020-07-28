Ruser_function_cmt <- function(phi, evt) {
  evt.elem <- function(k) {
    ix.i <- ix.save[[k]]
    evt.i <- evt[ix.i, -1]
    evid <- evt.i[, 2]
    list(
      obs_time = evt.i[evid == 0, 1],
      dose_time = evt.i[evid > 0, 1],
      dose = evt.i[evid > 0, 3],
      tinf = evt.i[evid > 0, 4]
    )
  }

  get_pars <- function(i) {
    lCL <- phi[i, 1]
    lV <- phi[i, 2]
    lKA <- phi[i, 3]

    CL <- exp(lCL)
    V <- exp(lV)
    KA <- exp(lKA)
    c(CL, V, KA, 0)
  }

  ix.save <- tapply(1:nrow(evt), evt[, 1], identity)
  ncmt <- 1
  oral <- TRUE
  infusion <- FALSE
  parameterization <- 1
  s <- lapply(1:nrow(phi), function(k) {
    ev <- evt.elem(k)
    params <- get_pars(k)
    lin_cmt(ev$obs_time, ev$dose_time, ev$dose, ev$tinf, params, oral, infusion, ncmt, parameterization)
  })

  rm(.Random.seed, envir = .GlobalEnv) # important! otherwise, arma uses the current .Random.seed to set seed for random number generation in arma randu and randn
  unlist(s)
}


Ruser_function_ode <- function(phi, evt) {
  # phi = fit$mpost_phi; evt=cfg$evt
  evt.elem <- function(k) {
    ix.i <- ix.save[[k]]
    evt[ix.i, -1]
  }

  get_pars <- function(i) {
    lCL <- phi[i, 1]
    lV <- phi[i, 2]
    lKA <- phi[i, 3]

    CL <- exp(lCL)
    V <- exp(lV)
    KA <- exp(lKA)
    c(KA = KA, KE = CL / V, V = V) # FIXME
  }

  ix.save <- tapply(1:nrow(evt), evt[, 1], identity)

  s <- lapply(1:nrow(phi), function(k) {
    events <- evt.elem(k)
    params <- get_pars(k)

    neq <- 2
    nlhs <- 0
    atol <- rtol <- 1e-8
    stiff <- 1
    transit_abs <- 0
    rc <- as.integer(0)
    time <- events[, 1]
    evid <- events[, 2]
    amt <- events[evid > 0, 3]
    ntime <- length(time)
    ret <- double(neq * ntime)
    lhs <- double(nlhs * ntime)
    inits <- double(neq) # FIXME

    x <- .C(
      "RxODE_mod_m1_ode_solver",
      as.integer(neq), as.double(params), as.double(time),
      as.integer(evid), length(time), inits, as.double(amt),
      ret, as.double(atol), as.double(rtol), as.integer(stiff),
      as.integer(transit_abs), as.integer(nlhs), lhs, rc
    )
    x <- matrix(x[[8]], ncol = neq, byrow = T)
    x[evid == 0, 2] / params["V"]
  })

  .GlobalEnv$.Random.seed <- NULL # important! otherwise, arma uses the current .Random.seed to set seed for random number generation in arma randu and randn
  unlist(s)
}
