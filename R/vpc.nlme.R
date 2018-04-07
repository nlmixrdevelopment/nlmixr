vpc.nlme = function (fit, nsim = 100, by = NULL, ...) 
{
    if (class(fit) == "nlmixr.ui.nlme") fit = as.nlme(fit)
    dat = getData(fit)

    if (!is.null(by)) {
        if (by %in% names(dat)) {
            dat$grp = eval(parse(text = paste0("dat$", by)))
        }
        else {
            msg = paste0(by, " not found in data")
            stop(msg)
        }
    }
    else dat$grp = T

    xd = subset(dat, EVID == 0)
    nsub = length(unique(xd$ID))
    ntim = dim(xd)[1]
    ord = rep(1:ntim, nsim)
    sim = rep(1:nsim, each = ntim)

    nlmeModList(fit$env)
    on.exit({
        nlmeModList(new.env(parent = emptyenv()))
    })
    ..ModList <- nlmeModList()
    options(warn = -1)
    s = lapply(1:nsim, sim.one, x = fit)
    xs = do.call("cbind", s)

    df = cbind(xd[ord, c("ID", "TIME", "grp")], DV = as.vector(xs), SIM = sim)
    if (!is.null(by)) {
        p = vpc::vpc(sim = df, obs = dat, strat = c("grp"), facet = "wrap", ...)
    } else {
        p = vpc::vpc(sim = df, obs = dat, ...)
    }
    print(p)
    invisible(df)
}

vpc.nlme(fit)
