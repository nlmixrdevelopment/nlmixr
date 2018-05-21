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
    ns <- loadNamespace("vpc");
    if (exists("vpc_vpc",ns)){
        vpcn <- "vpc_vpc"
    } else {
        vpcn <- "vpc"
    }
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    if (!is.null(by)) {
        call$strat <- c("grp")
        call$facet <- "wrap";
    }
    call <- call[names(call) %in% methods::formalArgs(getFromNamespace(vpcn,"vpc"))]
    p = do.call(getFromNamespace(vpcn,"vpc"), c(list(sim=df, obs=dat), call), envir = parent.frame(1))
    print(p)
    invisible(df)
}

