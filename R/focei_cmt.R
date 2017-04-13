## gnlmm.R: population PK/PD modeling library
##
## Copyright (C) 2014 - 2016  Wenping Wang
##
## This file is part of nlmixr.
##
## nlmixr is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## nlmixr is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with nlmixr.  If not, see <http:##www.gnu.org/licenses/>.

require(lbfgs)
require(lbfgsb3)
require(madness)
require(Rcpp)
#require(nlmixr)
require(parallel)
require(minqa)
require(Deriv)
lin_cmt_stan <- function(obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)
   .Call('lin_cmt_stan', obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)

llik_binomial <- function(y,n,params) {
	r = .Call('llik_binomial', y,n,params)
	r$J = diag(r$J)
	r
}
llik_normal <- function(y,params) .Call('llik_normal', y,params)
llik_beta <- function(y,params) .Call('llik_beta', y,params)
llik_neg_binomial <- function(y,params) .Call('llik_neg_binomial', y,params)
llik_poisson <- function(y,params) .Call('llik_poisson', y,params)
llik_betabinomial <- function(y,n,params) .Call('llik_betabinomial', y,n,params)
llik_student_t <- function(y,params) .Call('llik_student_t', y,params)



#-- for new gnlmm
getModelVars = function(blik, bpar, m1) {
	argsList = list(
		dpois = list(ix=2, dvdx=".arg1"),
		dbinom= list(ix=2:3, dvdx=".arg2"),
		dnorm = list(ix=2:3, dvdx=c(".arg1", ".arg2")),
		dbeta = list(ix=2:3, dvdx=c(".arg1", ".arg2")),
		dneg_binomial = list(ix=2:3, dvdx=c(".arg1", ".arg2")),
		dbetabinomial = list(ix=2:4, dvdx=c(".arg2", ".arg3")),
		dt = list(ix=2:4, dvdx=c(".arg1", ".arg2", ".arg3"))
	)

	blik.txt = deparse(blik)
	len = length(blik.txt)

	s = gsub("\\s+", "", blik.txt[len-1], perl=T)	#FIXME
	lp=regexpr("\\(", s)
	dist = substr(s, 1, lp-1)
	s = strsplit(substr(s, lp+1, 200), ",")[[1]]

	args = argsList[[dist]]
	args.ix = args$ix		#position of dist pars
	args.dvdx = args$dvdx	#args need dvdx
	narg = length(args.ix)
	blik.new.text = paste(c(
		blik.txt[2:(len-2)],
		paste0(".arg", 1:narg, "=", s[args.ix])
		), collapse="\n"
	)
	blik.new = parse(text=blik.new.text)
	
	dist.df = NULL
	if (dist=="dbinom") dist.df=s[2]	#binomial size
	if (dist=="dt") dist.df=s[2]		#t df

	#----------------------------
	f = deparse(blik)
	len = length(f)
	f = parse(text=f[c(-1,-len)])
	out = utils::getParseData(f)
	s = out$text[out$token %in% c("SYMBOL", "LEFT_ASSIGN", "EQ_ASSIGN")]
	ix = s %in% c("=", "<-")
	lhs = (1:length(ix))[ix] - 1
	ix[lhs] = TRUE
	rhs = setdiff(s, s[ix])

	ixLlik = len-1
	ss = sub("^\\s*\\w+\\(", "", deparse(blik)[len-1], perl=T)	#FIXME
	prob = strsplit(ss, ",")[[1]][3]	#FIXME

	states = m1$get.modelVars()$state
	state.llik = intersect(states, s[!ix])	#state used in llik

	f = deparse(bpar)
	len = length(f)
	f = parse(text=f[c(-1,-len)])
	out = utils::getParseData(f)
	s = out$text[out$token %in% c("SYMBOL", "LEFT_ASSIGN", "EQ_ASSIGN")]
	ix = s %in% c("=", "<-")
	lhs = (1:length(ix))[ix] - 1
	pars.llik=intersect(s[lhs], rhs)	#pars used in llik
	vars.par = s[lhs]					#vars def'ed in pars

	list(state.llik=state.llik, pars.llik=pars.llik, 
		 vars.par=vars.par, #prob=prob, ixLlik=ixLlik,
		 dist=dist, dist.df=dist.df, 
		 blik.new=blik.new, blik.new.text=blik.new.text,
		 args.dvdx=args.dvdx)
}


gnlmm2 <- function(llik, data, inits, syspar=NULL, 
	system=NULL, diag.xform=c("sqrt", "log", "identity"),
	..., control=list()) {
    #data
	if (is.null(data$ID)) stop('"ID" not found in data')
	if (is.null(data$EVID)) data$EVID = 0
	data.obs = subset(data, data$EVID == 0)
	data.sav = data
	names(data) <- tolower(names(data))		#needed in ev
    
    #model
    if (is.null(system)) {}
    else if (class(system) == "RxODE") {}
    else if (class(system) == "character") {
        obj <- basename(tempfile())
        system <- RxODE(model = system, modName = obj)
    }
    else {
        stop("invalid system input")
    }    

    #options
    con <- list(trace = 0, 
    	maxit = 100L, 
        atol.ode=1e-08, 
        rtol.ode=1e-08,
        reltol.inner = 1.0e-4,
        reltol.outer = 1.0e-3,
        optim.inner="lbfgs",
        optim.outer="newuoa",
        start.zero.inner=FALSE,
        mc.cores=1,
        nAQD=1,
        transit_abs=FALSE,
        cov=FALSE,
        eps=c(1e-8, 1e-3),		#finite difference step
        NOTRUN = F,
        DEBUG.INNER = F,
        rhobeg = .2,
        rhoend = 1e-3,
        iprint = 2,
        npt=NULL,
        do.optimHess=TRUE
        )
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

	square = function(x) x*x
	diag.xform <- match.arg(diag.xform)
	diag.xform.inv = c("sqrt"="square", "log"="exp", "identity"="identity")[diag.xform]

	#process inits
	lh = parseOM(inits$OMGA)
	nlh = sapply(lh, length)
	osplt = rep(1:length(lh), nlh)

	lini = list(inits$THTA, unlist(lh))
	nlini = sapply(lini, length)
	nsplt = rep(1:length(lini), nlini)

	om0 = genOM(lh)
	th0.om = lapply(1:length(lh), function(k)
	{
		m = genOM(lh[k])
		nr = nrow(m)
		mi = tryCatch(
			backsolve(chol(m), diag(nr)),
			error = function(e) {
				stop("OMEGA not positive-definite") 
			}
		)
		diag(mi) = eval(call(diag.xform, diag(mi)))
		mi[col(mi)>=row(mi)]
	})
	inits.vec = c(inits$THTA, unlist(th0.om), inits$SGMA)
	names(inits.vec) = NULL

	nTHTA = nlini[1]
	nETA  = nrow(om0)
	ID.all = unique(data[,"id"])
	ID.ord = order(ID.all)
	names(ID.ord) = ID.all
	nSUB  = length(ID.all)

	#gaussian quadrature nodes & wts
	nAQD = con$nAQD
	nw = gauss.quad(nAQD)
	mij = as.matrix(
		do.call("expand.grid", lapply(1:nETA, function(x) 1:nAQD))
	)
	nij = nrow(mij)

	
	#obj fn by AQD
	if (!is.null(syspar)) {
		bpar = body(syspar)
	}
	blik = body(llik)
	modVars = getModelVars(blik, bpar, system)
	

	#=== start of dvdx code
	proc.deriv = function() {
		getDeriv = function(pars) {
			npar = length(pars)
			s = deparse(bpar)
			for (i in 1:nETA) {
				s = gsub(sprintf("\\bETA\\[%d\\]", i), sprintf("ETA%d", i), s, perl=T)
			}
			s = gsub("initCondition", "#initCondition", s)
			len = length(s)
			s = s[-len]

			ix = matrix(unlist(expand.grid(1:npar, 1:nETA)), ncol=2)
			m = sapply(1:nrow(ix), function(k) {
				i = ix[k,1]
				j = ix[k,2]
				s1 = c("Deriv(~", s, pars[i], sprintf("}, \"ETA%d\")", j))
				a = paste(s1, collapse="\n")
				e = eval(parse(text=a))
				s = if (class(e)=="call") {
					s = deparse(e)
					gsub(sprintf("\\bETA%d\\b", j), sprintf("ETA[%d]", j), s, perl=T)
				} else "0"
				#parse(text=s)
				s
			})
			m
		}

		env = environment()
		dati = data.sav[data.sav$ID==ID.all[1], ]
		list2env(dati, env)
		THETA = tapply(inits.vec, nsplt, identity)[[1]]
		ETA <- madness(array(0, c(nETA,1))) 
		eval(bpar)

		px = as.list(env)
		madIx <- sapply(px, function(x) {
				if (class(x)=="madness") TRUE else FALSE
		})
		madVars = names(px)[madIx]

		pars = system$get.modelVars()$params
		px = as.list(env)[pars]
		madIx <- sapply(px, function(x) {
				if (class(x)=="madness") TRUE else FALSE
		})

		pars = pars[madIx]
		matode = getDeriv(pars)
		pars = setdiff(madVars, c(pars, "ETA"))
		matllk = getDeriv(pars)
		
		list(matode=matode, madIx=madIx, matllk=matllk, madllk=pars)	#deriv expr for ode pars
		                                                                #idx of ode pars that need deriv
		                                                                #llik pars that need deriv
		                                                                #deriv expr llik pars 
	}
	s = proc.deriv()

	# FIXME: chk par order of below 2 para against those in d(llik)/d(ETA) by chain rule
	#d(pars)/d(ETA) for ode
	madIx = s$madIx
	npar = sum(madIx)
	m = s$matode
	idx.dpde.ode = m!="0"
	expr.dpde.ode = lapply(m[idx.dpde.ode], function(k) parse(text=k))
	dpde.ode = matrix(0, npar, nETA)

	#d(pars)/d(ETA) for llik
	madVars.llk = s$madllk
	npar = length(madVars.llk)
	m = s$matllk
	idx.dpde.llk = m!="0"
	expr.dpde.llk = lapply(m[idx.dpde.llk], function(k) parse(text=k))
	dpde.llk =  matrix(0, npar, nETA)
	dimnames(dpde.llk) = list(madVars.llk, NULL)

	#d(args)/d(pars) for llik
	pars = c(modVars$state.llik, madVars.llk)
	npar = length(pars)
	lexpr = unlist(lapply(modVars$args.dvdx, function(arg) {
		s = sprintf("{\n%s\n%s", modVars$blik.new.text, arg)	#FIXME
		 
		expr.dldp = lapply(1:npar, function(k) {
			s1 = c("Deriv(~", s, sprintf("}, \"%s\")", pars[k]))
			a = paste(s1, collapse="\n")
			e = eval(parse(text=a))
			e
		})
		names(expr.dldp) = pars
		expr.dldp
	}))
	
	llik.narg = length(modVars$args.dvdx)	#llik.narg = # of args in density that need dvdx
											#may have efficiency gain to rm args that do not need dvdx
	llik.npar = npar
	m = matrix(lexpr, llik.npar, llik.narg)
	dadp.expr = c(t(m))						#chg the order of args & pars in llik; see dldp in ..fg()
	#=== ends of dvdx code 
	

	starts = matrix(0., nSUB, nETA)
	omga_save = NULL; update_starts=TRUE


	#algo starts
	obj.vec = function(th, noomga=FALSE) {
		th = th*inits.vec
		#if (con$DEBUG.INNER) print(th)

		if(noomga) th = c(th, omga_save)

		lth = tapply(th, nsplt, identity)
		THETA = lth[[1]]
		lD = tapply(lth[[2]], osplt, identity)
		Dinv.5 = genOMinv.5(lD)
		diag(Dinv.5) = eval(call(diag.xform.inv, diag(Dinv.5)))
		detDinv.5 = prod(diag(Dinv.5))

		ep = environment()

	llik2.subj = function(ix) {
		#ix=4; ETA=rep(0, nETA)
		dati = data.sav[data.sav$ID==ix, ]
		evi = dati[, c("TIME", "EVID", "AMT")]
		names(evi) = tolower(names(evi))
		ev = eventTable()
		ev$import.EventTable(evi)
		dati = dati[dati$EVID==0, ]


		# ETA => pars & stateVar => .args => llik
		# d(llik)/d(ETA) = d(llik)/d(args) * d(args)/d(ETA)
		# d(args)/d(ETA) = d(args)/d(pars) * d(pars)/d(ETA)
		# stateVar is parallel to pars when forming args, however, stateVar = f(pars(ETA))
		# hence, we need d(State)/d(ETA). d(State)/d(ETA) = d(State)/d(pars) * d(pars)/d(ETA)
		..g.fn = function(ETA) {
			env = environment()
			list2env(dati, env)
			eval(bpar)

			pars = system$get.modelVars()$params
			po = unlist(as.list(env)[pars])
			x = system$run(po, ev, initCondition)
			if(any(is.na(x))) {
				print(ID[1])
				print(po);
				print(head(x, 10))
				stop("NA in ODE solution")
			}

			#d(State)/d(ETA)
			whState = modVars$state.llik
			senState = paste0(whState, "_", pars[madIx]) 
			fxJ = list(fx=x[, whState], J=x[, senState])		#FIXME, t()

			dvdx = sapply(expr.dpde.ode, eval, envir=env)		#FIXME
			dpde.ode[idx.dpde.ode] = dvdx
			#d(State)/d(ETA) = d(State)/d(pars) * d(pars)/d(ETA)
			dvdxState = fxJ$J %*% dpde.ode						#FIXME

			#make state var & other symbols available for calc
			valState = fxJ$fx
			assign(whState, valState, envir=env)
			eval(modVars$blik.new)								#FIXME: here or at llik_binomial?

			#d(pars)/d(ETA)
			if (length(expr.dpde.llk)>0) {
			dvdx = sapply(expr.dpde.llk, eval, envir=env)		#FIXME
			dpde.llk[idx.dpde.llk] = dvdx
			} else {
			dpde.llk = NULL
			}

			#d(args)/d(pars) 
			ni = dim(x)[1]
			dadp = sapply(1:length(dadp.expr), function(k) {	#why lapply(m, eval) doesn't work?
				s = eval(dadp.expr[[k]])
				if (length(s)==1) rep(s, ni) else s
			})
			dim(dadp) = c(ni, llik.narg, llik.npar)				#FIXME

			#d(args)/d(ETA) = d(args)/d(pars) * d(pars)/d(ETA)
			dade = sapply(1:ni, function(k) {					#FIXME: vectorize?
				dpde = rbind(dvdxState[k,], dpde.llk)
				dadp[k,,] %*% dpde								#FIXME: need t()?
			})
			#dade = t(s)										#FIXME: t() can be rm'ed?
			dim(dade) = c(llik.narg, nETA, ni)
			

			#d(llik)/d(ETA) = d(llik)/d(args) * d(args)/d(ETA)
			if (modVars$dist=="dt") {
				if (length(.arg1)==1 && ni>1) .arg1 = rep(.arg1, ni)
				if (length(.arg2)==1 && ni>1) .arg2 = rep(.arg2, ni)
				if (length(.arg3)==1 && ni>1) .arg3 = rep(.arg3, ni)
				s = sapply(1:ni, function(k) {					#FIXME: vectorize?
					unlist(llik_student_t(DV[k], c(.arg1[k], .arg2[k], .arg3[k])))
				})
				s = list(fx=s[1,], J=t(s[-1,]))
			} else if (modVars$dist=="dbetabinomial") {
				if (length(.arg1)==1 && ni>1) .arg1 = rep(.arg1, ni)
				if (length(.arg2)==1 && ni>1) .arg2 = rep(.arg2, ni)
				if (length(.arg3)==1 && ni>1) .arg3 = rep(.arg3, ni)
				s = sapply(1:ni, function(k) {					#FIXME: vectorize?
					unlist(llik_betabinomial(DV[k], .arg1[k], c(.arg2[k], .arg3[k])))
				})
				s = list(fx=s[1,], J=t(s[-1,]))
			} else if (modVars$dist=="dneg_binomial") {
				if (length(.arg1)==1 && ni>1) .arg1 = rep(.arg1, ni)
				if (length(.arg2)==1 && ni>1) .arg2 = rep(.arg2, ni)
				s = sapply(1:ni, function(k) {					#FIXME: vectorize?
					unlist(llik_neg_binomial(DV[k], c(.arg1[k], .arg2[k])))
				})
				s = list(fx=s[1,], J=t(s[-1,]))
			} else if (modVars$dist=="dbeta") {
				if (length(.arg1)==1 && ni>1) .arg1 = rep(.arg1, ni)
				if (length(.arg2)==1 && ni>1) .arg2 = rep(.arg2, ni)
				s = sapply(1:ni, function(k) {					#FIXME: vectorize?
					unlist(llik_beta(DV[k], c(.arg1[k], .arg2[k])))
				})
				s = list(fx=s[1,], J=t(s[-1,]))
			} else if (modVars$dist=="dnorm") {
				if (length(.arg1)==1 && ni>1) .arg1 = rep(.arg1, ni)
				if (length(.arg2)==1 && ni>1) .arg2 = rep(.arg2, ni)
				s = sapply(1:ni, function(k) {					#FIXME: vectorize?
					if (.arg2[k]<0.0001) {
						#print("HAY!");print(.arg1[k]);print(.arg2[k])
					}
					unlist(llik_normal(DV[k], c(.arg1[k], .arg2[k])))
				})
				s = list(fx=s[1,], J=t(s[-1,]))
			} else if (modVars$dist=="dbinom") {
				if (length(.arg1)==1 && ni>1) .arg1 = rep(.arg1, ni)
				s = llik_binomial(DV, .arg1, c(.arg2))
			} else if (modVars$dist=="dpois") {
				s = llik_poisson(DV, c(.arg1))
			} else {
				stop("dist not supported")
			}

			dim(s$J) = c(ni, llik.narg)
			dlde = sapply(1:ni, function(k) {					#FIXME: vectorize?
				s$J[k,] %*% dade[,,k]							#FIXME: need t()?
			})		                                            
			s = rbind(s$fx, dlde)								#FIXME t()?
			s = apply(s, 1, sum)								#FIXME sum index; chg'ed w/ vec stan call

			#llik.dat = madness(val=matrix(s[1], 1, 1), dvdx=matrix(s[-1], 1, nETA))	
			#llik.eta = -crossprod(Dinv.5 %*% ETA)/2 -nETA/2*log(2*pi)+log(detDinv.5)
			#llik.eta = madness(val=val(llik.eta), dvdx=dvdx(llik.eta))
			llik.eta.val = -crossprod(Dinv.5 %*% ETA)/2 -nETA/2*log(2*pi)+log(detDinv.5)
			llik.eta.dvd = -t(ETA) %*% crossprod(Dinv.5)
			
			r = s[1] + c(llik.eta.val)
			attr(r, "dvdx") = s[-1] + c(llik.eta.dvd)
			r
		}
		fg <- function(par) {
			if(identical(par, pvd[[1]])) return(pvd)
			ym = ..g.fn(par)
			pvd <<- list(par, ym)
		}
		f = function(ETA)  -as.vector(fg(ETA)[[2]])
		g = function(ETA) -as.vector(attr(fg(ETA)[[2]],"dvdx"))
		
		pvd = NULL

		.wh = ID.ord[as.character(ix)]
		ETA.val = starts[.wh, ]
		..fit.inner = nlminb2(ETA.val, f, g, control=list(trace=FALSE, rel.tol=1e-4))
		#..fit.inner = lbfgs(f, g, ETA.val, invisible=T, ftol=1e-4) # epsilon=1e-3)
		if (con$do.optimHess) {
			..fit.inner$hessian = optimHess(..fit.inner$par, f, g)
		}
		if (con$DEBUG.INNER) {
			#print(..fit.inner$message)
		}

		#=========================================================
		if (con$do.optimHess) {
			Ginv.5 = tryCatch({
					.m <- chol(..fit.inner$hessian)
					backsolve(.m, diag(nETA))
				},
				error = function(e) {
					cat("Warning: Hessian not positive definite\n")
					print(..fit.inner$hessian)
					.m <- ..fit.inner$hessian
					#.m <- chol(.m+diag(nETA)*100)
					#.m[col(.m)!=row(.m)] = .001*.m[col(.m)!=row(.m)]
					.md = matrix(0, nETA, nETA)
					diag(.md) = abs(diag(.m))*1.1
					.m <- chol(.md)
					backsolve(.m, diag(nETA))
				}
			)
		} else {
			Ginv.5 = chol(..fit.inner$Hessian.inv)
		}
		det.Ginv.5 = prod(diag(Ginv.5))

		##-- AQD
		..lik.ij = lapply(1:nij, function(ix) {
			ij = mij[ix,]
			w = nw$weights[ij]
			z = nw$nodes[ij]
			a = ..fit.inner$par + sqrt(2)*Ginv.5%*%z
			f1 = exp(as.vector(..g.fn(a)))						#FIXME
			f2 = prod(w*exp(z^2))
			f1*f2
		})

		..lik = 2^(nETA/2)*det.Ginv.5*do.call("sum", ..lik.ij)
		c(-2 * log(..lik), .wh, ..fit.inner$par)
	}
	s = mclapply(ID.all, llik2.subj, mc.cores = con$mc.cores)	#FIXME
    m = matrix(unlist(s), ncol = 2 + nETA, byrow = T)
		
		if (update_starts) starts[m[, 2], ] <<- m[, 3:(2+nETA)]
		m[, 1]
	}

	nobjcall = 0
	obj = function(th, noomga=FALSE) {
		nobjcall <<- nobjcall+1
		s = obj.vec(th, noomga)
		r = sum(s)
		if (con$DEBUG.INNER) {
			print(rbind(
				c(nobjcall, r, th),
				c(nobjcall, r, th*inits.vec)
			))
		}
		attr(r, "subj") = s
		r
	}

	np = length(inits.vec)
	start = rep(1, np)
	args = list(start, obj, control=list(trace=con$trace, reltol=con$reltol.outer))
	
  if (!con$NOTRUN) {
	fit = if (con$optim.outer=="nmsimplex") do.call("nmsimplex", args)
	else if (con$optim.outer=="Nelder-Mead") {
		args = c(args, method=con$optim.outer)
		do.call("optim", args)
	}
	else if (con$optim.outer=="nlminb") {
		args = list(start, obj, control=list(trace=con$trace, rel.tol=con$reltol.outer))
		do.call("nlminb", args)
	}
	else {
		if (!is.null(con$npt)) npt=con$npt
		else npt = 2*np+1
		newuoa(start, obj, control=list(rhobeg=con$rhobeg, rhoend=con$rhoend, npt=npt, iprint=con$iprint))
	}
  } else fit=NULL
	
	fit = c(fit, obj=obj, list(ETA=starts, con=con, diag.xform=diag.xform, nsplt=nsplt, osplt=osplt, calls=list(data=data.sav, system=system, syspar=syspar)))
	fit$par.unscaled = fit$par*inits.vec
	attr(fit, "class") <- "gnlmm.fit"
	fit
}


#------------------
mat.indices = function(nETA){
	idx = do.call("rbind", 
		lapply(1:nETA, function(k) cbind(k:nETA, k)))
	H = matrix(1:(nETA^2), nETA, nETA)
	Hlo.idx = row(H)>=col(H)
	lo.idx = H[row(H)>col(H)]
	hi.idx = t(H)[row(H)>col(H)]

	list(idx=idx,			# (r, c) of lo-half 
	     Hlo.idx=Hlo.idx, 	# index of lo-half
	     lo.idx=lo.idx, 	# index of strict lo-half
	     hi.idx=hi.idx)		# index of strict hi-half
}


#' Print a gnlmm fit
#'
#' Print a generalized non-linear mixed effect model fit
#' 
#' @param x a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
print.focei.fit = function(fit) {
	attr(fit, "data") = NULL
	attr(fit, "ofv.FOCEi") = NULL
	fit$par = NULL
	fit$Hessian.inv = NULL
	print.default(fit)
}

#' Print a gnlmm fit
#'
#' Print a generalized non-linear mixed effect model fit
#' 
#' @param x a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
plot.focei.fit = function(fit) {
	ofv.FOCEi = attr(fit, "ofv.FOCEi")
	dat = attr(fit, "data")
	
	x = ofv.FOCEi(fit$par)
	df = cbind(dat, IPRED=unlist(
		lapply(attr(x,"subj"), function(s) attr(s,"fitted"))
	))

	require(ggplot2)
	p1 = ggplot(df, 
	  aes(x=IPRED, y=DV)) +
	  geom_point() +
	  geom_abline(slope=1, intercept=0, col="red")
	p2 = ggplot(df, 
	  aes(x=IPRED, y=DV-IPRED)) +
	  geom_point() +
	  geom_abline(slope=0, intercept=0, col="red")
	p3 = ggplot(df, 
	  aes(x=TIME, y=DV)) +
	  geom_point() +
	  geom_line(aes(x=TIME, y=IPRED), data=df, col="red") +
	  facet_wrap(~ID)

	list(p1, p2, p3)
}


proc.err = function(f) {
    s = unlist(lapply(attr(terms(f), "variables"), as.list))
    s = sapply(s, deparse)
    ix.add = match("add", s, nomatch = 0)
    ix.pro = match("prop", s, nomatch = 0)
    err.type = c("add", "prop", "combo")[(ix.add > 0) + 2 * 
        (ix.pro > 0)]
    sig.add = if (ix.add > 0) 
        as.numeric(s[ix.add + 1])
    else NULL
    sig.pro = if (ix.pro > 0) 
        as.numeric(s[ix.pro + 1])
    else NULL
    inits.err <- c(sig.add, sig.pro)
    if (any(is.na(inits.err) | inits.err <= 0)) 
        stop("error model misspecification")
    inits.err
}


#--------------------------------
require(RxODE)

ode = "
d/dt(depot) = -KA*depot;
d/dt(C2   ) = KA/V*depot - CL/V*C2;
d/dt(depot_KA) = (-KA) * (depot_KA) + -depot; 
d/dt(C2_KA) = (KA/V) * (depot_KA) + (-(CL/V)) * (C2_KA) + 1/V*depot; 
d/dt(depot_V) = (-KA) * (depot_V); 
d/dt(C2_V) = (KA/V) * (depot_V) + (-(CL/V)) * (C2_V) + -(KA/V^2*depot-CL/V^2*C2); 
d/dt(depot_CL) = (-KA) * (depot_CL); 
d/dt(C2_CL) = (KA/V) * (depot_CL) + (-(CL/V)) * (C2_CL) + -(1/V*C2); 
"

#FIXME
#m1 = RxODE(ode, modName="m1")

solveODE = function(parsODE, ev) {
	x = m1$run(parsODE, ev)
	list(fx=x[,"C2"], J=t(cbind(x[, c("C2_CL", "C2_V", "C2_KA")], 0)))
}



#' Fit a generalized nonlinear mixed-effect model
#'
#' Fit a generalized nonlinear mixed-effect model by adapative Gaussian quadrature (AQD)
#' 
#' @param llik log-likelihood function
#' @param data data to be fitted
#' @param inits initial values 
#' @param system an optional (compiled) RxODE object
#' @param syspar function: calculation of PK parameters
#' @param diag.xform transformation to diagonal elements of OMEGA during fitting
#' @param ... additional options
#' @param control additional optimization options
#' @return NULL
#' @details
#'    Fit a generalized nonlinear mixed-effect model by adapative Gaussian quadrature (AGQ)
#' 
#' @author Wenping Wang
#' @examples
#' \dontrun{
#' llik <- function()
#' {
#' 	if (group==1) lp = THETA[1]+THETA[2]*logtstd+ETA[1]
#' 	else          lp = THETA[3]+THETA[4]*logtstd+ETA[1]
#' 	lam = exp(lp)
#' 	dpois(y, lam, log=TRUE)
#' }
#' inits = list(THTA=c(1,1,1,1), OMGA=list(ETA[1]~1))
#' 
#' fit = gnlmm(llik, pump, inits, 
#' 	control=list(
#' 	    reltol.outer=1e-4,
#' 		optim.outer="nmsimplex",
#' 		nAQD=5
#' 	)
#' )
#' 
#' 
#' 
#' llik <- function()
#' {
#' 	lp = THETA[1]*x1+THETA[2]*x2+(x1+x2*THETA[3])*ETA[1]
#' 	p = pnorm(lp)
#' 	dbinom(x, m, p, log=TRUE)
#' }
#' inits = list(THTA=c(1,1,1), OMGA=list(ETA[1]~1))
#' 
#' gnlmm(llik, rats, inits, control=list(nAQD=7))
#' 
#' 
#' ode <- "
#' d/dt(depot) =-KA*depot;
#' d/dt(centr) = KA*depot - KE*centr;
#' "
#' sys1 = RxODE(ode)
#' 
#' pars <- function()
#' {
#' 	CL = exp(THETA[1] + ETA[1])#; if (CL>100) CL=100
#' 	KA = exp(THETA[2] + ETA[2])#; if (KA>20) KA=20
#' 	KE = exp(THETA[3])
#' 	V  = CL/KE
#' 	sig2 = exp(THETA[4])
#' }
#' llik <- function() {
#' 	pred = centr/V
#' 	dnorm(DV, pred, sd=sqrt(sig2), log=TRUE)
#' }
#' inits = list(THTA=c(-3.22, 0.47, -2.45, 0))
#' inits$OMGA=list(ETA[1]~.027, ETA[2]~.37)
#' #inits$OMGA=list(ETA[1]+ETA[2]~c(.027, .01, .37))
#' theo <- read.table("theo_md.txt", head=TRUE)
#' 
#' fit = gnlmm(llik, theo, inits, pars, sys1, 
#' 	control=list(trace=TRUE, nAQD=5))
#' 
#' cv = calcCov(fit)
#' cbind(fit$par[fit$nsplt==1], sqrt(diag(cv)))
#' 
#' }
focei.fit = function(
	data, 
	inits, 
	PKpars, 
	diag.xform=c("sqrt", "log", "identity"), 
	optim=c("newuoa", "lbfgsb3", "nlminb", "nelder-mead"), 
	model=list(), 
	control=list()
){
    #data = dat; PKpars=mypars; diag.xform="sqrt"; model=list(); control=list()

    #model options
	if(class(model)=="RxODE") {
		ODEmodel = TRUE
		mod = list(infusion=FALSE)
	}
	else if (is.list(model)) {
		ODEmodel = FALSE
		mod = list(
		ncmt=1,
		oral=T,
		tlag=F,
		infusion=F,
		parameterization=1
		)
		nmsC <- names(mod)
		mod[(namc <- names(model))] <- model
		if (length(noNms <- namc[!namc %in% nmsC])) 
			warning("unknown names in model: ", paste(noNms, collapse = ", "))
	}
	else {
		stop("wrong model input")
	}

    #FIXME
    #data
	if (is.null(data$ID)) stop('"ID" not found in data')
	if (is.null(data$EVID)) data$EVID = 0
	data.sav = data
    if (mod$infusion) {
        data <- fmt_infusion_data(data)
    }
    else {
        data$DUR <- -1
    }
    ds <- data[data$EVID > 0, c("ID", "TIME", "AMT", "DUR")]
    data$DUR <- NULL
    data <- data[data$EVID == 0, ]
	names(data.sav) <- tolower(names(data.sav))		#needed in ev
        
    #options
    con <- list(
        DEBUG.ODE = F,
        DEBUG = F, RESET.INITS.MAT = T, TRACE.INNER=F, TOL.INNER=1e-4,
        trace = 0, 
        atol.ode=1e-08, 
        rtol.ode=1e-08,
        reltol.outer = 1.0e-2,
        mc.cores=1,
        transit_abs=FALSE,
        NONMEM=TRUE,
        NOTRUN=FALSE,
        PRINT.PARS=FALSE,
        cov.method="hessian",
        factr=1e10,
        eps = c(sqrt(.Machine$double.eps),1e-6),
        mc.cores=1,
        rhobeg=.2, 
        rhoend=1e-2,
        npt=NULL
        )
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

	square = function(x) x*x
	diag.xform <- match.arg(diag.xform)
	diag.xform.inv = c("sqrt"="square", "log"="exp", "identity"="identity")[diag.xform]
	optim.method <- match.arg(optim)

	#process inits
	inits$SGMA = proc.err(inits$ERROR[[1]])
	
	lh = parseOM(inits$OMGA)
	nlh = sapply(lh, length)
	osplt = rep(1:length(lh), nlh)

	lini = list(inits$THTA, unlist(lh), inits$SGMA)	#CHECK!
	nlini = sapply(lini, length)
	nsplt = rep(1:length(lini), nlini)

	om0 = genOM(lh)
	th0.om = lapply(1:length(lh), function(k)
	{
		m = genOM(lh[k])
		nr = nrow(m)
		mi = tryCatch(
			backsolve(chol(m), diag(nr)),
			error = function(e) {
				stop("OMEGA not positive-definite") 
			}
		)
		diag(mi) = eval(call(diag.xform, diag(mi)))
		mi[col(mi)>=row(mi)]
	})
	inits.vec = c(inits$THTA, unlist(th0.om), sqrt(inits$SGMA))
	names(inits.vec) = NULL
    #print(inits.vec)

	nTHTA = nlini[1]
	nETA  = nrow(om0)
	ID.all = unique(data[,"ID"])	#FIXME
	ID.ord = order(ID.all)
	names(ID.ord) = ID.all
	nSUB  = length(ID.all)
	Hidx = mat.indices(nETA)
	

	#proc pars deriv
	pm <- list(
		c("CL", "V", "KA", "TLAG"),
		c("CL", "V", "CLD", "VT", "KA", "TLAG"),
		c("CL", "V", "CLD", "VT", "CLD2", "VT2", "KA", "TLAG"),
		c("KE", "V", "KA", "TLAG"),
		c("KE", "V", "K12", "K21", "KA", "TLAG"),
		c("KE", "V", "K12", "K21", "K13", "K31", "KA", "TLAG")
	)
	dim(pm)<-c(3,2)
	
	if (is.list(mod)) {
	ncmt=mod$ncmt; 
	oral=mod$oral; 
	parameterization=mod$parameterization
	pars = pm[[ncmt, parameterization]]
	npar = 2*(ncmt+oral)
	}

	# get symbolic derivatives
	bpar = body(PKpars)
	s = deparse(bpar)
	for (i in 1:nETA) {
		s = gsub(sprintf("\\bETA\\[%d\\]", i), sprintf("ETA%d", i), s, perl=T)
	}
	s = gsub("initCondition", "#initCondition", s)
	len = length(s)
	s = s[-len]

	ix = matrix(unlist(expand.grid(1:npar, 1:nETA)), ncol=2)
	m = sapply(1:nrow(ix), function(k) {
		i = ix[k,1]
		j = ix[k,2]
		s1 = c("Deriv(~", s, pars[i], sprintf("}, \"ETA%d\")", j))
		a = paste(s1, collapse="\n")
		e = eval(parse(text=a))
		s = if (class(e)=="call") {
			s = deparse(e)
			gsub(sprintf("\\bETA%d\\b", j), sprintf("ETA[%d]", j), s, perl=T)
		} else "0"
		#parse(text=s)
		s
	})

	idx.save = m!="0"
	expr.save = lapply(m[idx.save], function(k) parse(text=k))
	mat.save = matrix(0, npar, nETA)

	#algo starts	
	ofv.FOCEi.ind = function(pars, NONMEM=con$NONMEM) {

		if(con$PRINT.PARS) print(pars)
		pars = pars*inits.vec
		lth = tapply(pars, nsplt, identity)
		THETA = lth[[1]]
		sig2 = lth[[3]]^2	#FIXME
		lD = tapply(lth[[2]], osplt, identity)
		OMGAinv.5 = genOMinv.5(lD)
		diag(OMGAinv.5) = eval(call(diag.xform.inv, diag(OMGAinv.5)))
		log.det.OMGAinv.5 = sum(log(diag(OMGAinv.5)))
		OMGAinv   = crossprod(OMGAinv.5)

		
		env = environment()
		PKpars__ = PKpars
		environment(PKpars__) = env 


		#------------------------------
		llik.subj = mclapply(ID.all, function(subj) 
		{

			env = environment()
			
			if (ODEmodel) {
				dati = data.sav[data.sav$id==subj, ]
			}
			else {
				dati = data[data$ID==subj, ]
				dsi = ds[ds$ID==subj, ]
			}
			list2env(dati, env)
		   #print(ls(env=env)); str(ID)


			fg = function(ETA) {
				if(identical(ETA, pvd[[1]])) return(pvd)

				#if(nfcall==0) str(ID)
				nfcall <<- nfcall+1

				len = length(ETA)
				#xm  = madness(array(ETA, dim=c(len,1)))
				#pm = PKpars__(xm)
				list2env(list(ETA=ETA), env)
				pm = PKpars__(ETA)
				dvdx = sapply(expr.save, eval, envir=env)
				mat.save[idx.save] = dvdx


				#FIXME
				if (ODEmodel) {
					ev = eventTable()
					dati = data.sav[data.sav$id==subj, ]
					ev$import.EventTable(dati)
					assign("DV", dv[ev$get.obs.rec()], envir=env)
					parsODE = as.vector(val(pm))
					names(parsODE) = c("CL", "V", "KA", "TLAG")	#FIXME
					fxJ = solveODE(parsODE, ev)	#FIXME 
					if (con$DEBUG.ODE) print("i'm here :)")
				}
				else {
				   #fxJ = lin_cmt(TIME, dsi$TIME, dsi$AMT, dsi$DUR, val(pm), mod$oral, mod$infusion, mod$ncmt, mod$parameterization)
					fxJ = lin_cmt_stan(TIME, dsi$TIME, dsi$AMT, dsi$DUR, pm, mod$oral, mod$infusion, mod$ncmt, mod$parameterization)
				}				

				f = fxJ$fx
			   #fp = t(fxJ$J) %*% dvdx(pm)
				fp = fxJ$J[,1:npar] %*% mat.save		#FIXME: t() or not depending on stan version
				eps = DV - f
				c = 2*fp/f
				B = 2/(f^2*sig2)	

			   #lp = .5*apply(2*eps*fp/(f^2*sig2) + eps^2*2*fp/(f^3*sig2) - c, 2, sum) - OMGAinv %*% ETA
				lp = .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA

				#no need of log det(OMGAinv) when posthoc'ing; need it when Laplacian
				llik = -.5*sum(eps^2/(f^2*sig2) + log(f^2*sig2)) - .5*t(ETA) %*% OMGAinv %*% ETA

				pvd <<- list(ETA=ETA, llik=llik, lp=lp, f=f, fp=fp, eps=eps, B=B, c=c)
			}
			f = function(x) -fg(x)[[2]]
			g = function(x) -fg(x)[[3]]

			# pvd is a cache
			pvd = NULL; nfcall=0
			.wh = ID.ord[as.character(subj)]
			if(con$DEBUG>10 && .wh==1) print(inits.mat[.wh,])		#FIXME
			fit <- lbfgs(f, g, inits.mat[.wh,], invisible=1-con$TRACE.INNER, epsilon=con$TOL.INNER)
			if (is.na(fit$value)) {
				#cat("Warning: lbfgs step failed\n")
				fit <- optim(inits.mat[.wh,], f, control=list(trace=F, reltol=.001))
				#print(fit$convergence)
			}
			#if(con$RESET.INITS.MAT) inits.mat[.wh,] <<- fit$par

			#---------------------------------------
			# 2*pi in log(det(OMGA)) and log(-det(H)/(2*pi)) seems canceled
			# when doing laplacian, lbfgs result (value) is partial llik;
			# intermediates (fp, eps, c, B) can be used for a, H below
			# could use H from lbfgs update formula

			#cat("subj: ", subj, "\n");print(fit); print(pvd)
			pvd$a = with(pvd, if(NONMEM) fp else fp + eps*c)	# -a
			pvd$k = Hidx$idx[,1]; pvd$l = Hidx$idx[,2]

			H = matrix(0, nETA, nETA)
			H[Hidx$Hlo.idx] = with(pvd, 	# 1st fill the lower half
			  apply(a[,k]*B*a[,l] + c(-1, 1)[NONMEM+1] *c[,k]*c[,l], 2, sum))
			H[Hidx$hi.idx] = H[Hidx$lo.idx]

			H = -.5*H - OMGAinv
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
			#H.neg.5 = chol(-H)
			log.det.H.neg.5 = sum(log(diag(H.neg.5)))

			llik = pvd$llik + log.det.OMGAinv.5		#note no -1/2 with OMGAinv.5
			llik.lapl =llik - log.det.H.neg.5		#note no 1/2
			attr(llik.lapl, "fitted") = pvd$f
			attr(llik.lapl, "posthoc") = fit$par
			attr(llik.lapl, "wh") = .wh 
			
			llik.lapl

		}, mc.cores=con$mc.cores)
		
		m = t(sapply(llik.subj, function(x) {
			c(attr(x, "wh"), attr(x, "posthoc"))
		}))
		if(con$RESET.INITS.MAT) inits.mat[m[,1],] <<- m[,-1]
		
		
		llik.subj
	}

	ofv.FOCEi.vec = function(pars) {
		llik.subj = ofv.FOCEi.ind(pars)
		unlist(llik.subj)
	}
	
	ofv.FOCEi = function(pars) {
		llik.subj = ofv.FOCEi.ind(pars)
		llik = -2*do.call("sum", llik.subj)
		attr(llik, "subj") = llik.subj
		llik
	}


	eps = con$eps
	f = ofv.FOCEi
	g = function(x, ...) {
	  fx = f(x)
	  sapply(1:length(x), function(k) {
		xc = x
		dx = abs(xc[k]) * eps[1] + eps[2]
		xc[k] = xc[k] + dx
		(f(xc) - fx)/dx
	  })
	}


	inits.mat = matrix(0, nSUB, nETA)
	np = length(inits.vec) 
	start = rep(1, np)
	if (con$NOTRUN) {
		fit = list()
		fit$par = rep(1, length(inits.vec))
		fit$value = as.numeric(ofv.FOCEi(inits.vec))
		fit$convergence = -99
		fit$message = "notrun"
	}
	else {
		if (optim.method=="lbfgsb3") {
			fit = lbfgsb3(start, f, g, control=list(trace=T, factr = con$factr, pgtol = con$reltol.outer))
			fit$par = fit$prm
			fit$prm = NULL 
		}
		else if (optim.method=="nlminb") {
			fit = nlminb2(start, ofv.FOCEi, 
				control=list(trace=con$trace, rel.tol=con$reltol.outer))
			m = diag(inits.vec)
			if (con$cov.method=="hessian") { 
				fit$cov = 2*m %*% fit$Hessian.inv %*% m
				fit$Hessian.inv = NULL
			}
			else {
				R1 = optimHess(fit$par, ofv.FOCEi)
				fit$cov = m %*% (2*solve(R1)) %*% m
			}
			fit$se = sqrt(diag(fit$cov))
		}
		else if (optim.method=="newuoa") {
			if (!is.null(con$npt)) npt=con$npt
			else npt = 2*np+1
			fit = newuoa(start, ofv.FOCEi, 
				control=list(rhobeg=con$rhobeg, rhoend=con$rhoend, npt=npt, iprint=2*con$trace))
		}
		else
			fit = optim(start, ofv.FOCEi, control=list(trace=con$trace, reltol=con$reltol.outer))
	}
	fit$par.unscaled = fit$par*inits.vec
	attr(fit, "data") = data	#FIXME
	attr(fit, "ofv.FOCEi") = ofv.FOCEi	#FIXME
	attr(fit, "class") = "focei.fit"	
	fit
}
