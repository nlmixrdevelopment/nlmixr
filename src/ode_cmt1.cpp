#include "../inst/include/nlmixr_types.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

#include <vector>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math.hpp>
#include "PKPDLib_WW.h"

// Change to the Rcpp interface and let Rcpp handle the BEGIN_RCPP interface.
// It may have changed in the new RcppEigen...

//[[Rcpp::export]]
SEXP lin_cmt_stan(Eigen::Map<Eigen::VectorXd> obs_time,
		  Eigen::Map<Eigen::VectorXd> dose_time,
		  Eigen::Map<Eigen::VectorXd> dose,
		  Eigen::Map<Eigen::VectorXd> Tinf,
		  Eigen::Map<Eigen::VectorXd> params,
		  SEXP oralSEXP,
		  SEXP infusionSEXP,
		  SEXP ncmtSEXP,
		  SEXP parameterizationSEXP ) {
  const int oral = as<int>(oralSEXP);
  const int infusion = as<int>(infusionSEXP);
  const int ncmt = as<int>(ncmtSEXP);
  const int parameterization = as<int>(parameterizationSEXP);
  stan::math::lin_cmt_fun f(obs_time, dose_time, dose, Tinf, ncmt, oral, infusion, parameterization);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  // ncol = npar
  /// nrow = nobs
  stan::math::jacobian(f, params, fx, J);

  return Rcpp::List::create(Rcpp::Named("fx") = wrap(fx),
			    Rcpp::Named("J") = wrap(J));
}

extern void lin_cmt_stanC(double *obs_timeD, const int nobs, double *dose_timeD, const int ndose, double *doseD, double *TinfD,
			  double *paramsD, const int oral, const int infusion, const int ncmt, const int parameterization,
			  const int neta, double *fxD, double *dvdxD, double *fpD){
  Eigen::Map<Eigen::VectorXd> obs_time(obs_timeD, nobs);
  Eigen::Map<Eigen::VectorXd> dose_time(dose_timeD, ndose);
  Eigen::Map<Eigen::VectorXd> dose(doseD, ndose);
  Eigen::Map<Eigen::VectorXd> Tinf(TinfD, ndose);
  Eigen::Map<Eigen::VectorXd> params(paramsD, (int)(2*ncmt+2));
  stan::math::lin_cmt_fun f(obs_time, dose_time, dose, Tinf, ncmt, oral, infusion, parameterization);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  // Jacobian ncols=npars
  // nrows=nobs
  
  stan::math::jacobian(f, params, fx, J);
  std::copy(&fx[0],&fx[0]+nobs, fpD);
  // dvdx
  // ncol = netas
  // nrow = npars
  Eigen::Map<Eigen::MatrixXd> dvdx(dvdxD, 2*ncmt+2, neta);
  Eigen::Map<Eigen::MatrixXd> fp(fpD, nobs, neta);
  fp = J.transpose() * dvdx;
}


//===============================================================
struct binomial_llik {
  const Eigen::VectorXd y_, N_;
  binomial_llik(const Eigen::VectorXd& y, const Eigen::VectorXd& N) : y_(y), N_(N) { }
  
  template <typename T>
  Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    Eigen::Matrix<T, -1, 1> lp(y_.size());
    for (int n = 0; n < y_.size(); ++n)
      lp[n] = binomial_log(y_[n], N_[n], theta[n]);
    return lp;
  }
};

//[[Rcpp::export]]
SEXP llik_binomial_c(Eigen::Map<Eigen::VectorXd> y,
			  Eigen::Map<Eigen::VectorXd> N,
			  Eigen::Map<Eigen::VectorXd> params) {
    int i;
    for (i=0; i<params.size(); ++i) {
		if (params[i] > .99999) params[i] = .99999;
		if (params[i] < .00001) params[i] = .00001;
	}

	binomial_llik f(y, N);
	Eigen::VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
}


//===============================================================
struct poisson_llik {
  const Eigen::VectorXd y_;
  poisson_llik(const Eigen::VectorXd& y) : y_(y) { }
  
  template <typename T>
  Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    Eigen::Matrix<T, -1, 1> lp(y_.size());
    for (int n = 0; n < y_.size(); ++n)
      lp[n] = poisson_log(y_[n], theta[n]);
    return lp;
  }
};

//[[Rcpp::export]]
SEXP llik_poisson(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params) {
  poisson_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  return Rcpp::List::create(Rcpp::Named("fx") = fx,
			    Rcpp::Named("J") = J);
}


//===============================================================
struct normal_llik {
  const Eigen::VectorXd y_;
  normal_llik(const Eigen::VectorXd& y) : y_(y) { }

  template <typename T>
  Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    T mu = theta[0];
    T sigma = theta[1];
		
    if (sigma <= 0) {
      Rcpp::Rcout << "Warning: sigma <= 0" <<std::endl;
      sigma = 1.0e-12;
    }
		
    Eigen::Matrix<T, -1, 1> lp(y_.size());
    for (int n = 0; n < y_.size(); ++n)
      lp[n] = normal_log(y_[n], mu, sigma);
    return lp;
  }
};

//[[Rcpp::export]]
SEXP llik_normal(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params) {
  normal_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);
  
  return Rcpp::List::create(Rcpp::Named("fx") = fx,
			    Rcpp::Named("J") = J);
}


//===============================================================
struct betabinomial_llik {
	const Eigen::VectorXd y_, N_;
	betabinomial_llik(const Eigen::VectorXd& y, const Eigen::VectorXd& N) : y_(y), N_(N) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T alpha = theta[0];
		T beta  = theta[1];
		if (alpha <= 0) {
			Rcpp::Rcout << "Warning: alpha <= 0" <<std::endl;
			alpha = 1.0e-12;
		}
		if (beta <= 0) {
			Rcpp::Rcout << "Warning: beta <= 0" <<std::endl;
			beta = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = beta_binomial_log(y_[n], N_[n], alpha, beta);
		return lp;
	}
};

//[[Rcpp::export]]
SEXP llik_betabinomial(Eigen::Map<Eigen::VectorXd> y,
			      Eigen::Map<Eigen::VectorXd> N,
			      Eigen::Map<Eigen::VectorXd> params) {
  betabinomial_llik f(y, N);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);
	
  return Rcpp::List::create(Rcpp::Named("fx") = fx,
			    Rcpp::Named("J") = J);
}

//===============================================================
struct student_t_llik {
	const Eigen::VectorXd y_;
	student_t_llik(const Eigen::VectorXd& y) : y_(y) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T nu = theta[0];
		T mu = theta[1];
		T sigma = theta[2];

		if (nu <= 0) {
			Rcpp::Rcout << "Warning: nu <= 0" <<std::endl;
			nu = 1.0e-12;	//FIXME
		}
		if (sigma <= 0) {
			Rcpp::Rcout << "Warning: sigma <= 0" <<std::endl;
			sigma = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = student_t_log(y_[n], nu, mu, sigma);
		return lp;
	}
};

//[[Rcpp::export]]
SEXP llik_student_t(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params) {
  student_t_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  return Rcpp::List::create(Rcpp::Named("fx") = fx,
			    Rcpp::Named("J") = J);
}

//===============================================================
struct beta_llik {
	const Eigen::VectorXd y_;
	beta_llik(const Eigen::VectorXd& y) : y_(y) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T alpha = theta[0];
		T beta = theta[1];

		if (alpha <= 0) {
			Rcpp::Rcout << "Warning: alpha <= 0" <<std::endl;
			alpha = 1.0e-12;
		}
		if (beta <= 0) {
			Rcpp::Rcout << "Warning: beta <= 0" <<std::endl;
			beta = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = beta_log(y_[n], alpha, beta);
		return lp;
	}
};

//[[Rcpp::export]]
SEXP llik_beta(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params) {
  beta_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);
  return Rcpp::List::create(Rcpp::Named("fx") = fx,
			    Rcpp::Named("J") = J);
}


//===============================================================
struct neg_binomial_llik {
	const Eigen::VectorXd y_;
	neg_binomial_llik(const Eigen::VectorXd& y) : y_(y) { }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
		T alpha = theta[0];
		T beta = theta[1];

		if (alpha <= 0) {
			Rcpp::Rcout << "Warning: alpha <= 0" <<std::endl;
			alpha = 1.0e-12;
		}
		if (beta <= 0) {
			Rcpp::Rcout << "Warning: beta <= 0" <<std::endl;
			beta = 1.0e-12;
		}

		Eigen::Matrix<T, -1, 1> lp(y_.size());
		for (int n = 0; n < y_.size(); ++n)
		lp[n] = neg_binomial_log(y_[n], alpha, beta);
		return lp;
	}
};

//[[Rcpp::export]]
SEXP llik_neg_binomial(Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> params) {
  neg_binomial_llik f(y);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  stan::math::jacobian(f, params, fx, J);

  return Rcpp::List::create(Rcpp::Named("fx") = fx,
			    Rcpp::Named("J") = J);
}


// I'm not sure why the below is here...

#if 0

require(Rcpp)
dyn.load("nlmixr.so")
lin_cmt <- function(obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)
   .Call('lin_cmt_stan', obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)
a = lin_cmt(0:72*1.0, 0:1*24*1.0, rep(10.0,2), 0, c(.1, 1, .2, 0), 1, 0, 1, 1)
b = matrix(scan("1"),  nrow=73)
range(b[,1] - a$fx)
range(as.vector(b[,-1]) - a$J)


ode_sol <- function(inits,time,evid,amt,params,absolute_tolerance,relative_tolerance,nobs,wh)
   .Call('ode_sol', inits,time,evid,amt,params,absolute_tolerance,relative_tolerance,nobs,wh)

require(RxODE)
e = eventTable()
e$add.dosing(10, 2, 24)
e$add.sampling(0:72)
x = e$get.EventTable()
x$amt[is.na(x$amt)] = 0

ode_sol(c(0,0), x$time, x$evid, x$amt, c(.2, .1), 1e-8, 1e-8, 73, 2)

#endif

