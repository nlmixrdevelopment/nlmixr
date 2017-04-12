// includes from the plugin
#ifdef __STANDALONE__
#include <Eigen/Dense>
#else
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif
#endif

#include <stan/math.hpp>


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

RcppExport SEXP llik_binomial(SEXP ySEXP, SEXP NSEXP, SEXP paramsSEXP) {
BEGIN_RCPP

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type y_(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type N_(NSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params_(paramsSEXP);

    const VectorXd y(y_);
    const VectorXd N(N_);
    VectorXd params(params_);

    int i;
    for (i=0; i<params.size(); ++i) {
		if (params[i] > .99999) params[i] = .99999;
		if (params[i] < .00001) params[i] = .00001;
	}

	binomial_llik f(y, N);
    VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
END_RCPP
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

RcppExport SEXP llik_poisson(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type y_(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params_(paramsSEXP);

    const VectorXd y(y_);
    const VectorXd params(params_);

	poisson_llik f(y);
    VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
END_RCPP
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

RcppExport SEXP llik_normal(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type y_(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params_(paramsSEXP);

    const VectorXd y(y_);
    const VectorXd params(params_);

	normal_llik f(y);
    VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
END_RCPP
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

RcppExport SEXP llik_betabinomial(SEXP ySEXP, SEXP NSEXP, SEXP paramsSEXP) {
BEGIN_RCPP

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type y_(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type N_(NSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params_(paramsSEXP);

    const VectorXd y(y_);
    const VectorXd N(N_);
    const VectorXd params(params_);

	betabinomial_llik f(y, N);
    VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
END_RCPP
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

RcppExport SEXP llik_student_t(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type y_(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params_(paramsSEXP);

    const VectorXd y(y_);
    const VectorXd params(params_);

	student_t_llik f(y);
    VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
END_RCPP
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

RcppExport SEXP llik_beta(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type y_(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params_(paramsSEXP);

    const VectorXd y(y_);
    const VectorXd params(params_);

	beta_llik f(y);
    VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
END_RCPP
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

RcppExport SEXP llik_neg_binomial(SEXP ySEXP, SEXP paramsSEXP) {
BEGIN_RCPP

    using Eigen::VectorXd;
    Rcpp::traits::input_parameter< const VectorXd& >::type y_(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type params_(paramsSEXP);

    const VectorXd y(y_);
    const VectorXd params(params_);

	neg_binomial_llik f(y);
    VectorXd fx;
	Eigen::Matrix<double, -1, -1> J;
	stan::math::jacobian(f, params, fx, J);

    return Rcpp::List::create(Rcpp::Named("fx") = fx,
	                          Rcpp::Named("J") = J);
END_RCPP
}



#if 0

require(Rcpp)
dyn.load("ode_cmt1.dll")
lin_cmt <- function(obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)
   .Call('lin_cmt', obs_time,dose_time,dose,Tinf,params,oral,infusion,ncmt,parameterization)
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

