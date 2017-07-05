//chk consistenty of par order bewteen models
#ifndef __STAN__MATH__FUNCTIONS__PKPDLIB_HPP__
#define __STAN__MATH__FUNCTIONS__PKPDLIB_HPP__

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <vector>
#include <utility>
#include <cmath>
#include <iostream>


namespace stan {
  namespace math {

    using std::exp;
    using stan::math::exp;
    using std::sqrt;
    using stan::math::sqrt;
    using std::pow;
    using stan::math::pow;
    using std::acos;
    using stan::math::acos;
    using std::cos;
    using stan::math::cos;


int locate_dose_index(const Eigen::VectorXd& dose_time, const double obs_time){
  int m = 0;
  while(m < dose_time.size() && dose_time[m] <= obs_time) m++;
  return m-1;
}//subscript of dose


//micro to macro conversion
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 2>
micros2macros(const Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
              const int ncmt,
              const int oral,
              const int parameterization){

  Eigen::Matrix<T, Eigen::Dynamic, 2> g(ncmt,2);

  if (ncmt==1) {

    T k     ;
    T volume;

    if (parameterization == 1) {
      //parameterization: CL V
      k       = params[0] / params[1];    // ke = CL/V
      volume  = params[1];
    }
    else if (parameterization == 2) {
      //parameterization: k V
      k       = params[0];
      volume  = params[1];
    }
    else{
      throw std::invalid_argument("Wrong parameterization type.");
    }

    const T alpha = k;
    const T A     = 1.0 / volume;

    g(0,0) = alpha;
    g(0,1) = A;

    if (oral==1) {
      const T ka = params[2];
      g(0,1) = ka / (ka - alpha) * A;
    }
  }

  if (ncmt==2) {

    T k     ;
    T volume;
    T k12   ;
    T k21   ;

    if (parameterization == 1) {
      //parameterization: CL V Q V2
      k       = params[0] / params[1];    // ke = CL/V
      volume  = params[1];
      k12     = params[2] / params[1];    // k12 = Q/V
      k21     = params[2] / params[3];    // k21 = Q/V2
    }
    else if (parameterization == 2) {
      //parameterization: k V k12 k21
      k       = params[0];
      volume  = params[1];
      k12     = params[2];
      k21     = params[3];
    }
    else{
      throw std::invalid_argument("Wrong parameterization type.");
    }

    const T beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
    const T alpha = k21 * k / beta;

    const T A     = (alpha - k21) / (alpha - beta) / volume;
    const T B     = (beta - k21) / (beta - alpha) / volume;

    g(0,0) = alpha;
    g(1,0) = beta;
    g(0,1) = A;
    g(1,1) = B;

    if (oral==1) {
      const T ka = params[4];
      g(0,1) = ka / (ka - alpha) * A;
      g(1,1) = ka / (ka - beta) * B;
    }
  }

  if (ncmt==3) {

    T k     ;
    T volume;
    T k12   ;
    T k21   ;
    T k13   ;
    T k31   ;

    if (parameterization == 1) {
      //parameterization: CL V Q V2 Q2 V3
      k       = params[0] / params[1];    // ke = CL/V
      volume  = params[1];
      k12     = params[2] / params[1];    // k12 = Q/V
      k21     = params[2] / params[3];    // k21 = Q/V2
      k13     = params[4] / params[1];    // k12 = Q2/V
      k31     = params[4] / params[5];    // k21 = Q2/V3
    }
    else if (parameterization == 2) {
      //parameterization: k V k12 k21 k13 k31
      k       = params[0];
      volume  = params[1];
      k12     = params[2];
      k21     = params[3];
      k13     = params[4];
      k31     = params[5];
    }
    else{
      throw std::invalid_argument("Wrong parameterization type.");
    }

    const T a0      = k * k21 * k31;
    const T a1      = k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
    const T a2      = k + k12 + k13 + k21 + k31;

    const T p       = a1 - a2 * a2 / 3.0;
    const T q       = 2.0 * a2 * a2 * a2 / 27.0 - a1 * a2 /3.0 + a0;

    const T r1      = sqrt(-p * p * p / 27.0);
    const T r2      = 2 * pow(r1 , 1.0 / 3.0);

/*
    // coding errors by Fei
    const T theta   = acos(-q / r1 / 6.0);

    const T alpha   = -(cos(theta) * r2 - a2 / 3.0);
    const T beta    = -(cos(theta + 2.0 / 3.0) * r2 - a2 / 3.0);
    const T gamma   = -(cos(theta + 4.0 / 3.0) * r2 - a2 / 3.0);
*/
    const T theta   = acos(-q / (2.0 * r1)) / 3.0;

    const T alpha   = -(cos(theta) * r2 - a2 / 3.0);
    const T beta    = -(cos(theta + 2.0 / 3.0 * M_PI) * r2 - a2 / 3.0);
    const T gamma   = -(cos(theta + 4.0 / 3.0 * M_PI) * r2 - a2 / 3.0);

    const T A       = (k21 - alpha) * (k31 - alpha) / (alpha - beta) / (alpha - gamma) / volume;
    const T B       = (k21 - beta) * (k31 - beta) / (beta - alpha) / (beta - gamma) / volume;
    const T C       = (k21 - gamma) * (k31 - gamma) / (gamma - alpha) / (gamma - beta) / volume;

    g(0,0) = alpha;
    g(1,0) = beta;
    g(2,0) = gamma;
    g(0,1) = A;
    g(1,1) = B;
    g(2,1) = C;

    if (oral==1) {
      const T ka = params[6];
      g(0,1) = ka / (ka - alpha) * A;
      g(1,1) = ka / (ka - beta) * B;
      g(2,1) = ka / (ka - gamma) * C;
    }
  }

  return g;
}


template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
generic_cmt_interface(const Eigen::VectorXd& obs_time,
                      const Eigen::VectorXd& dose_time,
                      const Eigen::VectorXd& dose,
                      const Eigen::VectorXd& Tinf,
                      const Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
                      const int ncmt,
                      const int oral,
                      const int infusion,
                      const int parameterization){

  stan::math::check_greater_or_equal("generic_cmt_interface", "params.size()", params.size(), 2*ncmt + 2*oral);
  T ka   = params[2 * ncmt];
  T Tlag = params[2 * ncmt + 1];
  if (oral != 1) Tlag = 0.0;    //i.v.

  Eigen::Matrix<T, Eigen::Dynamic, 2> par(ncmt, 2);
  par = micros2macros(params, ncmt, oral, parameterization);

  Eigen::Matrix<T, Eigen::Dynamic, 1> g(obs_time.size());
  for(int i=0;i<obs_time.size();i++){
    g(i,0) = 0;
  }

  for(int b = 0; b < obs_time.size(); ++b){

    int m = locate_dose_index(dose_time, obs_time[b]);
    for(int l=0; l<=m; l++){    //superpostion

      T this_t = obs_time[b] - dose_time[l] - Tlag;
      if (this_t < 0) continue;

      T sum = 0.0;
      if (infusion > 0) {
        T t1 = this_t < Tinf[l] ? this_t : Tinf[l];        //during infusion
        T t2 = this_t > Tinf[l] ? this_t - Tinf[l] : 0.0;  // after infusion

        for (int i = 0; i < ncmt; i++)
          sum += par(i,1) / par(i,0) * (1 - exp(-par(i,0) * t1)) * exp(-par(i,0) * t2);

        g(b,0) += dose[l] / Tinf[l] * sum;
      }
      else {
        T res = oral == 1 ? exp(-ka * this_t) : 0.0;
        for (int i = 0; i < ncmt; i++)
          sum += par(i,1) * (exp(-par(i,0) * this_t) - res);

        g(b,0) += dose[l] * sum;
      }
    } // l
  } // b

  return g;
}


template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
linear_cmpt_1order_absor(const Eigen::VectorXd& obs_time,
                         const Eigen::VectorXd& dose_time,
                         const Eigen::VectorXd& dose,
                         const Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
                         const int ncmt,
                         const int parameterization){

  const int oral = 1;
  const int infusion = 0;
  const Eigen::VectorXd Tinf(1);

  return generic_cmt_interface(obs_time, dose_time, dose, Tinf, params, ncmt, oral, infusion, parameterization);
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
linear_cmpt_iv_bolus(const Eigen::VectorXd& obs_time,
                     const Eigen::VectorXd& dose_time,
                     const Eigen::VectorXd& dose,
                     const Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
                     const int ncmt,
                     const int parameterization){

  const int oral = 0;
  const int infusion = 0;
  const Eigen::VectorXd Tinf(1);

  return generic_cmt_interface(obs_time, dose_time, dose, Tinf, params, ncmt, oral, infusion, parameterization);
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
linear_cmpt_iv_infusion(const Eigen::VectorXd& obs_time,
                        const Eigen::VectorXd& dose_time,
                        const Eigen::VectorXd& dose,
                        const Eigen::VectorXd& Tinf,
                        const Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
                        const int ncmt,
                        const int parameterization){

  const int oral = 0;
  const int infusion = 1;

  return generic_cmt_interface(obs_time, dose_time, dose, Tinf, params, ncmt, oral, infusion, parameterization);
}


struct lin_cmt_fun {
	const Eigen::VectorXd obs_time_, dose_time_, dose_, Tinf_;
	const int ncmt_, oral_, infusion_, parameterization_;
	lin_cmt_fun(const Eigen::VectorXd& obs_time,
			const Eigen::VectorXd& dose_time,
			const Eigen::VectorXd& dose,
			const Eigen::VectorXd& Tinf,
			const int ncmt,
			const int oral,
			const int infusion,
			const int parameterization) :
			obs_time_(obs_time),
			dose_time_(dose_time),
			dose_(dose),
			Tinf_(Tinf),
			ncmt_(ncmt),
			oral_(oral),
			infusion_(infusion),
			parameterization_(parameterization)
			{ }

	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& params) const {
		return generic_cmt_interface(obs_time_, dose_time_, dose_, Tinf_, params, ncmt_, oral_, infusion_, parameterization_);
	}
};


  } // ns math

}// ns stan



#endif
