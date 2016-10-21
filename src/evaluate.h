// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2015  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#ifndef Rcpp_DE_evaluate_h_
#define Rcpp_DE_evaluate_h_

#include <Rcpp.h>

namespace Rcpp {

        double genrose(SEXP xs) {       // genrose function in C++
            Rcpp::NumericVector x(xs);
            int n = x.size();
            double sum = 1.0;
            for (int i=1; i<n; i++) {
                sum += 100*( ::pow(x[i-1]*x[i-1] - x[i], 2)) + (x[i] - 1)*(x[i] - 1);
            }
            return(sum);
        }

        double wild(SEXP xs) {          // wild function in C++
            Rcpp::NumericVector x(xs);
            int n = x.size();
            double sum = 0.0;
            for (int i=0; i<n; i++) {
                double xsq = x[i]*x[i];
                sum += 10 * ::sin(0.3 * x[i]) * ::sin(1.3 * xsq) + 0.00001 * xsq*xsq + 0.2 * x[i] + 80;
            }
            sum /= n;
            return(sum);
        }

        double rastrigin(SEXP xs) {     // rastrigin function in C++
            Rcpp::NumericVector x(xs);
            int n = x.size();
            double sum = 20.0;
            for (int i=0; i<n; i++) {
                sum += x[i]+2 - 10*::cos(M_2PI*x[i]);
            }
            return(sum);
        }

        class EvalBase {
        public:
            EvalBase() : neval(0) {};
            virtual double eval(SEXP par) = 0;
            unsigned long getNbEvals() { return neval; }
        protected:
            unsigned long int neval;
        };

        class EvalStandard : public EvalBase {
        public:
            EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
            double eval(SEXP par) {
                neval++;
                return defaultfun(par);
            }
        private:
            SEXP fcall, env;
            double defaultfun(SEXP par) {                       // essentialy same as the old evaluate
              //SEXP fn = ::Rf_lang3(fcall, par, R_DotsSymbol); // this could be done with Rcpp
                SEXP fn = ::Rf_lang2(fcall, par); // this could be done with Rcpp
                SEXP sexp_fvec = ::Rf_eval(fn, env);            // but is still a lot slower right now
                double f_result = REAL(sexp_fvec)[0];
                if (ISNAN(f_result))
                    ::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
                return(f_result);
            }
        };

        typedef double (*funcPtr)(SEXP);
        class EvalCompiled : public EvalBase {
        public:
            EvalCompiled( Rcpp::XPtr<funcPtr> xptr ) {
                funptr = *(xptr);
            };
            EvalCompiled( SEXP xps ) {
                Rcpp::XPtr<funcPtr> xptr(xps);
                funptr = *(xptr);
            };
            double eval(SEXP par) {
                neval++;
                return funptr(par);
            }
        private:
            funcPtr funptr;
        };

        RcppExport SEXP putFunPtrInXPtr(SEXP funname) {
            std::string fstr = Rcpp::as<std::string>(funname);
            if (fstr == "genrose")
                return(Rcpp::XPtr<funcPtr>(new funcPtr(&genrose)));
            else if (fstr == "wild")
                return(Rcpp::XPtr<funcPtr>(new funcPtr(&wild)));
            else
                return(Rcpp::XPtr<funcPtr>(new funcPtr(&rastrigin)));
        }


}

#endif
