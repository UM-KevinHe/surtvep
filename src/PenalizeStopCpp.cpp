#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_OPENMP
#define STRICT_R_HEADERS // needed on Windows, not macOS
#include <RcppArmadillo.h>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h> // for Rcpp::RcppArmadillo::sample
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// using namespace arma;
using namespace std;


// [[Rcpp::export]]
arma::uvec stra_sampling_cpp(int size, int prop){
  arma::uvec sample_id;
  
  int s2      = floor(size/prop*1.0);
  sample_id   = arma::randperm(size,s2);

  return sample_id;
}


List objfun_fixtra(const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta,
                   const arma::mat &Z_ti, const arma::vec &beta_ti, const bool &ti,
                   const unsigned int n_strata,
                   vector<arma::uvec> &idx_B_sp, vector<arma::uvec> &idx_fail,
                   vector<unsigned int> n_Z_strata,
                   vector<vector<unsigned int>> &idx_Z_strata,
                   vector<vector<unsigned int>> &istart,
                   vector<vector<unsigned int>> &iend,
                   const bool &parallel=false, const unsigned int &threads=1) {

  double norm_parm;
  vector<arma::vec>  hazard;
  if (ti) {
    norm_parm = max(arma::norm(theta, "inf"), arma::norm(beta_ti, "inf"));
  } else {
    norm_parm = arma::norm(theta, "inf");
  }
  double logplkd = 0.0;
  if (norm_parm < sqrt(arma::datum::eps)) { // theta and beta_ti are 0
    if (parallel) {
      for (unsigned int i = 0; i < n_strata; ++i) {
        arma::vec S0_fail(idx_fail[i].n_elem);
        omp_set_num_threads(threads);
        #pragma omp parallel
        {
          unsigned int id, j;
          id = omp_get_thread_num();
          double val_tmp = 0;
          for (j = istart[i][id]; j < iend[i][id]; ++j) {
            S0_fail(j) = n_Z_strata[i]-idx_fail[i](j);
            val_tmp += log(S0_fail(j));
          }
          #pragma omp atomic
          logplkd -= val_tmp;
        }
        hazard.push_back(1/S0_fail);
      }
    } else {
      for (unsigned int i = 0; i < n_strata; ++i) {
        arma::vec S0_fail(idx_fail[i].n_elem);
        for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
          S0_fail(j) = n_Z_strata[i]-idx_fail[i](j);
          logplkd -= log(S0_fail(j));
        }
        hazard.push_back(1/S0_fail);
      }
    }
  } else if (max(var(theta, 0, 1)) < sqrt(arma::datum::eps)) { // each row of theta is const
    for (unsigned int i = 0; i < n_strata; ++i) {
      arma::vec Z_tv_theta =
        Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * theta.col(0);
      arma::vec Z_ti_beta_ti;
      if (ti) {
        Z_ti_beta_ti =
          Z_ti.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * beta_ti;
      }
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
      arma::vec lincomb_fail(idx_fail[i].n_elem), S0_fail(idx_fail[i].n_elem);
      if (parallel) {
        omp_set_num_threads(threads);
        #pragma omp parallel
        {
          unsigned int id, j;
          id = omp_get_thread_num();
          double val_tmp = 0;
          if (ti) {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.subvec(idx_fail[i](j), n_Z_tv_theta-1) *
                accu(B_sp.row(j)) +
                Z_ti_beta_ti.subvec(idx_fail[i](j),n_Z_tv_theta-1);
              lincomb_fail(j) = lincomb(0);
              S0_fail(j) = sum(exp(lincomb));
              val_tmp += log(S0_fail(j));
            }
          } else {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.subvec(idx_fail[i](j), n_Z_tv_theta-1) *
                accu(B_sp.row(j));
              lincomb_fail(j) = lincomb(0);
              S0_fail(j) = sum(exp(lincomb));
              val_tmp += log(S0_fail(j));
            }
          }
          #pragma omp atomic
          logplkd -= val_tmp;
        }
      } else {
        if (ti) {
          for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
            arma::vec lincomb =
              Z_tv_theta.subvec(idx_fail[i](j),n_Z_tv_theta-1) *
              accu(B_sp.row(j)) +
              Z_ti_beta_ti.subvec(idx_fail[i](j),n_Z_tv_theta-1);
            lincomb_fail(j) = lincomb(0);
            S0_fail(j) = sum(exp(lincomb));
            logplkd -= log(S0_fail(j));
          }
        } else {
          for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
            arma::vec lincomb =
              Z_tv_theta.subvec(idx_fail[i](j), n_Z_tv_theta-1) *
              accu(B_sp.row(j));
            lincomb_fail(j) = lincomb(0);
            S0_fail(j) = sum(exp(lincomb));
            logplkd -= log(S0_fail(j));
          }
        }
      }
      logplkd += accu(lincomb_fail);
      hazard.push_back(1/S0_fail);
    }
  } else { // general theta
    for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_theta =
        Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * theta;
      arma::vec Z_ti_beta_ti;
      if (ti) {
        Z_ti_beta_ti =
          Z_ti.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * beta_ti;
      }
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
      arma::vec lincomb_fail(idx_fail[i].n_elem), S0_fail(idx_fail[i].n_elem);
      if (parallel) {
        omp_set_num_threads(threads);
        #pragma omp parallel
        {
          unsigned int id, j;
          id = omp_get_thread_num();
          double val_tmp = 0;
          if (ti) {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.rows(idx_fail[i](j), n_Z_tv_theta-1) *
                B_sp.row(j).t() +
                Z_ti_beta_ti.subvec(idx_fail[i](j),n_Z_tv_theta-1);
              lincomb_fail(j) = lincomb(0);
              S0_fail(j) = sum(exp(lincomb));
              val_tmp += log(S0_fail(j));
            }
          } else {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.rows(idx_fail[i](j), n_Z_tv_theta-1) *
                B_sp.row(j).t();
              lincomb_fail(j) = lincomb(0);
              S0_fail(j) = sum(exp(lincomb));
              val_tmp += log(S0_fail(j));
            }
          }
          #pragma omp atomic
          logplkd -= val_tmp;
        }
      } else {
        if (ti) {
          for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
            arma::vec lincomb =
              Z_tv_theta.rows(idx_fail[i](j), n_Z_tv_theta-1) * 
              B_sp.row(j).t() +
              Z_ti_beta_ti.subvec(idx_fail[i](j),n_Z_tv_theta-1);
            lincomb_fail(j) = lincomb(0);
            S0_fail(j) = sum(exp(lincomb));
            logplkd -= log(S0_fail(j));
          }
        } else {
          for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
            arma::vec lincomb =
              Z_tv_theta.rows(idx_fail[i](j), n_Z_tv_theta-1) * 
              B_sp.row(j).t();
            lincomb_fail(j) = lincomb(0);
            S0_fail(j) = sum(exp(lincomb));
            logplkd -= log(S0_fail(j));
          }
        }
      }
      logplkd += accu(lincomb_fail);
      hazard.push_back(1/S0_fail);
    }
  }
  logplkd /= Z_tv.n_rows;
  return List::create(_["logplkd"]=logplkd, _["hazard"]=hazard);
}


List obj_fixtra_bresties(const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta,
                         const arma::mat &Z_ti, const arma::vec &beta_ti, const bool &ti,
                         const unsigned int n_strata,
                         vector<arma::uvec> &idx_B_sp, vector<vector<arma::uvec>> &idx_fail,
                         vector<unsigned int> n_Z_strata,
                         vector<vector<unsigned int>> &idx_Z_strata,
                         vector<vector<unsigned int>> &istart,
                         vector<vector<unsigned int>> &iend,
                         const bool &parallel=false, const unsigned int &threads=1) {
  
  double norm_parm;
  vector<arma::vec>  hazard;
  if (ti) {
    norm_parm = max(arma::norm(theta, "inf"), arma::norm(beta_ti, "inf"));
  } else {
    norm_parm = arma::norm(theta, "inf");
  }
  double logplkd = 0.0;
  if (norm_parm < sqrt(arma::datum::eps)) { // theta and beta_ti are 0
    if (parallel) {
      for (unsigned int i = 0; i < n_strata; ++i) {
        arma::vec hazard_tmp(idx_fail[i].size());
        omp_set_num_threads(threads);
        #pragma omp parallel
        {
          unsigned int id, j;
          id = omp_get_thread_num();
          double val_tmp = 0;
          for (j = istart[i][id]; j < iend[i][id]; ++j) {
            double tmp = n_Z_strata[i]-idx_fail[i][j](0);
            hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
            val_tmp += idx_fail[i][j].n_elem*log(tmp);
          }
          #pragma omp atomic
          logplkd -= val_tmp;
        }
        hazard.push_back(hazard_tmp);
      }
    } else {
      for (unsigned int i = 0; i < n_strata; ++i) {
        arma::vec hazard_tmp(idx_fail[i].size());
        for (unsigned int j = 0; j < idx_fail[i].size(); ++j) {
          double tmp = n_Z_strata[i]-idx_fail[i][j](0);
          hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
          logplkd -= idx_fail[i][j].n_elem*log(tmp);
        }
        hazard.push_back(hazard_tmp);
      }
    }
  } else if (max(var(theta, 0, 1)) < sqrt(arma::datum::eps)) { // each row of theta is const
    for (unsigned int i = 0; i < n_strata; ++i) {
      arma::vec Z_tv_theta =
        Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * theta.col(0);
      arma::vec Z_ti_beta_ti;
      if (ti) {
        Z_ti_beta_ti =
          Z_ti.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * beta_ti;
      }
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
      arma::vec lincomb_fail(idx_fail[i].size()), hazard_tmp(idx_fail[i].size());
      if (parallel) {
        omp_set_num_threads(threads);
        #pragma omp parallel
        {
          unsigned int id, j;
          id = omp_get_thread_num();
          double val_tmp = 0;
          if (ti) {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.subvec(idx_fail[i][j](0), n_Z_tv_theta-1) *
                accu(B_sp.row(j)) +
                Z_ti_beta_ti.subvec(idx_fail[i][j](0),n_Z_tv_theta-1);
              lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
              double tmp = sum(exp(lincomb));
              hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
              val_tmp += idx_fail[i][j].n_elem*log(tmp);
            }
          } else {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.subvec(idx_fail[i][j](0), n_Z_tv_theta-1) *
                accu(B_sp.row(j));
              lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
              double tmp = sum(exp(lincomb));
              hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
              val_tmp += idx_fail[i][j].n_elem*log(tmp);
            }
          }
          #pragma omp atomic
          logplkd -= val_tmp;
        }
      } else {
        if (ti) {
          for (unsigned int j = 0; j < idx_fail[i].size(); ++j) {
            arma::vec lincomb =
              Z_tv_theta.subvec(idx_fail[i][j](0),n_Z_tv_theta-1) *
              accu(B_sp.row(j)) +
              Z_ti_beta_ti.subvec(idx_fail[i][j](0),n_Z_tv_theta-1);
            lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
            double tmp = sum(exp(lincomb));
            hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
            logplkd -= idx_fail[i][j].n_elem*log(tmp);
          }
        } else {
          for (unsigned int j = 0; j < idx_fail[i].size(); ++j) {
            arma::vec lincomb =
              Z_tv_theta.subvec(idx_fail[i][j](0), n_Z_tv_theta-1) *
              accu(B_sp.row(j));
            lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
            double tmp = sum(exp(lincomb));
            hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
            logplkd -= idx_fail[i][j].n_elem*log(tmp);
          }
        }
      }
      logplkd += accu(lincomb_fail);
      hazard.push_back(hazard_tmp);
    }
  } else { // general theta
    for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_theta =
        Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * theta;
      arma::vec Z_ti_beta_ti;
      if (ti) {
        Z_ti_beta_ti =
          Z_ti.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]) * beta_ti;
      }
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
      arma::vec lincomb_fail(idx_fail[i].size()), hazard_tmp(idx_fail[i].size());
      if (parallel) {
        omp_set_num_threads(threads);
        #pragma omp parallel
        {
          unsigned int id, j;
          id = omp_get_thread_num();
          double val_tmp = 0;
          if (ti) {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.rows(idx_fail[i][j](0), n_Z_tv_theta-1) *
                B_sp.row(j).t() +
                Z_ti_beta_ti.subvec(idx_fail[i][j](0),n_Z_tv_theta-1);
              lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
              double tmp = sum(exp(lincomb));
              hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
              val_tmp += idx_fail[i][j].n_elem*log(tmp);
            }
          } else {
            for (j = istart[i][id]; j < iend[i][id]; ++j) {
              arma::vec lincomb =
                Z_tv_theta.rows(idx_fail[i][j](0), n_Z_tv_theta-1) *
                B_sp.row(j).t();
              lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
              double tmp = sum(exp(lincomb));
              hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
              val_tmp += idx_fail[i][j].n_elem*log(tmp);
            }
          }
          #pragma omp atomic
          logplkd -= val_tmp;
        }
      } else {
        if (ti) {
          for (unsigned int j = 0; j < idx_fail[i].size(); ++j) {
            arma::vec lincomb =
              Z_tv_theta.rows(idx_fail[i][j](0), n_Z_tv_theta-1) * 
              B_sp.row(j).t() +
              Z_ti_beta_ti.subvec(idx_fail[i][j](0),n_Z_tv_theta-1);
            lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
            double tmp = sum(exp(lincomb));
            hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
            logplkd -= idx_fail[i][j].n_elem*log(tmp);
          }
        } else {
          for (unsigned int j = 0; j < idx_fail[i].size(); ++j) {
            arma::vec lincomb =
              Z_tv_theta.rows(idx_fail[i][j](0), n_Z_tv_theta-1) * 
              B_sp.row(j).t();
            lincomb_fail(j) = sum(lincomb.head(idx_fail[i][j].n_elem));
            double tmp = sum(exp(lincomb));
            hazard_tmp(j) = idx_fail[i][j].n_elem/tmp;
            logplkd -= idx_fail[i][j].n_elem*log(tmp);
          }
        }
      }
      logplkd += accu(lincomb_fail);
      hazard.push_back(hazard_tmp);
    }
  }
  logplkd /= Z_tv.n_rows;
  return List::create(_["logplkd"]=logplkd, _["hazard"]=hazard);
}

arma::mat spline_construct(const int knot,
                     const int p,
                     const std::string SplineType = "pspline"){

  arma::mat S_matrix     = arma::zeros<arma::mat>(p*knot, p*knot);

    arma::mat P_pre  = arma::zeros<arma::mat>(knot,knot);
    P_pre.diag().ones();
    P_pre      = diff(P_pre);
    arma::mat S_pre  = P_pre.t()*P_pre;

    S_matrix     = arma::zeros<arma::mat>(p*knot, p*knot);
    for (int i = 0; i < p; ++i)
    {
      S_matrix.submat(i*knot,i*knot,i*knot+knot-1, i*knot+knot-1) = S_pre;
    }

  return S_matrix;
}


arma::mat spline_construct2(const int knot,
                      const int p,
                      const std::string SplineType,
                      const arma::mat &SmoothMatrix){

  arma::mat S_matrix     = arma::zeros<arma::mat>(p*knot, p*knot);

  arma::mat S_pre    = SmoothMatrix;
  S_matrix     = arma::zeros<arma::mat>(p*knot, p*knot);
  for (int i = 0; i < p; ++i)
  {
    S_matrix.submat(i*knot,i*knot,i*knot+knot-1, i*knot+knot-1) = S_pre;
  }
  

  return S_matrix;
}

List stepinc_fixtra_spline(const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta,
                          const arma::mat &Z_ti, const arma::vec &beta_ti, 
                          arma::mat &S_matrix,
                          double &lambda_i,
                          arma::mat &lambda_i_mat,
                          const bool &difflambda,
                          const bool &ti,
                          const unsigned int n_strata,
                          vector<arma::uvec> &idx_B_sp, vector<arma::uvec> &idx_fail,
                          vector<vector<unsigned int>> &idx_Z_strata,
                          vector<vector<unsigned int>> &istart,
                          vector<vector<unsigned int>> &iend,
                          const std::string &method="Newton", const double &lambda=1e8,
                          const bool &parallel=false, const unsigned int &threads=1) {
  
  // int N = Z_tv.n_rows;

  arma::vec grad, grad_p; arma::mat info, info_p; // gradient and info matrix
  if (ti) {
    grad = arma::zeros<arma::vec> (theta.n_elem+beta_ti.n_elem);
    info = arma::zeros<arma::mat>(theta.n_elem+beta_ti.n_elem, theta.n_elem+beta_ti.n_elem);
  } else {
    grad = arma::zeros<arma::vec> (theta.n_elem);
    info = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  }
  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
    arma::mat Z_ti_strata; arma::vec Z_ti_beta_ti;
    if (ti) {
      Z_ti_strata = Z_ti.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
      Z_ti_beta_ti = Z_ti_strata * beta_ti;
    }
    arma::mat Z_tv_theta = Z_tv_strata * theta;
    arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
    unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
    if (parallel) {
      omp_set_num_threads(threads);
      #pragma omp parallel
      {
        unsigned int id, j;
        id = omp_get_thread_num();
        arma::vec grad_tmp(size(grad), arma::fill::zeros);
        arma::mat info_tmp(size(info), arma::fill::zeros);
        for (j = istart[i][id]; j < iend[i][id]; ++j) {
          arma::vec B_sp_tmp = B_sp.row(j).t();
          unsigned int arstart = idx_fail[i](j), arend = n_Z_tv_theta-1;
          arma::vec exp_lincomb = exp(Z_tv_theta.rows(arstart,arend) * B_sp_tmp);
          double S0 = accu(exp_lincomb);
          arma::mat Z_tv_exp =
            Z_tv_strata.rows(arstart,arend).each_col() % exp_lincomb;
          arma::vec S1_tv = sum(Z_tv_exp).t();
          arma::mat S2 = Z_tv_strata.rows(arstart,arend).t() * Z_tv_exp;
          grad_tmp += kron(Z_tv_strata.row(arstart).t()-S1_tv/S0, B_sp_tmp);
          info_tmp +=
            kron(S2/S0-S1_tv*S1_tv.t()/pow(S0,2),B_sp_tmp*B_sp_tmp.t());
        }
        #pragma omp critical (gradient)
          grad += grad_tmp;
        #pragma omp critical (information)
          info += info_tmp;
      }
    } else {
      for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
        arma::vec B_sp_tmp = B_sp.row(j).t();
        unsigned int arstart = idx_fail[i](j), arend = n_Z_tv_theta-1;
        arma::vec exp_lincomb = exp(Z_tv_theta.rows(arstart,arend) * B_sp_tmp);
        double S0 = accu(exp_lincomb);
        arma::mat Z_tv_exp =
          Z_tv_strata.rows(arstart,arend).each_col() % exp_lincomb;
        arma::vec S1_tv = sum(Z_tv_exp).t();
        arma::mat S2 = Z_tv_strata.rows(arstart,arend).t() * Z_tv_exp;
        grad += kron(Z_tv_strata.row(arstart).t()-S1_tv/S0, B_sp_tmp);
        info += kron(S2/S0-S1_tv*S1_tv.t()/pow(S0,2),B_sp_tmp*B_sp_tmp.t());
      }
    }
  }

  arma::vec step; // Newton step
  // info_p  = info + N*lambda_i*S_matrix;
  // grad_p  = grad - N*lambda_i*S_matrix*vectorise(theta, 1);

  if(difflambda == false){
    info_p  = info + lambda_i*S_matrix;
    grad_p  = grad - lambda_i*S_matrix*vectorise(theta.t(), 0);
  }
  else{
    info_p  = info + lambda_i_mat%S_matrix;
    grad_p  = grad - (lambda_i_mat%S_matrix)*vectorise(theta.t(), 0);
  }

  // add proximal
  if (method=="ProxN") {
    info_p.diag() += Z_tv.n_rows / lambda;
    step = solve(info_p, grad_p, arma::solve_opts::fast+arma::solve_opts::likely_sympd);
    info_p.diag() -= Z_tv.n_rows / lambda;
  } else if (method=="Newton") {
    step = solve(info_p, grad_p, arma::solve_opts::fast+arma::solve_opts::likely_sympd);
  }
  
  double inc = dot(grad_p, step) / Z_tv.n_rows; //increment, notice how to calculate in later part
  if (ti) {
    return List::create(_["step_tv"]=reshape(step.head(theta.n_elem),
                                      size(theta.t())).t(),
                        _["step_ti"]=step.tail(beta_ti.n_elem),
                        _["grad"]=grad, _["info"]=info, _["inc"]=inc);
  } else {
    return List::create(_["step"]=reshape(step, size(theta.t())).t(),
                        _["grad"]=grad, _["info"]=info, _["inc"]=inc);
  }
}

List TIC_J_penalized(const arma::mat &Z_tv, const arma::mat &B_spline, arma::mat &theta,
                      const unsigned int n_strata,
                      vector<arma::uvec> &idx_B_sp, vector<arma::uvec> &idx_fail,
                      vector<vector<unsigned int>> &idx_Z_strata,
                      const bool &TIC_prox,
                      double lambda_i,
                      arma::mat &lambda_i_mat,
                      const bool &difflambda,
                      arma::mat S_matrix) {

  // gradient and info matrix
  int N       = Z_tv.n_rows;
  arma::vec grad    = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_p  = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_tmp, grad_p_tmp;
  arma::mat info_J   = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  arma::mat info_J_p = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  arma::mat info_J_p_gic = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);


  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);

    arma::mat Z_tv_theta = Z_tv_strata * theta;
    arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
    unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;

    for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
      arma::vec B_sp_tmp = B_sp.row(j).t();
      unsigned int arstart = idx_fail[i](j), arend = n_Z_tv_theta-1;
      arma::vec exp_lincomb = exp(Z_tv_theta.rows(arstart,arend) * B_sp_tmp);
      double S0 = accu(exp_lincomb);
      arma::mat Z_tv_exp =
        Z_tv_strata.rows(arstart,arend).each_col() % exp_lincomb;
      arma::vec S1_tv    = sum(Z_tv_exp).t();
      grad_tmp     = kron(Z_tv_strata.row(arstart).t()-S1_tv/S0, B_sp_tmp);
      if(difflambda == false){
        grad_p_tmp   = grad_tmp - lambda_i*S_matrix*vectorise(theta.t(), 0)/N;
      }
      else{
        grad_p_tmp   = grad_tmp - lambda_i_mat%S_matrix*vectorise(theta.t(), 0)/N;
      }
      grad += grad_tmp;
      info_J       += grad_tmp*grad_tmp.t();
      info_J_p     += grad_p_tmp*grad_p_tmp.t();
      info_J_p_gic += grad_p_tmp*grad_tmp.t();
    }
  }

  return List::create(_["info_J"]=info_J,
                      _["info_J_p"]=info_J_p,
                      _["info_J_p_gic"]=info_J_p_gic);
}

List TIC_J_penalized_second(const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta,
                            const unsigned int n_strata,
                            vector<arma::uvec> &idx_B_sp, vector<arma::uvec> &idx_fail,
                            vector<vector<unsigned int>> &idx_Z_strata,
                            const bool &TIC_prox,
                            double lambda_i,
                            arma::mat &lambda_i_mat,
                            const bool &difflambda,
                            arma::mat S_matrix) {

  // gradient and info matrix
  int N     = Z_tv.n_rows;    //sample size 
  int p     = theta.n_rows; //dimension
  int K     = theta.n_cols; //number of knots   

  arma::vec grad    = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_p  = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_tmp, grad_p_tmp;

  arma::vec grad_tmp1 = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_tmp2 = arma::zeros<arma::vec> (theta.n_elem); 
  arma::mat info_J   = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  arma::mat info_J_p = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  arma::mat info_J_p_gic = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);

  arma::mat grad_part1 = arma::zeros<arma::mat>(p*K,N);
  arma::mat grad_part2 = arma::zeros<arma::mat>(p*K,N);
  arma::mat grad_all = arma::zeros<arma::mat>(p*K,N);

  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
    arma::mat Z_tv_theta = Z_tv_strata * theta;
    arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
    unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
    for (unsigned int j = 0; j < idx_fail[i].n_elem; ++j) {
      arma::vec B_sp_tmp = B_sp.row(j).t();
      unsigned int arstart = idx_fail[i](j), arend = n_Z_tv_theta-1;
      arma::vec exp_lincomb = exp(Z_tv_theta.rows(arstart,arend) * B_sp_tmp);
      double S0 = accu(exp_lincomb);
      arma::mat Z_tv_exp =
        Z_tv_strata.rows(arstart,arend).each_col() % exp_lincomb;
      arma::vec S1_tv = sum(Z_tv_exp).t();
      grad_tmp1 = kron(Z_tv_strata.row(arstart).t()-S1_tv/S0, B_sp_tmp);
      grad += grad_tmp1;
      grad_part1.col(arstart) = grad_tmp1;
    }
  }

  //calculate S0 and S1:
  NumericVector S0_second_all;
  vector<arma::vec>  S1_second_all;
  for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
      arma::mat Z_tv_theta = Z_tv_strata * theta;
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;

      for (unsigned int k = 0; k < idx_fail[i].size(); ++k){
        arma::vec B_sp_second_tmp = B_sp.row(k).t();
        unsigned int arstart_second = idx_fail[i](k), arend_second = n_Z_tv_theta-1;
        arma::vec exp_lincomb_second = exp(Z_tv_theta.rows(arstart_second,arend_second) * B_sp_second_tmp);
        double S0_second = accu(exp_lincomb_second);
        S0_second_all.push_back(S0_second);  
        arma::mat Z_tv_exp_second =
            Z_tv_strata.rows(arstart_second,arend_second).each_col() % exp_lincomb_second;
        arma::vec S1_tv_second = sum(Z_tv_exp_second).t();
        S1_second_all.push_back(S1_tv_second);
      }
  }



  for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
      arma::mat Z_tv_theta = Z_tv_strata * theta;
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      // unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;

      for (unsigned int j = 0; j < Z_tv_strata.n_rows; ++j) {
        grad_tmp2 =  arma::zeros<arma::vec> (theta.n_elem);
        for (unsigned int k = 0; k < idx_fail[i].size(); ++k) { // note the k = i here.
          if(j<idx_fail[i](k)) continue;
          arma::vec B_sp_second_tmp = B_sp.row(k).t();
          //unsigned int arstart_second = idx_fail[i](k), arend_second = n_Z_tv_theta-1;
          //vec exp_lincomb_second = exp(Z_tv_theta.rows(arstart_second,arend_second) * B_sp_second_tmp);
          double S0_second = S0_second_all[i*idx_fail[i].size()+k];     
          arma::vec S1_tv_second = S1_second_all[i*idx_fail[i].size()+k];
          double dlambda = 1/S0_second;
          arma::vec exp_zi_beta = exp(Z_tv_theta.row(j) * B_sp_second_tmp);
          grad_tmp2 += kron(dlambda*exp_zi_beta[0]*(Z_tv_strata.row(j).t() - S1_tv_second/S0_second), B_sp_second_tmp);//*exp_zi_beta[0];
        }

        grad_part2.col(j) = grad_tmp2;
      }
  }

  for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);

      for (unsigned int j = 0; j < Z_tv_strata.n_rows; ++j) {
        grad_all.col(j) = grad_part1.col(j) - grad_part2.col(j);

        arma::vec grad_tmp = grad_all.col(j);
        info_J       += grad_tmp * grad_tmp.t();

        if(difflambda == false){
          grad_p_tmp   = grad_tmp - lambda_i*S_matrix*vectorise(theta.t(), 0)/N;
        }
        else{
          grad_p_tmp   = grad_tmp - lambda_i_mat%S_matrix*vectorise(theta.t(), 0)/N;
        }        

        info_J_p     += grad_p_tmp*grad_p_tmp.t();
        info_J_p_gic += grad_p_tmp*grad_tmp.t();
      }
  }

  return List::create(_["info_J"]=info_J,
                      _["info_J_p"]=info_J_p,
                      _["info_J_p_gic"]=info_J_p_gic);
}


List TIC_J_penalized_second_bresties(const arma::mat &Z_tv, const arma::mat &B_spline, arma::mat &theta,
                                    const unsigned int n_strata,
                                    vector<arma::uvec> &idx_B_sp, 
                                    vector<vector<arma::uvec>> &idx_fail,
                                    vector<vector<unsigned int>> &idx_Z_strata,
                                    const bool &TIC_prox,
                                    double lambda_i,
                                    arma::mat &lambda_i_mat,
                                    const bool &difflambda,
                                    arma::mat S_matrix) {

  // gradient and info matrix
  int N     = Z_tv.n_rows;    //sample size 
  int p     = theta.n_rows; //dimension
  int K     = theta.n_cols; //number of knots   

  arma::vec grad    = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_p  = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_tmp, grad_p_tmp;

  arma::vec grad_tmp1 = arma::zeros<arma::vec> (theta.n_elem);
  arma::vec grad_tmp2 = arma::zeros<arma::vec> (theta.n_elem); 
  arma::mat info_J   = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  arma::mat info_J_p = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  arma::mat info_J_p_gic = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);

  arma::mat grad_part1 = arma::zeros<arma::mat>(p*K,N);
  arma::mat grad_part2 = arma::zeros<arma::mat>(p*K,N);
  arma::mat grad_all = arma::zeros<arma::mat>(p*K,N);

  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
    arma::mat Z_tv_theta = Z_tv_strata * theta;
    arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
    unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
    for (unsigned int j = 0; j < idx_fail[i].size(); ++j) {
      arma::vec B_sp_tmp = B_sp.row(j).t();
      unsigned int arstart = idx_fail[i][j](0), arend = n_Z_tv_theta-1;
          // nar = idx_fail[i][j].n_elem;
      arma::vec exp_lincomb = exp(Z_tv_theta.rows(arstart,arend) * B_sp_tmp);
      double S0 = accu(exp_lincomb);
      arma::mat Z_tv_exp =
        Z_tv_strata.rows(arstart,arend).each_col() % exp_lincomb;
      arma::vec S1_tv    = sum(Z_tv_exp).t();
      for(unsigned int k = 0; k < idx_fail[i][j].n_elem; ++k){
        grad_tmp = kron((Z_tv_strata.row(idx_fail[i][j](k))).t()-S1_tv/S0, B_sp_tmp);
        grad_part1.col(idx_fail[i][j](k)) = grad_tmp;
      }
    }
  }

  //calculate S0 and S1:
  NumericVector S0_second_all, dlambda_all;
  vector<arma::vec>  S1_second_all;
  for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
      arma::mat Z_tv_theta = Z_tv_strata * theta;
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;

      for (unsigned int k = 0; k < idx_fail[i].size(); ++k){
        arma::vec B_sp_second_tmp = B_sp.row(k).t();
        unsigned int arstart_second = idx_fail[i][k](0), arend_second = n_Z_tv_theta-1,
            nar_second = idx_fail[i][k].n_elem;
        arma::vec exp_lincomb_second = exp(Z_tv_theta.rows(arstart_second,arend_second) * B_sp_second_tmp);
        double S0_second = accu(exp_lincomb_second);
        S0_second_all.push_back(S0_second); 
        dlambda_all.push_back(nar_second/S0_second);  
        arma::mat Z_tv_exp_second =
            Z_tv_strata.rows(arstart_second,arend_second).each_col() % exp_lincomb_second;
        arma::vec S1_tv_second = sum(Z_tv_exp_second).t();
        S1_second_all.push_back(S1_tv_second);
      }
  }


  for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
      arma::mat Z_tv_theta = Z_tv_strata * theta;
      arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
      // unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;

      for (unsigned int j = 0; j < Z_tv_strata.n_rows; ++j) {
        grad_tmp2 =  arma::zeros<arma::vec> (theta.n_elem);
        for (unsigned int k = 0; k < idx_fail[i].size(); ++k) { // note the k = i here.
          if(j<idx_fail[i][k](0)) {
            continue;
          }
          arma::vec B_sp_second_tmp = B_sp.row(k).t();
          // unsigned int arstart_second = idx_fail[i][k](0), arend_second = n_Z_tv_theta-1,
          //      nar_second = idx_fail[i][k].n_elem;
          // arma::vec exp_lincomb_second = exp(Z_tv_theta.rows(arstart_second,arend_second) * B_sp_second_tmp);
          // arma::mat Z_tv_exp_second =
          //   Z_tv_strata.rows(arstart_second,arend_second).each_col() % exp_lincomb_second;
          double S0_second = S0_second_all[i*idx_fail[i].size() + k];
          arma::vec S1_tv_second = S1_second_all[i*idx_fail[i].size() + k];
          double dlambda = dlambda_all[i*idx_fail[i].size() + k];
          arma::vec exp_zi_beta = exp(Z_tv_theta.row(j) * B_sp_second_tmp);
          grad_tmp2 += kron(dlambda*exp_zi_beta[0]*(Z_tv_strata.row(j).t() - S1_tv_second/S0_second), B_sp_second_tmp);//*exp_zi_beta[0];
        }
        grad_part2.col(j) = grad_tmp2;
      }
  }

  for (unsigned int i = 0; i < n_strata; ++i) {
      arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);

      for (unsigned int j = 0; j < Z_tv_strata.n_rows; ++j) {
        grad_all.col(j) = grad_part1.col(j) - grad_part2.col(j);

        arma::vec grad_tmp = grad_all.col(j);
        info_J       += grad_tmp * grad_tmp.t();

        if(difflambda == false){
          grad_p_tmp   = grad_tmp - lambda_i*S_matrix*vectorise(theta.t(), 0)/N;
        }
        else{
          grad_p_tmp   = grad_tmp - lambda_i_mat%S_matrix*vectorise(theta.t(), 0)/N;
        }        

        info_J_p     += grad_p_tmp*grad_p_tmp.t();
        info_J_p_gic += grad_p_tmp*grad_tmp.t();
      }
  }

  return List::create(_["info_J"]=info_J,
                      _["info_J_p"]=info_J_p,
                      _["info_J_p_gic"]=info_J_p_gic);
}


List spline_udpate(const arma::mat &Z_tv, const arma::mat &B_spline, arma::mat &theta, 
                  const arma::mat &Z_ti, const arma::vec &beta_ti, 
                  arma::mat &S_matrix,
                  double &lambda_i,
                  arma::mat &lambda_i_mat, arma::mat &lambda_S_matrix, 
                  const bool &difflambda,
                  const bool &ti,
                  const unsigned int n_strata,
                  vector<arma::uvec> &idx_B_sp, vector<arma::uvec> &idx_fail,
                  vector<unsigned int> n_Z_strata,
                  vector<vector<unsigned int>> &idx_Z_strata,
                  vector<vector<unsigned int>> &istart,
                  vector<vector<unsigned int>> &iend,
                  const std::string &method="Newton", const double &lambda=1e8, const double &factor = 1.0,
                  const bool &parallel=false, const unsigned int &threads=1,
                  const unsigned int &iter_max=20,
                  const double &tol=1e-10, 
                  const double &s=1e-2, const double &t=0.6,
                  const std::string &btr="dynamic",
                  const std::string &stop="incre",
                  const bool &TIC_prox = false,
                  const bool &fixedstep = false,
                  const bool &ICLastOnly = false) {

  int N     = Z_tv.n_rows;    //sample size 

  unsigned int iter = 0, btr_max = 1000, btr_ct = 0, iter_NRstop = 0;
  bool NRstop = false;
  double crit = 1.0, v = 1.0, logplkd_init = 0, logplkd, diff_logplkd, 
    inc, rhs_btr = 0;
  arma::vec crit123(3);
  crit123.fill(1.0);
  List objfun_list, update_list;
  NumericVector logplkd_vec;
  objfun_list = objfun_fixtra(Z_tv, B_spline, theta, Z_ti, beta_ti, ti, n_strata,
                            idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata,
                            istart, iend, parallel, threads);
  logplkd = objfun_list["logplkd"];
  logplkd -= lambda_i*as_scalar(vectorise(theta,1)*S_matrix*vectorise(theta,1).t())/(N*1.0);
  List theta_list = List::create(theta), grad_list;
    
  NumericVector AIC_all, TIC_all, TIC2_all, GIC_all;
  arma::mat VarianceMatrix;

  while (iter < iter_max) {
    ++iter;
    update_list = stepinc_fixtra_spline(Z_tv, B_spline, theta, Z_ti, beta_ti, 
                                        S_matrix, lambda_i, lambda_i_mat, difflambda,
                                        ti, n_strata,
                                        idx_B_sp, idx_fail, idx_Z_strata, istart, iend,
                                        method, lambda*pow(factor,iter-1), parallel, threads);
    v = 1.0; // reset step size
    arma::mat step = update_list["step"];
    inc = update_list["inc"];
    if (btr=="none") {
      theta += step;
      crit = inc / 2;
      objfun_list = objfun_fixtra(Z_tv, B_spline, theta, Z_ti, beta_ti, ti,
                                  n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                  idx_Z_strata, istart, iend, parallel, threads);
      logplkd = objfun_list["logplkd"];
      logplkd -= lambda_i*as_scalar(vectorise(theta,1)*S_matrix*vectorise(theta,1).t())/(N*1.0);
    } else {
      arma::mat theta_tmp = theta + step;
      objfun_list = objfun_fixtra(Z_tv, B_spline, theta_tmp, Z_ti, beta_ti, ti,
                              n_strata, idx_B_sp, idx_fail, n_Z_strata,
                              idx_Z_strata, istart, iend, parallel,
                              threads);
      double logplkd_tmp = objfun_list["logplkd"];
      logplkd_tmp -= lambda_i*as_scalar(vectorise(theta_tmp,1)*S_matrix*vectorise(theta_tmp,1).t())/(N*1.0);
      diff_logplkd = logplkd_tmp - logplkd;
      if (btr=="dynamic")      rhs_btr = inc;
      else if (btr=="static")  rhs_btr = 1.0;
      //btr_ct = 0;
      while (diff_logplkd < s * v * rhs_btr && btr_ct < btr_max) {
        ++btr_ct;
        v *= t;
        theta_tmp   = theta + v * step;
        objfun_list = objfun_fixtra(Z_tv, B_spline, theta_tmp, Z_ti, beta_ti, ti,
                                    n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                    idx_Z_strata, istart, iend, parallel,
                                    threads);
        double logplkd_tmp = objfun_list["logplkd"];
        logplkd_tmp -= lambda_i*as_scalar(vectorise(theta_tmp,1)*S_matrix*vectorise(theta_tmp,1).t())/(N*1.0);
        diff_logplkd  = logplkd_tmp - logplkd;
      }
      theta = theta_tmp;
      if (iter==1) logplkd_init = logplkd;
      if (stop=="incre")
        crit = inc / 2;
      else if (stop=="relch")
        crit = abs(diff_logplkd/(diff_logplkd+logplkd));
      else if (stop=="ratch")
        crit = abs(diff_logplkd/(diff_logplkd+logplkd-logplkd_init));
      else if (stop=="all"){
        crit123(0) = abs(diff_logplkd/(diff_logplkd+logplkd-logplkd_init));
        arma::vec grad = update_list["grad"];
        crit123(1) = abs(accu(grad));
        if(iter > 1){
          arma::mat theta_tmp2 = theta_list[iter-2];
          double diff_theta = abs(accu(theta_tmp2-theta));
          crit123(2) = diff_theta;
        }
        crit = max(crit123);
        Rcout<<"crit = "<<crit<<endl;
      }
      logplkd += diff_logplkd;
    }
    Rcout << "Iter " << iter << ": Obj fun = " << setprecision(7) << fixed << 
     logplkd << "; Stopping crit = " << setprecision(7) << scientific << 
       crit << ";" << endl;
    logplkd_vec.push_back(logplkd);
    theta_list.push_back(theta);

    grad_list.push_back(update_list["grad"]);
    //force the algorithms to run specified steps
    if(!fixedstep) {
      if(crit < tol)
        break;
    }
    //record NR's stopping step:
    if(NRstop==false){
      if(crit < tol){
        iter_NRstop = iter;
        NRstop=true;
      }
    }


  }

  objfun_list = objfun_fixtra(Z_tv, B_spline, theta, Z_ti, beta_ti, ti,
                              n_strata, idx_B_sp, idx_fail, n_Z_strata,
                              idx_Z_strata, istart, iend, parallel, threads);
  logplkd = objfun_list["logplkd"];

  if(ICLastOnly == true){
    List J_tmp;
    J_tmp = TIC_J_penalized_second(Z_tv, B_spline, theta, n_strata, idx_B_sp, idx_fail, idx_Z_strata, TIC_prox, 
                            lambda_i, lambda_i_mat, difflambda, S_matrix);
    arma::mat J    = J_tmp["info_J"];
    arma::mat J_p  = J_tmp["info_J_p"];
    arma::mat J_p_GIC = J_tmp["info_J_p_gic"];

    double df;
    arma::mat info = update_list["info"];
    arma::mat info_lambda = info + lambda_S_matrix;
    arma::mat info_lambda_inv;
    info_lambda_inv = inv(info_lambda);
    //AIC:
    df = trace(info * info_lambda_inv);
    AIC_all.push_back(-2*logplkd*N + 2*df);
    //TIC
    df = trace(info*info_lambda_inv*J*info_lambda_inv);
    TIC_all.push_back(-2*logplkd*N + 2*df);
    //TIC2:
    df = trace(info*info_lambda_inv*J_p*info_lambda_inv);
    TIC2_all.push_back(-2*logplkd*N + 2*df);
    //GIC:
    df = trace(info_lambda_inv*J_p_GIC);
    GIC_all.push_back(-2*logplkd*N + 2*df);

    VarianceMatrix = info_lambda_inv*J_p*info_lambda_inv;
  }
  else{
    arma::mat info = update_list["info"];
    VarianceMatrix = info;
  }  

  arma::mat info = update_list["info"];
  arma::vec grad = update_list["grad"];

  return List::create(_["theta"]=theta,
                      _["logplkd"]=logplkd, 
                      _["info"]=info, 
                      _["grad"]=grad,
                      _["grad_list"]=grad_list,
                      _["logplkd_vec"]=logplkd_vec,
                      _["AIC_all"]=AIC_all,
                      _["TIC_all"]=TIC_all,
                      _["TIC2_all"]=TIC2_all,
                      _["GIC_all"]=GIC_all,
                      _["theta_list"]=theta_list,
                      _["iter_NRstop"]=iter_NRstop,
                      _["VarianceMatrix"]=VarianceMatrix);
}



// [[Rcpp::export]]
List surtiver_fixtra_fit_penalizestop(const arma::vec &event, const IntegerVector &count_strata,
                             const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta_init,
                             const arma::mat &Z_ti, const arma::vec &beta_ti_init,
                             const arma::vec &lambda_spline,
                             const arma::mat &SmoothMatrix,
                             const arma::vec &effectsize,
                             const std::string &SplineType = "pspline",
                             const std::string &method="Newton",
                             const double lambda=1e8, const double &factor=1.0,
                             const bool &parallel=false, const unsigned int &threads=1,
                             const double &tol=1e-10, const unsigned int &iter_max=20,
                             const double &s=1e-2, const double &t=0.6,
                             const std::string &btr="dynamic",
                             const std::string &stop="incre",
                             const bool &TIC_prox = false,
                             const bool &fixedstep = true,
                             const bool &difflambda = false,
                             const bool &ICLastOnly = false) {

  bool ti = arma::norm(Z_ti, "inf") > sqrt(arma::datum::eps);
  IntegerVector cumsum_strata = cumsum(count_strata);
  unsigned int n_strata = cumsum_strata.length();
  cumsum_strata.push_front(0);

  vector<arma::uvec> idx_fail, idx_B_sp;
  vector<vector<unsigned int>> idx_Z_strata;
  // each element of idx_Z_strata contains start/end row indices of Z_strata
  vector<unsigned int> n_fail, n_Z_strata;
  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::uvec idx_fail_tmp =
      find(event.rows(cumsum_strata[i], cumsum_strata[i+1]-1)==1);
    n_fail.push_back(idx_fail_tmp.n_elem);
    vector<unsigned int> idx_Z_strata_tmp;
    idx_Z_strata_tmp.push_back(cumsum_strata[i]+idx_fail_tmp(0));
    idx_Z_strata_tmp.push_back(cumsum_strata[i+1]-1);
    idx_Z_strata.push_back(idx_Z_strata_tmp);
    n_Z_strata.push_back(cumsum_strata[i+1]-cumsum_strata[i]-
      idx_fail_tmp(0));
    idx_B_sp.push_back(cumsum_strata[i]+idx_fail_tmp);
    idx_fail_tmp -= idx_fail_tmp(0);
    idx_fail.push_back(idx_fail_tmp);
  }

  // istart and iend for each thread when parallel=true
  vector<arma::vec>  cumsum_ar;
  vector<vector<unsigned int>> istart, iend;
  if (parallel) {
    for (unsigned int i = 0; i < n_strata; ++i) {
      double scale_fac = as_scalar(idx_fail[i].tail(1));
      cumsum_ar.push_back(
        (double)n_Z_strata[i] / scale_fac * arma::regspace(1,n_fail[i]) -
          arma::cumsum(arma::conv_to<arma::vec> ::from(idx_fail[i])/scale_fac));
      vector<unsigned int> istart_tmp, iend_tmp;
      for (unsigned int id = 0; id < threads; ++id) {
        istart_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
          cumsum_ar[i](n_fail[i]-1)/(double)threads*id, 1)));
        iend_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
          cumsum_ar[i](n_fail[i]-1)/(double)threads*(id+1), 1)));
        if (id == threads-1) {
          iend_tmp.pop_back();
          iend_tmp.push_back(n_fail[i]);
        }
      }
      istart.push_back(istart_tmp);
      iend.push_back(iend_tmp);
    }
  }

  arma::mat theta = theta_init; arma::vec beta_ti = beta_ti_init;
  // double crit = 1.0, v = 1.0, logplkd_init = 0, logplkd, diff_logplkd, 
  //   inc, rhs_btr = 0;
  double logplkd;
  List objfun_list, update_list;
  NumericVector logplkd_vec, iter_NR_all;
  objfun_list = objfun_fixtra(Z_tv, B_spline, theta, Z_ti, beta_ti, ti, n_strata,
                              idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata,
                              istart, iend, parallel, threads);
  logplkd = objfun_list["logplkd"];
  List theta_list;
  
  //for tic:
  int N     = Z_tv.n_rows;    //sample size 
  int p     = theta_init.n_rows; //dimension
  int K     = theta_init.n_cols; //number of knots   
  
  List SplineUdpate, J_tmp;

  int n_lambda         = lambda_spline.n_elem;
  arma::cube theta_all(p, K, n_lambda);   

  arma::vec AIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec BIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC2 = arma::zeros<arma::vec> (n_lambda);
  arma::vec GIC  = arma::zeros<arma::vec> (n_lambda);

  arma::vec AIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec BIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC2_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec GIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  NumericVector AIC_all, TIC_all, TIC2_all, GIC_all;
  //mat likelihood_mat = arma::zeros<arma::mat>(iter_max, n_lambda);    for fixed step usage

  arma::mat S_matrix;
  if(SplineType == "pspline") {
    S_matrix        = spline_construct(K, p, SplineType);  
  }
  else{
    S_matrix        = spline_construct2(K, p, SplineType, SmoothMatrix);  
  }

  arma::vec lambda_effct = arma::zeros<arma::vec> (p*K);
  arma::mat lambda_i_mat;
  //List VarianceMatrix, VarianceMatrix2, VarianceMatrix_I, VarianceMatrix_J;
  arma::mat VarianceMatrix;
  arma::mat info;

  for (int i = 0; i < n_lambda; ++i)
  {
    //generate a lambda matrix if different lambda is used for different covariate:
    if(difflambda){

      for (int j = 0; j < p; ++j)
      {
        arma::vec sublambda = arma::zeros<arma::vec> (K);
        for (int jj = 0; jj < K; ++jj)
        {
          sublambda(jj) = lambda_spline(i) * sqrt(N /  effectsize(j));
        }
        lambda_effct.subvec(j*K, ((j+1)*K-1)) = sublambda;
      }      
    }
    lambda_i_mat = repmat(lambda_effct,1,p*K);
    arma::mat lambda_S_matrix;
    arma::mat theta_ilambda  = theta_init;
    double lambda_i = lambda_spline(i);
    if(difflambda == false){
      lambda_S_matrix = lambda_i*S_matrix;
    }
    else{
      lambda_S_matrix = lambda_i_mat%S_matrix;
    }

    SplineUdpate    = spline_udpate(Z_tv, B_spline, theta_ilambda,
                                    Z_ti, beta_ti, 
                                    S_matrix, lambda_i, lambda_i_mat, lambda_S_matrix, difflambda,
                                    ti, n_strata, idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata, istart, iend,
                                    method, lambda, factor,
                                    parallel, threads, iter_max, tol, s, t, btr, stop, TIC_prox, fixedstep, ICLastOnly);

    arma::mat info_tmp              = SplineUdpate["info"];
    double logplkd          = SplineUdpate["logplkd"];
    AIC_all                 = SplineUdpate["AIC_all"];
    TIC_all                 = SplineUdpate["TIC_all"];
    TIC2_all                = SplineUdpate["TIC2_all"];
    GIC_all                 = SplineUdpate["GIC_all"];
    theta_list              = SplineUdpate["theta_list"];
    int iter_tmp            = SplineUdpate["iter_NRstop"];
    arma::mat VarianceMatrix_tmp        = SplineUdpate["VarianceMatrix"];

    theta_all.slice(i)    = theta_ilambda;
    logplkd_vec.push_back(logplkd);
    iter_NR_all.push_back(iter_tmp);
    VarianceMatrix = VarianceMatrix_tmp;
    info           = info_tmp;
    //Variance matrix:
    // arma::mat VarianceMatrix_tmp = info_lambda_inv*J_p*info_lambda_inv;
    // arma::mat VarianceMatrix_tmp2 = info_lambda_inv*J*info_lambda_inv;     
    // VarianceMatrix.push_back(VarianceMatrix_tmp);
    // VarianceMatrix2.push_back(VarianceMatrix_tmp2);
    // VarianceMatrix_I.push_back(info_lambda_inv);
    // VarianceMatrix_J.push_back(inv(J_p));
    Rcout<<fixed<<"current lambda done: "<< lambda_spline(i) <<endl;

  }

  List grad_list            = SplineUdpate["grad_list"]; 

  return List::create(_["theta"]=theta,
                      _["logplkd"]=logplkd, 
                      _["theta_all"]=theta_all,
                      _["theta_list"]=theta_list,
                      _["AIC_all"]=AIC_all,
                      _["TIC_all"]=TIC_all,
                      _["TIC2_all"]=TIC2_all,
                      _["GIC_all"]=GIC_all,
                      _["logplkd_vec"]=logplkd_vec,
                      _["SplineType"]=SplineType,
                      _["iter_NR_all"]=iter_NR_all,
                      _["VarianceMatrix"]=VarianceMatrix,
                      _["info"] = info,
                      _["grad_list"]=grad_list);
}






List stepinc_fixtra_spline_bresties(const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta,
                                    const arma::mat &Z_ti, const arma::vec &beta_ti, 
                                    arma::mat &S_matrix,
                                    double &lambda_i,
                                    arma::mat &lambda_i_mat,
                                    const bool &difflambda,
                                    const bool &ti,
                                    const unsigned int n_strata,
                                    vector<arma::uvec> &idx_B_sp, 
                                    vector<vector<arma::uvec>> &idx_fail,
                                    vector<vector<unsigned int>> &idx_Z_strata,
                                    vector<vector<unsigned int>> &istart,
                                    vector<vector<unsigned int>> &iend,
                                    const std::string &method="Newton", const double &lambda=1e8,
                                    const bool &parallel=false, const unsigned int &threads=1) {
  
  // int N = Z_tv.n_rows;

  arma::vec grad, grad_p; arma::mat info, info_p; // gradient and info matrix
  if (ti) {
    grad = arma::zeros<arma::vec> (theta.n_elem+beta_ti.n_elem);
    info = arma::zeros<arma::mat>(theta.n_elem+beta_ti.n_elem, theta.n_elem+beta_ti.n_elem);
  } else {
    grad = arma::zeros<arma::vec> (theta.n_elem);
    info = arma::zeros<arma::mat>(theta.n_elem, theta.n_elem);
  }
  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::mat Z_tv_strata = Z_tv.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
    arma::mat Z_ti_strata; arma::vec Z_ti_beta_ti;
    if (ti) {
      Z_ti_strata = Z_ti.rows(idx_Z_strata[i][0], idx_Z_strata[i][1]);
      Z_ti_beta_ti = Z_ti_strata * beta_ti;
    }
    arma::mat Z_tv_theta = Z_tv_strata * theta;
    arma::mat B_sp = B_spline.rows(idx_B_sp[i]);
    unsigned int n_Z_tv_theta = Z_tv_theta.n_rows;
    if (parallel) {
      omp_set_num_threads(threads);
      #pragma omp parallel
      {
        unsigned int id, j;
        id = omp_get_thread_num();
        arma::vec grad_tmp(size(grad), arma::fill::zeros);
        arma::mat info_tmp(size(info), arma::fill::zeros);
        for (j = istart[i][id]; j < iend[i][id]; ++j) {
          arma::vec B_sp_tmp = B_sp.row(j).t();
          unsigned int arstart = idx_fail[i][j](0), arend = n_Z_tv_theta-1, 
            nar = idx_fail[i][j].n_elem;
          arma::vec exp_lincomb = exp(Z_tv_theta.rows(arstart,arend) * B_sp_tmp);
          double S0 = accu(exp_lincomb);
          arma::mat Z_tv_exp =
            Z_tv_strata.rows(arstart,arend).each_col() % exp_lincomb;
          arma::vec S1_tv = sum(Z_tv_exp).t();
          arma::mat S2 = Z_tv_strata.rows(arstart,arend).t() * Z_tv_exp;
          grad_tmp += kron(sum(Z_tv_strata.rows(idx_fail[i][j])).t()-nar*S1_tv/S0, B_sp_tmp);
          info_tmp += nar*
            kron(S2/S0-S1_tv*S1_tv.t()/pow(S0,2),B_sp_tmp*B_sp_tmp.t());
        }
        #pragma omp critical (gradient)
          grad += grad_tmp;
        #pragma omp critical (information)
          info += info_tmp;
      }
    } else {
      for (unsigned int j = 0; j < idx_fail[i].size(); ++j) {
        arma::vec B_sp_tmp = B_sp.row(j).t();
        unsigned int arstart = idx_fail[i][j](0), arend = n_Z_tv_theta-1,
            nar = idx_fail[i][j].n_elem;
        arma::vec exp_lincomb = exp(Z_tv_theta.rows(arstart,arend) * B_sp_tmp);
        double S0 = accu(exp_lincomb);
        arma::mat Z_tv_exp =
          Z_tv_strata.rows(arstart,arend).each_col() % exp_lincomb;
        arma::vec S1_tv = sum(Z_tv_exp).t();
        arma::mat S2 = Z_tv_strata.rows(arstart,arend).t() * Z_tv_exp;
        grad += kron(sum(Z_tv_strata.rows(idx_fail[i][j])).t()-nar*S1_tv/S0, B_sp_tmp);
        info += nar*
          (kron(S2/S0-S1_tv*S1_tv.t()/pow(S0,2),B_sp_tmp*B_sp_tmp.t()));
      }
    }
  }

  arma::vec step; // Newton step
  // info_p  = info + N*lambda_i*S_matrix;
  // grad_p  = grad - N*lambda_i*S_matrix*vectorise(theta, 1);

  if(difflambda == false){
    info_p  = info + lambda_i*S_matrix;
    grad_p  = grad - lambda_i*S_matrix*vectorise(theta.t(), 0);
  }
  else{
    info_p  = info + lambda_i_mat%S_matrix;
    grad_p  = grad - (lambda_i_mat%S_matrix)*vectorise(theta.t(), 0);
  }

  step = solve(info_p, grad_p, arma::solve_opts::fast+arma::solve_opts::likely_sympd);
  
  double inc = dot(grad_p, step) / Z_tv.n_rows; //increment, notice how to calculate in later part
  if (ti) {
    return List::create(_["step_tv"]=reshape(step.head(theta.n_elem),
                                      size(theta.t())).t(),
                        _["step_ti"]=step.tail(beta_ti.n_elem),
                        _["grad"]=grad, _["info"]=info, _["inc"]=inc);
  } else {
    return List::create(_["step"]=reshape(step, size(theta.t())).t(),
                        _["grad"]=grad, _["info"]=info, _["inc"]=inc);
  }
}


List spline_udpate_bresties(const arma::mat &Z_tv, const arma::vec &time, 
                  const arma::mat &B_spline, 
                  arma::mat &theta,
                  const arma::mat &Z_ti, const arma::vec &beta_ti, 
                  arma::mat &S_matrix,
                  double &lambda_i,
                  arma::mat &lambda_i_mat, arma::mat &lambda_S_matrix, 
                  const bool &difflambda,
                  const bool &ti,
                  const unsigned int n_strata,
                  vector<arma::uvec> &idx_B_sp, 
                  vector<vector<arma::uvec>> &idx_fail,
                  vector<unsigned int> n_Z_strata,
                  vector<vector<unsigned int>> &idx_Z_strata,
                  vector<vector<unsigned int>> &istart,
                  vector<vector<unsigned int>> &iend,
                  const std::string &method="Newton", const double &lambda=1e8, const double &factor = 1.0,
                  const bool &parallel=false, const unsigned int &threads=1,
                  const unsigned int &iter_max=20,
                  const double &tol=1e-10, 
                  const double &s=1e-2, const double &t=0.6,
                  const std::string &btr="dynamic",
                  const std::string &stop="incre",
                  const bool &TIC_prox = false,
                  const bool &fixedstep = true,
                  const bool &ICLastOnly = false) {
    
    int N     = Z_tv.n_rows;    //sample size 

    unsigned int iter = 0, btr_max = 100, btr_ct = 0;
    double crit = 1.0, v = 1.0, logplkd_init = 0, logplkd, diff_logplkd, 
      inc, rhs_btr = 0;
    arma::vec crit123(3);
    crit123.fill(1.0);
    List objfun_list, update_list;
    NumericVector logplkd_vec;
    objfun_list = obj_fixtra_bresties(Z_tv, B_spline, theta, Z_ti, beta_ti, 
                                      ti, n_strata, idx_B_sp, idx_fail, 
                                      n_Z_strata, idx_Z_strata,
                                      istart, iend, parallel, threads);
    logplkd = objfun_list["logplkd"];
    logplkd -= lambda_i*as_scalar(vectorise(theta,1)*S_matrix*vectorise(theta,1).t())/(N*1.0);

    List theta_list = List::create(theta);
    
    NumericVector AIC_all, TIC_all, TIC2_all, GIC_all;
    NumericVector AIC_trace, TIC_trace, TIC2_trace, GIC_trace;
    arma::mat VarianceMatrix;

    //while (iter < iter_max) {    
    while (iter < iter_max) {
      ++iter;
      update_list = stepinc_fixtra_spline_bresties( Z_tv, B_spline, theta, Z_ti, beta_ti, 
                                                    S_matrix, lambda_i, lambda_i_mat, difflambda,
                                                    ti, n_strata,
                                                    idx_B_sp, idx_fail, idx_Z_strata, istart, iend,
                                                    method, lambda*pow(factor,iter-1), parallel, threads);
      v = 1.0; // reset step size
      arma::mat step = update_list["step"];
      inc = update_list["inc"];
      if (btr=="none") {
        theta += step;
        crit = inc / 2;
        objfun_list = obj_fixtra_bresties(Z_tv, B_spline, theta, Z_ti, beta_ti, ti,
                                          n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                          idx_Z_strata, istart, iend, parallel, threads);
        logplkd = objfun_list["logplkd"];
        logplkd -= lambda_i*as_scalar(vectorise(theta,1)*S_matrix*vectorise(theta,1).t())/(N*1.0);
      } else {
        arma::mat theta_tmp = theta + step;
        objfun_list = obj_fixtra_bresties(Z_tv, B_spline, theta_tmp, Z_ti, beta_ti, ti,
                                          n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                          idx_Z_strata, istart, iend, parallel, threads);
        double logplkd_tmp = objfun_list["logplkd"];
        diff_logplkd = logplkd_tmp - logplkd;
        if (btr=="dynamic")      rhs_btr = inc;
        else if (btr=="static")  rhs_btr = 1.0;
        while (diff_logplkd < s * v * rhs_btr && btr_ct < btr_max) {
          ++btr_ct;
          v *= t;
          theta_tmp = theta + v * step;
          objfun_list = obj_fixtra_bresties(Z_tv, B_spline, theta_tmp, Z_ti, beta_ti, ti,
                                            n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                            idx_Z_strata, istart, iend, parallel, threads);
          double logplkd_tmp = objfun_list["logplkd"];
          logplkd_tmp -= lambda_i*as_scalar(vectorise(theta_tmp,1)*S_matrix*vectorise(theta_tmp,1).t())/(N*1.0);
          diff_logplkd  = logplkd_tmp - logplkd;
        }
        theta = theta_tmp;
        if (iter==1) logplkd_init = logplkd;
        if (stop=="incre")
          crit = inc / 2;
        else if (stop=="relch")
          crit = abs(diff_logplkd/(diff_logplkd+logplkd));
        else if (stop=="ratch")
          crit = abs(diff_logplkd/(diff_logplkd+logplkd-logplkd_init));
        else if (stop=="all"){
          crit123(0) = abs(diff_logplkd/(diff_logplkd+logplkd-logplkd_init));
          arma::vec grad = update_list["grad"];
          crit123(1) = abs(accu(grad));
          if(iter > 1){
            arma::mat theta_tmp2 = theta_list[iter-2];
            double diff_theta = abs(accu(theta_tmp2-theta));
            crit123(2) = diff_theta;
          }
          crit = max(crit123);
          Rcout<<"crit = "<<crit<<endl;
        }        
        logplkd += diff_logplkd;
      }
      Rcout << "Iter " << iter << ": Obj fun = " << setprecision(7) << fixed << 
       logplkd << "; Stopping crit = " << setprecision(7) << scientific << 
         crit << ";" << endl;
      logplkd_vec.push_back(logplkd);
      theta_list.push_back(theta);

      if(!fixedstep) {
        if(crit < tol)
          break;
      }

      // if(penalizestop == true){

      //   List J_tmp;
      //   J_tmp = TIC_J_penalized_second_bresties(Z_tv, B_spline, theta, n_strata, idx_B_sp, idx_fail, idx_Z_strata, TIC_prox, 
      //                           lambda_i, lambda_i_mat, difflambda, S_matrix);
      //   arma::mat J    = J_tmp["info_J"];
      //   arma::mat J_p  = J_tmp["info_J_p"];
      //   arma::mat J_p_GIC = J_tmp["info_J_p_gic"];

      //   double df;
      //   arma::mat info = update_list["info"];
      //   arma::mat info_lambda = info + lambda_S_matrix;
      //   arma::mat info_lambda_inv;
      //   info_lambda_inv = inv(info_lambda);

      //   //AIC:
      //   df = trace(info * info_lambda_inv);
      //   AIC_all.push_back(-2*logplkd*N + 2*df);
      //   //TIC
      //   df = trace(info*info_lambda_inv*J*info_lambda_inv);
      //   TIC_all.push_back(-2*logplkd*N + 2*df);
      //   //TIC2:
      //   df = trace(info*info_lambda_inv*J_p*info_lambda_inv);
      //   TIC2_all.push_back(-2*logplkd*N + 2*df);
      //   //GIC:
      //   df = trace(info_lambda_inv*J_p_GIC);
      //   GIC_all.push_back(-2*logplkd*N + 2*df);

      //   VarianceMatrix = info_lambda_inv*J_p*info_lambda_inv;
      // }

    }

    objfun_list = obj_fixtra_bresties(Z_tv, B_spline, theta, Z_ti, beta_ti, ti,
                                          n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                          idx_Z_strata, istart, iend, parallel, threads);
    logplkd = objfun_list["logplkd"];

    if(ICLastOnly == true){
      List J_tmp;
      J_tmp = TIC_J_penalized_second_bresties(Z_tv, B_spline, theta, n_strata, idx_B_sp, idx_fail, idx_Z_strata, TIC_prox, 
                              lambda_i, lambda_i_mat, difflambda, S_matrix);
      arma::mat J    = J_tmp["info_J"];
      arma::mat J_p  = J_tmp["info_J_p"];
      arma::mat J_p_GIC = J_tmp["info_J_p_gic"];

      double df;
      arma::mat info = update_list["info"];
      arma::mat info_lambda = info + lambda_S_matrix;
      arma::mat info_lambda_inv;
      if(TIC_prox == false){
        info_lambda_inv = inv(info_lambda);
      }
      else{
        info_lambda.diag() += N / 1e9;
        info_lambda_inv = inv(info_lambda);
        info_lambda.diag() -= N / 1e9;
      }
      //AIC:
      df = trace(info * info_lambda_inv);
      AIC_all.push_back(-2*logplkd*N + 2*df);
      AIC_trace.push_back(2*df);
      //TIC
      df = trace(info*info_lambda_inv*J*info_lambda_inv);
      TIC_all.push_back(-2*logplkd*N + 2*df);
      TIC_trace.push_back(2*df);
      //TIC2:
      df = trace(info*info_lambda_inv*J_p*info_lambda_inv);
      TIC2_all.push_back(-2*logplkd*N + 2*df);
      TIC2_trace.push_back(2*df);
      //GIC:
      df = trace(info_lambda_inv*J_p_GIC);
      GIC_all.push_back(-2*logplkd*N + 2*df);
      GIC_trace.push_back(2*df);

      VarianceMatrix = info_lambda_inv*J_p*info_lambda_inv;
    }
    else{
      arma::mat info = update_list["info"];
      VarianceMatrix = info;
    }    

    arma::mat info = update_list["info"];
    arma::vec grad = update_list["grad"];

    return List::create(_["theta"]=theta,
                        _["logplkd"]=logplkd, 
                        _["info"]=info, 
                        _["grad"]=grad,
                        _["logplkd_vec"]=logplkd_vec,
                        _["AIC_all"]=AIC_all,
                        _["TIC_all"]=TIC_all,
                        _["TIC2_all"]=TIC2_all,
                        _["GIC_all"]=GIC_all,
                        _["AIC_trace"]=AIC_trace,
                        _["TIC_trace"]=TIC_trace,
                        _["TIC2_trace"]=TIC2_trace,
                        _["GIC_trace"]=GIC_trace,
                        _["theta_list"]=theta_list,
                        _["VarianceMatrix"]=VarianceMatrix);
}


// [[Rcpp::export]]
List surtiver_fixtra_fit_penalizestop_bresties(const arma::vec &event, const arma::vec &time, 
                                                const IntegerVector &count_strata,
                                                const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta_init,
                                                const arma::mat &Z_ti, const arma::vec &beta_ti_init,
                                                const arma::vec &lambda_spline,
                                                const arma::mat &SmoothMatrix,
                                                const arma::vec &effectsize,
                                                const std::string &SplineType = "pspline",
                                                const std::string &method="Newton",
                                                const double &lambda=1e8, const double &factor=1.0,
                                                const bool &parallel=false, const unsigned int &threads=1,
                                                const double &tol=1e-10, const unsigned int &iter_max=20,
                                                const double &s=1e-2, const double &t=0.6,
                                                const std::string &btr="dynamic",
                                                const std::string &stop="incre",
                                                const bool &TIC_prox = false,
                                                const bool &fixedstep = true,
                                                const bool &difflambda = false,
                                                const bool &ICLastOnly = false) {

  bool ti = arma::norm(Z_ti, "inf") > sqrt(arma::datum::eps);
  IntegerVector cumsum_strata = cumsum(count_strata);
  unsigned int n_strata = cumsum_strata.length();
  cumsum_strata.push_front(0);

  vector<vector<arma::uvec>> idx_fail;
  vector<vector<unsigned int>> idx_fail_1st, idx_Z_strata;
  vector<unsigned int> n_fail_time, n_Z_strata;
  // each element of idx_Z_strata contains start/end row indices of Z_strata
  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::vec event_tmp = event.rows(cumsum_strata[i], cumsum_strata[i+1]-1),
      time_tmp = time.rows(cumsum_strata[i], cumsum_strata[i+1]-1);
    arma::uvec idx_fail_tmp = find(event_tmp==1);
    vector<unsigned int> idx_Z_strata_tmp;
    idx_Z_strata_tmp.push_back(cumsum_strata[i]+idx_fail_tmp(0));
    idx_Z_strata_tmp.push_back(cumsum_strata[i+1]-1);
    idx_Z_strata.push_back(idx_Z_strata_tmp);
    n_Z_strata.push_back(cumsum_strata[i+1]-cumsum_strata[i]-
      idx_fail_tmp(0));
    arma::vec time_fail_tmp = time_tmp.elem(idx_fail_tmp);
    idx_fail_tmp -= idx_fail_tmp(0);
    arma::vec uniq_t = unique(time_fail_tmp);
    n_fail_time.push_back(uniq_t.n_elem);
    vector<arma::uvec> idx_fail_tmp_tmp;
    vector<unsigned int> idx_fail_1st_tmp;
    for (arma::vec::iterator j = uniq_t.begin(); j < uniq_t.end(); ++j) {
      arma::uvec tmp = idx_fail_tmp.elem(find(time_fail_tmp==*j));
      idx_fail_tmp_tmp.push_back(tmp);
      idx_fail_1st_tmp.push_back(tmp[0]);
    }
    idx_fail.push_back(idx_fail_tmp_tmp);
    idx_fail_1st.push_back(idx_fail_1st_tmp);
  }
  IntegerVector n_failtime = wrap(n_fail_time);
  IntegerVector cumsum_failtime = cumsum(n_failtime);
  cumsum_failtime.push_front(0);
  vector<arma::uvec> idx_B_sp;
  for (unsigned int i = 0; i < n_strata; ++i) {
    idx_B_sp.push_back(arma::regspace<arma::uvec>(cumsum_failtime[i],
                                      cumsum_failtime[i+1]-1));
  }
  // istart and iend for each thread when parallel=true
  vector<arma::vec>  cumsum_ar;
  vector<vector<unsigned int>> istart, iend;
  if (parallel) {
    for (unsigned int i = 0; i < n_strata; ++i) {
      double scale_fac = as_scalar(idx_fail[i].back().tail(1));
      cumsum_ar.push_back(
        (double)n_Z_strata[i] / scale_fac * arma::regspace(1, idx_fail[i].size()) -
          arma::cumsum(arma::conv_to<arma::vec> ::from(idx_fail_1st[i])/scale_fac));
      vector<unsigned int> istart_tmp, iend_tmp;
      for (unsigned int id = 0; id < threads; ++id) {
        istart_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
          cumsum_ar[i](idx_fail[i].size()-1)/(double)threads*id, 1)));
        iend_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
          cumsum_ar[i](idx_fail[i].size()-1)/(double)threads*(id+1), 1)));
        if (id == threads-1) {
          iend_tmp.pop_back();
          iend_tmp.push_back(idx_fail[i].size());
        }
      }
      istart.push_back(istart_tmp);
      iend.push_back(iend_tmp);
    }
  }


  // iterative algorithm with backtracking line search
  //unsigned int iter = 0, btr_max = 1000, btr_ct = 0;      //set as 100 instead of 1000 to improve speed
  arma::mat theta = theta_init; arma::vec beta_ti = beta_ti_init;
  List theta_list = List::create(theta), beta_ti_list = List::create(beta_ti);

  double logplkd;
  List objfun_list, update_list;
  NumericVector logplkd_vec;
  objfun_list = obj_fixtra_bresties(Z_tv, B_spline, theta, Z_ti, beta_ti, 
                                    ti, n_strata, idx_B_sp, idx_fail, 
                                    n_Z_strata, idx_Z_strata,
                                    istart, iend, parallel, threads);
  logplkd = objfun_list["logplkd"];
  
  //for tic:
  int N     = Z_tv.n_rows;    //sample size 
  int p     = theta_init.n_rows; //dimension
  int K     = theta_init.n_cols; //number of knots   
  
  List SplineUdpate, J_tmp;

  int n_lambda         = lambda_spline.n_elem;
  arma::cube theta_all(p, K, n_lambda);   

  arma::vec AIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec BIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC2 = arma::zeros<arma::vec> (n_lambda);
  arma::vec GIC  = arma::zeros<arma::vec> (n_lambda);

  arma::vec AIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec BIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC2_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec GIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  NumericVector AIC_all, TIC_all, TIC2_all, GIC_all;
  NumericVector AIC_trace, TIC_trace, TIC2_trace, GIC_trace;
  //mat likelihood_mat = arma::zeros<arma::mat>(iter_max, n_lambda);    for fixed step usage

  arma::mat S_matrix;
  if(SplineType == "pspline") {
    S_matrix        = spline_construct(K, p, SplineType);  
  }
  else{
    S_matrix        = spline_construct2(K, p, SplineType, SmoothMatrix);  
  }

  arma::vec lambda_effct = arma::zeros<arma::vec> (p*K);
  arma::mat lambda_i_mat;
  // List VarianceMatrix, VarianceMatrix2, VarianceMatrix_I, VarianceMatrix_J;
  arma::mat VarianceMatrix;
  arma::mat info;

  for (int i = 0; i < n_lambda; ++i)
  {
    //generate a lambda matrix if different lambda is used for different covariate:
    if(difflambda){
      for (int j = 0; j < p; ++j)
      {
        arma::vec sublambda = arma::zeros<arma::vec> (K);
        for (int jj = 0; jj < K; ++jj)
        {
          sublambda(jj) = lambda_spline(i) * sqrt(N /  effectsize(j));
        }
        lambda_effct.subvec(j*K, ((j+1)*K-1)) = sublambda;
      }      
    }
    lambda_i_mat = repmat(lambda_effct,1,p*K);
    arma::mat lambda_S_matrix;
    arma::mat theta_ilambda  = theta_init;
    double lambda_i = lambda_spline(i);
    if(difflambda == false){
      lambda_S_matrix = lambda_i*S_matrix;
    }
    else{
      lambda_S_matrix = lambda_i_mat%S_matrix;
    }

    SplineUdpate    = spline_udpate_bresties(Z_tv, time, B_spline, theta_ilambda,
                                            Z_ti, beta_ti, 
                                            S_matrix, lambda_i, lambda_i_mat, lambda_S_matrix, difflambda,
                                            ti, n_strata, idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata, istart, iend,
                                            method, lambda, factor,
                                            parallel, threads, iter_max, tol, s, t, btr, stop, TIC_prox, fixedstep, ICLastOnly);

    arma::mat info_tmp              = SplineUdpate["info"];
    arma::vec grad              = SplineUdpate["grad"];
    double logplkd        = SplineUdpate["logplkd"];
    AIC_all               = SplineUdpate["AIC_all"];
    TIC_all               = SplineUdpate["TIC_all"];
    TIC2_all               = SplineUdpate["TIC2_all"];
    GIC_all               = SplineUdpate["GIC_all"];
    theta_list            = SplineUdpate["theta_list"];
    AIC_trace             = SplineUdpate["AIC_trace"];
    TIC_trace             = SplineUdpate["TIC_trace"];
    TIC2_trace             = SplineUdpate["TIC2_trace"];
    GIC_trace             = SplineUdpate["GIC_trace"];
    arma::mat VarianceMatrix_tmp        = SplineUdpate["VarianceMatrix"];

    theta_all.slice(i)    = theta_ilambda;
    logplkd_vec.push_back(logplkd);
    VarianceMatrix = VarianceMatrix_tmp;
    info = info_tmp;
    // Rcout<<fixed<<"current lambda done: "<< lambda_spline(i) <<endl;


  }

  return List::create(_["theta"]=theta,
                      _["logplkd"]=logplkd, 
                      _["theta_all"]=theta_all,
                      _["theta_list"]=theta_list,
                      _["AIC_all"]=AIC_all,
                      _["TIC_all"]=TIC_all,
                      _["TIC2_all"]=TIC2_all,
                      _["GIC_all"]=GIC_all,
                      _["AIC_trace"]=AIC_trace,
                      _["TIC_trace"]=TIC_trace,
                      _["TIC2_trace"]=TIC2_trace,
                      _["GIC_trace"]=GIC_trace,
                      _["logplkd_vec"]=logplkd_vec,
                      _["SplineType"]=SplineType,
                      _["VarianceMatrix"]=VarianceMatrix,
                      _["info"]=info);
}




// [[Rcpp::export]]
List LogPartialTest(const arma::vec &event, const arma::mat &Z_tv, const arma::mat &B_spline,  
                    const arma::vec &event_test, const arma::mat &Z_tv_test, const arma::mat &B_spline_test,
                    const List &theta_list,
                    const bool &parallel=false, const unsigned int &threads=1,
                    const bool TestAll = true){

    // int p     = Z_tv.n_cols; //dimension
    // int K     = B_spline.n_cols; //number of knots  

    //for training data: //////////////////////////////////////////////////////////////
    int N    = Z_tv.n_rows;
    IntegerVector count_strata = IntegerVector::create(N);
    IntegerVector cumsum_strata = cumsum(count_strata);
    unsigned int n_strata = cumsum_strata.length();
    cumsum_strata.push_front(0);

    vector<arma::uvec> idx_fail, idx_B_sp;
    vector<vector<unsigned int>> idx_Z_strata;
    // each element of idx_Z_strata contains start/end row indices of Z_strata
    vector<unsigned int> n_fail, n_Z_strata;
    for (unsigned int i = 0; i < n_strata; ++i) {
      arma::uvec idx_fail_tmp =
        find(event.rows(cumsum_strata[i], cumsum_strata[i+1]-1)==1);
      n_fail.push_back(idx_fail_tmp.n_elem);
      vector<unsigned int> idx_Z_strata_tmp;
      idx_Z_strata_tmp.push_back(cumsum_strata[i]+idx_fail_tmp(0));
      idx_Z_strata_tmp.push_back(cumsum_strata[i+1]-1);
      idx_Z_strata.push_back(idx_Z_strata_tmp);
      n_Z_strata.push_back(cumsum_strata[i+1]-cumsum_strata[i]-
        idx_fail_tmp(0));
      idx_B_sp.push_back(cumsum_strata[i]+idx_fail_tmp);
      idx_fail_tmp -= idx_fail_tmp(0);
      idx_fail.push_back(idx_fail_tmp);
    }

    // istart and iend for each thread when parallel=true
    vector<arma::vec>  cumsum_ar;
    vector<vector<unsigned int>> istart, iend;
    if (parallel) {
      for (unsigned int i = 0; i < n_strata; ++i) {
        double scale_fac = as_scalar(idx_fail[i].tail(1));
        cumsum_ar.push_back(
          (double)n_Z_strata[i] / scale_fac * arma::regspace(1,n_fail[i]) -
            arma::cumsum(arma::conv_to<arma::vec> ::from(idx_fail[i])/scale_fac));
        vector<unsigned int> istart_tmp, iend_tmp;
        for (unsigned int id = 0; id < threads; ++id) {
          istart_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
            cumsum_ar[i](n_fail[i]-1)/(double)threads*id, 1)));
          iend_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
            cumsum_ar[i](n_fail[i]-1)/(double)threads*(id+1), 1)));
          if (id == threads-1) {
            iend_tmp.pop_back();
            iend_tmp.push_back(n_fail[i]);
          }
        }
        istart.push_back(istart_tmp);
        iend.push_back(iend_tmp);
      }
    }


    //for testing data:    ///////////////////////////////////////////////////////////
    int N2    = Z_tv_test.n_rows;
    IntegerVector cumsum_strata_test = IntegerVector::create(N2);
    unsigned int n_strata_test = cumsum_strata_test.length();
    cumsum_strata_test.push_front(0);

    vector<arma::uvec> idx_fail_test, idx_B_sp_test;
    vector<vector<unsigned int>> idx_Z_strata_test;
    vector<unsigned int> n_fail_test, n_Z_strata_test;
    for (unsigned int i = 0; i < n_strata_test; ++i) {
      arma::uvec idx_fail_tmp =
        find(event_test.rows(cumsum_strata_test[i], cumsum_strata_test[i+1]-1)==1);
      n_fail_test.push_back(idx_fail_tmp.n_elem);
      vector<unsigned int> idx_Z_strata_tmp;
      idx_Z_strata_tmp.push_back(cumsum_strata_test[i]+idx_fail_tmp(0));
      idx_Z_strata_tmp.push_back(cumsum_strata_test[i+1]-1);
      idx_Z_strata_test.push_back(idx_Z_strata_tmp);
      n_Z_strata_test.push_back(cumsum_strata_test[i+1]-cumsum_strata_test[i]-
        idx_fail_tmp(0));
      idx_B_sp_test.push_back(cumsum_strata_test[i]+idx_fail_tmp);
      idx_fail_tmp -= idx_fail_tmp(0);
      idx_fail_test.push_back(idx_fail_tmp);
    }

    // istart and iend for each thread when parallel=true
    vector<arma::vec>  cumsum_ar_test;
    vector<vector<unsigned int>> istart_test, iend_test;
    if (parallel) {
      for (unsigned int i = 0; i < n_strata_test; ++i) {
        double scale_fac = as_scalar(idx_fail_test[i].tail(1));
        cumsum_ar_test.push_back(
          (double)n_Z_strata_test[i] / scale_fac * arma::regspace(1,n_fail_test[i]) -
            arma::cumsum(arma::conv_to<arma::vec> ::from(idx_fail_test[i])/scale_fac));
        vector<unsigned int> istart_tmp, iend_tmp;
        for (unsigned int id = 0; id < threads; ++id) {
          istart_tmp.push_back(as_scalar(find(cumsum_ar_test[i] >=
            cumsum_ar_test[i](n_fail_test[i]-1)/(double)threads*id, 1)));
          iend_tmp.push_back(as_scalar(find(cumsum_ar_test[i] >=
            cumsum_ar_test[i](n_fail_test[i]-1)/(double)threads*(id+1), 1)));
          if (id == threads-1) {
            iend_tmp.pop_back();
            iend_tmp.push_back(n_fail_test[i]);
          }
        }
        istart_test.push_back(istart_tmp);
        iend_test.push_back(iend_tmp);
      }
    }

    List objfun_list, objfun_list_test, update_list;

    int theta_length = theta_list.length();
    NumericVector logplkd_vec, logplkd_vec_test;
    double logplkd, logplkd_test;

    arma::mat Z_ti = arma::zeros<arma::mat>(1,1);
    arma::mat Z_ti_test = arma::zeros<arma::mat>(1,1);
    arma::vec beta_ti = arma::zeros<arma::vec> (1);
    bool ti = false;

    for (int i = 0; i < theta_length; ++i){
      arma::mat theta_tmp = theta_list[i];
      objfun_list         = objfun_fixtra(Z_tv, B_spline, theta_tmp, Z_ti, beta_ti, ti, n_strata,
                                          idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata, 
                                          istart, iend, parallel, threads);
      objfun_list_test    = objfun_fixtra(Z_tv_test, B_spline_test, theta_tmp, Z_ti, beta_ti, ti, n_strata_test, 
                                          idx_B_sp_test, idx_fail_test, n_Z_strata_test, idx_Z_strata_test,
                                          istart_test, iend_test, parallel, threads);

      logplkd = objfun_list["logplkd"];
      logplkd_test = objfun_list_test["logplkd"];

      logplkd_vec.push_back(logplkd);
      logplkd_vec_test.push_back(logplkd_test);

      if(TestAll == false){
        if(i>0){
          if(logplkd_test < logplkd_vec_test[i-1]){
            break;
          }
        }
      }
    }

    List result;

    result["likelihood_all"] = logplkd_vec;
    result["likelihood_all_test"] = logplkd_vec_test;

    return result;
}



//Variance Matrix Calculation
// [[Rcpp::export]]
List VarianceMatrixCalculate(const arma::vec &event, const arma::mat &Z_tv, const arma::mat &B_spline,
                            arma::mat &theta,
                            double lambda_i,
                            const arma::mat &SmoothMatrix,
                            const std::string &SplineType = "smooth-spline",
                            const std::string &method="Newton",
                            const double &lambda=1e8, const double &factor = 1.0,
                            const bool &parallel=false, const unsigned int &threads=1){


    int p     = Z_tv.n_cols; //dimension
    int K     = B_spline.n_cols; //number of knots  
    int N     = Z_tv.n_rows;
    arma::mat Z_ti = arma::zeros<arma::mat>(N,p);
    arma::vec beta_ti = arma::zeros<arma::vec> (p);
    bool difflambda = false, ti = false, TIC_prox = false;

    IntegerVector count_strata = IntegerVector::create(N);
    IntegerVector cumsum_strata = cumsum(count_strata);
    unsigned int n_strata = cumsum_strata.length();
    cumsum_strata.push_front(0);

    vector<arma::uvec> idx_fail, idx_B_sp;
    vector<vector<unsigned int>> idx_Z_strata;
    // each element of idx_Z_strata contains start/end row indices of Z_strata
    vector<unsigned int> n_fail, n_Z_strata;
    for (unsigned int i = 0; i < n_strata; ++i) {
      arma::uvec idx_fail_tmp =
        find(event.rows(cumsum_strata[i], cumsum_strata[i+1]-1)==1);
      n_fail.push_back(idx_fail_tmp.n_elem);
      vector<unsigned int> idx_Z_strata_tmp;
      idx_Z_strata_tmp.push_back(cumsum_strata[i]+idx_fail_tmp(0));
      idx_Z_strata_tmp.push_back(cumsum_strata[i+1]-1);
      idx_Z_strata.push_back(idx_Z_strata_tmp);
      n_Z_strata.push_back(cumsum_strata[i+1]-cumsum_strata[i]-
        idx_fail_tmp(0));
      idx_B_sp.push_back(cumsum_strata[i]+idx_fail_tmp);
      idx_fail_tmp -= idx_fail_tmp(0);
      idx_fail.push_back(idx_fail_tmp);
    }

    // istart and iend for each thread when parallel=true
    vector<arma::vec>  cumsum_ar;
    vector<vector<unsigned int>> istart, iend;
    if (parallel) {
      for (unsigned int i = 0; i < n_strata; ++i) {
        double scale_fac = as_scalar(idx_fail[i].tail(1));
        cumsum_ar.push_back(
          (double)n_Z_strata[i] / scale_fac * arma::regspace(1,n_fail[i]) -
            arma::cumsum(arma::conv_to<arma::vec> ::from(idx_fail[i])/scale_fac));
        vector<unsigned int> istart_tmp, iend_tmp;
        for (unsigned int id = 0; id < threads; ++id) {
          istart_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
            cumsum_ar[i](n_fail[i]-1)/(double)threads*id, 1)));
          iend_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
            cumsum_ar[i](n_fail[i]-1)/(double)threads*(id+1), 1)));
          if (id == threads-1) {
            iend_tmp.pop_back();
            iend_tmp.push_back(n_fail[i]);
          }
        }
        istart.push_back(istart_tmp);
        iend.push_back(iend_tmp);
      }
    }

    arma::vec lambda_effct = arma::zeros<arma::vec> (p*K);
    arma::mat lambda_i_mat = repmat(lambda_effct,1,p*K);

    List SplineUdpate;
    arma::mat S_matrix;
    if(SplineType == "pspline") {
      S_matrix        = spline_construct(K, p, SplineType);  
    }
    else{
      S_matrix        = spline_construct2(K, p, SplineType, SmoothMatrix);  
    }    

    SplineUdpate    = stepinc_fixtra_spline(Z_tv, B_spline, theta, Z_ti, beta_ti, 
                                            S_matrix, lambda_i, lambda_i_mat, difflambda,
                                            ti, n_strata,
                                            idx_B_sp, idx_fail, idx_Z_strata, istart, iend,
                                            method, lambda, parallel, threads);

    arma::mat info = SplineUdpate["info"];

    List J_tmp;
    J_tmp = TIC_J_penalized_second(Z_tv, B_spline, theta, n_strata, idx_B_sp, idx_fail, idx_Z_strata, TIC_prox, 
                                  lambda_i, lambda_i_mat, difflambda, S_matrix);

    arma::mat J    = J_tmp["info_J"];
    arma::mat J_p  = J_tmp["info_J_p"];
    arma::mat J_p_GIC = J_tmp["info_J_p_gic"];

    arma::mat lambda_S_matrix = lambda_i*S_matrix;
    arma::mat info_lambda = info + lambda_S_matrix;

    arma::mat info_lambda_inv;
    if (method == "ProxN")
    {
      info_lambda.diag() += Z_tv.n_rows / 1e9;
      info_lambda_inv = inv(info_lambda);
      info_lambda.diag() -= Z_tv.n_rows / 1e9;
    }
    else {
      info_lambda_inv = inv(info_lambda);
    }

    arma::mat VarianceMatrix = info_lambda_inv*J_p*info_lambda_inv;
    arma::mat VarianceMatrix2 = info_lambda_inv*J*info_lambda_inv;     
    arma::mat VarianceMatrix_I = info_lambda_inv;
    // arma::mat VarianceMatrix_J = inv(J_p);
    arma::mat VarianceMatrix_gray = info_lambda_inv*info*info_lambda_inv;


    return List::create(_["theta"]=theta,
                        _["VarianceMatrix"]=VarianceMatrix,
                        _["VarianceMatrix2"]=VarianceMatrix2,
                        _["VarianceMatrix_I"]=VarianceMatrix_I,
                        // _["VarianceMatrix_J"]=VarianceMatrix_J,
                        _["VarianceMatrix_gray"]=VarianceMatrix_gray,
                        _["SplineType"]=SplineType,
                        _["lambda_i"]=lambda_i);
}

// Variance Matrix Calculation
// [[Rcpp::export]]
List VarianceMatrixCalculate_bresties(const arma::vec &event, const arma::vec &time, 
                                      const arma::mat &Z_tv, 
                                      const arma::mat &B_spline,
                                      arma::mat &theta,
                                      double lambda_i,
                                      const arma::mat &SmoothMatrix,
                                      const std::string &SplineType = "smooth-spline",
                                      const std::string &method="Newton",
                                      const double &lambda=1e8, const double &factor = 1.0,
                                      const bool &parallel=false, const unsigned int &threads=1){


    int p     = Z_tv.n_cols; //dimension
    int K     = B_spline.n_cols; //number of knots  
    int N     = Z_tv.n_rows;
    arma::mat Z_ti = arma::zeros<arma::mat>(N,p);
    arma::vec beta_ti = arma::zeros<arma::vec> (p);
    bool difflambda = false, ti = false, TIC_prox = false;

    IntegerVector count_strata = IntegerVector::create(N);
    IntegerVector cumsum_strata = cumsum(count_strata);
    unsigned int n_strata = cumsum_strata.length();
    cumsum_strata.push_front(0);

    vector<vector<arma::uvec>> idx_fail;
    vector<vector<unsigned int>> idx_fail_1st, idx_Z_strata;
    vector<unsigned int> n_fail_time, n_Z_strata;
    // each element of idx_Z_strata contains start/end row indices of Z_strata
    for (unsigned int i = 0; i < n_strata; ++i) {
      arma::vec event_tmp = event.rows(cumsum_strata[i], cumsum_strata[i+1]-1),
      time_tmp = time.rows(cumsum_strata[i], cumsum_strata[i+1]-1);
      arma::uvec idx_fail_tmp = find(event_tmp==1);
      vector<unsigned int> idx_Z_strata_tmp;
      idx_Z_strata_tmp.push_back(cumsum_strata[i]+idx_fail_tmp(0));
      idx_Z_strata_tmp.push_back(cumsum_strata[i+1]-1);
      idx_Z_strata.push_back(idx_Z_strata_tmp);
      n_Z_strata.push_back(cumsum_strata[i+1]-cumsum_strata[i]-
        idx_fail_tmp(0));
      arma::vec time_fail_tmp = time_tmp.elem(idx_fail_tmp);
      idx_fail_tmp -= idx_fail_tmp(0);
      arma::vec uniq_t = unique(time_fail_tmp);
      n_fail_time.push_back(uniq_t.n_elem);
      vector<arma::uvec> idx_fail_tmp_tmp;
      vector<unsigned int> idx_fail_1st_tmp;
      for (arma::vec::iterator j = uniq_t.begin(); j < uniq_t.end(); ++j) {
        arma::uvec tmp = idx_fail_tmp.elem(find(time_fail_tmp==*j));
        idx_fail_tmp_tmp.push_back(tmp);
        idx_fail_1st_tmp.push_back(tmp[0]);
      }
      idx_fail.push_back(idx_fail_tmp_tmp);
      idx_fail_1st.push_back(idx_fail_1st_tmp);
    }
    IntegerVector n_failtime = wrap(n_fail_time);
    IntegerVector cumsum_failtime = cumsum(n_failtime);
    cumsum_failtime.push_front(0);
    vector<arma::uvec> idx_B_sp;
    for (unsigned int i = 0; i < n_strata; ++i) {
      idx_B_sp.push_back(arma::regspace<arma::uvec>(cumsum_failtime[i],
                                        cumsum_failtime[i+1]-1));
    }
    // istart and iend for each thread when parallel=true
    vector<arma::vec>  cumsum_ar;
    vector<vector<unsigned int>> istart, iend;
    if (parallel) {
      for (unsigned int i = 0; i < n_strata; ++i) {
        double scale_fac = as_scalar(idx_fail[i].back().tail(1));
        cumsum_ar.push_back(
          (double)n_Z_strata[i] / scale_fac * arma::regspace(1, idx_fail[i].size()) -
            arma::cumsum(arma::conv_to<arma::vec> ::from(idx_fail_1st[i])/scale_fac));
        vector<unsigned int> istart_tmp, iend_tmp;
        for (unsigned int id = 0; id < threads; ++id) {
          istart_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
            cumsum_ar[i](idx_fail[i].size()-1)/(double)threads*id, 1)));
          iend_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
            cumsum_ar[i](idx_fail[i].size()-1)/(double)threads*(id+1), 1)));
          if (id == threads-1) {
            iend_tmp.pop_back();
            iend_tmp.push_back(idx_fail[i].size());
          }
        }
        istart.push_back(istart_tmp);
        iend.push_back(iend_tmp);
      }
    }

    arma::vec lambda_effct = arma::zeros<arma::vec> (p*K);
    arma::mat lambda_i_mat = repmat(lambda_effct,1,p*K);

    List SplineUdpate;
    arma::mat S_matrix;
    if(SplineType == "pspline") {
      S_matrix        = spline_construct(K, p, SplineType);  
    }
    else{
      S_matrix        = spline_construct2(K, p, SplineType, SmoothMatrix);  
    }    

    SplineUdpate    = stepinc_fixtra_spline_bresties(Z_tv, B_spline, theta, Z_ti, beta_ti, 
                                                    S_matrix, lambda_i, lambda_i_mat, difflambda,
                                                    ti, n_strata,
                                                    idx_B_sp, idx_fail, idx_Z_strata, istart, iend,
                                                    method, lambda, parallel, threads);

    arma::mat info = SplineUdpate["info"];

    List J_tmp;
    J_tmp = TIC_J_penalized_second_bresties(Z_tv, B_spline, theta, n_strata, idx_B_sp, idx_fail, idx_Z_strata, TIC_prox, 
                                            lambda_i, lambda_i_mat, difflambda, S_matrix);

    arma::mat J       = J_tmp["info_J"];
    arma::mat J_p     = J_tmp["info_J_p"];
    arma::mat J_p_GIC = J_tmp["info_J_p_gic"];

    arma::mat lambda_S_matrix = lambda_i*S_matrix;
    arma::mat info_lambda = info + lambda_S_matrix;

    arma::mat info_lambda_inv;
    if (method == "ProxN")
    {
      info_lambda.diag() += Z_tv.n_rows / 1e9;
      info_lambda_inv = inv(info_lambda);
      info_lambda.diag() -= Z_tv.n_rows / 1e9;
    }
    else {
      info_lambda_inv = inv(info_lambda);
    }

    arma::mat VarianceMatrix = info_lambda_inv*J_p*info_lambda_inv;
    arma::mat VarianceMatrix2 = info_lambda_inv*J*info_lambda_inv;     
    arma::mat VarianceMatrix_I = info_lambda_inv;
    // arma::mat VarianceMatrix_J = inv(J_p);
    arma::mat VarianceMatrix_gray = info_lambda_inv*info*info_lambda_inv;

    return List::create(_["theta"]=theta,
                        _["VarianceMatrix"]=VarianceMatrix,
                        _["VarianceMatrix2"]=VarianceMatrix2,
                        _["VarianceMatrix_I"]=VarianceMatrix_I,
                        // _["VarianceMatrix_J"]=VarianceMatrix_J,
                        _["VarianceMatrix_gray"]=VarianceMatrix_gray,
                        _["SplineType"]=SplineType,
                        _["lambda_i"]=lambda_i);
}



// this function is for lambda from large to small
List spline_udpate_lambdafromlarge(const arma::mat &Z_tv, const arma::mat &B_spline, arma::mat &theta, 
                  const arma::mat &Z_ti, const arma::vec &beta_ti, 
                  arma::mat &S_matrix,
                  double &lambda_i,
                  arma::mat &lambda_i_mat, arma::mat &lambda_S_matrix, 
                  const bool &difflambda,
                  const bool &ti,
                  const unsigned int n_strata,
                  vector<arma::uvec> &idx_B_sp, vector<arma::uvec> &idx_fail,
                  vector<unsigned int> n_Z_strata,
                  vector<vector<unsigned int>> &idx_Z_strata,
                  vector<vector<unsigned int>> &istart,
                  vector<vector<unsigned int>> &iend,
                  const std::string &method="Newton", const double &lambda=1e8, const double &factor = 1.0,
                  const bool &parallel=false, const unsigned int &threads=1,
                  const unsigned int &iter_max=20,
                  const double &tol=1e-10, 
                  const double &s=1e-2, const double &t=0.6,
                  const std::string &btr="dynamic",
                  const std::string &stop="incre",
                  const bool &TIC_prox = false,
                  const bool &fixedstep = false,
                  const bool &ICLastOnly = false) {

  int N     = Z_tv.n_rows;    //sample size 
  int p     = Z_tv.n_cols;
  int K     = theta.n_cols;

  unsigned int iter = 0, btr_max = 1000, btr_ct = 0, iter_NRstop = 0;
  bool NRstop = false;
  double crit = 1.0, v = 1.0, logplkd_init = 0, logplkd, diff_logplkd, 
    inc, rhs_btr = 0;
  List objfun_list, update_list;
  NumericVector logplkd_vec;
  arma::mat theta_init_allzero = arma::zeros<arma::mat>(p, K);
  objfun_list = objfun_fixtra(Z_tv, B_spline, theta_init_allzero, Z_ti, beta_ti, ti, n_strata,
                            idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata,
                            istart, iend, parallel, threads);
  logplkd = objfun_list["logplkd"];
  logplkd -= lambda_i*as_scalar(vectorise(theta,1)*S_matrix*vectorise(theta,1).t())/(N*1.0);
  List theta_list = List::create(theta);
    
  NumericVector AIC_all, TIC_all, TIC2_all, GIC_all;
  arma::mat VarianceMatrix;

  while (iter < iter_max) {
    ++iter;
    update_list = stepinc_fixtra_spline(Z_tv, B_spline, theta, Z_ti, beta_ti, 
                                        S_matrix, lambda_i, lambda_i_mat, difflambda,
                                        ti, n_strata,
                                        idx_B_sp, idx_fail, idx_Z_strata, istart, iend,
                                        method, lambda*pow(factor,iter-1), parallel, threads);
    v = 1.0; // reset step size
    arma::mat step = update_list["step"];
    inc = update_list["inc"];
    if (btr=="none") {
      theta += step;
      crit = inc / 2;
      objfun_list = objfun_fixtra(Z_tv, B_spline, theta, Z_ti, beta_ti, ti,
                                  n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                  idx_Z_strata, istart, iend, parallel, threads);
      logplkd = objfun_list["logplkd"];
      logplkd -= lambda_i*as_scalar(vectorise(theta,1)*S_matrix*vectorise(theta,1).t())/(N*1.0);
    } else {
      arma::mat theta_tmp = theta + step;
      objfun_list = objfun_fixtra(Z_tv, B_spline, theta_tmp, Z_ti, beta_ti, ti,
                              n_strata, idx_B_sp, idx_fail, n_Z_strata,
                              idx_Z_strata, istart, iend, parallel,
                              threads);
      double logplkd_tmp = objfun_list["logplkd"];
      logplkd_tmp -= lambda_i*as_scalar(vectorise(theta_tmp,1)*S_matrix*vectorise(theta_tmp,1).t())/(N*1.0);
      diff_logplkd = logplkd_tmp - logplkd;
      if (btr=="dynamic")      rhs_btr = inc;
      else if (btr=="static")  rhs_btr = 1.0;
      //btr_ct = 0;
      while (diff_logplkd < s * v * rhs_btr && btr_ct < btr_max) {
        ++btr_ct;
        v *= t;
        theta_tmp   = theta + v * step;
        objfun_list = objfun_fixtra(Z_tv, B_spline, theta_tmp, Z_ti, beta_ti, ti,
                                    n_strata, idx_B_sp, idx_fail, n_Z_strata,
                                    idx_Z_strata, istart, iend, parallel,
                                    threads);
        double logplkd_tmp = objfun_list["logplkd"];
        logplkd_tmp -= lambda_i*as_scalar(vectorise(theta_tmp,1)*S_matrix*vectorise(theta_tmp,1).t())/(N*1.0);
        diff_logplkd  = logplkd_tmp - logplkd;
      }
      theta = theta_tmp;
      //if (iter==1) logplkd_init = logplkd;
      if (stop=="incre")
        crit = inc / 2;
      else if (stop=="relch")
        crit = abs(diff_logplkd/(diff_logplkd+logplkd));
      else if (stop=="ratch")
        crit = abs(diff_logplkd/(diff_logplkd+logplkd-logplkd_init));
      logplkd += diff_logplkd;
    }
    Rcout << "Iter " << iter << ": Obj fun = " << setprecision(7) << fixed << 
     logplkd << "; Stopping crit = " << setprecision(7) << scientific << 
       crit << ";" << endl;
    logplkd_vec.push_back(logplkd);
    theta_list.push_back(theta);

    //force the algorithms to run specified steps
    if(!fixedstep) {
      if(crit < tol)
        break;
    }
    //record NR's stopping step:
    if(NRstop==false){
      if(crit < tol){
        iter_NRstop = iter;
        NRstop=true;
      }
    }

  }

  objfun_list = objfun_fixtra(Z_tv, B_spline, theta, Z_ti, beta_ti, ti,
                              n_strata, idx_B_sp, idx_fail, n_Z_strata,
                              idx_Z_strata, istart, iend, parallel, threads);
  logplkd = objfun_list["logplkd"];
  if(ICLastOnly == true){
    List J_tmp;
    J_tmp = TIC_J_penalized_second(Z_tv, B_spline, theta, n_strata, idx_B_sp, idx_fail, idx_Z_strata, TIC_prox, 
                            lambda_i, lambda_i_mat, difflambda, S_matrix);
    arma::mat J    = J_tmp["info_J"];
    arma::mat J_p  = J_tmp["info_J_p"];
    arma::mat J_p_GIC = J_tmp["info_J_p_gic"];

    double df;
    arma::mat info = update_list["info"];
    arma::mat info_lambda = info + lambda_S_matrix;
    arma::mat info_lambda_inv;
    info_lambda_inv = inv(info_lambda);
    //AIC:
    df = trace(info * info_lambda_inv);
    AIC_all.push_back(-2*logplkd*N + 2*df);
    //TIC
    df = trace(info*info_lambda_inv*J*info_lambda_inv);
    TIC_all.push_back(-2*logplkd*N + 2*df);
    //TIC2:
    df = trace(info*info_lambda_inv*J_p*info_lambda_inv);
    TIC2_all.push_back(-2*logplkd*N + 2*df);
    //GIC:
    df = trace(info_lambda_inv*J_p_GIC);
    GIC_all.push_back(-2*logplkd*N + 2*df);

    VarianceMatrix = info_lambda_inv*J_p*info_lambda_inv;
  }
  else{
    arma::mat info = update_list["info"];
    VarianceMatrix = info;
  }  

  arma::mat info = update_list["info"];
  arma::vec grad = update_list["grad"];

  return List::create(_["theta"]=theta,
                      _["logplkd"]=logplkd, 
                      _["info"]=info, 
                      _["grad"]=grad,
                      _["logplkd_vec"]=logplkd_vec,
                      _["AIC_all"]=AIC_all,
                      _["TIC_all"]=TIC_all,
                      _["TIC2_all"]=TIC2_all,
                      _["GIC_all"]=GIC_all,
                      _["theta_list"]=theta_list,
                      _["iter_NRstop"]=iter_NRstop,
                      _["VarianceMatrix"]=VarianceMatrix);
}



// this function is for lambda from large to small
// [[Rcpp::export]]
List surtiver_fixtra_fit_penalizestop_lambdafromlarge(const arma::vec &event, const IntegerVector &count_strata,
                             const arma::mat &Z_tv, const arma::mat &B_spline, const arma::mat &theta_init,
                             const arma::mat &Z_ti, const arma::vec &beta_ti_init,
                             const arma::vec &lambda_spline,
                             const arma::mat &SmoothMatrix,
                             const arma::vec &effectsize,
                             const std::string &SplineType = "pspline",
                             const std::string &method="Newton",
                             const double lambda=1e8, const double &factor=1.0,
                             const bool &parallel=false, const unsigned int &threads=1,
                             const double &tol=1e-10, const unsigned int &iter_max=20,
                             const double &s=1e-2, const double &t=0.6,
                             const std::string &btr="dynamic",
                             const std::string &stop="incre",
                             const bool &TIC_prox = false,
                             const bool &fixedstep = true,
                             const bool &difflambda = false,
                             const bool &ICLastOnly = false) {

  bool ti = arma::norm(Z_ti, "inf") > sqrt(arma::datum::eps);
  IntegerVector cumsum_strata = cumsum(count_strata);
  unsigned int n_strata = cumsum_strata.length();
  cumsum_strata.push_front(0);

  vector<arma::uvec> idx_fail, idx_B_sp;
  vector<vector<unsigned int>> idx_Z_strata;
  // each element of idx_Z_strata contains start/end row indices of Z_strata
  vector<unsigned int> n_fail, n_Z_strata;
  for (unsigned int i = 0; i < n_strata; ++i) {
    arma::uvec idx_fail_tmp =
      find(event.rows(cumsum_strata[i], cumsum_strata[i+1]-1)==1);
    n_fail.push_back(idx_fail_tmp.n_elem);
    vector<unsigned int> idx_Z_strata_tmp;
    idx_Z_strata_tmp.push_back(cumsum_strata[i]+idx_fail_tmp(0));
    idx_Z_strata_tmp.push_back(cumsum_strata[i+1]-1);
    idx_Z_strata.push_back(idx_Z_strata_tmp);
    n_Z_strata.push_back(cumsum_strata[i+1]-cumsum_strata[i]-
      idx_fail_tmp(0));
    idx_B_sp.push_back(cumsum_strata[i]+idx_fail_tmp);
    idx_fail_tmp -= idx_fail_tmp(0);
    idx_fail.push_back(idx_fail_tmp);
  }

  // istart and iend for each thread when parallel=true
  vector<arma::vec>  cumsum_ar;
  vector<vector<unsigned int>> istart, iend;
  if (parallel) {
    for (unsigned int i = 0; i < n_strata; ++i) {
      double scale_fac = as_scalar(idx_fail[i].tail(1));
      cumsum_ar.push_back(
        (double)n_Z_strata[i] / scale_fac * arma::regspace(1,n_fail[i]) -
          arma::cumsum(arma::conv_to<arma::vec> ::from(idx_fail[i])/scale_fac));
      vector<unsigned int> istart_tmp, iend_tmp;
      for (unsigned int id = 0; id < threads; ++id) {
        istart_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
          cumsum_ar[i](n_fail[i]-1)/(double)threads*id, 1)));
        iend_tmp.push_back(as_scalar(find(cumsum_ar[i] >=
          cumsum_ar[i](n_fail[i]-1)/(double)threads*(id+1), 1)));
        if (id == threads-1) {
          iend_tmp.pop_back();
          iend_tmp.push_back(n_fail[i]);
        }
      }
      istart.push_back(istart_tmp);
      iend.push_back(iend_tmp);
    }
  }

  arma::mat theta = theta_init; arma::vec beta_ti = beta_ti_init;
  // double crit = 1.0, v = 1.0, logplkd_init = 0, logplkd, diff_logplkd, 
  //   inc, rhs_btr = 0;
  double logplkd;
  List objfun_list, update_list;
  NumericVector logplkd_vec, iter_NR_all;
  objfun_list = objfun_fixtra(Z_tv, B_spline, theta, Z_ti, beta_ti, ti, n_strata,
                              idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata,
                              istart, iend, parallel, threads);
  logplkd = objfun_list["logplkd"];
  List theta_list;
  
  //for tic:
  int N     = Z_tv.n_rows;    //sample size 
  int p     = theta_init.n_rows; //dimension
  int K     = theta_init.n_cols; //number of knots   
  
  List SplineUdpate, J_tmp;

  int n_lambda         = lambda_spline.n_elem;
  arma::cube theta_all(p, K, n_lambda);   

  arma::vec AIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec BIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC  = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC2 = arma::zeros<arma::vec> (n_lambda);
  arma::vec GIC  = arma::zeros<arma::vec> (n_lambda);

  arma::vec AIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec BIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec TIC2_secondterm = arma::zeros<arma::vec> (n_lambda);
  arma::vec GIC_secondterm = arma::zeros<arma::vec> (n_lambda);
  NumericVector AIC_all, TIC_all, TIC2_all, GIC_all;
  //mat likelihood_mat = arma::zeros<arma::mat>(iter_max, n_lambda);    for fixed step usage

  arma::mat S_matrix;
  if(SplineType == "pspline") {
    S_matrix        = spline_construct(K, p, SplineType);  
  }
  else{
    S_matrix        = spline_construct2(K, p, SplineType, SmoothMatrix);  
  }

  arma::vec lambda_effct = arma::zeros<arma::vec> (p*K);
  arma::mat lambda_i_mat;
  //List VarianceMatrix, VarianceMatrix2, VarianceMatrix_I, VarianceMatrix_J;
  arma::mat VarianceMatrix;

  for (int i = 0; i < n_lambda; ++i)
  {
    //generate a lambda matrix if different lambda is used for different covariate:
    if(difflambda){

      for (int j = 0; j < p; ++j)
      {
        arma::vec sublambda = arma::zeros<arma::vec> (K);
        for (int jj = 0; jj < K; ++jj)
        {
          sublambda(jj) = lambda_spline(i) * sqrt(N /  effectsize(j));
        }
        lambda_effct.subvec(j*K, ((j+1)*K-1)) = sublambda;
      }      
    }
    lambda_i_mat = repmat(lambda_effct,1,p*K);
    arma::mat lambda_S_matrix;

    //use the theta estimated from the last lambda
    arma::mat theta_ilambda;
    if(i>0){
      arma::mat theta_ilambda  = theta_all.slice(i-1);
    }
    else{
      theta_ilambda = theta_init;
    }
    double lambda_i = lambda_spline(i);
    if(difflambda == false){
      lambda_S_matrix = lambda_i*S_matrix;
    }
    else{
      lambda_S_matrix = lambda_i_mat%S_matrix;
    }

    SplineUdpate    = spline_udpate_lambdafromlarge(Z_tv, B_spline, theta_ilambda,
                                    Z_ti, beta_ti, 
                                    S_matrix, lambda_i, lambda_i_mat, lambda_S_matrix, difflambda,
                                    ti, n_strata, idx_B_sp, idx_fail, n_Z_strata, idx_Z_strata, istart, iend,
                                    method, lambda, factor,
                                    parallel, threads, iter_max, tol, s, t, btr, stop, TIC_prox, fixedstep, ICLastOnly);

    //mat info              = SplineUdpate["info"];
    //vec grad              = SplineUdpate["grad"];
    double logplkd        = SplineUdpate["logplkd"];
    AIC_all               = SplineUdpate["AIC_all"];
    TIC_all               = SplineUdpate["TIC_all"];
    TIC2_all               = SplineUdpate["TIC2_all"];
    GIC_all               = SplineUdpate["GIC_all"];
    theta_list            = SplineUdpate["theta_list"];
    int iter_tmp          = SplineUdpate["iter_NRstop"];
    arma::mat VarianceMatrix_tmp        = SplineUdpate["VarianceMatrix"];

    theta_all.slice(i)    = theta_ilambda;
    logplkd_vec.push_back(logplkd);
    iter_NR_all.push_back(iter_tmp);
    VarianceMatrix = VarianceMatrix_tmp;

    //Variance matrix:
    // arma::mat VarianceMatrix_tmp = info_lambda_inv*J_p*info_lambda_inv;
    // arma::mat VarianceMatrix_tmp2 = info_lambda_inv*J*info_lambda_inv;     
    // VarianceMatrix.push_back(VarianceMatrix_tmp);
    // VarianceMatrix2.push_back(VarianceMatrix_tmp2);
    // VarianceMatrix_I.push_back(info_lambda_inv);
    // VarianceMatrix_J.push_back(inv(J_p));
    Rcout<<fixed<<"current lambda done: "<< lambda_spline(i) <<endl;

  }

  return List::create(_["theta"]=theta,
                      _["logplkd"]=logplkd, 
                      _["theta_all"]=theta_all,
                      _["theta_list"]=theta_list,
                      _["AIC_all"]=AIC_all,
                      _["TIC_all"]=TIC_all,
                      _["TIC2_all"]=TIC2_all,
                      _["GIC_all"]=GIC_all,
                      _["logplkd_vec"]=logplkd_vec,
                      _["SplineType"]=SplineType,
                      _["iter_NR_all"]=iter_NR_all,
                      _["VarianceMatrix"]=VarianceMatrix);
}


// main version (ignore above)
// [[Rcpp::export]]    
List Lambda_estimate_ties2(int knot,
                           arma::colvec &delta,
                           arma::mat &z,
                           arma::mat &b_spline, 
                           arma::mat &theta, 
                           arma::vec &tieseq){
  
  int n          = delta.n_elem;
  int ntie       = tieseq.n_elem;
  //int p          = theta.n_rows;
  
  arma::vec S0         = arma::zeros<arma::vec>(ntie);
  arma::vec lambda     = arma::zeros<arma::vec>(ntie);
  
  arma::mat z_temp;
  arma::mat zbeta_mat;
  
  arma::vec cumtie     = cumsum(tieseq);
  arma::vec delta_temp;
  
  int       index;
  for (int i = 0; i < ntie; ++i)
  {
    if(i == 0){
      index = i;
    }
    else{
      index = cumtie(i-1);
    }
    z_temp     = z.rows(index, n-1);
    zbeta_mat  = exp(z_temp * theta * b_spline.row(i).t());
    S0(i)      = accu(zbeta_mat);
    delta_temp = delta.subvec(index,cumtie(i)-1);
    lambda(i)  = accu(delta_temp)/S0(i);
  }
  
  
  
  List result;
  result["S0"]        = S0;
  result["lambda"]    = lambda;
  result["z_temp"]    = z_temp;
  result["zbeta_mat"] = zbeta_mat;
  result["delta_temp"]= delta_temp;
  result["cumtie"]    = cumtie;
  result["ntie"]      = ntie;
  result["index"]     = index;
  return result;
}



