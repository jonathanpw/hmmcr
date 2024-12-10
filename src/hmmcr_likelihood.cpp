#include <RcppArmadillo.h>
#include "hmmcr_response.h"
#include "hmmcr_bookkeeping.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>


arma::mat P_state(const arma::vec& state_from, const arma::vec& state_to, const arma::mat& Pi) {

  // Find index to positive elements in state_from and state_to
  arma::uvec indices_from = arma::find(state_from > 0);
  arma::uvec indices_to = arma::find(state_to > 0);

  // Create sub-matrix
  arma::mat Pi_state = Pi(indices_from, indices_to);

  return Pi_state;
}


arma::mat D_state(const arma::vec& state_to, const arma::mat& Di) {

  // Find index to positive elements in state_to
  arma::uvec indices = arma::find(state_to > 0);

  // Create sub-matrix
  arma::mat Di_state = Di(indices, indices);

  return Di_state;
}


arma::vec f_init(const arma::vec& init, const arma::vec& state) {
  arma::vec f_init;

  for (int i = 0; i < init.size(); i++) {
    if (state[i] > 0) {
      f_init.resize(f_init.size() + 1);
      f_init(f_init.size() - 1) = init(i);
    }
  }
  return f_init;
}


double log_fn_i_fun(const Para&       para,
                    arma::vec&        f,
                    const arma::vec&  y,
                    const arma::mat&  x,
                    const arma::mat&  z,
                    const arma::mat&  S,
                    const std::string model) {

  double log_norm = 0;
  bool stop = false;

  // Note: Number of states is defined by length of "a"
  int nr_state = para.getInit().n_elem; // gud: move this out or put it into para

  // Initialize Pi, D, and other temporary variables
  // arma::mat D(nr_state, nr_state, arma::fill::zeros);
  arma::mat D(nr_state, nr_state, arma::fill::eye);

  arma::mat Pi(nr_state, nr_state, arma::fill::eye);

  arma::mat tmp_P, tmp_D;
  arma::vec v;
  double norm_v;

  // iterate over each time point, except the first
  for (int k = 1; k < y.size(); k++) {

    // the transition probability matrix
    compute_Pi_matrix(Pi, z.row(k).t(), para);

    // response function
    compute_D_matrix(D, y, x, para, k, model);

    // compute PiDi product
    tmp_P = P_state(S.row(k - 1).t(), S.row(k).t(), Pi);
    tmp_D = D_state(S.row(k).t(), D);

    // compute v and norm of v
    v = (tmp_P * tmp_D).t() * f;
    norm_v = arma::norm(v);
    
    // check if the norm of v is zero
    if (norm_v == 0) {
      stop = true;
      break;
    } else {
      // update f and log_norm
      f = v / norm_v;
      log_norm += std::log(norm_v);
    }
  }
  
  // if stop is true, return -DBL_MAX
  if (stop) {
    return -INFINITY;
  } else {
    // otherwise, return the logarithm of the sum of f and log_norm
    return std::log(arma::sum(f)) + log_norm;
  }
}


// TODO: Add "model" to log_fn_fun, log_fn_i_fun, D_state and Pi matrix functions.
double log_fn_fun(const arma::vec& id,
                  const Para&      para,
                  const arma::vec& y,
                  const arma::mat& x,
                  const arma::mat& z,
                  const arma::mat& S,
                  const std::string     model) {

  arma::vec y_i(y.n_elem);
  arma::mat x_i(x.n_rows, x.n_cols);
  arma::mat z_i(z.n_rows, z.n_cols);
  arma::mat S_i(S.n_rows, S.n_cols);
  arma::vec f_i(para.getInit().n_elem); // gud: change to nr_states

  int N = id.n_elem;
  int id_previous = id(0);
  int id_now = -1;

  int index_start = 0;
  int index_end = -1;

  double log_fn = 0;
  double log_fn_i = 0;

  // Should this start at i = 1 or i = 0?
  for (int i = 1; i < N; i++) {

    id_now = id(i);

    if (id_now != id_previous || i == N - 1) {
      

      index_end = i - 1;

      if (i == N - 1) index_end = i;

      y_i = y.subvec(index_start, index_end);
      x_i = x.rows(index_start, index_end);
      z_i = z.rows(index_start, index_end);
      S_i = S.rows(index_start, index_end);

      f_i = f_init(para.getInit(), S_i.row(0).t());
      log_fn_i = log_fn_i_fun(para, f_i, y_i, x_i, z_i, S_i, model);
      
      if (!std::isfinite(log_fn_i)) {
        log_fn = -INFINITY;
        break;
      } else {
        log_fn = log_fn + log_fn_i;
      }
      index_start = i;
    }
    id_previous = id_now;
  }
  return(log_fn);
}

double log_mvnorm(const arma::vec& x, const arma::vec& mean, const arma::mat& sigma) {
  arma::vec d = x - mean;
  double log_det_sigma;
  double sign;

  arma::log_det(log_det_sigma, sign, sigma);
  arma::mat inv_sigma = arma::inv_sympd(sigma);

  double log_density = -0.5 * arma::dot(d, inv_sigma * d) - 0.5 * log_det_sigma - (x.n_elem / 2.0) * log(2.0 * M_PI);

  return log_density;
}


double log_post(const arma::vec& prior_par,
                const Para&      para,
                const arma::vec& y,
                const arma::mat& X,
                const arma::mat& Z,
                const arma::mat& S,
                const arma::vec& id,
                std::string            model) {

  double log_fn = -INFINITY;
  double log_prior_dens = 0;

  
  // Model 0: Three hidden state (peace, conflict, escalation) with NB-INGARCH response function ---


  // if (a(0) <= a(1) && a(1) <= a(2) && beta(0, 0) <= beta(0, 1)) {
  // gud: move this part to the prior object! Make prior with model as input!

    // compute log-density
    log_fn = log_fn_fun(id, para, y, X, Z, S, model);

    // only compute log-prior density if log_fn is "finite"
    if (std::isfinite(log_fn)) {

      // TODO: make this a function and add model as an option
      arma::vec prior_mean = prior_par.rows(0, para.getPara().n_elem-1);
      arma::vec prior_sd = prior_par.rows( para.getPara().n_elem, para.getPara().n_elem*2-1);

      log_prior_dens = log_mvnorm(para.getPara(), prior_mean, arma::diagmat(prior_sd));
    } 
  // ---------------------------------------------------------------------------
  
  if (std::isfinite(log_fn) && std::isfinite(log_prior_dens)) {
    return log_fn + log_prior_dens; 
  } else {
    return -INFINITY;
  }
}
