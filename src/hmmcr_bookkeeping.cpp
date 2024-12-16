// Para.h
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

#include "hmmcr_bookkeeping.h"

Para::Para(const arma::vec& par_init,
           const arma::vec& par_index,
           const arma::vec& par_unique,
           const arma::mat& Z,
           const arma::mat& X,
           const int nr_state,
           const std::string model)
  : par_index(par_index),
    par_unique(par_unique),
    nr_state(nr_state),
    nr_col_Z(Z.n_cols),
    nr_col_X(X.n_cols),
    nr_col_theta(1),
    m_par(par_init.n_elem),
    init(nr_state),
    zeta(nr_state * (nr_state - 1), nr_col_Z),
    beta(nr_state, nr_col_X)
{
  m_par.zeros(par_init.n_elem);
  init.ones(nr_state);
  zeta.zeros(nr_state * (nr_state - 1), nr_col_Z);
  beta.zeros(nr_state, nr_col_X);

  if (model == "ingarch") {
    nr_col_theta = 2;
  }

  theta.zeros(nr_state, nr_col_theta);

  setPara(par_init);
}

void Para::setPara(const arma::vec& par) {

  m_par = par;

  // indexes
  std::vector<size_t> k = {1, 0, 0, 0};
  std::vector<size_t> l = {0, 0, 0, 0};

  for (size_t i = 0; i < par_unique.size(); ++i) {

    int par_index_i = par_index[i];
    int par_unique_i = par_unique[i];

    // get values or set to default
    double value = (par_unique_i > 0) ? m_par(par_unique_i - 1) : 0.0;

    switch (par_index_i) {

    case 0: { // init
      init(k[0]++) = exp(value);
      break;
    }
    case 1: { // zeta
      zeta(k[1], l[1]++) = value;
      if (l[1] >= nr_col_Z) {
        l[1] = 0;
        k[1]++;
      }
      break;
    }
    case 2: { // beta
      beta(k[2], l[2]++) = value;
      if (l[2] >= nr_col_X) {
        l[2] = 0;
        k[2]++;
      }
      break;
    }
    case 3: { // theta
      theta(k[3], l[3]++) = exp(value);
      if (l[3] >= nr_col_theta) {
        l[3] = 0;
        k[3]++;
      }
      break;
    }
    default:
      throw std::runtime_error("Invalid value in par_index.");
    }
  }
}

arma::vec Para::getPara() const {
  return m_par;
}

arma::vec Para::getInit() const {
  return init / arma::sum(init);
}

arma::mat Para::getZeta() const {
  return zeta;
}

arma::rowvec Para::getZeta(size_t i, size_t j) const {
  if (i == j) {
    throw std::invalid_argument("Diagonal elements are not part of zeta.");
  }
  if (i > nr_state || j > nr_state) {
    throw std::out_of_range("State indices out of range.");
  }

  size_t index = (i - 1) * (nr_state - 1) + (j - 1) - (i <= j);

  std::cout << index << std::endl;

  if (index >= zeta.n_rows) {
    throw std::out_of_range("Computed zeta index is out of range.");
  }
  return zeta.row(index);
}

arma::mat Para::getBeta() const {
  return beta;
}

arma::rowvec Para::getBeta(size_t state) const {
  if (state >= beta.n_rows) {
    throw std::out_of_range("State index out of range in beta.");
  }
  return beta.row(state);
}

arma::mat Para::getTheta() const {
  return theta;
}

arma::rowvec Para::getTheta(size_t state) const {
  if (state >= theta.n_rows) {
    throw std::out_of_range("State index out of range in theta.");
  }
  return theta.row(state);
}


