#pragma once

#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

class Para {
private:
  arma::vec m_par;
  arma::vec par_index;
  arma::vec par_unique;
  arma::vec init;
  arma::mat zeta;
  arma::mat beta;
  arma::mat theta;

  int nr_state;
  int nr_col_Z;
  int nr_col_X;
  int nr_col_theta;

public:
  // Constructor declaration
  Para(const arma::vec& par_init,
       const arma::vec& par_index,
       const arma::vec& par_unique,
       const arma::mat& Z,
       const arma::mat& X,
       const int nr_state,
       const std::string model);

  // Method declarations
  void setPara(const arma::vec& new_par);

  // Getters
  arma::vec getPara() const;
  arma::vec getInit() const;
  arma::mat getZeta() const;
  arma::rowvec getZeta(size_t i, size_t j) const;
  arma::mat getBeta() const;
  arma::rowvec getBeta(size_t state) const;
  arma::mat getTheta() const;
  arma::rowvec getTheta(size_t state) const;
};
