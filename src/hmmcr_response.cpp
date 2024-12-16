#include <RcppArmadillo.h>
#include "hmmcr_bookkeeping.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>


double subset_mean(const arma::vec& b, int k, int l) {
  return (k >= 0 && l < b.size() && k <= l) ? arma::mean(b.subvec(k, l)) : 0;
}

void compute_Pi_matrix(arma::mat& Pi, const arma::vec& z, const Para& para) {

  arma::mat Q = arma::zeros<arma::mat>(Pi.n_rows, Pi.n_cols);
  arma::vec q = para.getZeta()*z;

  int l = 0;

  for (int i = 0; i < Pi.n_rows; i++) {
    for (int j = 0; j < Pi.n_cols; j++) {
      Q(i, j) = (i == j) ? 1 : exp(q[l++]);
    }
  }

  arma::vec Q_rowSum = sum(Q, 1);

  for (int i = 0; i < Pi.n_rows; i++) {
    Pi.row(i) = Q.row(i) / Q_rowSum(i);
  }
}

void compute_D_matrix(arma::mat &D, const arma::vec& y, const arma::mat& x, const Para& para, const int k, const std::string model) {

  // Make D a clean diagonal matrix
  D.eye();

  // S state NB-INGARCH model ------------------------------------------------------------------
  if (model == "ingarch") {
    if (k > 3) {

      for (int i = 0; i < D.n_rows; i++) {

        double size_i = para.getTheta(i)(0);
        double c_i = para.getTheta(i)(1);

        if (i > 0) {
          size_i = size_i + exp(arma::dot(x.row(k), para.getBeta(i))) * subset_mean(y, k - 4, k - 1);
        }

        D(i, i) = R::dnbinom(y(k), size_i, c_i / (1 + c_i), false);
      }
    }
  }
  // -----------------------------------------------------------------------------------------------


  // Poisson ---------------------------------------------------------------------------------------
  if (model == "poisson") {

    // std::cout << "poisson" << std::endl;

    for (int i = 0; i < D.n_rows; i++) {

      // std::cout << exp(arma::dot(x.row(k), para.getBeta(i))) << std::endl;

      D(i, i) = R::dpois(y(k), exp(arma::dot(x.row(k), para.getBeta(i))), false);
    }
  }
  // -----------------------------------------------------------------------------------------------


  // Gaussian --------------------------------------------------------------------------------------
  if (model == "gaussian") {
    for (int i = 0; i < D.n_rows; i++) {
      double sigma_i = para.getTheta(i)(0);
      D(i, i) = R::dnorm(y(k), arma::dot(x.row(k), para.getBeta(i)), sigma_i, false);
    }
  }
  // -----------------------------------------------------------------------------------------------


  // log-Gaussian ----------------------------------------------------------------------------------
  if (model == "log-gaussian") {
    for (int i = 0; i < D.n_rows; i++) {
      double sigma_i = para.getTheta(i)(0);
      D(i, i) = R::dlnorm(y(k), arma::dot(x.row(k), para.getBeta(i)), sigma_i, false);
    }
  }
  // -----------------------------------------------------------------------------------------------


  // NB --------------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------

}

















