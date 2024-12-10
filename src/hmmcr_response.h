#pragma once

#include <RcppArmadillo.h>
#include "hmmcr_bookkeeping.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

void compute_Pi_matrix(arma::mat& Pi, const arma::vec& z, const Para& para);


void compute_D_matrix(arma::mat &D, 
                      const arma::vec& y, 
											const arma::mat& x, 
											const Para& para, 
											const int k, 
											const std::string model);
















