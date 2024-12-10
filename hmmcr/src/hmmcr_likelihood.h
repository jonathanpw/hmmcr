#pragma once

#include <RcppArmadillo.h>
#include "hmmcr_response.h"
#include "hmmcr_bookkeeping.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

double log_post(const arma::vec& prior_par,
                const Para&      para,
                const arma::vec& y,
                const arma::mat& X,
                const arma::mat& Z,
                const arma::mat& S,
                const arma::vec& id,
                std::string model);
