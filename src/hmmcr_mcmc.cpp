#include <RcppArmadillo.h>
#include <vector>
#include <stdexcept>

using namespace Rcpp;

#include "hmmcr_likelihood.h"
#include "hmmcr_response.h"
#include "hmmcr_bookkeeping.h"


void debugPrint(const Para& manager) {
  std::cout << "=== Debug Output for ParManager ===" << std::endl;

  // Init
  std::cout << "Init vector" << std::endl;
  std::cout << manager.getInit() << std::endl;

  // Beta matrix
  const arma::mat& beta = manager.getBeta();
  std::cout << "Beta Matrix (" << beta.n_rows << " x " << beta.n_cols << "):" << std::endl;
  for (size_t i = 0; i < beta.n_rows; ++i) {
    for (size_t j = 0; j < beta.n_cols; ++j) {
      std::cout << std::setw(8) << beta(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Zeta matrix
  const arma::mat& zeta = manager.getZeta();
  std::cout << "Zeta Matrix (" << zeta.n_rows << " x " << zeta.n_cols << "):" << std::endl;
  for (size_t i = 0; i < zeta.n_rows; ++i) {
    for (size_t j = 0; j < zeta.n_cols; ++j) {
      std::cout << std::setw(8) << zeta(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Theta matrix
  const arma::mat& theta = manager.getTheta();
  std::cout << "Theta Matrix (" << theta.n_rows << " x " << theta.n_cols << "):" << std::endl;
  for (size_t i = 0; i < theta.n_rows; ++i) {
    for (size_t j = 0; j < theta.n_cols; ++j) {
      std::cout << std::setw(8) << theta(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "===================================" << std::endl;
}


// [[Rcpp::export]]
double log_Gaussian( arma::vec y, arma::vec pars) {
  return -arma::sum(pow(y-pars(0),2)/(2*exp(pars(1))) + .5*pars(1));
}



// [[Rcpp::export]]
Rcpp::List mcmc_routine( const arma::vec par_index,
                         const arma::vec par_init,
                         const arma::vec par_prior,
                         const arma::vec par_unique,
                         const arma::vec y,
                         const arma::mat X,
                         const arma::mat Z,
                         const arma::mat S,
                         const arma::vec id,
                         const std::string model,
                         const int nr_states,
                         const Rcpp::List groups,
                         const int steps,
                         const int burnin){


  // Helper lambda function
  auto mod = [](int a, int n) -> int { return a - floor(a/n)*n;	};

  const int n_pars = par_init.n_elem;
  arma::vec pars( n_pars, arma::fill::zeros);
  arma::vec proposal( n_pars, arma::fill::zeros);
  pars = par_init;
  arma::mat chain( steps, n_pars, arma::fill::zeros);
  chain.row(0) = pars.t();

  const int n_groups = groups.size();
  std::vector<arma::mat> pcov;
  for(int j = 0; j < n_groups; j++){
    arma::uvec ind_j = groups(j);
    pcov.push_back(arma::eye( ind_j.n_rows, ind_j.n_rows));
  }

  arma::vec pscale( n_groups, arma::fill::value(.01));
  arma::mat temp_chain( 1000, n_pars, arma::fill::zeros);

  // gudmund: New par object
  Para para(par_init, par_index, par_unique, Z, X, nr_states, model);
  para.setPara(par_init);

  // debugPrint(para); // gud

  // Evaluate the log_post of the initial pars
  double targ_prev = log_post(par_prior,
                              para,
                              y,
                              X,
                              Z,
                              S,
                              id,
                              model);

  if(!std::isfinite(targ_prev)){
    std::cout << targ_prev << "\n";
    std::cout << "Infinite log-posterior; change initial parameters" << "\n";
    return Rcpp::List::create();
  }

  arma::vec accept( n_groups, arma::fill::zeros);

  // Begin the MCMC algorithm --------------------------------------------------
  for(unsigned int ttt = 1; ttt < steps; ttt++){
    for(int j = 0; j < n_groups; j++){
      arma::uvec ind_j = groups(j);

      // Propose an update
      // proposal.col(0) = pars;
      proposal = pars; // gudmund:
      double pscale_j = pscale(j);
      arma::mat pcov_j = pcov[j];
      proposal.rows(ind_j) = arma::mvnrnd( pars.rows(ind_j), pscale_j*pcov_j);

      // Compute the log density for the proposal

      // start_time = std::chrono::steady_clock::now(); // gud: Starttiming

      para.setPara(proposal); // gud
      double targ = log_post(par_prior,
                             para,
                             y,
                             X,
                             Z,
                             S,
                             id,
                             model);

      // targ_debug = targ; // gud
      // end_time = std::chrono::steady_clock::now(); // gud

      // std::cout << targ << "\n"; // gud

      // Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!std::isfinite(targ)){
          std::cout << "bad proposal" << "\n";
          proposal.col(0) = pars;
          proposal.rows(ind_j) = arma::mvnrnd( pars.rows(ind_j), pscale_j*pcov_j);

          para.setPara(proposal); // gud
          targ = log_post(par_prior,
                          para,
                          y,
                          X,
                          Z,
                          S,
                          id,
                          model);
        }
      }

      // Evaluate the Metropolis-Hastings ratio
      if( targ - targ_prev > log(arma::randu()) ){
        targ_prev = targ;
        pars.rows(ind_j) = proposal(ind_j);
        accept(j) += 1;
      }
      chain.submat( arma::uvec {ttt}, ind_j) = pars.rows(ind_j).t();

      // Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        // During the burnin period update the proposal covariance in each step
        // to capture the relationships within the parameters vectors for each
        // transition.  This helps with mixing.
        if(ttt == 100)  pscale(j) = 1;

        if(998 < ttt){
          temp_chain = chain.rows(ttt-999,ttt);
          arma::uvec nonduplicate = find_unique(temp_chain.col(ind_j(0)));
          pcov[j] = cov(temp_chain( nonduplicate, ind_j));
        }
        if(pcov[j].n_rows==1) pcov[j] = arma::eye( ind_j.n_rows, ind_j.n_rows);

        // Tune the proposal covariance for each transition to achieve
        // reasonable acceptance ratios.
        if(mod(ttt,30) == 0){
          if(mod(ttt,480) == 0){
            accept(j) = 0;

          } else if( accept(j) / mod(ttt,480) < .4 ){
            pscale(j) = pow(.75,2)*pscale(j);

          } else if( accept(j) / mod(ttt,480) > .5 ){
            pscale(j) = pow(1.25,2)*pscale(j);
          }
        }
      }
      // -----------------------------------------------------------------------
    }


    // Restart the acceptance ratio at burnin.
    if(ttt == burnin)  accept.zeros();

    if(mod(ttt, 1000) == 0) {
      std::cout << ttt << "\n";
      // std::cout << targ_debug << "\n";
      //
      // std::chrono::duration<double> elapsed_time = end_time - start_time;
      // std::cout << "Iteration " << ttt
      //           << " completed in " << elapsed_time.count() << " seconds.\n";
      //
      // start_time = std::chrono::steady_clock::now();
    }
  }
  // ---------------------------------------------------------------------------

  debugPrint(para); // gud

  std::cout << accept/(steps-burnin) << "\n";
  return Rcpp::List::create( chain, pscale, pcov);
}
