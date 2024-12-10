# rm(list = ls())
#
# 
# seed = 1
# set.seed(seed)

library(Rcpp, quietly=T)
library(data.table, quietly=T)
Rcpp::sourceCpp("src/hmmcr_mcmc.cpp")


.hmmcr_convert_to_state_matrix <- function(s) {

  tmp <- function(s, max_length) {
    digits <- as.numeric(strsplit(s, "::")[[1]])
    result <- numeric(max_length)
    result[seq_along(digits)] <- digits
    return(result)
  }

  max_length <- max(sapply(strsplit(s, "::"), length))

  return(t(sapply(s, tmp, max_length = max_length)))
  # return(as.matrix(sapply(s, tmp, max_length = max_length)))
}

# .hmmcr_make_par_index <- function(nr_states, nr_zeta, nr_beta, nr_theta) {
#   c(rep(0, nr_states - 1), rep(1, nr_zeta), rep(2, nr_beta), rep(3, nr_theta))
# }
# 
# .hmmcr_make_par_init <- function(par_index) {
#   return(rnorm(n = length(par_index), mean = 0, sd = 20))
# }

.hmmcr_make_par <- function(nr_states, col_z, col_x, model) {

  # init
  par_index <- rep(0, nr_states - 1)
  par_init <- rnorm(nr_states - 1, mean = 0, sd = 1)
  par_unique <- 1:(nr_states - 1)
  
  if (model == "ingarch") {

    # zeta
    nr_zeta <- col_z*(nr_states)*(nr_states - 1)
    par_index <- c(par_index, rep(1, nr_zeta))
    par_init <- c(par_init, rep(0, nr_zeta))
    par_unique <- c(par_unique, max(par_unique) + 1:nr_zeta)

    # beta
    nr_beta <- col_x*nr_states
    par_index <- c(par_index, rep(2, nr_beta))
    par_init <- c(par_init, rep(0, col_x*(nr_states - 1))) # rnorm(nr_beta, mean = 0, sd = 20))
    par_unique <- c(par_unique, rep(-1, col_x))
    par_unique <- c(par_unique, max(par_unique) + 1:(col_x*(nr_states - 1))) # + 1:nr_beta)

    # theta
    nr_theta <- 2*nr_states # nr_states + 1
    par_index <- c(par_index, rep(3, nr_theta))
    par_init <- c(par_init, rnorm(nr_states + 1, mean = 0, sd = 1)) # rnorm(nr_theta, mean = 0, sd = 20))
    par_unique <- c(par_unique, max(par_unique) + c(rbind(1:nr_states, nr_states + 1)))
  }
  
  if (model == "poisson") {
    
    # zeta
    nr_zeta <- col_z*(nr_states)*(nr_states - 1)
    par_index <- c(par_index, rep(1, nr_zeta))
    par_init <- c(par_init, rep(0.0, nr_zeta))
    par_unique <- c(par_unique, max(par_unique) + 1:nr_zeta)
    
    # beta
    nr_beta <- col_x*nr_states 
    par_index <- c(par_index, rep(2, nr_beta))
    par_init <- c(par_init, rep(0.0, nr_beta))
    # par_init <- c(par_init, 5, 0, 5, 0, 5, 0)
    par_unique <- c(par_unique, max(par_unique) + 1:nr_beta) 

  }
  
  if (model == "gaussian") {
    
    # zeta
    nr_zeta <- col_z*(nr_states)*(nr_states - 1)
    par_index <- c(par_index, rep(1, nr_zeta))
    par_init <- c(par_init, rep(0.0, nr_zeta))
    par_unique <- c(par_unique, max(par_unique) + 1:nr_zeta)
    
    # beta
    nr_beta <- col_x*nr_states 
    par_index <- c(par_index, rep(2, nr_beta))
    par_init <- c(par_init, rep(0.0, nr_beta))
    par_unique <- c(par_unique, max(par_unique) + 1:nr_beta) 
    
    # theta
    nr_theta <- nr_states
    par_index <- c(par_index, rep(3, nr_theta))
    par_init <- c(par_init, rnorm(nr_states, mean = 0, sd = 1))
    par_unique <- c(par_unique, max(par_unique) + rbind(1:nr_states))
    
  }
  
  par_index_unique <- par_index[(par_unique > 0) & !duplicated(par_unique)]
  par_groups <-list(which(par_index_unique!=1)-1, c(which(par_index_unique==1))-1)
  

  return(list(par_index = par_index, par_init = par_init, par_unique = par_unique, par_groups = par_groups))
}

.hmmcr_make_prior_par <- function(par_index, model) {
  # prior_par <- NULL
  # if (model == "normal-ingarch") {
  prior_par <- c(rep(0, length(par_index)), rep(20, length(par_index)))
  # }
  return(prior_par)
}

.hmmcr_normalise_column <- function(x) {

  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- sd(x, na.rm = TRUE)

  if (x_sd != 0) {
    return((x - x_mean)/x_sd)
  } else {
    return(x)
  }
}

hmmcr <- function(formula_response,
                  formula_state,
                  index,
                  group,
                  model = "normal-ingarch",
                  burnin = 100,
                  steps = 1000,
                  par = list(init = NULL, index = NULL, unique = NULL, lower = NULL, upper = NULL, prior = NULL, group = NULL),
                  # prior = list(par, index), 
                  normalize = TRUE,
                  data = NULL) {

  # TODO: 
  # 1) introduce a par_lower and par_upper? 
  # 2) introducing constraints to parameters (via prior) is complicated? 
  # 3) make a prior class and in addition to model add prior as to select type of prior
  # 4) add model to Para or Prior class and add para conditions here. 
  # 5) add prior_par and prior_par_index
  # 6) change to model = "normal-ingarch" to "normal::ingarch::3" where the last number is the size of rolling window, then we need a new hyper_par and hyper_par_index thing for this
  # 7) check std::isfinite() and perhaps let log_post and log_fn_fun and log_fn_i_fun return std::infinite?
  # 8) rename log_fn_fun and log_fn_i_fun to log_f and log_f_i or log_density
  # 9) add check of input 
  # 10) make par_unique start at 0 
  
  mf_response <- model.frame(formula_response, data)
  mf_state <- model.frame(formula_state, data)


  # check stuff ----------------------------------------------------------------
  if (is.character(index)) {
    # data_index <- NULL 
    # try(eval(parse(text = paste("data_index <-", index))))
    stopifnot(index %in% colnames(data))
    data_index <- data[[index]]
  } else {
    data_index <- index
  }

  if (is.character(group)) {
    stopifnot(group %in% colnames(data))
    data_group <- data[[group]]
  } else {
    data_group <- group
  }
  # ----------------------------------------------------------------------------

  # make sorting order and sort index and group accordingly
  data_order <- order(data_group, data_index)
  data_index <- data_index[data_order]
  data_group <- data_group[data_order]



  # Extract data for response formula ------------------------------------------
  data_y <- model.response(mf_response)
  data_X <- model.matrix(attr(mf_response, "terms"), mf_response)

  data_X <- apply(data_X, 2, .hmmcr_normalise_column)

  data_y <- data_y[data_order]
  data_X <- data_X[data_order, ]
  # ----------------------------------------------------------------------------



  # Extract data for state formula ----------------------------------------------------------------------------
  data_S <- .hmmcr_convert_to_state_matrix(model.response(mf_state))
  colnames(data_S) <- paste("S", 1:ncol(data_S))
  data_Z <- model.matrix(attr(mf_state, "terms"), mf_state)

  data_S <- as.matrix(data_S[data_order, ])
  data_Z <- as.matrix(data_Z[data_order, ])
  # ----------------------------------------------------------------------------



  # par and prior --------------------------------------------------------------
  n <- length(data_y)
  p <- dim(data_X)[2]
  q <- dim(data_Z)[2]
  nr_states <- dim(data_S)[2]

  par_tmp <- .hmmcr_make_par(nr_states = nr_states, col_z = dim(data_Z)[2], col_x = dim(data_X)[2], model = model)

  par_index <- par$index
  if (is.null(par$index)) {
    par_index <- par_tmp$par_index
  }

  par_init <- par$init
  if (is.null(par$init)) {
    par_init <- par_tmp$par_init
  }

  par_unique <- par$unique
  if (is.null(par_unique)) {
    par_unique <- par_tmp$par_unique
  }

  par_groups <- par$group
  if (is.null(par_groups)) {
    par_groups <- par_tmp$par_groups
  }

  # prior
  par_prior <- .hmmcr_make_prior_par(par_init, model = model)
  # ----------------------------------------------------------------------------
  # par_init <- par_init*0
  
  # print(par_init)
  # print(par_prior)
  # print(cbind(par_index, par_unique))
  
  # print(summary(data_y))
  # print(summary(data_X))
  # print(summary(data_Z))
  
  # model_int <- 0
  # if (model == "poisson") {model_int <- 1}
  # if (model == "gaussian") {model_int <- 2}
  
  
  
  # mcmc -----------------------------------------------------------------------
  mcmc_out = mcmc_routine(par_init = par_init,
                          par_index = par_index,
                          par_prior = par_prior,
                          par_unique = par_unique,
                          y = as.matrix(data_y),
                          X = as.matrix(data_X),
                          Z = as.matrix(data_Z),
                          S = as.matrix(data_S),
                          id = as.matrix(data_group),
                          model = model,
                          nr_states = nr_states,
                          groups = par_groups,
                          steps = steps,
                          burnin = burnin)
  # ----------------------------------------------------------------------------
  
  # summary and post-processing ------------------------------------------------
  
  
  # ----------------------------------------------------------------------------
  
}


