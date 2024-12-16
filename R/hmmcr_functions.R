library(Rcpp, quietly = TRUE)
library(data.table, quietly = TRUE)

.hmmcr_convert_to_state_matrix <- function(s) {
  tmp <- function(s, max_length) {
    digits <- as.numeric(strsplit(s, "::")[[1]])
    result <- numeric(max_length)
    result[digits] <- digits
    return(result)
  }

  max_length <- length(unique(unlist(strsplit(s, "::"))))

  return(t(sapply(s, tmp, max_length = max_length)))
}

.hmmcr_make_par <- function(nr_states, col_z, col_x, model) {
  # init
  par_index <- rep(0, nr_states - 1)
  par_init <- rep(0.0, nr_states - 1)
  par_unique <- 1:(nr_states - 1)

  # zeta
  nr_zeta <- col_z * (nr_states) * (nr_states - 1)
  par_index <- c(par_index, rep(1, nr_zeta))
  par_init <- c(par_init, rep(0.0, nr_zeta))
  par_unique <- c(par_unique, max(par_unique) + 1:nr_zeta)

  # beta
  nr_beta <- col_x * nr_states
  par_index <- c(par_index, rep(2, nr_beta))
  par_init <- c(par_init, rep(0.0, nr_beta))
  par_unique <- c(par_unique, max(par_unique) + 1:nr_beta)

  if (model == "ingarch") {
    # beta
    nr_beta <- col_x * nr_states
    par_index <- c(par_index, rep(2, nr_beta))
    par_init <- c(par_init, rep(0, col_x * (nr_states - 1)))
    par_unique <- c(par_unique, rep(-1, col_x))
    par_unique <- c(par_unique, max(par_unique) + 1:(col_x * (nr_states - 1)))

    # theta
    nr_theta <- 2 * nr_states # nr_states + 1
    par_index <- c(par_index, rep(3, nr_theta))
    par_init <- c(par_init, rnorm(nr_states + 1, mean = 0, sd = 1))
    par_unique <- c(par_unique, max(par_unique) + c(rbind(1:nr_states, nr_states + 1)))
  }

  if (model == "poisson") {
    # no additional parameters
  }

  if (model == "gaussian" || model == "log-gaussian") {
    nr_theta <- nr_states
    par_index <- c(par_index, rep(3, nr_theta))
    par_init <- c(par_init, rep(1, nr_theta))
    par_unique <- c(par_unique, max(par_unique) + rbind(1:nr_states))
  }

  return(list(par_index = par_index, par_init = par_init, par_unique = par_unique))
}

.hmmcr_make_prior_par <- function(par_index) {
  return(c(rep(0, length(par_index)), rep(20, length(par_index))))
}

.hmmcr_normalise_column <- function(x) {
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- sd(x, na.rm = TRUE)

  if (x_sd != 0) {
    return((x - x_mean) / x_sd)
  } else {
    return(x)
  }
}

.hmmcr_check_par <- function(par_init,
                             par_index,
                             par_unique,
                             par_prior,
                             nr_states,
                             col_z,
                             col_x,
                             model,
                             debug) {
  col_theta <- 1 + (model == "ingarch")

  tmp_init <- rep(NA, nr_states - 1)
  tmp_zeta <- matrix(NA, nrow = nr_states * (nr_states - 1), ncol = col_z)
  tmp_beta <- matrix(NA, nrow = nr_states, ncol = col_x)
  tmp_theta <- matrix(NA, nrow = nr_states, ncol = col_theta)

  tmp_par_index <- par_index[(par_unique > 0) & !duplicated(par_unique)]
  tmp_par_index[length(tmp_par_index)] <- 4

  k1 <- 1
  k2 <- 1
  k3 <- 1
  l2 <- 1
  l3 <- 1
  k4 <- 1
  l4 <- 1

  for (i in seq_along(par_unique)) {
    par_index_i <- par_index[i]
    par_unique_i <- par_unique[i]

    if (par_index_i == 0) {
      tmp_init[k1] <- ifelse(par_unique_i > 0, par_init[par_unique_i], 0)
      k1 <- k1 + 1
    }

    if (par_index_i == 1) {
      tmp_zeta[k2, l2] <- ifelse(par_unique_i > 0, par_init[par_unique_i], 0)
      l2 <- l2 + 1
      if (l2 > col_z) {
        l2 <- 1
        k2 <- k2 + 1
      }
    }

    if (par_index_i == 2) {
      tmp_beta[k3, l3] <- ifelse(par_unique_i > 0, par_init[par_unique_i], 0)
      l3 <- l3 + 1
      if (l3 > col_x) {
        l3 <- 1
        k3 <- k3 + 1
      }
    }

    if (par_index_i == 3) {
      tmp_theta[k4, l4] <- ifelse(par_unique_i > 0, par_init[par_unique_i], 0)
      l4 <- l4 + 1
      if (l4 > col_theta) {
        l4 <- 1
        k4 <- k4 + 1
      }
    }
  }

  if (debug) {
    cat("\nModel:", paste0(toupper(substr(model, 1, 1)), substr(model, 2, nchar(model))), "\n")

    cat("\nInit vector: \n ")
    tmp_init <- c(1, exp(tmp_init))
    print(tmp_init / sum(tmp_init))

    cat(paste0("\nZeta Matrix (", nrow(tmp_zeta), ", ", ncol(tmp_zeta), "): \n"))
    print(tmp_zeta)

    cat(paste0("\nBeta Matrix (", nrow(tmp_beta), ", ", ncol(tmp_beta), "): \n"))
    print(tmp_beta)

    cat(paste0("\nTheta Matrix (", nrow(tmp_theta), ", ", ncol(tmp_theta), "): \n"))
    print(tmp_theta)
  }

  if (anyNA(c(1, exp(tmp_init)))) {
    stop("Error: Initial parameters for the initial probabilities are not properly defined.")
  }

  if (anyNA(tmp_zeta)) {
    stop("Error: Initial parameters for the Zeta matrix are not properly defined.")
  }

  if (anyNA(tmp_beta)) {
    stop("Error: Initial parameters for the Beta matrix are not properly defined.")
  }

  if (model != "poisson") {
    if (anyNA(tmp_theta)) {
      stop("Error: Initial parameters for the Theta matrix are not properly defined.")
    }
  }
}


if (FALSE) {
  n <- 100
  m <- 3
  formula_response <- y ~ x1 * x2 + x3
  formula_state <- state ~ x1 + x2
  y <- rbinom(n * m, 0, 1)
  state <- rep("1::2::3", n * m)
  x1 <- rnorm(n * m)
  x2 <- rnorm(n * m)
  x3 <- rnorm(n * m)
  z <- rep(1, n)
  index <- 1:n
  group <- rep("4", n)
  burnin <- 100
  steps <- 1000
  model <- "poisson"
  normalize <- TRUE
}

hmmcr <- function(formula_response,
                  formula_state,
                  index = NULL,
                  group = NULL,
                  model = "poissin",
                  burnin = 1000,
                  steps = 1000,
                  par = list(init = NULL, index = NULL, unique = NULL, prior = NULL),
                  normalize = TRUE,
                  data = NULL,
                  debug = FALSE) {
  # TODO:
  # 3) make a prior class and in addition to model add prior as to select type
  #    of prior
  # 4) add model to Para or Prior class and add para conditions here.
  # 6) change to model = "normal-ingarch" to "normal::ingarch::3" where the
  #    last number is the size of rolling window, then we need a new hyper_par
  #    and hyper_par_index thing for this
  # 8) rename log_fn_fun and log_fn_i_fun to log_f and log_f_i or log_density
  # 10) make par_unique start at 0


  # check arguments ------------------------------------------------------------
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("Error: 'data' must be a data.frame.")
    }
  }

  if (is.null(index)) {
    data_index <- NULL
  } else {
    if (length(index) == 1 && all(is.character(index))) {
      stopifnot(index %in% colnames(data))
      data_index <- data[[index]]
    } else {
      data_index <- index
    }
  }

  if (is.null(group)) {
    data_group <- NULL
  } else {
    if (length(group) == 1 && is.character(group)) {
      stopifnot(group %in% colnames(data))
      data_group <- data[[group]]
    } else {
      data_group <- group
    }
  }

  valid_models <- c("ingarch", "poisson", "gaussian", "log-gaussian", "nbinomial")
  if (!model %in% valid_models) {
    stop(
      sprintf(
        "Invalid model '%s'. Valid models are: %s",
        model,
        paste(valid_models, collapse = ", ")
      )
    )
  }

  if (!(is.finite(burnin) && burnin == floor(burnin) && burnin >= 1)) {
    stop("Error: 'burnin' must be a positive integer.")
  }

  if (!(is.finite(steps) && steps == floor(steps) && steps >= 1)) {
    stop("Error: 'steps' must be a positive integer.")
  }

  if (!is.logical(normalize)) {
    stop("Error: 'normalize' must be TRUE/FALSE.")
  }
  # ----------------------------------------------------------------------------


  # data from formula ----------------------------------------------------------
  mf_response <- model.frame(formula_response, data, na.action = na.fail)
  mf_state <- model.frame(formula_state, data, na.action = na.fail)

  # response
  data_y <- model.response(mf_response)
  data_x <- model.matrix(attr(mf_response, "terms"), mf_response)
  data_x <- apply(data_x, 2, .hmmcr_normalise_column)

  # state
  data_s <- .hmmcr_convert_to_state_matrix(model.response(mf_state))
  colnames(data_s) <- paste("S", seq_len(ncol(data_s)))
  data_z <- model.matrix(attr(mf_state, "terms"), mf_state)
  # ----------------------------------------------------------------------------


  # check index and group  -----------------------------------------------------
  if (!all(is.numeric(data_y))) {
    stop(paste("Error: The response", names(mf_response)[1], "variable must be numeric."))
  }

  if (nrow(data_z) != nrow(data_x)) {
    stop("Error: The number of samples in the state and response formulas must be identical.")
  }

  if (is.null(data_index)) {
    data_index <- seq_len(length(data_y))
  }

  if (is.null(data_group)) {
    data_group <- rep(1, length(data_y))
  }

  if (!all(is.numeric(data_group))) {
    data_group <- as.numeric(factor(data_group))
    warning("Warning: 'group' is converted to integers.")
  } else {
    if (!all(data_group == floor(data_group))) {
      data_group <- as.numeric(factor(data_group))
      warning("Warning: 'group' is converted to integers.")
    }
  }
  # ----------------------------------------------------------------------------


  # make sorting order and sort index and group accordingly
  data_order <- order(data_group, data_index)
  data_index <- data_index[data_order]
  data_group <- data_group[data_order]

  # check that the difference between each index for each group is the same
  data_group_index <- split(data_index, data_group)
  for (g in names(data_group_index)) {
    if (length(unique(diff(data_group_index[[g]]))) != 1) {
      stop("Error: Differences in 'index' are not constant for each group.")
    }
  }

  # extract data for response formula
  data_y <- data_y[data_order]
  data_x <- data_x[data_order, ]

  # extract data for state formula
  data_s <- as.matrix(data_s[data_order, ], drop = FALSE)
  data_z <- as.matrix(data_z[data_order, ])
  # ----------------------------------------------------------------------------



  # construct par and prior ----------------------------------------------------
  # browser()
  nr_states <- sum(unique(c(data_s)) > 0)

  par_tmp <- .hmmcr_make_par(
    nr_states = nr_states,
    col_z = dim(data_z)[2],
    col_x = dim(data_x)[2],
    model = model
  )

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
    par_index_unique <- par_index[(par_unique > 0) & !duplicated(par_unique)]
    par_groups <- list(which(par_index_unique != 1) - 1, c(which(par_index_unique == 1)) - 1)
  }

  # construct prior
  par_prior <- .hmmcr_make_prior_par(par_init)

  # check of par arguments
  if (length(par_index) != length(par_unique)) {
    stop("Error: The length of 'par_index' must be the same as length of 'par_unique'.")
  }

  if (2 * length(par_init) != length(par_prior)) {
    stop("Error: The length of 'par_prior' must be 2 times the length of 'par_init'.")
  }

  if (!all(is.numeric(par_init))) {
    stop("Error: 'par_init' must be numeric.")
  }

  if (!all(is.numeric(par_prior))) {
    stop("Error: 'par_prior' must be numeric.")
  }

  if (!(is.numeric(par_index) && all(par_index == floor(par_index)) && all(par_index > -1))) {
    stop("Error: Elements of 'par_index' must all positive integers.")
  }

  if (!(is.numeric(par_unique) && all(par_unique == floor(par_unique)) && all(par_unique > -2))) {
    stop("Error: Elements of 'par_unique' must either be a positive integer or -1.")
  }
  # ----------------------------------------------------------------------------

  # Debug
  .hmmcr_check_par(par_init,
    par_index,
    par_unique,
    par_prior,
    nr_states,
    col_z = dim(data_z)[2],
    col_x = dim(data_x)[2],
    model,
    debug
  )


  cat("\nMCMC Routine - Start \n")

  # the mcmc routine -----------------------------------------------------------
  mcmc_out <- mcmc_routine(
    par_init = par_init,
    par_index = par_index,
    par_prior = par_prior,
    par_unique = par_unique,
    y = data_y,
    X = data_x,
    Z = data_z,
    S = data_s,
    id = data_group,
    model = model,
    nr_states = nr_states,
    groups = par_groups,
    steps = steps,
    burnin = burnin
  )
  cat("\nMCMC Routine - End \n")
  # ----------------------------------------------------------------------------


  # summary and post-processing ------------------------------------------------
  # ----------------------------------------------------------------------------

  return(list(output = mcmc_out))
}
