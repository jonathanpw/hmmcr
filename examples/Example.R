library(hmmcr)

seed = 1
set.seed(seed)


# simulate and plot data
n <- 900
x <- rnorm(n)

beta <- matrix( c(-10, .1,
                  1,  .1,
                  2,  .2), ncol=2, byrow=T)


zeta <- c(-3, -2, -1, -1, -4, -2)

q <- zeta
Q <- matrix(c(1,        exp(q[1]), exp(q[2]),
              exp(q[3]),         1, exp(q[4]),
              exp(q[5]), exp(q[6]),         1), ncol=3, byrow=T)
P <- Q / rowSums(Q)


s <- rep( NA, n); s[1] <- sample( x=1:3, size=1, prob=c(.8, .15, .05))
y <- rep( NA, n); y[1] <- rpois( n=1, lambda=exp(c(1,x[1])%*%beta[s[1],]))

for( i in 2:n){
  s[i] <- sample( x=1:3, size=1, prob=P[s[i-1],])
  y[i] <- rpois( n=1, lambda=exp(c(1,x[i])%*%beta[s[i],]))
}

plot(y)

# setup and run hmmcr
data_1 <- data.frame(
  y = y,
  state = "1::2::3",
  x = x,
  z = 1,
  index = 1:n,
  group = 1
)

burnin <- 5000
steps <- 50000
n_post <- 10000

mcmc_out <- hmmcr(formula_response = y ~ x,
                  formula_state = state ~ -1 + z,
                  index = "index",
                  group = "group",
                  burnin = burnin,
                  steps = steps,
                  data = data_1,
                  model = "poisson")

index_post <- (steps - n_post + 1):steps

init <- c(1, exp(colMeans(mcmc_out$output[[1]][index_post, 1:2])))
print(init/sum(init))

zeta <- colMeans(mcmc_out$output[[1]][index_post, 3:8])
print(zeta)

beta <- colMeans(mcmc_out$output[[1]][index_post, -(1:8)])
print(beta)

save_object <- mcmc_out$output
save(save_object, file = paste0("mcmc_pars_", seed, ".rda"))
trace_plot(mcmc_out$output[[1]], n_post, steps, burnin)

print("See the generated file trace_plot.pdf in the working directory")


summary.hmmcr <- function(object, ...) {
  if (!inherits(object, "hmmcr")) {
    stop("The object is not of class 'hmmcr'")
  }

  index_post <- (object$steps - n_post + 1):object$steps

  init <- c(1, exp(colMeans(object$output[[1]][index_post, 1:2])))
  zeta <- colMeans(object$output[[1]][index_post, 3:8])
  beta <- colMeans(object$output[[1]][index_post, -(1:8)])

  # Organize details of summary
  summary_list <- list(
    init, zeta, beta)

  class(summary_list) <- "summary.hmmcr"
  return(summary_list)
}

# Print function for hmmcr objects
print.hmmcr <- function(x, ...) {
  if (!inherits(x, "hmmcr")) {
    stop("The object is not of class 'hmmcr'")
  }

  cat("Hidden Markov Model - MCMC Results:\n\n")
  cat("Output Summary:\n")

  # Use the summary method for concise details
  summary_out <- summary(x)

  for (param in names(summary_out$parameters)) {
    cat("\nParameter:", param, "\n")
    param_summary <- summary_out$parameters[[param]]
    cat("Mean:", param_summary$mean, "\n")
    cat("SD:", param_summary$sd, "\n")
    cat("Median:", param_summary$median, "\n")
    cat("2.5% Quantile:", param_summary$q2.5, "\n")
    cat("97.5% Quantile:", param_summary$q97.5, "\n")
  }

  cat("\nMCMC Steps:", summary_out$total_steps, "\nBurn-In Period:", summary_out$burnin, "\n")
}

# You can add these functions to your R environment for use with objects returned by the hmmcr function
