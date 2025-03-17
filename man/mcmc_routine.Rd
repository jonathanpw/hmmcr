\name{mcmc_routine}
\alias{mcmc_routine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The mcmc routine for sampling the parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mcmc_routine( data, init_par, prior_par, par_index, steps, burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A list containing the following components:}

  \item{y}{A vector of fatality counts corresponding to each time index and each
  country.}

  \item{X}{A matrix with the same number of rows as columns of the y argument,
  where the columns correspond to explanatory variables included in the
  autoregressive response model.}

  \item{Z}{A matrix with the same number of rows as columns of the y argument,
  where the columns correspond to explanatory variables included in the
  transition probability matrix.}

  \item{S}{A matrix with the same number of rows as columns of the y argument,
  and number of columns corresponding to the number of states to be specified in
  the HMM.  Each row, corresponding to some time index of some country,
  indicates the possible states for that time index and country with the value 1
  in the columns assocaited with the possible states; value 0 otherwise.}

  \item{id}{A vector the same length as the y argument, each component
  corresponding to a country index value.}

  \item{model}{Only model 0 is available for now.}

  \item{par_index}{A vector of indicies for each of the HMM parameters.
  Indices must be specified in the order: ``pi_logit'' (initial state probability
  parameters), ``zeta'' (transition rate parameters), ``beta'' (response
  function coefficient parameters), ``log_a_vec'' (baseline response function
  parameters), and ``log_c'' (negative binomial probability parameter).}

  \item{init_par}{A vector specifying the initial values for each HMM parameter,
  where the order of the parameters is as specified in the par_index argument.}

  \item{prior_par}{A vector twice the length of par_index the, the first half of
  the components specifying the prior mean for each HMM parameter, and the
  second half of the components specifying the prior standard deviation for each
  HMM parameter.  The order of the parameters is as specified in the par_index
  argument.}

  \item{groups}{A list, each component of which is a vector of indices to be
  updated together in a Metropolis-within-Gibbs fashion.  The number of list
  components determines the number of Metropolis-within-Gibbs steps within each
  step of the MCMC algorithm.}

  \item{steps}{Number of steps to run the MCMC algorithm for.}

  \item{burnin}{Number of steps to adaptively learn the covariance matrix used
  for the Gaussian random walk proposals.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Example 1 --------------------------------------------------------------------


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
y <- rep( NA, n); y[1] <- rpois( n=1, lambda=exp(c(1,x[1])\%*\%beta[s[1],]))

for( i in 2:n){
  s[i] <- sample( x=1:3, size=1, prob=P[s[i-1],])
  y[i] <- rpois( n=1, lambda=exp(c(1,x[i])\%*\%beta[s[i],]))
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

index_post <- (steps - burnin - n_post + 1):(steps - burnin)

init <- c(1, exp(colMeans(mcmc_out$output[[1]][index_post, 1:2])))
print(init/sum(init))

zeta <- colMeans(mcmc_out$output[[1]][index_post, 3:8])
print(zeta)

beta <- colMeans(mcmc_out$output[[1]][index_post, -(1:8)])
print(beta)

save(mcmc_out$output, file = paste0("mcmc_pars_", seed, ".rda"))
trace_plot(mcmc_out$output[[1]], n_post, steps, burnin)

print("See the generated file trace_plot.pdf in the working directory")

# ------------------------------------------------------------------------------

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
