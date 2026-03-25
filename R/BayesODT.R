#' Fit a Bayesian Hierarchical ordinal diagnostic test model
#'
#' Fit the Bayesian Hierarchical ordinal diagnostic test model for multirater ordinal data using MCMC.
#' The model is a two-stage latent variable model for ordinal ratings and allows the true binary disease status to be either observed or unobserved.
#'
#' @param w An \eqn{m \times n} matrix of ordinal ratings with values in \code{1, \dots, K}. Missing ratings should be coded as \code{NA}.
#' @param D An optional binary vector of length \code{m}. If supplied, it represents the true disease status for each subject; if \code{NULL}, the disease status is treated as unknown.
#' @param inits An optional list of initial values. If \code{NULL}, default initial values are generated internally.
#' @param prior An optional list of prior hyperparameters. If \code{NULL}, default prior values are used.
#' @param n.chains Number of MCMC chains.
#' @param n.iter Total number of MCMC iterations per chain.
#' @param n.burnin Number of burn-in iterations per chain.
#' @param n.thin Thinning interval.
#'
#' @details
#' Let \eqn{W_{ij}} denote the ordinal rating assigned by rater \eqn{j} to subject \eqn{i}, and let \eqn{u_i} denote a latent disease severity.
#' The model links disease status and ordinal ratings through a probit formulation, with rater-specific bias parameters \eqn{a_j}, magnifier parameters \eqn{b_j}, and common cutpoints.
#' See Kim, Lin, and Nelson (2021) for model details. Missing ratings are allowed.
#'
#' @return An object of class \code{"BayesODT"} containing posterior samples from each chain, the input data, MCMC settings, dimension information, parameter labels, and the matched function call.
#'
#' @references
#' Kim C, Lin X, Nelson KP (2021).
#' \emph{Measuring rater bias in diagnostic tests with ordinal ratings}.
#' \emph{Statistics in Medicine}, \strong{40}(18), 4014--4033.
#' doi:10.1002/sim.9011.
#'
#' @examples
#' set.seed(123)
#' w <- matrix(sample(1:3, 24, replace = TRUE), nrow = 8, ncol = 3)
#' fit <- BayesODT(w, n.chains = 1, n.iter = 20, n.burnin = 10, n.thin = 2)
#' fit
#'
#' @export

BayesODT <- function(w, D = NULL,
                     inits = NULL, prior = NULL,
                     n.chains = 1, n.iter = 4000, n.burnin = 2000, n.thin = 5) {
  ####### argument Validation ######
  .validate_mcmc(n.iter, n.burnin, n.thin)
  .validate_w(w)
  .validate_D(D, w)
  .validate_init_chains(inits, w, n.chains)
  .validate_prior(prior)

  ######## fill default value in init & prior
  init <- .fill_init(inits, w, D, n.chains)
  prior <- .fill_prior(prior)

  ####### Preprocessing ########
  K <- length(unique(w[!is.na(w)]))
  m <- nrow(w)
  n <- ncol(w)

  mu_a <- prior$mu_a
  tau_a <- prior$tau_a
  mu_b <- prior$mu_b
  tau_b <- prior$tau_b
  delta <- prior$delta
  gamma <- prior$gamma
  tau <- prior$tau
  mu_u <- prior$mu_u
  tau_u <- prior$tau_u

  # observed rating positions by category
  W_list <- vector("list", K)
  for (k in seq_len(K)) {
    W <- which(w == k, arr.ind = TRUE)
    i <- W[, 1]
    j <- W[, 2]
    idx <- (j - 1L) * m + i # linear index for faster computation
    W_list[[k]] <- list(i = i, j = j, idx = idx)
  }

  # missing rating positions
  W_miss <- which(is.na(w), arr.ind = TRUE)
  if (nrow(W_miss) > 0L) {
    i_miss <- W_miss[ ,1]
    j_miss <- W_miss[ ,2]
    idx_miss <- (j_miss - 1L) * m + i_miss
    W_miss <- list(i = i_miss, j = j_miss, idx = idx_miss)
  } else {
    W_miss <- NULL
  }

  # Precompute indices for D = 0 and D = 1
  if (!is.null(D)) {
    idx1 <- which(D == 1L)
    idx0 <- which(D == 0L)
  } else {
    idx1 <- idx0 <- NULL
  }

  # Run the Rcpp MCMC for each chain
  fit_list <- vector("list", n.chains)

  for (i in seq_len(n.chains)) {
    a <- init[[i]]$a
    b <- init[[i]]$b
    u <- init[[i]]$u
    z <- init[[i]]$z
    z0 <- init[[i]]$z0
    A <- c(init[[i]]$A, Inf)

    state <- list(z = z, z0 = z0, u = u, a = a, b = b, A = A,
                  tau_a = tau_a, mu_b = mu_b, tau_b = tau_b)
    fit_list[[i]] <- MCMC_loop(state = state,
                               D = D,
                               mu_a = mu_a, tau_a = tau_a, mu_b = mu_b, tau_b = tau_b, mu_u = mu_u, tau_u = tau_u, gamma = gamma, delta = delta, tau = tau,
                               K = K,
                               W_list = W_list, W_miss = W_miss, idx1 = idx1, idx0 = idx0,
                               n_iter = n.iter, n_burnin = n.burnin, n_thin = n.thin)
  }
  out <- list(chains = fit_list, D = D, w = w)
  out$mcmc_info <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.save = fit_list[[1]]$num)
  out$dim_info <- list(m = m, n = n, K = K)
  out$param_names <- list(A = paste0("A[", seq_len(K - 1), "]"),
                          a = paste0("a[", seq_len(n), "]"),
                          tau_a = "tau_a",
                          b = paste0("b[", seq_len(n), "]"),
                          mu_b = "mu_b",
                          tau_b = "tau_b",
                          u = paste0("u[", seq_len(m), "]"))
  out$call <- match.call() # Store the function call
  class(out) <- "BayesODT"
  out

}
