#' Fit the Bayesian ordinal diagnostic-trial model
#'
#' Fits the Bayesian ordinal diagnostic-trial model for ordinal rating data, with optional known disease status.
#' Posterior samples are obtained by MCMC.
#' The function returns chain-specific posterior samples for the cutpoints, rater-specific parameters, latent disease severities, and hyperparameters, along with the last state of each chain.
#'
#' @param w A numeric matrix of ordinal ratings. Observed entries must be consecutive integers `1, 2, ..., K`. Missing values (`NA`) are allowed.
#' @param D Optional binary vector of true disease status, coded as `0` and `1`. If `NULL`, disease status is treated as unknown.
#' @param inits Optional initial values. For a single chain, this should be a named list containing any subset of `z`, `z0`, `a`, `b`, `u`, and `A`.
#'   For multiple chains, this can be either a single named list used as a template for all chains, or an unnamed list of chain-specific named lists.
#' @param prior Optional named list of prior hyperparameters. Allowed elements are `mu_a`, `tau_a`, `mu_b`, `tau_b`, `delta`, `gamma`, `tau`, `mu_u`, and `tau_u`. Missing elements are filled with default values.
#' @param n.chains Number of MCMC chains.
#' @param n.iter Total number of MCMC iterations per chain.
#' @param n.burnin Number of burn-in iterations per chain.
#' @param n.thin Thinning interval.
#'
#' @return An object of class `"BayesODT"`, which is a list containing:
#' \describe{
#'   \item{chains}{A list of chain-specific posterior samples.}
#'   \item{D}{The supplied disease status vector.}
#'   \item{w}{The supplied rating matrix.}
#'   \item{mcmc_info}{A list containing MCMC settings and the number of saved draws.}
#'   \item{dim_info}{A list containing `m`, `n`, and `K`.}
#'   \item{param_names}{Parameter names used internally.}
#'   \item{call}{The matched function call.}
#' }
#'
#' Each element of `chains` contains posterior samples for `A`, `a`, `tau_a`, `b`, `mu_b`, `tau_b`, and `u`, together with `num` and `state_last`.
#' @details
#' The model is fit by Gibbs sampling with latent-variable augmentation.
#' When `D` is provided, the known disease-status version of the model is used.
#' When `D = NULL`, the unknown disease-status version is used.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' w <- matrix(sample(1:5, 50, replace = TRUE), nrow = 10, ncol = 5)
#' D <- sample(0:1, 10, replace = TRUE)
#'
#' fit <- BayesODT(w, D, n.chains = 2, n.iter = 2000, n.burnin = 1000, n.thin = 5)
#'
#' fit$mcmc_info
#' names(fit$chains[[1]])
#' }
#'
#' @export

BayesODT<- function(w, D = NULL,
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


  # observed rating positons by category
  W_list <- vector("list", K)
  for (k in seq_len(K)) {
    W <- which(w == k, arr.ind = TRUE)
    i <- W[, 1]
    j <- W[, 2]
    idx <- (j - 1L) * m + i
    W_list[[k]] <- list(i = i, j = j, idx = idx)
  }

  # missing positions
  W_miss <- which(is.na(w), arr.ind = TRUE)
  if (nrow(W_miss) > 0L) {
    i_miss <- W_miss[ ,1]
    j_miss <- W_miss[ ,2]
    idx_miss <- (j_miss - 1L) * m + i_miss # z의 선형 인덱스 (컬럼 메이저)
    W_miss <- list(i = i_miss, j = j_miss, idx = idx_miss)
  } else {
    W_miss <- NULL
  }

  # D가 0, 1일 때의 index 미리 만들기
  if (!is.null(D)) {
    idx1 <- which(D == 1L)
    idx0 <- which(D == 0L)
  } else {
    idx1 <- idx0 <- NULL
  }

  #### 체인당 rcpp MCMC 작동 ####
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
  OUT <- list(chains = fit_list, D = D, w = w)
  OUT$mcmc_info <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.save = fit_list[[1]]$num)
  OUT$dim_info <- list(m = m, n = n, K = K)
  OUT$param_names <- list(A = paste0("A[", seq_len(K - 1), "]"),
                          a = paste0("a[", seq_len(n), "]"),
                          tau_a = "tau_a",
                          b = paste0("b[", seq_len(n), "]"),
                          mu_b = "mu_b",
                          tau_b = "tau_b",
                          u = paste0("u[", seq_len(m), "]"))
  OUT$call <- match.call()
  class(OUT) <- "BayesODT"
  OUT

}
