#' Simulate multirater ordinal diagnostic test data
#'
#' Generate one simulated dataset from a multirater ordinal diagnostic model.
#' The function returns an ordinal rating matrix, the true binary disease status,
#' and the true AUC computed from the generated data.
#'
#' @param m Integer. Number of items (subjects).
#' @param n Integer. Number of raters.
#' @param K Integer. Number of ordinal categories.
#' @param alpha Numeric vector of length \code{K - 1}. Ordered cutpoints for the
#'   ordinal rating model.
#' @param mu_a Numeric. Mean of the rater bias parameters \eqn{a_j}.
#' @param var_a Numeric. Variance of the rater bias parameters \eqn{a_j}.
#' @param mu_b Numeric. Mean of the rater magnifier parameters \eqn{b_j}.
#' @param var_b Numeric. Variance of the rater magnifier parameters \eqn{b_j}.
#' @param u_mean1 Numeric. Mean of the first mixture component for the latent
#'   disease severity \eqn{u_i}.
#' @param u_sd1 Numeric. Standard deviation of the first mixture component for
#'   \eqn{u_i}.
#' @param u_mean2 Numeric. Mean of the second mixture component for the latent
#'   disease severity \eqn{u_i}.
#' @param u_sd2 Numeric. Standard deviation of the second mixture component for
#'   \eqn{u_i}.
#' @param u_weight Numeric between 0 and 1. Mixing weight for the first mixture
#'   component of \eqn{u_i}.
#'
#' @details
#' The latent disease severity \eqn{u_i} is generated from a two-component
#' normal mixture distribution. The true disease status \eqn{D_i} is then
#' sampled from a Bernoulli distribution with success probability
#' \eqn{\Phi(u_i)}, where \eqn{\Phi} is the standard normal cumulative
#' distribution function.
#'
#' For each rater \eqn{j} and item \eqn{i}, the ordinal rating \eqn{W_{ij}} is
#' generated from category probabilities
#' \deqn{
#' P(W_{ij} = k) =
#' \Phi(\alpha_k - (a_j + b_j u_i)) -
#' \Phi(\alpha_{k-1} - (a_j + b_j u_i)),
#' }
#' for \eqn{k = 1, \dots, K}.
#'
#' The returned \code{true_AUC} is computed from the generated diseased and
#' nondiseased rating distributions.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{W}{An \eqn{m \times n} matrix of ordinal ratings.}
#'   \item{D}{A binary vector of length \eqn{m} containing the true disease status.}
#'   \item{true_AUC}{The true AUC computed from the generated data.}
#' }
#'
#' @examples
#' sim <- simulate_data()
#' str(sim$W)
#' table(sim$D)
#' sim$true_AUC
#'
#' @export
simulate_data <- function(m = 90, n = 50, K = 4,
                          alpha = c(-1, 1, 2),
                          mu_a = 0, var_a = 0.5,
                          mu_b = 1.7, var_b = 0.5,
                          u_mean1 = -1.5, u_sd1 = 0.5,
                          u_mean2 = 1.5, u_sd2 = 0.5,
                          u_weight = 0.9) {
  # ------- input check -----------
  if (!is.numeric(m) || length(m) != 1L || m < 1) stop("`m` must be a positive integer.", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || n < 1) stop("`n` must be a positive integer.", call. = FALSE)
  if (!is.numeric(K) || length(K) != 1L || K < 3) stop("`K` must be an integer >= 3.", call. = FALSE)
  if (length(alpha) != (K - 1L)) stop("`cutpoints` must have length `K - 1`.", call. = FALSE)
  if (any(diff(alpha) <= 0)) stop("`cutpoints` must be strictly increasing.", call. = FALSE)
  if (u_weight <= 0 || u_weight >= 1) stop("`u_weight` must be between 0 and 1.", call. = FALSE)
  if (var_a <= 0 || var_b <= 0 || u_sd1 <= 0 || u_sd2 <= 0) stop("Variances/SDs must be positive.", call. = FALSE)

  m <- as.integer(m)
  n <- as.integer(n)
  K <- as.integer(K)

  alpha_vec <- c(-Inf, alpha, Inf)

  aj_vec <- rnorm(n, mean = mu_a, sd = sqrt(var_a))
  aj_vec <- aj_vec - mean(aj_vec)
  bj_vec <- rnorm(n, mean = mu_b, sd = sqrt(var_b))

  # generate latent disease severity u_i
  u_true = rep(NA, m)
  D_true = rep(NA, m)

  for (i in 1:m) {
    t1 = runif(1)
    if (t1 < u_weight) {
      u_true[i] = rnorm(1, u_mean1, u_sd1) # d.true[i] = 0
    } else {
      u_true[i] = rnorm(1, u_mean2, u_sd2) # d.true[i] = 1
    }
  }

  ui_vec <- u_true
  u_prob_vec = pnorm(ui_vec)

  for (i in 1:m) {
    D_true[i] = rbinom(1,1,u_prob_vec[i])
  }
  Di_vec <- D_true

  prob_vec = rep(0, K)
  wij_vec  = NULL
  for (j in 1:n) {
    for (i in 1:m) {
      prob_vec = rep(0, K)
      for (k in 1:K) {
        prob_vec[k] = pnorm(alpha_vec[k+1]- (aj_vec[j] + bj_vec[j]*ui_vec[i])) - pnorm(alpha_vec[k]  - (aj_vec[j] + bj_vec[j]*ui_vec[i]))
      }
      value = rmultinom(1, size = 1, prob=prob_vec)
      for (q in 1:K) {
        if (value[q,1] == 1) {
          wij = q
        }
      }
      wij_vec = c(wij_vec, wij)
    }
  }

  cattable=table(wij_vec)
  round(prop.table(cattable)*50,0)

  fullrater = rep(seq(1,n), m)
  fullrater = sort(fullrater)
  fullitem  = rep(seq(1,m), n)
  fulldivec = rep(Di_vec, n)


  # Simulated Binary Dataset contains rater, item and classification and true disease status

  simdataset=data.frame(fullrater, fullitem, wij_vec, fulldivec)
  colnames(simdataset) = c("rater","item","wij","di")

  y0.table <- prop.table(table(factor(subset(simdataset$wij, simdataset$di==0), levels = 1:K)))
  y1.table <- prop.table(table(factor(subset(simdataset$wij, simdataset$di==1), levels = 1:K)))

  true_AUC <- 0
  for (k in 1:(K-1)) {
    true_AUC <- true_AUC + y0.table[k] * sum(y1.table[(k+1):K])
  }
  true_AUC <- true_AUC + 0.5 * sum(y0.table * y1.table)
  true_AUC <- unname(true_AUC)

  D <- simdataset$di[1:m]


  W <- matrix(ncol=n,nrow=m)
  ii <- 0
  for(i in unique(simdataset$item)){
    jj <- 0
    ii <- ii + 1
    for(j in unique(simdataset$rater)){
      jj <- jj + 1
      W[ii,jj] <- simdataset$wij[which(simdataset$item==i & simdataset$rater==j)]
    }
  }

  list(W = W,
       D = D,
       true_AUC = true_AUC)
}


#' Goodness-of-Fit Test
#'
#' Check goodness of fit of the fitted Bayesian hierarchical ordinal diagnostic
#' test model using posterior predictive checking.
#'
#' This function computes posterior predictive probabilities for each ordinal
#' rating category and evaluates a posterior predictive P-value based on the
#' chi-square-type discrepancy statistic described in Kim, Lin, and Nelson (2021).
#'
#' @param fit A fitted model object of class \code{"BayesODT"} returned by
#'   \code{\link{BayesODT}}.
#' @param seed Optional integer seed for generating replicated data.
#'
#' @return A list of class \code{"GoF.BayesODT"} containing:
#' \describe{
#'   \item{\code{Estimated.Prob}}{Posterior mean of model-based probabilities
#'   for each ordinal rating category.}
#'   \item{\code{Estimated.Prob.ci}}{Pointwise 95\% posterior credible intervals
#'   for model-based probabilities.}
#'   \item{\code{Empirical.Prob}}{Empirical proportions of observed ordinal
#'   rating categories. Missing ratings are excluded.}
#'   \item{\code{Observed.Count}}{Observed counts for each ordinal rating category.}
#'   \item{\code{Expected.Count}}{Posterior mean expected counts for each ordinal
#'   rating category.}
#'   \item{\code{Observed.T}}{Discrepancy statistics computed from the observed data.}
#'   \item{\code{Replicated.T}}{Discrepancy statistics computed from posterior
#'   predictive replicated data.}
#'   \item{\code{Posterior.P.value}}{Posterior predictive P-value.}
#' }
#'
#' @details
#' For each posterior draw \eqn{\theta^{(r)}}, the function computes
#' \eqn{p_k^{(r)} = P(W_{ij}=k \mid \theta^{(r)})}, generates one replicated
#' dataset \eqn{W^{rep,r}}, and compares
#' \deqn{
#' T(W^{rep,r}, \theta^{(r)}) \ge T(W, \theta^{(r)}),
#' }
#' where
#' \deqn{
#' T(W, \theta) =
#' \sum_{k=1}^{K}
#' \frac{(N_k - p_k N)^2}{p_k N}.
#' }
#'
#' The posterior predictive P-value is estimated by the proportion of posterior
#' draws for which the replicated discrepancy statistic is greater than or equal
#' to the observed discrepancy statistic.
#'
#' @examples
#' set.seed(123)
#' w <- matrix(sample(1:3, 24, replace = TRUE), nrow = 8, ncol = 3)
#' fit <- BayesODT(w, n.chains = 1, n.iter = 20, n.burnin = 10, n.thin = 2)
#' GoF(fit)
#'
#' @export
GoF <- function(fit, seed = NULL) {

  if (!inherits(fit, "BayesODT")) {
    stop("'fit' must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (is.null(fit$chains) || length(fit$chains) == 0L) {
    stop("'fit$chains' is missing or empty.", call. = FALSE)
  }

  if (is.null(fit$w)) {
    stop("'fit$w' is missing.", call. = FALSE)
  }

  w <- fit$w

  if (!is.matrix(w)) {
    stop("'fit$w' must be a matrix.", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("'seed' must be a single numeric value.", call. = FALSE)
    }

    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }

    set.seed(seed)

    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  if (!is.null(fit$dim_info)) {
    m <- fit$dim_info$m
    n <- fit$dim_info$n
    K <- fit$dim_info$K
  } else {
    m <- nrow(w)
    n <- ncol(w)
    K <- max(w, na.rm = TRUE)
  }

  obs_pos <- which(!is.na(w), arr.ind = TRUE)

  if (nrow(obs_pos) == 0L) {
    stop("'fit$w' contains no observed ratings.", call. = FALSE)
  }

  i_obs <- obs_pos[, 1L]
  j_obs <- obs_pos[, 2L]

  w_obs <- as.integer(w[!is.na(w)])
  N <- length(w_obs)

  observed_count <- tabulate(w_obs, nbins = K)
  empirical_prob <- observed_count / N

  prob_list <- vector("list", length(fit$chains))
  T_obs_list <- vector("list", length(fit$chains))
  T_rep_list <- vector("list", length(fit$chains))

  for (ch in seq_along(fit$chains)) {

    chain <- fit$chains[[ch]]

    if (is.null(chain$A) || is.null(chain$a) ||
        is.null(chain$b) || is.null(chain$u)) {
      stop(
        "Each chain must contain posterior samples named 'A', 'a', 'b', and 'u'.",
        call. = FALSE
      )
    }

    A_mat <- as.matrix(chain$A)
    a_mat <- as.matrix(chain$a)
    b_mat <- as.matrix(chain$b)
    u_mat <- as.matrix(chain$u)

    n_save <- nrow(a_mat)

    if (nrow(A_mat) != n_save ||
        nrow(b_mat) != n_save ||
        nrow(u_mat) != n_save) {
      stop(
        "Posterior sample sizes are inconsistent within chain ",
        ch,
        ".",
        call. = FALSE
      )
    }

    if (ncol(a_mat) != n) {
      stop(
        "The number of columns in 'a' does not match the number of raters.",
        call. = FALSE
      )
    }

    if (ncol(b_mat) != n) {
      stop(
        "The number of columns in 'b' does not match the number of raters.",
        call. = FALSE
      )
    }

    if (ncol(u_mat) != m) {
      stop(
        "The number of columns in 'u' does not match the number of subjects.",
        call. = FALSE
      )
    }

    if (ncol(A_mat) == K - 1L) {
      alpha_mat <- cbind(A_mat, Inf)
    } else if (ncol(A_mat) == K) {
      alpha_mat <- A_mat
      alpha_mat[, K] <- Inf
    } else {
      stop(
        "The number of columns in 'A' must be either K - 1 or K.",
        call. = FALSE
      )
    }

    prob_mat <- matrix(NA_real_, nrow = n_save, ncol = K)
    T_obs <- numeric(n_save)
    T_rep <- numeric(n_save)

    for (s in seq_len(n_save)) {

      a <- a_mat[s, ]
      b <- b_mat[s, ]
      u <- u_mat[s, ]
      alpha <- alpha_mat[s, ]

      eta_obs <- a[j_obs] + b[j_obs] * u[i_obs]

      cum_prob <- matrix(NA_real_, nrow = N, ncol = K)

      for (k in seq_len(K)) {
        cum_prob[, k] <- stats::pnorm(alpha[k] - eta_obs)
      }

      cell_prob <- cum_prob
      cell_prob[, 1L] <- cum_prob[, 1L]

      if (K >= 2L) {
        for (k in 2L:K) {
          cell_prob[, k] <- cum_prob[, k] - cum_prob[, k - 1L]
        }
      }

      cell_prob[cell_prob < 0] <- 0
      cell_prob[cell_prob > 1] <- 1

      row_sum <- rowSums(cell_prob)

      if (any(row_sum <= 0)) {
        stop(
          "Invalid posterior predictive probabilities were produced.",
          call. = FALSE
        )
      }

      cell_prob <- cell_prob / row_sum

      p_k <- colMeans(cell_prob)
      prob_mat[s, ] <- p_k

      expected_count <- p_k * N
      expected_count <- pmax(expected_count, .Machine$double.eps)

      T_obs[s] <- sum((observed_count - expected_count)^2 / expected_count)

      if (K == 1L) {
        w_rep <- rep.int(1L, N)
      } else {
        u_rand <- stats::runif(N)
        cum_cell_prob <- t(apply(cell_prob[, -K, drop = FALSE], 1L, cumsum))
        w_rep <- rowSums(u_rand > cum_cell_prob) + 1L
      }

      replicated_count <- tabulate(w_rep, nbins = K)
      T_rep[s] <- sum((replicated_count - expected_count)^2 / expected_count)
    }

    prob_list[[ch]] <- prob_mat
    T_obs_list[[ch]] <- T_obs
    T_rep_list[[ch]] <- T_rep
  }

  prob_all <- do.call(rbind, prob_list)
  T_obs_all <- unlist(T_obs_list, use.names = FALSE)
  T_rep_all <- unlist(T_rep_list, use.names = FALSE)

  col_names <- paste0("Pr(W=", seq_len(K), ")")
  colnames(prob_all) <- col_names

  estimated_prob <- colMeans(prob_all)
  names(estimated_prob) <- col_names

  estimated_prob_ci <- apply(
    prob_all,
    2L,
    stats::quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )

  rownames(estimated_prob_ci) <- c("2.5%", "97.5%")
  colnames(estimated_prob_ci) <- col_names

  names(empirical_prob) <- col_names
  names(observed_count) <- col_names

  expected_count <- estimated_prob * N
  names(expected_count) <- col_names

  posterior_p_value <- mean(T_rep_all >= T_obs_all, na.rm = TRUE)

  out <- list(
    Estimated.Prob = estimated_prob,
    Estimated.Prob.ci = estimated_prob_ci,
    Empirical.Prob = empirical_prob,
    Observed.Count = observed_count,
    Expected.Count = expected_count,
    Observed.T = T_obs_all,
    Replicated.T = T_rep_all,
    Posterior.P.value = posterior_p_value,
    N = N,
    K = K,
    call = match.call()
  )

  class(out) <- "GoF.BayesODT"
  out
}


#' Print Goodness-of-Fit Test Results
#'
#' Print method for objects of class \code{"GoF.BayesODT"}.
#'
#' @param x An object of class \code{"GoF.BayesODT"} returned by \code{\link{GoF}}.
#' @param digits Number of digits to print.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.GoF.BayesODT <- function(x, digits = 4, ...) {

  est <- x$Estimated.Prob
  ci <- x$Estimated.Prob.ci
  emp <- x$Empirical.Prob
  obs_count <- x$Observed.Count
  exp_count <- x$Expected.Count

  K <- length(est)

  tab <- data.frame(
    Category = seq_len(K),
    Observed.Count = as.numeric(obs_count),
    Expected.Count = as.numeric(exp_count),
    Estimated.Prob = as.numeric(est),
    Lower = as.numeric(ci[1L, ]),
    Upper = as.numeric(ci[2L, ]),
    Empirical.Prob = as.numeric(emp),
    Difference = as.numeric(est - emp)
  )

  cat("\nGoodness-of-Fit Test for BayesODT\n")
  cat("---------------------------------\n")
  cat("Method: Posterior predictive check\n")
  cat("Discrepancy: sum_k (N_k - p_k * N)^2 / (p_k * N)\n\n")

  cat("Posterior predictive P-value:",
      formatC(x$Posterior.P.value, digits = digits, format = "f"), "\n")
  cat("Number of observed ratings:", x$N, "\n")
  cat("Number of posterior draws:", length(x$Observed.T), "\n\n")

  print(tab, digits = digits, row.names = FALSE)

  cat("\nInterpretation:\n")
  cat("  Values of the posterior predictive P-value near 0 or 1 may indicate\n")
  cat("  lack of fit. Values away from 0 and 1 suggest that the observed rating\n")
  cat("  distribution is not unusual under the fitted posterior predictive model.\n\n")

  invisible(x)
}
