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
simulate_data <- function(m = 50, n = 10, K = 4,
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
