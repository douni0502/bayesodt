#' Generic function for ROC curve extraction
#' Extracts ROC curve information from a fitted object.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class `"BayesODT_roc"`.
#' @export
bayes_roc <- function(object, ...) {
  UseMethod("bayes_roc")
}


#' Compute ROC curves for a fitted BayesODT object
#'
#' Extract ROC curve summaries from a fitted `"BayesODT"` object.
#' Posterior draws from all saved MCMC chains are combined before posterior-based ROC summaries are computed.
#'
#' @param object An object of class `"BayesODT"`.
#' @param type Type of ROC curve to compute. One of `"individual"`, `"ROC1"`, `"ROC2"`, `"smoothed"`, `"pROC"`, or `"smoothed.pROC"`.
#' @param rater Integer rater index used when `type = "individual"`.
#' @param h Numeric vector of threshold values used when `type = "smoothed"`.
#' @param probs Numeric vector of probabilities used to compute credible intervals for posterior-based AUC estimates.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class `"BayesODT_roc"` containing ROC coordinates, AUC summaries, and related metadata.
#'
#' @details
#' For posterior-model-based ROC types, saved posterior draws are first combined across chains.
#' The `"individual"`, `"ROC1"`, and `"ROC2"` summaries are obtained by posterior averaging, whereas `"smoothed"` uses posterior mean parameter estimates on a fine threshold grid.
#'
#' The empirical ROC types (`"pROC"` and `"smoothed.pROC"`) are computed directly from the observed ratings and therefore require observed disease status `D`.
#'
#' @export
bayes_roc.BayesODT <- function(object,
                               type = c("individual", "ROC1", "ROC2",
                                        "smoothed", "pROC", "smoothed.pROC"),
                               rater = 1L,
                               h = seq(-5, 5, length.out = 200),
                               probs = c(0.025, 0.975),
                               ...) {
  type <- match.arg(type)

  .check_BayesODT_for_roc(object)

  if (type %in% c("individual", "ROC1", "ROC2", "smoothed")) {
    post <- .combine_roc_posterior_samples(object)

    out <- switch(type,
                  individual = .make_individual_roc(post, rater = rater, probs = probs),
                  ROC1 = .make_ROC1(post, probs = probs),
                  ROC2 = .make_ROC2(post, probs = probs),
                  smoothed = .make_smoothed_roc(post, h = h, probs = probs)
    )
  } else if (type == "pROC") {
    out <- .make_empirical_roc(object, smooth = FALSE)
  } else if (type == "smoothed.pROC") {
    out <- .make_empirical_roc(object, smooth = TRUE)
  }

  out$call <- match.call()
  class(out) <- c("BayesODT_roc", class(out))
  out
}

# Compute a credible interval from posterior AUC samples
.compute_auc_ci <- function(auc_samples, probs = c(0.025, 0.975)) {
  if (length(probs) != 2L || anyNA(probs) || any(probs <= 0 | probs >= 1) || probs[1] >= probs[2]) {
    stop("`probs` must be a numeric vector of length 2 with 0 < probs[1] < probs[2] < 1.", call. = FALSE)
  }
  stats::quantile(auc_samples, probs = probs, names = FALSE, na.rm = TRUE)
}

# Validate that a fitted BayesODT object contains the components required for ROC computation
.check_BayesODT_for_roc <- function(object) {
  if (!inherits(object, "BayesODT")) stop("`object` must inherit from class 'BayesODT'.", call. = FALSE)
  if (is.null(object$chains) || !is.list(object$chains) || length(object$chains) == 0L) stop("The fitted object has no valid `chains` component.", call. = FALSE)

  # Required posterior samples used in ROC calculations
  ch1 <- object$chains[[1]]
  needed <- c("A", "a", "b", "u")
  miss <- needed[vapply(needed, function(x) is.null(ch1[[x]]), logical(1))]
  if (length(miss) > 0L) stop(sprintf("The fitted object is missing required chain component(s): %s", paste(miss, collapse = ", ")), call. = FALSE)

  if (is.null(object$w)) stop("The fitted object is missing `w`.", call. = FALSE)
}


# Combine posterior samples with the same name across MCMC chains
.combine_chain_samples <- function(chains, name) {
  vals <- lapply(chains, function(ch) ch[[name]])
  vals <- vals[!vapply(vals, is.null, logical(1))]

  # If the component is absent from all chains, return NULL
  if (length(vals) == 0L) {
    return(NULL)
  }
  # Ensure all extracted values are matrices
  vals <- lapply(vals, function(x) {
    if (!is.matrix(x)) as.matrix(x) else x
  })
  # Combine posterior draws from all chains into one matrix
  do.call(rbind, vals)
}

# Combine posterior samples required for ROC computation
.combine_roc_posterior_samples <- function(object) {
  A <- .combine_chain_samples(object$chains, "A")
  a <- .combine_chain_samples(object$chains, "a")
  b <- .combine_chain_samples(object$chains, "b")
  u <- .combine_chain_samples(object$chains, "u")

  if (is.null(A) || is.null(a) || is.null(b) || is.null(u)) stop("Failed to collect posterior samples needed for ROC computation.", call. = FALSE)
  if (!all(c(nrow(A), nrow(a), nrow(b), nrow(u)) == nrow(A))) stop("Posterior sample sizes do not match across parameters.", call. = FALSE)

  # Return combined posterior samples and their dimensions
  list(A = A, a = a, b = b, u = u,
       num = nrow(A), n = ncol(a), m = ncol(u), K = ncol(A))
}


# Convert cumulative probabilities P(W <= k) to category probabilities P(W = k)
.cumprob_to_mass <- function(cdf_leq) {
  cdf_leq <- as.numeric(cdf_leq)
  K_cut <- length(cdf_leq) # K_cut is the number of cumulative cutpoints
  K_cat <- K_cut + 1L # K_cat is the number of ordinal categories

  out <- numeric(K_cat)

  out[1] <- cdf_leq[1]
  if (K_cut >= 2L) {
    for (k in 2:K_cut) {
      out[k] <- cdf_leq[k] - cdf_leq[k - 1L]
    }
  }
  out[K_cat] <- 1 - cdf_leq[K_cut]

  out
}


# Compute the ordinal AUC from category probabilities
.compute_auc <- function(p0, p1) {
  p0 <- as.numeric(p0) # p0[k] = P(W = k | D = 0)
  p1 <- as.numeric(p1) # p1[k] = P(W = k | D = 1)

  if (length(p0) != length(p1)) stop("`p0` and `p1` must have the same length.", call. = FALSE)

  K_cat <- length(p0)
  auc <- 0

  for (k in seq_len(K_cat - 1L)) {
    auc <- auc + p0[k] * sum(p1[(k + 1L):K_cat])
  }
  auc <- auc + 0.5 * sum(p0 * p1)

  auc
}


# Compute posterior ROC components used by individual, ROC1, and ROC2
.compute_roc_components <- function(post) {
  num <- post$num
  n <- post$n
  m <- post$m
  K_cut <- post$K
  K_cat <- K_cut + 1L

  # ROC coordinates and category probabilities
  sp <- array(NA_real_, dim = c(n, K_cut, num))
  se <- array(NA_real_, dim = c(n, K_cut, num))
  sp.mean <- matrix(NA_real_, nrow = num, ncol = K_cut)
  se.mean <- matrix(NA_real_, nrow = num, ncol = K_cut)

  p0.ind <- array(NA_real_, dim = c(n, K_cat, num))
  p1.ind <- array(NA_real_, dim = c(n, K_cat, num))
  p0.mean <- matrix(NA_real_, nrow = num, ncol = K_cat)
  p1.mean <- matrix(NA_real_, nrow = num, ncol = K_cat)

  for (s in seq_len(num)) {
    A_s <- post$A[s, ]
    a_s <- post$a[s, ]
    b_s <- post$b[s, ]
    u_s <- post$u[s, ]

    wd <- pnorm(u_s)
    wn <- 1 - wd
    sum_wd <- sum(wd)
    sum_wn <- sum(wn)

    # individual raters
    for (j in seq_len(n)) {
      eta_j <- a_s[j] + b_s[j] * u_s
      cdf0_j <- numeric(K_cut)
      cdf1_j <- numeric(K_cut)

      for (k in seq_len(K_cut)) {
        cdf_leq <- pnorm(A_s[k] - eta_j)
        p_pos   <- 1 - cdf_leq

        cdf0_j[k] <- sum(cdf_leq * wn) / sum_wn
        cdf1_j[k] <- sum(cdf_leq * wd) / sum_wd
        sp[j, k, s] <- sum(p_pos * wn) / sum_wn
        se[j, k, s] <- sum(p_pos * wd) / sum_wd
      }

      p0.ind[j, , s] <- .cumprob_to_mass(cdf0_j)
      p1.ind[j, , s] <- .cumprob_to_mass(cdf1_j)
    }

    # Mean rater
    eta_bar <- mean(a_s) + mean(b_s) * u_s
    cdf0_bar <- numeric(K_cut)
    cdf1_bar <- numeric(K_cut)

    for (k in seq_len(K_cut)) {
      cdf_leq <- pnorm(A_s[k] - eta_bar)
      p_pos   <- 1 - cdf_leq

      cdf0_bar[k] <- sum(cdf_leq * wn) / sum_wn
      cdf1_bar[k] <- sum(cdf_leq * wd) / sum_wd
      sp.mean[s, k] <- sum(p_pos * wn) / sum_wn
      se.mean[s, k] <- sum(p_pos * wd) / sum_wd
    }

    p0.mean[s, ] <- .cumprob_to_mass(cdf0_bar)
    p1.mean[s, ] <- .cumprob_to_mass(cdf1_bar)
  }

  list(sp = sp, se = se, sp.mean = sp.mean, se.mean = se.mean,
       p0.ind = p0.ind, p1.ind = p1.ind,
       p0.mean = p0.mean, p1.mean = p1.mean,
       n = n, m = m, K = K_cut, K_cat = K_cat, num = num)
}


# Finalize ROC coordinates by adding endpoints and sorting the points
.finish_roc_curve <- function(fpr, tpr) {
  fpr <- c(0, fpr, 1)
  tpr <- c(0, tpr, 1)

  ord <- order(fpr, tpr)

  list(fpr = fpr[ord], tpr = tpr[ord])
}


# Compute an ROC summary for a single rater
.make_individual_roc <- function(post, rater = 1L, probs = c(0.025, 0.975)) {
  comp <- .compute_roc_components(post)

  if (length(rater) != 1L || is.na(rater) || rater %% 1 != 0 || rater < 1L || rater > comp$n) {
    stop(sprintf("`rater` must be an integer between 1 and %d.", comp$n), call. = FALSE)
  }
  rater <- as.integer(rater)

  # Posterior-averaged ROC coordinates
  fpr <- vapply(seq_len(comp$K), function(k) mean(comp$sp[rater, k, ]), numeric(1))
  tpr <- vapply(seq_len(comp$K), function(k) mean(comp$se[rater, k, ]), numeric(1))

  auc_samples <- vapply(seq_len(comp$num), function(s) .compute_auc(comp$p0.ind[rater, , s], comp$p1.ind[rater, , s]), numeric(1))

  # Point-estimate AUC from posterior mean category probabilities
  p0_mat <- matrix(comp$p0.ind[rater, , ], nrow = comp$K_cat, ncol = comp$num)
  p1_mat <- matrix(comp$p1.ind[rater, , ], nrow = comp$K_cat, ncol = comp$num)

  p0_hat <- rowMeans(p0_mat)
  p1_hat <- rowMeans(p1_mat)

  out <- .finish_roc_curve(fpr, tpr)
  out$cutpoint_index <- seq_len(comp$K)
  out$auc <- .compute_auc(p0_hat, p1_hat)
  out$auc_samples <- auc_samples
  out$auc_ci <- .compute_auc_ci(auc_samples, probs = probs)
  out$type <- "individual"
  out$rater <- rater
  out
}


# Compute ROC1, the overall ROC obtained by averaging over raters
.make_ROC1 <- function(post, probs = c(0.025, 0.975)) {
  comp <- .compute_roc_components(post)

  # Posterior-averaged ROC1 coordinates
  fpr <- vapply(seq_len(comp$K), function(k) mean(comp$sp[, k, ]), numeric(1))
  tpr <- vapply(seq_len(comp$K), function(k) mean(comp$se[, k, ]), numeric(1))

  auc_samples <- vapply(seq_len(comp$num),
                        function(s) {
                          fpr_s <- vapply(seq_len(comp$K), function(k) mean(comp$sp[, k, s]), numeric(1))
                          tpr_s <- vapply(seq_len(comp$K), function(k) mean(comp$se[, k, s]), numeric(1))

                          cdf0_s <- 1 - fpr_s
                          cdf1_s <- 1 - tpr_s

                          p0_s <- .cumprob_to_mass(cdf0_s)
                          p1_s <- .cumprob_to_mass(cdf1_s)

                          .compute_auc(p0_s, p1_s)
                          },
                        numeric(1))

  # Point-estimate AUC from posterior-averaged ROC1 probabilities
  cdf0_hat <- 1 - fpr
  cdf1_hat <- 1 - tpr
  p0_hat <- .cumprob_to_mass(cdf0_hat)
  p1_hat <- .cumprob_to_mass(cdf1_hat)

  out <- .finish_roc_curve(fpr, tpr)
  out$cutpoint_index <- seq_len(comp$K)
  out$auc <- .compute_auc(p0_hat, p1_hat)
  out$auc_samples <- auc_samples
  out$auc_ci <- .compute_auc_ci(auc_samples, probs = probs)
  out$type <- "ROC1"
  out
}

# Compute ROC2, the ROC curve based on the posterior mean rater
.make_ROC2 <- function(post, probs = c(0.025, 0.975)) {
  comp <- .compute_roc_components(post)

  # Posterior-averaged ROC2 coordinates
  fpr <- vapply(seq_len(comp$K), function(k) mean(comp$sp.mean[, k]), numeric(1))
  tpr <- vapply(seq_len(comp$K), function(k) mean(comp$se.mean[, k]), numeric(1))

  auc_samples <- vapply(seq_len(comp$num), function(s) .compute_auc(comp$p0.mean[s, ], comp$p1.mean[s, ]), numeric(1))

  # Point-estimate AUC from posterior mean category probabilities
  p0_hat <- colMeans(comp$p0.mean)
  p1_hat <- colMeans(comp$p1.mean)

  out <- .finish_roc_curve(fpr, tpr)
  out$cutpoint_index <- seq_len(comp$K)
  out$auc <- .compute_auc(p0_hat, p1_hat)
  out$auc_samples <- auc_samples
  out$auc_ci <- .compute_auc_ci(auc_samples, probs = probs)
  out$type <- "ROC2"
  out
}


# Compute a smoothed ROC curve from posterior mean parameters
.make_smoothed_roc <- function(post,
                               h = seq(-10, 10, by = 0.1),
                               probs = c(0.025, 0.975)) {
  h <- as.numeric(h)

  a_hat <- mean(post$a)
  b_hat <- mean(post$b)
  u_hat <- colMeans(post$u)

  # Weights for D = 0 and D = 1
  wd <- pnorm(u_hat)
  wn <- 1 - wd
  sum_wd <- sum(wd)
  sum_wn <- sum(wn)

  # Smoothed ROC coordinates
  fpr <- numeric(length(h))
  tpr <- numeric(length(h))

  for (g in seq_along(h)) {
    cdf_leq <- pnorm(h[g] - (a_hat + b_hat * u_hat))
    p_pos <- 1 - cdf_leq

    fpr[g] <- sum(p_pos * wn) / sum_wn
    tpr[g] <- sum(p_pos * wd) / sum_wd
  }

  ord <- order(fpr, tpr)
  fpr <- fpr[ord]
  tpr <- tpr[ord]
  h_ord <- h[ord]

  # AUC from trapezoidal integration
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)

  list(fpr = fpr, tpr = tpr,
       grid = h_ord,
       auc = auc,
       auc_samples = NA_real_, auc_ci = NA_real_,
       type = "smoothed")
}


# Compute an empirical ROC curve using observed ratings
.make_empirical_roc <- function(object, smooth = FALSE) {
  if (is.null(object$D)) stop("Empirical ROC requires observed `D`, but `object$D` is NULL.", call. = FALSE)
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package `pROC` is required for empirical ROC.", call. = FALSE)

  y_mat <- matrix(rep(object$D, ncol(object$w)), nrow = nrow(object$w), ncol = ncol(object$w))

  ok <- !is.na(object$w)
  response  <- as.vector(y_mat[ok])
  predictor <- as.vector(object$w[ok])

  # fit the empirical ROC curve
  roc_obj <- pROC::roc(response = response, predictor = predictor, smooth = smooth, quiet = TRUE)

  thresholds <- if (!is.null(roc_obj$thresholds)) roc_obj$thresholds else NULL

  list(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities, thresholds = thresholds,
       auc = as.numeric(roc_obj$auc),
       type = if (smooth) "smoothed.pROC" else "pROC")
}

#' Print a BayesODT ROC object
#' Displays a concise summary of an object returned by [bayes_roc()].
#' @param x An object of class `"BayesODT_roc"`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The input object `x`, invisibly.
#' @export
print.BayesODT_roc <- function(x, ...) {
  if (!inherits(x, "BayesODT_roc")) stop("`x` must inherit from class 'BayesODT_roc'.", call. = FALSE)

  cat("ROC object for BayesODT\n")
  cat("Type: ", x$type, "\n", sep = "")
  cat("Number of points: ", length(x$fpr), "\n", sep = "")

  if (!is.null(x$rater)) {
    cat("Rater: ", x$rater, "\n", sep = "")
  }

  if (!is.null(x$auc) && !is.na(x$auc)) {
    cat("AUC: ", round(x$auc, 4), "\n", sep = "")
  }

  if (!is.null(x$auc_ci) && all(!is.na(x$auc_ci))) {
    cat("AUC credible interval: [",
        round(x$auc_ci[1], 4), ", ",
        round(x$auc_ci[2], 4), "]\n", sep = "")
  }

  invisible(x)
}
