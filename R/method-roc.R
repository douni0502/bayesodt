#' Compute ROC curves from a fitted BayesODT object
#'
#' Generic function for extracting ROC curve information from fitted objects.
#'
#' @param object A fitted object.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object representing ROC curve information.
#' @export
bayes_roc <- function(object) {
  UseMethod("bayes_roc")
}


#' Compute ROC curves for BayesODT objects
#'
#' Extract ROC curve information from a fitted `"BayesODT"` object.
#' Posterior draws are combined across all saved MCMC chains before the
#' ROC summaries are computed.
#'
#' @param object An object of class `"BayesODT"`.
#' @param type Type of ROC curve to compute. One of `"individual"`,
#'   `"ROC1"`, `"ROC2"`, `"smoothed"`, `"pROC"`, or `"smoothed.pROC"`.
#' @param rater Integer rater index used when `type = "individual"`.
#' @param h Numeric vector of threshold values used when
#'   `type = "smoothed"`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class `"BayesODT_roc"` containing false-positive
#'   rates, true-positive rates, AUC, and related ROC information.
#'
#' @details
#' For posterior-model-based ROC types (`"individual"`, `"ROC1"`,
#' `"ROC2"`, and `"smoothed"`), saved posterior draws are first combined
#' across chains and then averaged to obtain ROC summaries.
#'
#' The empirical ROC types (`"pROC"` and `"smoothed.pROC"`) are computed
#' directly from the observed ratings and therefore require observed
#' disease status `D`.
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

    out <- switch(
      type,
      individual = .make_individual_roc(post, rater = rater, probs = probs),
      ROC1       = .make_ROC1(post, probs = probs),
      ROC2       = .make_ROC2(post, probs = probs),
      smoothed   = .make_smoothed_roc(post, h = h, probs = probs)
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

# CI 계산 함수
.compute_auc_ci <- function(auc_samples, probs = c(0.025, 0.975)) {
  if (length(probs) != 2L || anyNA(probs) || probs[1] >= probs[2]) {
    stop("`probs` must be a numeric vector of length 2 with probs[1] < probs[2].",
         call. = FALSE)
  }

stats::quantile(auc_samples, probs = probs, names = FALSE, na.rm = TRUE)
}

# validation function
.check_BayesODT_for_roc <- function(object) {
  if (!inherits(object, "BayesODT")) {
    stop("`object` must inherit from class 'BayesODT'.", call. = FALSE)
  }

  if (is.null(object$chains) || !is.list(object$chains) || length(object$chains) == 0L) {
    stop("The fitted object has no valid `chains` component.", call. = FALSE)
  }

  ch1 <- object$chains[[1]]
  needed <- c("A", "a", "b", "u")
  miss <- needed[vapply(needed, function(x) is.null(ch1[[x]]), logical(1))]

  if (length(miss) > 0L) {
    stop(sprintf(
      "The fitted object is missing required chain component(s): %s",
      paste(miss, collapse = ", ")
    ), call. = FALSE)
  }

  if (is.null(object$w)) {
    stop("The fitted object is missing `w`.", call. = FALSE)
  }
}


# helper: combine posterior draws across chains
.combine_chain_samples <- function(chains, name) {
  vals <- lapply(chains, function(ch) ch[[name]])
  vals <- vals[!vapply(vals, is.null, logical(1))]

  if (length(vals) == 0L) {
    return(NULL)
  }

  vals <- lapply(vals, function(x) {
    if (!is.matrix(x)) as.matrix(x) else x
  })

  do.call(rbind, vals)
}


# collect posterior draws needed for ROC calculation
.combine_roc_posterior_samples <- function(object) {
  A <- .combine_chain_samples(object$chains, "A")  # cutpoints
  a <- .combine_chain_samples(object$chains, "a")  # rater bias
  b <- .combine_chain_samples(object$chains, "b")  # rater magnifier
  u <- .combine_chain_samples(object$chains, "u")  # latent disease severity

  if (is.null(A) || is.null(a) || is.null(b) || is.null(u)) {
    stop("Failed to collect posterior samples needed for ROC computation.",
         call. = FALSE)
  }

  if (!all(c(nrow(A), nrow(a), nrow(b), nrow(u)) == nrow(A))) {
    stop("Posterior sample sizes do not match across parameters.",
         call. = FALSE)
  }

  list(
    A = A,
    a = a,
    b = b,
    u = u,
    num = nrow(A),
    n = ncol(a),
    m = ncol(u),
    K = ncol(A)
  )
}


# author-style AUC from discrete ROC points
.compute_auc_author <- function(fpr, tpr) {
  fpr <- as.numeric(fpr)
  tpr <- as.numeric(tpr)

  # simulation code connects points from (1,1) -> ... -> (0,0)
  ord <- order(fpr, decreasing = TRUE)
  fpr <- fpr[ord]
  tpr <- tpr[ord]

  fpr_full <- c(1, fpr, 0)
  tpr_full <- c(1, tpr, 0)

  sum((head(fpr_full, -1) - tail(fpr_full, -1)) *
        (head(tpr_full, -1) + tail(tpr_full, -1)) / 2)
}


# common ROC computation function
.compute_roc_components <- function(post) {
  num <- post$num
  n <- post$n
  m <- post$m
  K_cut <- post$K

  sp <- array(NA_real_, dim = c(n, K_cut, num))
  se <- array(NA_real_, dim = c(n, K_cut, num))
  sp.mean <- matrix(NA_real_, nrow = num, ncol = K_cut)
  se.mean <- matrix(NA_real_, nrow = num, ncol = K_cut)

  for (s in seq_len(num)) {
    A_s <- post$A[s, ]   # cutpoints
    a_s <- post$a[s, ]   # rater bias
    b_s <- post$b[s, ]   # rater magnifier
    u_s <- post$u[s, ]   # latent disease severity

    wd <- pnorm(u_s)
    wn <- 1 - wd

    sum_wd <- sum(wd)
    sum_wn <- sum(wn)

    # individual raters
    for (j in seq_len(n)) {
      eta_j <- a_s[j] + b_s[j] * u_s

      for (k in seq_len(K_cut)) {
        p_pos <- 1 - pnorm(A_s[k] - eta_j)  # P(W > k | u)

        sp[j, k, s] <- sum(p_pos * wn) / sum_wn   # FPR
        se[j, k, s] <- sum(p_pos * wd) / sum_wd   # TPR
      }
    }

    # mean rater
    a_bar <- mean(a_s)
    b_bar <- mean(b_s)
    eta_bar <- a_bar + b_bar * u_s

    for (k in seq_len(K_cut)) {
      p_pos <- 1 - pnorm(A_s[k] - eta_bar)

      sp.mean[s, k] <- sum(p_pos * wn) / sum_wn
      se.mean[s, k] <- sum(p_pos * wd) / sum_wd
    }
  }

  list(
    sp = sp,
    se = se,
    sp.mean = sp.mean,
    se.mean = se.mean,
    n = n,
    m = m,
    K = K_cut,
    num = num
  )
}


# add ROC endpoints
.finish_roc_curve <- function(fpr, tpr) {
  fpr <- c(0, fpr, 1)
  tpr <- c(0, tpr, 1)

  ord <- order(fpr, tpr)

  list(
    fpr = fpr[ord],
    tpr = tpr[ord]
  )
}


.make_individual_roc <- function(post, rater = 1L, probs = c(0.025, 0.975)) {
  comp <- .compute_roc_components(post)

  if (length(rater) != 1L || is.na(rater) || rater < 1L || rater > comp$n) {
    stop(sprintf("`rater` must be an integer between 1 and %d.", comp$n),
         call. = FALSE)
  }

  fpr <- vapply(seq_len(comp$K), function(k) mean(comp$sp[rater, k, ]), numeric(1))
  tpr <- vapply(seq_len(comp$K), function(k) mean(comp$se[rater, k, ]), numeric(1))

  auc_samples <- vapply(
    seq_len(comp$num),
    function(s) .compute_auc_author(comp$sp[rater, , s], comp$se[rater, , s]),
    numeric(1)
  )

  out <- .finish_roc_curve(fpr, tpr)
  out$cutpoint_index <- seq_len(comp$K)
  out$auc <- .compute_auc_author(fpr, tpr)
  out$auc_samples <- auc_samples
  out$auc_ci <- .compute_auc_ci(auc_samples, probs = probs)
  out$type <- "individual"
  out$rater <- rater
  out
}


.make_ROC1 <- function(post, probs = c(0.025, 0.975)) {
  comp <- .compute_roc_components(post)

  fpr <- vapply(seq_len(comp$K), function(k) mean(comp$sp[, k, ]), numeric(1))
  tpr <- vapply(seq_len(comp$K), function(k) mean(comp$se[, k, ]), numeric(1))

  auc_samples <- c(
    sapply(seq_len(comp$num), function(s) {
      sapply(seq_len(comp$n), function(j) {
        .compute_auc_author(comp$sp[j, , s], comp$se[j, , s])
      })
    })
  )

  out <- .finish_roc_curve(fpr, tpr)
  out$cutpoint_index <- seq_len(comp$K)
  out$auc <- .compute_auc_author(fpr, tpr)
  out$auc_samples <- auc_samples
  out$auc_ci <- .compute_auc_ci(auc_samples, probs = probs)
  out$type <- "ROC1"
  out
}


.make_ROC2 <- function(post, probs = c(0.025, 0.975)) {
  comp <- .compute_roc_components(post)

  fpr <- vapply(seq_len(comp$K), function(k) mean(comp$sp.mean[, k]), numeric(1))
  tpr <- vapply(seq_len(comp$K), function(k) mean(comp$se.mean[, k]), numeric(1))

  auc_samples <- vapply(
    seq_len(comp$num),
    function(s) .compute_auc_author(comp$sp.mean[s, ], comp$se.mean[s, ]),
    numeric(1)
  )

  out <- .finish_roc_curve(fpr, tpr)
  out$cutpoint_index <- seq_len(comp$K)
  out$auc <- .compute_auc_author(fpr, tpr)
  out$auc_samples <- auc_samples
  out$auc_ci <- .compute_auc_ci(auc_samples, probs = probs)
  out$type <- "ROC2"
  out
}


# smoothed ROC
.make_smoothed_roc <- function(post,
                               h = seq(-5, 5, length.out = 200),
                               probs = c(0.025, 0.975)) {
  h <- as.numeric(h)
  num <- post$num

  fpr_mat <- matrix(NA_real_, nrow = num, ncol = length(h))
  tpr_mat <- matrix(NA_real_, nrow = num, ncol = length(h))

  for (s in seq_len(num)) {
    a_s <- post$a[s, ]
    b_s <- post$b[s, ]
    u_s <- post$u[s, ]

    wd <- pnorm(u_s)
    wn <- 1 - wd

    a_bar <- mean(a_s)
    b_bar <- mean(b_s)

    for (g in seq_along(h)) {
      p_pos <- 1 - pnorm(h[g] - (a_bar + b_bar * u_s))
      fpr_mat[s, g] <- sum(p_pos * wn) / sum(wn)
      tpr_mat[s, g] <- sum(p_pos * wd) / sum(wd)
    }
  }

  fpr <- colMeans(fpr_mat)
  tpr <- colMeans(tpr_mat)

  auc_samples <- vapply(
    seq_len(num),
    function(s) .compute_auc_author(fpr_mat[s, ], tpr_mat[s, ]),
    numeric(1)
  )

  ord <- order(fpr)

  list(
    fpr = fpr[ord],
    tpr = tpr[ord],
    grid = h[ord],
    auc = .compute_auc_author(fpr, tpr),
    auc_samples = auc_samples,
    auc_ci = .compute_auc_ci(auc_samples, probs = probs),
    type = "smoothed"
  )
}


#' Print a BayesODT ROC object
#'
#' Displays basic information about an object returned by `bayes_roc()`.
#'
#' @param x An object of class `"BayesODT_roc"`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The input object `x`, invisibly.
#' @export
print.BayesODT_roc <- function(x, ...) {
  if (!inherits(x, "BayesODT_roc")) {
    stop("Object is not of class 'BayesODT_roc'.")
  }

  cat("ROC object for BayesODT\n")
  cat("Type: ", x$type, "\n", sep = "")
  cat("Number of points: ", length(x$fpr), "\n", sep = "")

  if (!is.null(x$rater)) {
    cat("Rater: ", x$rater, "\n", sep = "")
  }

  if (!is.null(x$auc)) {
    cat("AUC: ", round(x$auc, 4), "\n", sep = "")
  }

  if (!is.null(x$auc_ci)) {
    cat("AUC credible interval: [",
        round(x$auc_ci[1], 4), ", ",
        round(x$auc_ci[2], 4), "]\n", sep = "")
  }

  invisible(x)
}

.make_empirical_roc <- function(object, smooth = FALSE) {
  if (is.null(object$D)) {
    stop("Empirical ROC requires observed `D`, but `object$D` is NULL.",
         call. = FALSE)
  }

  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package `pROC` is required for empirical ROC.", call. = FALSE)
  }

  y_mat <- matrix(rep(object$D, ncol(object$w)),
                  nrow = nrow(object$w),
                  ncol = ncol(object$w))

  ok <- !is.na(object$w)

  response  <- as.vector(y_mat[ok])
  predictor <- as.vector(object$w[ok])

  roc_obj <- pROC::roc(
    response = response,
    predictor = predictor,
    smooth = smooth,
    quiet = TRUE
  )

  thresholds <- if (!is.null(roc_obj$thresholds)) roc_obj$thresholds else NULL

  list(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    thresholds = thresholds,
    auc = as.numeric(roc_obj$auc),
    type = if (smooth) "smoothed.pROC" else "pROC"
  )
}
