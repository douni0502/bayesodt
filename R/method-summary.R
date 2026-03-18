#' Summarize posterior samples from a BayesODT fit
#'
#' Computes posterior summaries for parameters stored in a fitted
#' `"BayesODT"` object by combining saved draws across MCMC chains.
#' For each parameter, the posterior mean, standard deviation, and
#' equal-tailed credible interval are reported.
#'
#' @param object An object of class `"BayesODT"`.
#' @param probs A numeric vector of length 2 giving the lower and upper
#'   posterior probabilities for the credible interval.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class `"summary.BayesODT"`, which is a list containing
#'   posterior summaries for available parameters, together with basic
#'   model information.
#'
#' @details
#' Posterior draws are combined across all chains before summaries are computed.
#' Matrix-valued parameters are summarized column-wise. For interpretability,
#' the posterior draws of `a` are centered by subtracting the grand mean of all
#' saved `a` draws, and the same shift is applied to the posterior draws of `A`.
#'
#' @export
summary.BayesODT <- function(object,
                             probs = c(0.025, 0.975),
                             ...) {
  if (!inherits(object, "BayesODT")) {
    stop("Object is not of class 'BayesODT'.", call. = FALSE)
  }
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs)) {
    stop("`probs` must be a numeric vector of length 2.", call. = FALSE)
  }
  if (any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must satisfy 0 <= probs[1] < probs[2] <= 1.", call. = FALSE)
  }

  .summ_mat <- function(samples, probs) {
    if (!is.matrix(samples)) {
      samples <- as.matrix(samples)
    }
    out <- cbind(
      mean = colMeans(samples),
      sd = apply(samples, 2, stats::sd),
      t(apply(samples, 2, stats::quantile, probs = probs, names = FALSE))
    )
    colnames(out)[3:4] <- paste0(probs * 100, "%")
    out
  }

  .summ_vec <- function(samples, probs) {
    samples <- as.numeric(samples)
    out <- c(
      mean = mean(samples),
      sd = stats::sd(samples),
      stats::quantile(samples, probs = probs, names = FALSE)
    )
    names(out)[3:4] <- paste0(probs * 100, "%")
    out
  }

  .combine_chain_samples <- function(chains, name) {
    vals <- lapply(chains, function(ch) ch[[name]])
    vals <- vals[!vapply(vals, is.null, logical(1))]

    if (length(vals) == 0L) {
      return(NULL)
    }

    first <- vals[[1]]

    if (is.matrix(first) || is.vector(first)) {
      vals <- lapply(vals, function(x) {
        if (!is.matrix(x)) as.matrix(x) else x
      })
      return(do.call(rbind, vals))
    }

    NULL
  }

  chains <- object$chains
  if (is.null(chains) || !is.list(chains) || length(chains) == 0L) {
    stop("`object$chains` is missing or invalid.", call. = FALSE)
  }

  # combine posterior samples across chains
  A_samples     <- .combine_chain_samples(chains, "A")
  a_samples     <- .combine_chain_samples(chains, "a")
  tau_a_samples <- .combine_chain_samples(chains, "tau_a")
  b_samples     <- .combine_chain_samples(chains, "b")
  mu_b_samples  <- .combine_chain_samples(chains, "mu_b")
  tau_b_samples <- .combine_chain_samples(chains, "tau_b")
  u_samples     <- .combine_chain_samples(chains, "u")

  # apply the same global shift used in the original summary code:
  # a <- a - mean(a), A <- A - mean(a)
  a_shift <- NULL
  if (!is.null(a_samples)) {
    a_shift <- mean(a_samples)
    a_samples <- a_samples - a_shift
  }
  if (!is.null(A_samples) && !is.null(a_shift)) {
    A_samples <- A_samples - a_shift
  }

  # missing-rating information
  if (!is.null(object$w)) {
    n_total <- length(object$w)
    n_missing <- sum(is.na(object$w))
    n_observed <- n_total - n_missing
    prop_missing <- n_missing / n_total

    missing_info <- list(
      n_total = n_total,
      n_observed = n_observed,
      n_missing = n_missing,
      prop_missing = prop_missing
    )
  } else {
    missing_info <- NULL
  }

  result <- list(
    call = object$call,
    mcmc_info = object$mcmc_info,
    dims = object$dim_info,
    status = if (is.null(object$D)) "unknown" else "known",
    probs = probs,
    missing_info = missing_info
  )

  if (!is.null(A_samples)) {
    out_A <- .summ_mat(A_samples, probs)
    if (!is.null(object$param_names$A)) {
      rownames(out_A) <- object$param_names$A
    } else {
      rownames(out_A) <- paste0("A[", seq_len(ncol(A_samples)), "]")
    }
    result$A <- out_A
  }

  if (!is.null(a_samples)) {
    out_a <- .summ_mat(a_samples, probs)
    if (!is.null(object$param_names$a)) {
      rownames(out_a) <- object$param_names$a
    } else {
      rownames(out_a) <- paste0("a[", seq_len(ncol(a_samples)), "]")
    }
    result$a <- out_a
  }

  if (!is.null(b_samples)) {
    out_b <- .summ_mat(b_samples, probs)
    if (!is.null(object$param_names$b)) {
      rownames(out_b) <- object$param_names$b
    } else {
      rownames(out_b) <- paste0("b[", seq_len(ncol(b_samples)), "]")
    }
    result$b <- out_b
  }

  if (!is.null(u_samples)) {
    out_u <- .summ_mat(u_samples, probs)
    if (!is.null(object$param_names$u)) {
      rownames(out_u) <- object$param_names$u
    } else {
      rownames(out_u) <- paste0("u[", seq_len(ncol(u_samples)), "]")
    }
    result$u <- out_u
  }

  if (!is.null(tau_a_samples)) result$tau_a <- .summ_vec(tau_a_samples, probs)
  if (!is.null(mu_b_samples)) result$mu_b <- .summ_vec(mu_b_samples, probs)
  if (!is.null(tau_b_samples)) result$tau_b <- .summ_vec(tau_b_samples, probs)

  class(result) <- "summary.BayesODT"
  result
}


#' Print a summary of a BayesODT fit
#'
#' Displays posterior summaries from a `"summary.BayesODT"` object,
#' including basic model information and posterior summaries for
#' available parameters.
#'
#' @param x An object of class `"summary.BayesODT"`.
#' @param digits Number of digits to print.
#' @param max_items Maximum number of rows to print for each matrix-valued
#'   parameter block.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The input object `x`, invisibly.
#'
#' @export
print.summary.BayesODT <- function(x, digits = 3, max_items = 6, ...) {
  if (!inherits(x, "summary.BayesODT")) stop("Object is not of class 'summary.BayesODT'.")
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) stop("`digits` must be a single non-negative number.")
  if (!is.numeric(max_items) || length(max_items) != 1L || is.na(max_items) || max_items < 1) stop("`max_items` must be a single positive number.")

  max_items <- as.integer(max_items)

  cat("\nBayesODT summary\n")
  cat("----------------\n")

  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  if (!is.null(x$dims)) {
    m <- x$dims$m
    n <- x$dims$n
    K <- x$dims$K

    cat("Data:\n")
    cat("  Patients (m):   ", m, "\n", sep = "")
    cat("  Raters (n):     ", n, "\n", sep = "")
    cat("  Categories (K): ", K, "\n", sep = "")
  }

  if (!is.null(x$missing_info)) {
    cat("Missing ratings:\n")
    cat("  Total ratings:    ", x$missing_info$n_total, "\n", sep = "")
    cat("  Observed ratings: ", x$missing_info$n_observed, "\n", sep = "")
    cat("  Missing ratings:  ", x$missing_info$n_missing, "\n", sep = "")
    cat("  Missing rate:     ", round(100 * x$missing_info$prop_missing, digits), "%\n", sep = "")
  }

  if (!is.null(x$status)) {
    cat("Disease status: ", if (identical(x$status, "known")) "Known" else "Unknown", "\n", sep = "")
  }

  if (!is.null(x$mcmc_info)) {
    cat("\nMCMC:\n")
    if (!is.null(x$mcmc_info$n.chains)) cat("  Number of chains:      ", x$mcmc_info$n.chains, "\n", sep = "")
    if (!is.null(x$mcmc_info$n.iter)) cat("  Iterations per chain:  ", x$mcmc_info$n.iter, "\n", sep = "")
    if (!is.null(x$mcmc_info$n.burnin)) cat("  Burn-in per chain:     ", x$mcmc_info$n.burnin, "\n", sep = "")
    if (!is.null(x$mcmc_info$n.thin)) cat("  Thinning interval:     ", x$mcmc_info$n.thin, "\n", sep = "")
    if (!is.null(x$mcmc_info$n.save)) cat("  Saved draws per chain: ", x$mcmc_info$n.save, "\n", sep = "")
  }

  if (!is.null(x$probs)) {
    cat("\nCredible interval: ", x$probs[1] * 100, "%, ", x$probs[2] * 100, "%\n", sep = "")
  }

  .print_block <- function(mat, title, accessor = NULL) {
    if (is.null(mat)) return(invisible(NULL))
    if (!is.matrix(mat)) mat <- as.matrix(mat)

    cat("\n", title, "\n", sep = "")
    cat(strrep("-", nchar(title)), "\n", sep = "")

    nr <- nrow(mat)
    show_n <- min(nr, max_items)

    print(utils::head(round(mat, digits = digits), show_n))

    if (nr > show_n) {
      cat("... (", nr - show_n, " more rows)\n", sep = "")
      if (!is.null(accessor)) {
        cat("Use `summary(fit)$", accessor, "` to access the full table.\n", sep = "")
      }
    }

    invisible(NULL)
  }

  .print_hyper <- function(vec, title, accessor = NULL) {
    if (is.null(vec)) {
      return(invisible(NULL))
    }

    cat("\n", title, "\n", sep = "")
    cat(strrep("-", nchar(title)), "\n", sep = "")
    print(round(vec, digits = digits))

    invisible(NULL)
  }

  .print_block(x$A, "A: cutpoints", "A")
  .print_block(x$a, "a: rater-specific bias parameters", "a")
  .print_block(x$b, "b: rater-specific magnifier parameters", "b")
  .print_block(x$u, "u: latent disease severity parameters", "u")

  .print_hyper(x$tau_a, "tau_a: precision for a", "tau_a")
  .print_hyper(x$mu_b,  "mu_b: population mean for b", "mu_b")
  .print_hyper(x$tau_b, "tau_b: precision for b", "tau_b")

  cat("\n")
  invisible(x)
}
