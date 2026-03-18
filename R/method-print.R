#' Print a BayesODT object
#'
#' Displays basic information about a fitted `"BayesODT"` object,
#' including data dimensions, missing-rating information, disease-status
#' setting, and MCMC information.
#'
#' @param x An object of class `"BayesODT"`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The input object `x`, invisibly.
#' @export
print.BayesODT <- function(x, ...) {
  if (!inherits(x, "BayesODT")) {
    stop("Object is not of class 'BayesODT'.")
  }

  m <- if (!is.null(x$dim_info$m)) {
    x$dim_info$m
  } else if (!is.null(x$w)) {
    nrow(x$w)
  } else {
    NA_integer_
  }

  n <- if (!is.null(x$dim_info$n)) {
    x$dim_info$n
  } else if (!is.null(x$w)) {
    ncol(x$w)
  } else {
    NA_integer_
  }

  K <- if (!is.null(x$dim_info$K)) {
    x$dim_info$K
  } else {
    NA_integer_
  }

  if (!is.null(x$w)) {
    n_total <- length(x$w)
    n_missing <- sum(is.na(x$w))
    n_observed <- n_total - n_missing
    prop_missing <- n_missing / n_total
  } else {
    n_total <- NA_integer_
    n_missing <- NA_integer_
    n_observed <- NA_integer_
    prop_missing <- NA_real_
  }

  cat("\nBayesODT fit\n")
  cat("-----------\n")

  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  cat("Data:\n")
  cat("  Patients (m):   ", m, "\n", sep = "")
  cat("  Raters (n):     ", n, "\n", sep = "")
  cat("  Categories (K): ", K, "\n", sep = "")

  if (!is.na(n_total)) {
    cat("  Total ratings:  ", n_total, "\n", sep = "")
    cat("  Observed:       ", n_observed, "\n", sep = "")
    cat("  Missing:        ", n_missing, "\n", sep = "")
    cat("  Missing rate:   ", round(100 * prop_missing, 2), "%\n", sep = "")
  }

  cat("  Disease status: ", if (is.null(x$D)) "Unknown" else "Known", "\n", sep = "")

  if (!is.null(x$mcmc_info)) {
    cat("\nMCMC:\n")
    if (!is.null(x$mcmc_info$n.chains)) {
      cat("  Chains:         ", x$mcmc_info$n.chains, "\n", sep = "")
    }
    if (!is.null(x$mcmc_info$n.iter)) {
      cat("  Iterations:     ", x$mcmc_info$n.iter, "\n", sep = "")
    }
    if (!is.null(x$mcmc_info$n.burnin)) {
      cat("  Burn-in:        ", x$mcmc_info$n.burnin, "\n", sep = "")
    }
    if (!is.null(x$mcmc_info$n.thin)) {
      cat("  Thin:           ", x$mcmc_info$n.thin, "\n", sep = "")
    }
    if (!is.null(x$mcmc_info$n.save)) {
      cat("  Saved draws:    ", x$mcmc_info$n.save, " per chain\n", sep = "")
    }
  }

  if (!is.null(x$chains) && is.list(x$chains) && length(x$chains) > 0) {
    nm <- names(x$chains[[1]])
    if (!is.null(nm)) {
      cat("\nStored objects in each chain:\n")
      cat("  ", paste(nm, collapse = ", "), "\n", sep = "")
    }
  }

  cat("\nUse summary(x) for posterior summaries.\n")
  cat("Use plot(x) for posterior plots.\n\n")

  invisible(x)
}
