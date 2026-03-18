#' Coerce objects to coda mcmc objects
#'
#' Generic function for converting objects to \code{"mcmc"} objects.
#'
#' @param x An R object.
#' @param ... Further arguments passed to methods.
#'
#' @return An object, typically of class \code{"mcmc"}.
#' @export
as.mcmc <- function(x, ...) {
  UseMethod("as.mcmc")
}

#' Coerce a BayesODT object to coda mcmc objects
#' Convert selected posterior samples from a `"BayesODT"` object to a `coda::mcmc` object or a `coda::mcmc.list` object.
#' @param x A fitted object of class `"BayesODT"`.
#' @param param Character string specifying which parameter to extract. One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s).
#'   If `index = NULL`, all components of the selected parameter are returned.
#'   For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`, only `index = 1` or `index = NULL` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object. If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object. If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#'
#' @return
#' If `x` contains one chain, returns an object of class `"mcmc"`.
#' If `x` contains multiple chains, returns an object of class `"mcmc.list"`.
#'
#' @export
as.mcmc.BayesODT <- function(x,
                             param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                             index = 1L,
                             start = NULL,
                             thin = NULL) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `as.mcmc.BayesODT()`. Please install it first.", call. = FALSE)
  }
  if (is.null(x$chains) || !is.list(x$chains) || length(x$chains) == 0L) {
    stop("The fitted object has no valid `chains` component.", call. = FALSE)
  }

  param <- match.arg(param)

  if (!is.null(index)) {
    if (!is.numeric(index) || anyNA(index) || any(index < 1) || any(index != as.integer(index))) {
      stop("`index` must be NULL or a vector of positive integers.", call. = FALSE)
    }
    index <- as.integer(index)
  }

  if (is.null(start)) {
    start <- 1L
  }

  if (is.null(thin)) {
    if (!is.null(x$mcmc_info) && !is.null(x$mcmc_info$n.thin)) {
      thin <- x$mcmc_info$n.thin
    } else {
      thin <- 1L
    }
  }

  mcmc_list <- lapply(seq_along(x$chains), function(i) {
    samples <- x$chains[[i]][[param]]

    if (is.null(samples)) {
      stop(sprintf("Chain %d does not contain samples for `%s`.", i, param),
           call. = FALSE)
    }

    # scalar/vector hyperparameters -> 1-column matrix
    if (!is.matrix(samples)) {
      samples <- as.matrix(samples)
    }

    # index = NULL 이면 전체 열 선택
    index_use <- if (is.null(index)) seq_len(ncol(samples)) else index

    if (any(index_use > ncol(samples))) {
      stop(sprintf(
        "`index` contains values larger than the number of available `%s` parameters (%d) in chain %d.",
        param, ncol(samples), i
      ), call. = FALSE)
    }

    samples_sub <- samples[, index_use, drop = FALSE]

    cn <- NULL
    if (!is.null(x$param_names) && !is.null(x$param_names[[param]])) {
      pnames <- x$param_names[[param]]

      if (length(pnames) == 1L && ncol(samples) == 1L) {
        cn <- pnames
      } else {
        cn <- pnames[index_use]
      }
    } else {
      if (ncol(samples) == 1L && length(index_use) == 1L && identical(index_use, 1L)) {
        cn <- param
      } else {
        cn <- paste0(param, "[", index_use, "]")
      }
    }

    colnames(samples_sub) <- cn

    coda::mcmc(samples_sub, start = start, thin = thin)
  })

  if (length(mcmc_list) == 1L) {
    return(mcmc_list[[1]])
  }

  coda::mcmc.list(mcmc_list)
}








#' Trace plot for BayesODT posterior samples
#' Plot an MCMC trace plot for selected posterior samples from a 'BayesODT' object using package 'coda'.
#' @param x A fitted object of class `BayesODT"`.
#' @param param Character string specifying which parameter to plot. One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s). For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`, only `index = 1` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object. If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object. If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#'
#' @return Invisibly returns the created `"mcmc"` object or `"mcmc.list"` object.
#'
#' @export
plot_trace <- function(x,
                       param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                       index = 1L,
                       start = NULL,
                       thin = NULL,
                       legend_pos = "topright",
                       ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `plot_trace()`. Please install it first.", call. = FALSE)
  }

  param <- match.arg(param)

  if (param == "A") {
    mcmc_obj <- as.mcmc(x, param = "A", index = NULL, start = start, thin = thin)

    .make_iter <- function(obj) {
      mcpar <- attr(obj, "mcpar")
      seq.int(from = mcpar[1], by = mcpar[3], length.out = nrow(as.matrix(obj)))
    }

    if (inherits(mcmc_obj, "mcmc")) {
      mat <- as.matrix(mcmc_obj)
      iter <- .make_iter(mcmc_obj)

      matplot(
        x = iter,
        y = mat,
        type = "l",
        lty = 1,
        lwd = 1,
        xlab = "Iterations",
        ylab = "Value",
        main = "Trace of A",
        ...
      )
      legend(
        legend_pos,
        legend = colnames(mat),
        lty = 1,
        lwd = 1,
        col = seq_len(ncol(mat)),
        bty = "n"
      )
    } else {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)

      nch <- length(mcmc_obj)
      graphics::par(mfrow = c(nch, 1))

      for (i in seq_len(nch)) {
        mat <- as.matrix(mcmc_obj[[i]])
        iter <- .make_iter(mcmc_obj[[i]])

        matplot(
          x = iter,
          y = mat,
          type = "l",
          lty = 1,
          lwd = 1,
          xlab = "Iterations",
          ylab = "Value",
          main = paste("Trace of A (chain", i, ")"),
          ...
        )
        legend(
          legend_pos,
          legend = colnames(mat),
          lty = 1,
          lwd = 1,
          col = seq_len(ncol(mat)),
          bty = "n"
        )
      }
    }

    return(invisible(mcmc_obj))
  }

  mcmc_obj <- as.mcmc(x, param = param, index = index, start = start, thin = thin)
  coda::traceplot(mcmc_obj)
  invisible(mcmc_obj)
}
#' Effective sample size for BayesODT posterior samples
#'
#' Compute the effective sample size (ESS) for selected posterior samples
#' from a `"BayesODT"` object using \pkg{coda}.
#'
#' @param x A fitted object of class `"BayesODT"`.
#' @param param Character string specifying which parameter to evaluate.
#'   One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s).
#'   For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`,
#'   only `index = 1` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object.
#'   If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object.
#'   If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#' @param ... Not used.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{parameter}{Parameter name.}
#'   \item{ESS}{Estimated effective sample size.}
#' }
#'
#' @export
effective_size <- function(x,
                           param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                           index = 1L,
                           start = NULL,
                           thin = NULL,
                           ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `effective_size()`. Please install it first.",
         call. = FALSE)
  }

  mcmc_obj <- as.mcmc(
    x,
    param = param,
    index = index,
    start = start,
    thin = thin
  )

  ess <- coda::effectiveSize(mcmc_obj)

  data.frame(
    parameter = names(ess),
    ESS = as.numeric(ess),
    row.names = NULL
  )
}

#' Geweke diagnostic for BayesODT posterior samples
#'
#' Compute Geweke convergence diagnostics for selected posterior samples
#' from a `"BayesODT"` object using \pkg{coda}.
#'
#' @param x A fitted object of class `"BayesODT"`.
#' @param param Character string specifying which parameter to evaluate.
#'   One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s).
#'   For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`,
#'   only `index = 1` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object.
#'   If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object.
#'   If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#' @param frac1 Fraction of the beginning of the chain to use.
#' @param frac2 Fraction of the end of the chain to use.
#' @param ... Not used.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{chain}{Chain index.}
#'   \item{parameter}{Parameter name.}
#'   \item{z_score}{Geweke Z-score.}
#' }
#'
#' @export
geweke_diag <- function(x,
                        param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                        index = 1L,
                        start = NULL,
                        thin = NULL,
                        frac1 = 0.1,
                        frac2 = 0.5,
                        ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `geweke_diag()`. Please install it first.",
         call. = FALSE)
  }

  if (!is.numeric(frac1) || length(frac1) != 1L || is.na(frac1) ||
      frac1 <= 0 || frac1 >= 1) {
    stop("`frac1` must be a single number in (0, 1).", call. = FALSE)
  }

  if (!is.numeric(frac2) || length(frac2) != 1L || is.na(frac2) ||
      frac2 <= 0 || frac2 >= 1) {
    stop("`frac2` must be a single number in (0, 1).", call. = FALSE)
  }

  if (frac1 + frac2 >= 1) {
    stop("`frac1 + frac2` must be less than 1.", call. = FALSE)
  }

  mcmc_obj <- as.mcmc(
    x,
    param = param,
    index = index,
    start = start,
    thin = thin
  )

  # single chain
  if (inherits(mcmc_obj, "mcmc")) {
    z <- coda::geweke.diag(mcmc_obj, frac1 = frac1, frac2 = frac2)$z

    return(data.frame(
      chain = 1L,
      parameter = names(z),
      z_score = as.numeric(z),
      row.names = NULL
    ))
  }

  # multiple chains
  if (inherits(mcmc_obj, "mcmc.list")) {
    out_list <- lapply(seq_along(mcmc_obj), function(i) {
      z <- coda::geweke.diag(mcmc_obj[[i]], frac1 = frac1, frac2 = frac2)$z

      data.frame(
        chain = i,
        parameter = names(z),
        z_score = as.numeric(z),
        row.names = NULL
      )
    })

    return(do.call(rbind, out_list))
  }

  stop("`as.mcmc(x, ...)` did not return a valid `mcmc` or `mcmc.list` object.",
       call. = FALSE)
}

#' Gelman-Rubin diagnostic for BayesODT posterior samples
#'
#' Compute Gelman-Rubin convergence diagnostics for selected posterior samples
#' from a `"BayesODT"` object using \pkg{coda}.
#'
#' @param x A fitted object of class `"BayesODT"`.
#' @param param Character string specifying which parameter to evaluate.
#'   One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s).
#'   For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`,
#'   only `index = 1` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object.
#'   If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object.
#'   If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#' @param autoburnin Logical; passed to `coda::gelman.diag()`.
#' @param multivariate Logical; passed to `coda::gelman.diag()`.
#' @param transform Logical; passed to `coda::gelman.diag()`.
#' @param confidence Confidence level for the upper confidence limit.
#' @param ... Not used.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{parameter}{Parameter name.}
#'   \item{psrf}{Point estimate of the potential scale reduction factor.}
#'   \item{upper_ci}{Upper confidence limit of the PSRF.}
#' }
#'
#' @export
gelman_diag <- function(x,
                        param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                        index = 1L,
                        start = NULL,
                        thin = NULL,
                        autoburnin = FALSE,
                        multivariate = FALSE,
                        transform = FALSE,
                        confidence = 0.95,
                        ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `gelman_diag()`. Please install it first.",
         call. = FALSE)
  }

  mcmc_obj <- as.mcmc(
    x,
    param = param,
    index = index,
    start = start,
    thin = thin
  )

  if (!inherits(mcmc_obj, "mcmc.list")) {
    stop("`gelman_diag()` requires multiple chains. Fit the model with `n.chains >= 2`.",
         call. = FALSE)
  }

  if (length(mcmc_obj) < 2L) {
    stop("`gelman_diag()` requires at least two chains.", call. = FALSE)
  }

  gd <- coda::gelman.diag(
    mcmc_obj,
    confidence = confidence,
    transform = transform,
    autoburnin = autoburnin,
    multivariate = multivariate
  )

  out <- data.frame(
    parameter = rownames(gd$psrf),
    psrf = gd$psrf[, 1],
    upper_ci = gd$psrf[, 2],
    row.names = NULL
  )

  if (isTRUE(multivariate) && !is.null(gd$mpsrf)) {
    attr(out, "mpsrf") <- gd$mpsrf
  }

  out
}


#' Heidelberger-Welch diagnostic for BayesODT posterior samples
#'
#' Compute Heidelberger-Welch convergence diagnostics for selected posterior
#' samples from a `"BayesODT"` object using \pkg{coda}.
#'
#' @param x A fitted object of class `"BayesODT"`.
#' @param param Character string specifying which parameter to evaluate.
#'   One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s).
#'   For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`,
#'   only `index = 1` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object.
#'   If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object.
#'   If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#' @param eps Target value for the ratio of half-width to the sample mean.
#' @param pvalue Significance level used in the stationarity test.
#' @param ... Not used.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{chain}{Chain index.}
#'   \item{parameter}{Parameter name.}
#'   \item{stationarity_test}{Result of the stationarity test.}
#'   \item{start_iteration}{Starting iteration used after truncation.}
#'   \item{stationarity_pvalue}{P-value of the stationarity test.}
#'   \item{halfwidth_test}{Result of the half-width test.}
#'   \item{mean_estimate}{Estimated posterior mean.}
#'   \item{halfwidth}{Estimated half-width.}
#' }
#'
#' @export
heidel_diag <- function(x,
                        param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                        index = 1L,
                        start = NULL,
                        thin = NULL,
                        eps = 0.1,
                        pvalue = 0.05,
                        ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `heidel_diag()`. Please install it first.",
         call. = FALSE)
  }

  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0) {
    stop("`eps` must be a single positive number.", call. = FALSE)
  }

  if (!is.numeric(pvalue) || length(pvalue) != 1L || is.na(pvalue) ||
      pvalue <= 0 || pvalue >= 1) {
    stop("`pvalue` must be a single number in (0, 1).", call. = FALSE)
  }

  mcmc_obj <- as.mcmc(
    x,
    param = param,
    index = index,
    start = start,
    thin = thin
  )

  .format_heidel_output <- function(out, chain_id) {
    out_mat <- unclass(out)

    out_df <- data.frame(
      chain = chain_id,
      parameter = rownames(out_mat),
      out_mat,
      row.names = NULL,
      check.names = FALSE
    )

    names(out_df) <- c(
      "chain",
      "parameter",
      "stationarity_test",
      "start_iteration",
      "stationarity_pvalue",
      "halfwidth_test",
      "mean_estimate",
      "halfwidth"
    )

    out_df
  }

  # single chain
  if (inherits(mcmc_obj, "mcmc")) {
    out <- coda::heidel.diag(mcmc_obj, eps = eps, pvalue = pvalue)
    return(.format_heidel_output(out, chain_id = 1L))
  }

  # multiple chains
  if (inherits(mcmc_obj, "mcmc.list")) {
    out_list <- lapply(seq_along(mcmc_obj), function(i) {
      out <- coda::heidel.diag(mcmc_obj[[i]], eps = eps, pvalue = pvalue)
      .format_heidel_output(out, chain_id = i)
    })

    return(do.call(rbind, out_list))
  }

  stop("`as.mcmc(x, ...)` did not return a valid `mcmc` or `mcmc.list` object.",
       call. = FALSE)
}

#' Density plot for BayesODT posterior samples
#'
#' Plot posterior density estimates for selected posterior samples from a
#' `"BayesODT"` object using \pkg{coda}.
#'
#' @param x A fitted object of class `"BayesODT"`.
#' @param param Character string specifying which parameter to plot.
#'   One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s).
#'   For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`,
#'   only `index = 1` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object.
#'   If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object.
#'   If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#' @param ask Logical; if `TRUE`, ask before moving to a new figure.
#' @param ... Further graphical arguments passed to `coda::densplot()`.
#'
#' @return Invisibly returns the created `"mcmc"` or `"mcmc.list"` object.
#' @export
plot_density <- function(x,
                         param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                         index = 1L,
                         start = NULL,
                         thin = NULL,
                         ask = FALSE,
                         ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `plot_density()`. Please install it first.",
         call. = FALSE)
  }

  mcmc_obj <- as.mcmc(
    x,
    param = param,
    index = index,
    start = start,
    thin = thin
  )

  p <- coda::nvar(mcmc_obj)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (isTRUE(ask)) {
    old_ask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(old_ask), add = TRUE)
  }

  if (p == 1L) {
    coda::densplot(mcmc_obj, ...)
  } else {
    nr <- ceiling(sqrt(p))
    nc <- ceiling(p / nr)
    graphics::par(mfrow = c(nr, nc))

    for (j in seq_len(p)) {
      coda::densplot(
        mcmc_obj[, j],
        main = paste("Density of", coda::varnames(mcmc_obj)[j]),
        ...
      )
    }
  }

  invisible(mcmc_obj)
}


#' Autocorrelation plot for BayesODT posterior samples
#'
#' Plot autocorrelation functions for selected posterior samples from a
#' `"BayesODT"` object using \pkg{coda}.
#'
#' @param x A fitted object of class `"BayesODT"`.
#' @param param Character string specifying which parameter to plot.
#'   One of `"A"`, `"a"`, `"tau_a"`, `"b"`, `"mu_b"`, `"tau_b"`, or `"u"`.
#' @param index Integer index (or vector of indices) of the selected component(s).
#'   For scalar hyperparameters such as `"tau_a"`, `"mu_b"`, and `"tau_b"`,
#'   only `index = 1` is valid.
#' @param start Optional starting iteration for the `coda::mcmc` object.
#'   If `NULL`, uses 1.
#' @param thin Optional thinning interval for the `coda::mcmc` object.
#'   If `NULL`, uses `x$mcmc_info$n.thin` when available, otherwise 1.
#' @param ... Further graphical arguments passed to `coda::autocorr.plot()`.
#'
#' @return Invisibly returns the created `"mcmc"` or `"mcmc.list"` object.
#' @export
plot_autocorr <- function(x,
                          param = c("A", "a", "tau_a", "b", "mu_b", "tau_b", "u"),
                          index = 1L,
                          start = NULL,
                          thin = NULL,
                          ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for `plot_autocorr()`. Please install it first.",
         call. = FALSE)
  }

  mcmc_obj <- as.mcmc(
    x,
    param = param,
    index = index,
    start = start,
    thin = thin
  )

  coda::autocorr.plot(mcmc_obj, ...)
  invisible(mcmc_obj)
}
