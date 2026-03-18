#' Plot method for BayesODT objects
#'
#' Plot posterior summaries of model parameters or ROC curves from a fitted
#' `"BayesODT"` object.
#'
#' @param x A fitted object of class `"BayesODT"`.
#' @param which Type of plot. One of `"ABU"`, `"ROC"`, or `"Udist"`.
#' @param type For `which = "ROC"`, the ROC type. One of `"individual"`, `"ROC1"`, `"ROC2"`, `"smoothed"`, `"pROC"`, or `"smoothed.pROC"`.
#' @param rater Integer rater index, or an integer vector of rater indices, used when `type = "individual"`. At most 5 raters can be plotted at once.
#' @param h Numeric vector of thresholds used when `type = "smoothed"`.
#' @param probs Length-2 numeric vector of posterior probabilities used for credible intervals in `which = "ABU"`.
#'
#' @return Invisibly returns the plotted object.
#' @export
plot.BayesODT <- function(x,
                          which = c("ABU", "ROC", "Udist"),
                          type = c("individual", "ROC1", "ROC2",
                                   "smoothed", "pROC", "smoothed.pROC"),
                          rater = 1L,
                          h = seq(-5, 5, length.out = 200),
                          probs = c(0.025, 0.975)){
  which <- match.arg(which)

  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (which == "ABU") {
    .plot_ABU_BayesODT(x, probs = probs)
    return(invisible(x))
  }

  if (which == "Udist") {
    .plot_Udist_BayesODT(x)
    return(invisible(x))
  }

  type <- match.arg(type)

  if (type == "individual" && length(rater) > 1L) {
    if (length(rater) > 5L) {
      stop("At most 5 raters can be plotted at once for `type = \"individual\"`.",
           call. = FALSE)
    }
    if (!is.numeric(rater) || anyNA(rater) || any(rater < 1) ||
        any(rater != as.integer(rater))) {
      stop("`rater` must be a vector of positive integer rater indices.",
           call. = FALSE)
    }

    rater <- as.integer(rater)
    cols <- seq_along(rater)
    pchs <- seq_along(rater)

    roc_list <- lapply(rater, function(j) {
      bayes_roc(x, type = "individual", rater = j, h = h)
    })

    main_txt <- paste0("Estimated Individual ROC (Raters ",
                       paste(rater, collapse = ", "), ")")

    plot(roc_list[[1]], col = cols[1], pch = pchs[1], main = main_txt)
    if (length(roc_list) > 1L) {
      for (i in 2:length(roc_list)) {
        plot(roc_list[[i]], add = TRUE, col = cols[i], pch = pchs[i])
      }
    }

    graphics::legend("bottomright",
                     legend = paste("Rater", rater),
                     col = cols,
                     lty = 1,
                     pch = pchs,
                     bty = "n")

    return(invisible(roc_list))
  }

  roc_obj <- bayes_roc(x, type = type, rater = rater, h = h)
  plot(roc_obj)
  invisible(roc_obj)
}

#' Plot method for BayesODT ROC objects
#'
#' Plot ROC curves from an object returned by `bayes_roc()`.
#'
#' @param x An object of class `"BayesODT_roc"`.
#' @param add Logical; if `TRUE`, add the ROC curve to an existing plot.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Plot title.
#' @param col Color for points/lines.
#' @param lwd Line width.
#' @param lty Line type.
#' @param pch Plotting symbol.
#'
#' @return Invisibly returns `x`.
#' @export
plot.BayesODT_roc <- function(x,
                              add = FALSE,
                              xlab = "False Positive Rate",
                              ylab = "True Positive Rate",
                              main = NULL,
                              col = 1,
                              lwd = 1,
                              lty = 1,
                              pch = 1) {
  if (!inherits(x, "BayesODT_roc")) {
    stop("`x` must be an object of class 'BayesODT_roc'.", call. = FALSE)
  }

  if (is.null(main)) {
    main <- switch(
      x$type,
      individual = sprintf("Estimated Individual ROC (Rater %d)", x$rater),
      ROC1 = "Estimated ROC1",
      ROC2 = "Estimated ROC2",
      smoothed = "Smoothed ROC",
      pROC = "Empirical ROC",
      smoothed.pROC = "Smoothed Empirical ROC",
      "ROC"
    )
  }

  is_smooth <- x$type %in% c("smoothed", "pROC", "smoothed.pROC")

  if (!add) {
    if (is_smooth) {
      graphics::plot(x$fpr, x$tpr,
                     type = "l",
                     xlim = c(0, 1), ylim = c(0, 1),
                     xlab = xlab, ylab = ylab,
                     main = main,
                     col = col, lwd = lwd, lty = lty)
    } else {
      graphics::plot(x$fpr, x$tpr,
                     type = "n",
                     xlim = c(0, 1), ylim = c(0, 1),
                     xlab = xlab, ylab = ylab,
                     main = main)
      graphics::lines(x$fpr, x$tpr, col = col, lwd = lwd, lty = lty)
      graphics::points(x$fpr, x$tpr, col = col, pch = pch)
    }
    graphics::abline(0, 1, lty = 2, col = "gray")
  } else {
    if (is_smooth) {
      graphics::lines(x$fpr, x$tpr, col = col, lwd = lwd, lty = lty)
    } else {
      graphics::lines(x$fpr, x$tpr, col = col, lwd = lwd, lty = lty)
      graphics::points(x$fpr, x$tpr, col = col, pch = pch)
    }
  }

  invisible(x)
}

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


# internal helper for ABU plot
.plot_ABU_BayesODT <- function(fit,
                               probs = c(0.025, 0.975)) {
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs)) {
    stop("`probs` must be a numeric vector of length 2.", call. = FALSE)
  }
  if (any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must satisfy 0 <= probs[1] < probs[2] <= 1.", call. = FALSE)
  }

  if (is.null(fit$chains) || !is.list(fit$chains) || length(fit$chains) == 0L) {
    stop("The fitted object has no valid `chains` component.", call. = FALSE)
  }

  a.data <- .combine_chain_samples(fit$chains, "a")
  b.data <- .combine_chain_samples(fit$chains, "b")
  u.data <- .combine_chain_samples(fit$chains, "u")

  needed <- c(a = !is.null(a.data), b = !is.null(b.data), u = !is.null(u.data))
  if (!all(needed)) {
    miss <- names(needed)[!needed]
    stop(sprintf("The fitted object is missing required component(s): %s",
                 paste(miss, collapse = ", ")),
         call. = FALSE)
  }

  n <- ncol(a.data)
  m <- ncol(u.data)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = c(3, 1), mar = c(5, 5.5, 2, 1))

  # a_j : row-wise centering
  a.centered <- sweep(a.data, 1, rowMeans(a.data), FUN = "-")
  a.mean <- colMeans(a.centered)
  a.lo <- apply(a.centered, 2, stats::quantile, probs = probs[1], names = FALSE)
  a.hi <- apply(a.centered, 2, stats::quantile, probs = probs[2], names = FALSE)

  graphics::plot(seq_len(n), a.mean,
                 ylim = c(min(a.lo) - 0.1, max(a.hi) + 0.1),
                 ylab = expression(a[j]),
                 xlab = "rater",
                 cex.lab = 1.2)
  graphics::abline(h = 0, col = "gray")
  graphics::arrows(seq_len(n), a.lo, seq_len(n), a.hi,
                   angle = 90, code = 3, length = 0.05)

  # b_j
  b.mean <- colMeans(b.data)
  b.lo <- apply(b.data, 2, stats::quantile, probs = probs[1], names = FALSE)
  b.hi <- apply(b.data, 2, stats::quantile, probs = probs[2], names = FALSE)

  graphics::plot(seq_len(n), b.mean,
                 ylim = c(min(b.lo) - 0.1, max(b.hi) + 0.1),
                 ylab = expression(b[j]),
                 xlab = "rater",
                 cex.lab = 1.2)
  graphics::abline(h = mean(b.data), col = "gray")
  graphics::arrows(seq_len(n), b.lo, seq_len(n), b.hi,
                   angle = 90, code = 3, length = 0.05)

  # u_i
  u.mean <- colMeans(u.data)
  u.lo <- apply(u.data, 2, stats::quantile, probs = probs[1], names = FALSE)
  u.hi <- apply(u.data, 2, stats::quantile, probs = probs[2], names = FALSE)

  graphics::plot(seq_len(m), u.mean,
                 ylim = c(min(u.lo) - 0.1, max(u.hi) + 0.1),
                 ylab = expression(u[i]),
                 xlab = "subject",
                 cex.lab = 1.2)
  graphics::arrows(seq_len(m), u.lo, seq_len(m), u.hi,
                   angle = 90, code = 3, length = 0.03)

  if (!is.null(fit$D)) {
    idx1 <- which(fit$D == 1)
    if (length(idx1) > 0L) {
      graphics::points(idx1, u.mean[idx1], pch = 16, col = "red")
    }
  }

  invisible(fit)
}

.plot_Udist_BayesODT <- function(fit,
                                 adjust = 1,
                                 xlab = "u",
                                 ylab = "Density",
                                 main = "Estimation of the distribution of patient latent disease severity u") {
  if (!inherits(fit, "BayesODT")) {
    stop("`fit` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (is.null(fit$chains) || !is.list(fit$chains) || length(fit$chains) == 0L) {
    stop("The fitted object has no valid `chains` component.", call. = FALSE)
  }

  u.data <- .combine_chain_samples(fit$chains, "u")
  if (is.null(u.data)) {
    stop("The fitted object is missing required posterior sample `u`.",
         call. = FALSE)
  }

  if (is.null(fit$D)) {
    stop("The fitted object is missing required component `D`.", call. = FALSE)
  }

  if (!is.numeric(fit$D) && !is.logical(fit$D)) {
    stop("`fit$D` must be a binary numeric/logical vector.", call. = FALSE)
  }

  D <- as.integer(fit$D)
  if (anyNA(D) || !all(D %in% c(0L, 1L))) {
    stop("`fit$D` must contain only 0/1 values.", call. = FALSE)
  }

  if (ncol(u.data) != length(D)) {
    stop("Length of `fit$D` must match the number of subjects in posterior sample `u`.",
         call. = FALSE)
  }

  u.mean <- colMeans(u.data)

  idx0 <- which(D == 0L)
  idx1 <- which(D == 1L)

  if (length(idx0) < 2L || length(idx1) < 2L) {
    stop("Need at least two subjects in each disease group to estimate densities.",
         call. = FALSE)
  }

  den0 <- stats::density(u.mean[idx0], adjust = adjust)
  den1 <- stats::density(u.mean[idx1], adjust = adjust)
  denA <- stats::density(u.mean, adjust = adjust)

  xlim <- range(c(den0$x, den1$x, denA$x))
  ylim <- c(0, max(c(den0$y, den1$y, denA$y)))

  graphics::plot(den1,
                 type = "l",
                 xlim = xlim,
                 ylim = ylim,
                 xlab = xlab,
                 ylab = ylab,
                 main = main,
                 lty = 1)
  graphics::lines(den0, lty = 2)
  graphics::lines(denA, lty = 1, col = "gray70")

  graphics::legend("topright",
                   legend = c("D=1", "D=0", "All"),
                   lty = c(1, 2, 1),
                   col = c("black", "black", "gray70"),
                   bty = "n")

  invisible(fit)
}

