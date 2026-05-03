#' Group raters based on posterior summaries of bias and magnifier
#'
#' Group raters using rule-based classification based on posterior summaries
#' of rater-specific bias \eqn{a_j} and magnifier \eqn{b_j}.
#'
#' By default, raters are divided into four groups using the posterior means
#' of centered \eqn{a_j} and \eqn{b_j}.
#'
#' Group codes are defined as:
#' \describe{
#'   \item{1}{low bias / low magnifier}
#'   \item{2}{low bias / high magnifier}
#'   \item{3}{high bias / low magnifier}
#'   \item{4}{high bias / high magnifier}
#' }
#'
#' @param x A fitted object of class \code{"BayesODT"}.
#' @param method Grouping rule. Currently only \code{"quadrant"} is supported.
#' @param a_cut Cutoff for bias. If \code{NULL}, uses 0.
#' @param b_cut Cutoff for magnifier. If \code{NULL}, uses the mean of posterior
#'   mean magnifiers.
#' @param probs Length-2 numeric vector of posterior probabilities used for
#'   credible intervals in the returned summary.
#' @param labels Character vector of length 4 giving group labels in the order:
#'   \code{c("low bias / low magnifier", "low bias / high magnifier",
#'   "high bias / low magnifier", "high bias / high magnifier")}.
#' @param include_summary Logical; if \code{TRUE}, include posterior summary
#'   columns for \eqn{a_j} and \eqn{b_j}.
#' @param ... Not used.
#'
#' @return A data frame with rater group assignments.
#' @export
group_raters <- function(x,
                         probs = c(0.025, 0.975),
                         ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs)) {
    stop("`probs` must be a numeric vector of length 2.", call. = FALSE)
  }

  if (any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must satisfy 0 <= probs[1] < probs[2] <= 1.", call. = FALSE)
  }

  draws <- .extract_ab_draws(x)
  A <- draws$A
  B <- draws$B

  n <- ncol(A)

  a_mean <- colMeans(A)
  a_lo <- apply(A, 2, stats::quantile, probs = probs[1], names = FALSE)
  a_hi <- apply(A, 2, stats::quantile, probs = probs[2], names = FALSE)

  b_mean <- colMeans(B)
  b_lo <- apply(B, 2, stats::quantile, probs = probs[1], names = FALSE)
  b_hi <- apply(B, 2, stats::quantile, probs = probs[2], names = FALSE)

  a_cut <- 0
  b_cut <- mean(b_mean)

  bias_high <- a_mean >= a_cut
  mag_high  <- b_mean >= b_cut

  group_code <- integer(n)

  group_code[!bias_high & !mag_high] <- 1L  # low bias / low magnifier
  group_code[!bias_high &  mag_high] <- 2L  # low bias / high magnifier
  group_code[ bias_high & !mag_high] <- 3L  # high bias / low magnifier
  group_code[ bias_high &  mag_high] <- 4L  # high bias / high magnifier


  group_labels = c("low bias / low magnifier",
                   "low bias / high magnifier",
                   "high bias / low magnifier",
                   "high bias / high magnifier")
  out <- data.frame(
    rater = seq_len(n),
    group = group_code,
    group_label = group_labels[group_code],
    stringsAsFactors = FALSE
  )

  out$a_mean <- a_mean
  out$a_L <- a_lo
  out$a_U <- a_hi
  out$b_mean <- b_mean
  out$b_L <- b_lo
  out$b_U <- b_hi



  group_legend <- data.frame(
    group = 1:4,
    meaning = group_labels,
    stringsAsFactors = FALSE
  )

  attr(out, "a_cut") <- a_cut
  attr(out, "b_cut") <- b_cut
  attr(out, "group_legend") <- group_legend

  class(out) <- c("BayesODT_rater_group", class(out))

  out
}

#' Plot raters in the bias-magnifier plane
#'
#' Plot posterior summaries of rater-specific bias \eqn{a_j} and magnifier
#' \eqn{b_j} in a two-dimensional map. Each point corresponds to one rater.
#'
#' The \eqn{a_j} values are row-wise centered so that the plot represents
#' relative rater bias.
#'
#' If \code{group = TRUE}, raters are colored according to the group assignments
#' returned by \code{group_raters()}.
#'
#' @param x A fitted object of class \code{"BayesODT"}.
#' @param probs Length-2 numeric vector of posterior probabilities used for
#'   credible intervals.
#' @param show_labels Logical; if \code{TRUE}, label each point by rater index.
#' @param show_intervals Logical; if \code{TRUE}, draw horizontal and vertical
#'   credible intervals for \eqn{a_j} and \eqn{b_j}.
#' @param group Logical; if \code{TRUE}, color points by rater group.
#'
#' @return Invisibly returns a data frame of posterior summaries used in the plot.
#' @export
plot_rater_map <- function(x,
                           probs = c(0.025, 0.975),
                           show_intervals = FALSE,
                           group = FALSE) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs)) {
    stop("`probs` must be a numeric vector of length 2.", call. = FALSE)
  }

  if (any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must satisfy 0 <= probs[1] < probs[2] <= 1.", call. = FALSE)
  }

  if (!is.logical(group) || length(group) != 1L || is.na(group)) {
    stop("`group` must be TRUE or FALSE.", call. = FALSE)
  }

  draws <- .extract_ab_draws(x)
  A <- draws$A
  B <- draws$B

  n <- ncol(A)


  a_mean <- colMeans(A)
  a_L <- apply(A, 2, stats::quantile, probs = probs[1], names = FALSE)
  a_U <- apply(A, 2, stats::quantile, probs = probs[2], names = FALSE)

  b_mean <- colMeans(B)
  b_L <- apply(B, 2, stats::quantile, probs = probs[1], names = FALSE)
  b_U <- apply(B, 2, stats::quantile, probs = probs[2], names = FALSE)

  if (isTRUE(show_intervals)) {
    xlim <- range(c(a_L, a_U), finite = TRUE)
    ylim <- range(c(b_L, b_U), finite = TRUE)
  } else {
    xlim <- range(a_mean, finite = TRUE)
    ylim <- range(b_mean, finite = TRUE)
  }

  point_col <- rep(1, n)
  group_info <- NULL

  group_labels <- c(
    "low bias / low magnifier",
    "low bias / high magnifier",
    "high bias / low magnifier",
    "high bias / high magnifier"
  )

  group_cols <- c("gray40", "dodgerblue3", "tomato", "darkgreen")

  if (isTRUE(group)) {
    group_info <- group_raters(x, probs = probs)

    if (!("group" %in% names(group_info))) {
      stop("`group_raters()` must return a column named `group`.", call. = FALSE)
    }

    if (!all(group_info$group %in% 1:4)) {
      stop("`group_raters()` must return group codes 1, 2, 3, or 4.", call. = FALSE)
    }

    point_col <- group_cols[group_info$group]
  }

  graphics::plot(
    a_mean, b_mean,
    xlim = xlim,
    ylim = ylim,
    xlab = expression(a[j] ~ "(bias)"),
    ylab = expression(b[j] ~ "(magnifier)"),
    main = "Scatter Plot of Rater",
    pch = 19,
    col = point_col
  )

  graphics::abline(v = 0, lty = 2, col = "gray")
  graphics::abline(h = mean(b_mean), lty = 2, col = "gray")

  if (isTRUE(show_intervals)) {
    graphics::segments(
      x0 = a_L, y0 = b_mean,
      x1 = a_U, y1 = b_mean,
      col = point_col
    )

    graphics::segments(
      x0 = a_mean, y0 = b_L,
      x1 = a_mean, y1 = b_U,
      col = point_col
    )
  }


  graphics::text(a_mean, b_mean, labels = seq_len(n),
                 pos = 3, cex = 0.8, col = point_col)

  if (isTRUE(group)) {
    graphics::legend(
      "topright",
      legend = paste0(1:4, ": ", group_labels),
      col = group_cols,
      pch = 19,
      bty = "n"
    )
  }

  out <- data.frame(
    rater = seq_len(n),
    a_mean = a_mean,
    a_L = a_L,
    a_U = a_U,
    b_mean = b_mean,
    b_L = b_L,
    b_U = b_U
  )

  if (isTRUE(group)) {
    out$group <- group_info$group
    out$group_label <- group_labels[group_info$group]
  }

  invisible(out)
}

.extract_ab_draws <- function(x) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  # multi-chain preferred
  if (!is.null(x$chains) && is.list(x$chains) && length(x$chains) > 0L) {
    a_list <- lapply(seq_along(x$chains), function(i) {
      samp <- x$chains[[i]]$a
      if (is.null(samp)) {
        stop(sprintf("Chain %d does not contain posterior samples for `a`.", i),
             call. = FALSE)
      }
      if (!is.matrix(samp)) samp <- as.matrix(samp)
      samp
    })

    b_list <- lapply(seq_along(x$chains), function(i) {
      samp <- x$chains[[i]]$b
      if (is.null(samp)) {
        stop(sprintf("Chain %d does not contain posterior samples for `b`.", i),
             call. = FALSE)
      }
      if (!is.matrix(samp)) samp <- as.matrix(samp)
      samp
    })

    # check dimensions across chains
    a_ncol <- vapply(a_list, ncol, integer(1))
    b_ncol <- vapply(b_list, ncol, integer(1))

    if (length(unique(a_ncol)) != 1L) {
      stop("All chains must have the same number of raters in `a`.", call. = FALSE)
    }
    if (length(unique(b_ncol)) != 1L) {
      stop("All chains must have the same number of raters in `b`.", call. = FALSE)
    }
    if (a_ncol[1] != b_ncol[1]) {
      stop("Posterior samples for `a` and `b` must have the same number of raters.",
           call. = FALSE)
    }

    A <- do.call(rbind, a_list)
    B <- do.call(rbind, b_list)

    return(list(A = A, B = B))
  }

  # fallback: old single-chain storage
  if (is.null(x$A) || is.null(x$B)) {
    stop("The fitted object must contain posterior samples for both bias (`a`) and magnifier (`b`).",
         call. = FALSE)
  }

  A <- x$A
  B <- x$B

  if (!is.matrix(A)) A <- as.matrix(A)
  if (!is.matrix(B)) B <- as.matrix(B)

  if (ncol(A) != ncol(B)) {
    stop("`A` and `B` must have the same number of raters (columns).",
         call. = FALSE)
  }

  list(A = A, B = B)
}
