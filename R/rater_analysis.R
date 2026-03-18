#' Group raters based on posterior summaries of bias and magnifier
#'
#' Group raters using rule-based classification based on posterior summaries
#' of rater-specific bias \eqn{a_j} and magnifier \eqn{b_j}.
#'
#' By default, raters are divided into four groups using the posterior means
#' of centered \eqn{a_j} and \eqn{b_j}, relative to user-specified cutoffs.
#'
#' @param x A fitted object of class \code{"BayesODT"}.
#' @param method Grouping rule. Currently only \code{"quadrant"} is supported.
#' @param a_cut Cutoff for bias. If \code{NULL}, uses 0.
#' @param b_cut Cutoff for magnifier. If \code{NULL}, uses the median of posterior
#'   mean magnifiers.
#' @param probs Length-2 numeric vector of posterior probabilities used for
#'   credible intervals in the returned summary.
#' @param labels Character vector of length 4 giving group labels in the order:
#'   \code{c("low bias / low magnifier", "high bias / low magnifier",
#'   "low bias / high magnifier", "high bias / high magnifier")}.
#' @param include_summary Logical; if \code{TRUE}, include posterior summary
#'   columns for \eqn{a_j} and \eqn{b_j}.
#' @param ... Not used.
#'
#' @return A data frame with rater group assignments.
#' @export
group_raters <- function(x,
                         method = c("quadrant"),
                         a_cut = NULL,
                         b_cut = NULL,
                         probs = c(0.025, 0.975),
                         labels = c("low bias / low magnifier",
                                    "high bias / low magnifier",
                                    "low bias / high magnifier",
                                    "high bias / high magnifier"),
                         include_summary = TRUE,
                         ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  method <- match.arg(method)

  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs)) {
    stop("`probs` must be a numeric vector of length 2.", call. = FALSE)
  }

  if (any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must satisfy 0 <= probs[1] < probs[2] <= 1.", call. = FALSE)
  }

  if (!is.character(labels) || length(labels) != 4L) {
    stop("`labels` must be a character vector of length 4.", call. = FALSE)
  }

  draws <- .extract_ab_draws(x)
  A <- draws$A
  B <- draws$B

  n <- ncol(A)

  A_centered <- sweep(A, 1, rowMeans(A), FUN = "-")

  a_mean <- colMeans(A_centered)
  a_lo <- apply(A_centered, 2, stats::quantile, probs = probs[1], names = FALSE)
  a_hi <- apply(A_centered, 2, stats::quantile, probs = probs[2], names = FALSE)

  b_mean <- colMeans(B)
  b_lo <- apply(B, 2, stats::quantile, probs = probs[1], names = FALSE)
  b_hi <- apply(B, 2, stats::quantile, probs = probs[2], names = FALSE)

  if (is.null(a_cut)) a_cut <- 0
  if (is.null(b_cut)) b_cut <- stats::median(b_mean)

  if (!is.numeric(a_cut) || length(a_cut) != 1L || is.na(a_cut)) {
    stop("`a_cut` must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(b_cut) || length(b_cut) != 1L || is.na(b_cut)) {
    stop("`b_cut` must be a single numeric value.", call. = FALSE)
  }

  bias_high <- a_mean >= a_cut
  mag_high  <- b_mean >= b_cut

  group <- character(n)
  group[!bias_high & !mag_high] <- labels[1]
  group[ bias_high & !mag_high] <- labels[2]
  group[!bias_high &  mag_high] <- labels[3]
  group[ bias_high &  mag_high] <- labels[4]

  out <- data.frame(
    rater = seq_len(n),
    group = group,
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_summary)) {
    out$a_mean <- a_mean
    out$a_lo <- a_lo
    out$a_hi <- a_hi
    out$b_mean <- b_mean
    out$b_lo <- b_lo
    out$b_hi <- b_hi
  }

  attr(out, "a_cut") <- a_cut
  attr(out, "b_cut") <- b_cut
  attr(out, "method") <- method
  attr(out, "labels") <- labels

  out
}

#' Cluster raters based on posterior summaries of bias and magnifier
#'
#' Cluster raters using hierarchical clustering on posterior means of centered
#' \eqn{a_j} and \eqn{b_j}. If \code{k} is not supplied, the number of clusters
#' is chosen automatically by maximizing the average silhouette width.
#'
#' @param x A fitted object of class \code{"BayesODT"}.
#' @param k Number of clusters. If \code{NULL}, choose automatically.
#' @param k_max Maximum number of clusters considered when \code{k = NULL}.
#' @param scale Logical; if \code{TRUE}, standardize \eqn{a_j} and \eqn{b_j}
#'   before clustering.
#' @param hclust_method Linkage method passed to \code{stats::hclust()}.
#' @param probs Length-2 numeric vector of posterior probabilities used for
#'   credible intervals in the returned summary.
#' @param ... Not used.
#'
#' @return A data frame with cluster assignments and posterior summaries.
#'   The chosen clustering object is stored in \code{attr(, "fit")}. If
#'   \code{k = NULL}, the selected number of clusters is stored in
#'   \code{attr(, "k_selected")} and the silhouette scores are stored in
#'   \code{attr(, "silhouette_scores")}.
#' @export
cluster_raters <- function(x,
                           k = NULL,
                           k_max = 6L,
                           scale = TRUE,
                           hclust_method = "complete",
                           probs = c(0.025, 0.975),
                           ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required for automatic cluster selection in `cluster_raters()`. Please install it first.",
         call. = FALSE)
  }

  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs)) {
    stop("`probs` must be a numeric vector of length 2.", call. = FALSE)
  }

  if (any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must satisfy 0 <= probs[1] < probs[2] <= 1.", call. = FALSE)
  }

  if (!is.null(k)) {
    if (!is.numeric(k) || length(k) != 1L || is.na(k) ||
        k < 1 || k != as.integer(k)) {
      stop("`k` must be NULL or a positive integer.", call. = FALSE)
    }
    k <- as.integer(k)
  }

  if (!is.numeric(k_max) || length(k_max) != 1L || is.na(k_max) ||
      k_max < 2 || k_max != as.integer(k_max)) {
    stop("`k_max` must be an integer greater than or equal to 2.", call. = FALSE)
  }
  k_max <- as.integer(k_max)

  draws <- .extract_ab_draws(x)
  A <- draws$A
  B <- draws$B

  n <- ncol(A)
  if (n < 2L) {
    stop("At least two raters are required for clustering.", call. = FALSE)
  }

  if (!is.null(k) && k > n) {
    stop("`k` cannot be larger than the number of raters.", call. = FALSE)
  }

  A_centered <- sweep(A, 1, rowMeans(A), FUN = "-")

  a_mean <- colMeans(A_centered)
  a_lo <- apply(A_centered, 2, stats::quantile, probs = probs[1], names = FALSE)
  a_hi <- apply(A_centered, 2, stats::quantile, probs = probs[2], names = FALSE)

  b_mean <- colMeans(B)
  b_lo <- apply(B, 2, stats::quantile, probs = probs[1], names = FALSE)
  b_hi <- apply(B, 2, stats::quantile, probs = probs[2], names = FALSE)

  X <- cbind(a_mean = a_mean, b_mean = b_mean)
  X_clust <- if (isTRUE(scale)) base::scale(X) else X

  d <- stats::dist(X_clust)
  hc <- stats::hclust(d, method = hclust_method)

  silhouette_scores <- NULL

  if (is.null(k)) {
    k_candidates <- 2L:min(k_max, n - 1L)

    if (length(k_candidates) == 0L) {
      k <- 1L
    } else {
      silhouette_scores <- vapply(k_candidates, function(kk) {
        cl <- stats::cutree(hc, k = kk)
        sil <- cluster::silhouette(cl, dmatrix = as.matrix(d))
        mean(sil[, "sil_width"])
      }, numeric(1))

      k <- k_candidates[which.max(silhouette_scores)]
      names(silhouette_scores) <- paste0("k=", k_candidates)
    }
  }

  cluster_id <- stats::cutree(hc, k = k)

  out <- data.frame(
    rater = seq_len(n),
    cluster = factor(cluster_id, levels = sort(unique(cluster_id))),
    a_mean = a_mean,
    a_lo = a_lo,
    a_hi = a_hi,
    b_mean = b_mean,
    b_lo = b_lo,
    b_hi = b_hi,
    stringsAsFactors = FALSE
  )

  attr(out, "method") <- "hclust"
  attr(out, "k_selected") <- k
  attr(out, "k_input") <- k
  attr(out, "k_max") <- k_max
  attr(out, "scale") <- scale
  attr(out, "fit") <- hc
  attr(out, "silhouette_scores") <- silhouette_scores

  out
}

#' Plot raters in the bias-magnifier plane
#'
#' Plot posterior summaries of rater-specific bias \eqn{a_j} and magnifier
#' \eqn{b_j} in a two-dimensional map. Each point corresponds to one rater.
#'
#' The \eqn{a_j} values are row-wise centered, following the convention used
#' in \code{plot_ABU()}.
#'
#' @param x A fitted object of class \code{"BayesODT"}.
#' @param probs Length-2 numeric vector of posterior probabilities used for
#'   credible intervals.
#' @param show_labels Logical; if \code{TRUE}, label each point by rater index.
#' @param show_intervals Logical; if \code{TRUE}, draw horizontal and vertical
#'   credible intervals for \eqn{a_j} and \eqn{b_j}.
#' @param color_by Coloring rule. One of \code{"none"}, \code{"group"},
#'   or \code{"cluster"}.
#' @param a_cut Cutoff for bias grouping. Used only when
#'   \code{color_by = "group"}. If \code{NULL}, uses 0.
#' @param b_cut Cutoff for magnifier grouping. Used only when
#'   \code{color_by = "group"}. If \code{NULL}, uses the median of posterior
#'   mean magnifiers.
#' @param group_labels Character vector of length 4 giving group labels in the
#'   order used by \code{group_raters()}.
#' @param group_cols Character vector of length 4 specifying colors for groups.
#' @param cluster_k Number of clusters used when \code{color_by = "cluster"}.
#'   If \code{NULL}, the number of clusters is chosen automatically.
#' @param cluster_k_max Maximum number of clusters considered when
#'   \code{cluster_k = NULL}.
#' @param cluster_scale Logical; if \code{TRUE}, standardize variables before clustering.
#' @param cluster_cols Optional vector of colors for clusters. If \code{NULL},
#'   colors are generated automatically.
#' @param legend_pos Position for the legend.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Plot title.
#' @param pch Plotting symbol.
#' @param cex Point size.
#' @param ... Further graphical arguments passed to \code{plot()}.
#'
#' @return Invisibly returns a data frame of posterior summaries used in the plot.
#' @export
plot_rater_map <- function(x,
                           probs = c(0.025, 0.975),
                           show_labels = TRUE,
                           show_intervals = FALSE,
                           color_by = c("none", "group", "cluster"),
                           a_cut = NULL,
                           b_cut = NULL,
                           group_labels = c("low bias / low magnifier",
                                            "high bias / low magnifier",
                                            "low bias / high magnifier",
                                            "high bias / high magnifier"),
                           group_cols = c("gray40", "tomato", "dodgerblue3", "darkgreen"),
                           cluster_k = NULL,
                           cluster_k_max = 6L,
                           cluster_scale = TRUE,
                           cluster_cols = NULL,
                           legend_pos = "topright",
                           xlab = expression(a[j] ~ "(bias)"),
                           ylab = expression(b[j] ~ "(magnifier)"),
                           main = "Scatter plot of Rater",
                           pch = 19,
                           cex = 1,
                           ...) {
  if (!inherits(x, "BayesODT")) {
    stop("`x` must be an object of class 'BayesODT'.", call. = FALSE)
  }

  color_by <- match.arg(color_by)

  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs)) {
    stop("`probs` must be a numeric vector of length 2.", call. = FALSE)
  }

  if (any(probs < 0 | probs > 1) || probs[1] >= probs[2]) {
    stop("`probs` must satisfy 0 <= probs[1] < probs[2] <= 1.", call. = FALSE)
  }

  if (!is.character(group_labels) || length(group_labels) != 4L) {
    stop("`group_labels` must be a character vector of length 4.", call. = FALSE)
  }

  if (!is.character(group_cols) || length(group_cols) != 4L) {
    stop("`group_cols` must be a character vector of length 4.", call. = FALSE)
  }

  draws <- .extract_ab_draws(x)
  A <- draws$A
  B <- draws$B

  n <- ncol(A)

  A_centered <- sweep(A, 1, rowMeans(A), FUN = "-")

  a_mean <- colMeans(A_centered)
  a_lo <- apply(A_centered, 2, stats::quantile, probs = probs[1], names = FALSE)
  a_hi <- apply(A_centered, 2, stats::quantile, probs = probs[2], names = FALSE)

  b_mean <- colMeans(B)
  b_lo <- apply(B, 2, stats::quantile, probs = probs[1], names = FALSE)
  b_hi <- apply(B, 2, stats::quantile, probs = probs[2], names = FALSE)

  xlim <- range(c(a_lo, a_hi))
  ylim <- range(c(b_lo, b_hi))

  point_col <- rep(1, n)
  legend_labels <- NULL
  legend_cols <- NULL
  extra_info <- NULL

  if (color_by == "group") {
    extra_info <- group_raters(
      x,
      a_cut = a_cut,
      b_cut = b_cut,
      probs = probs,
      labels = group_labels,
      include_summary = FALSE
    )

    f <- factor(extra_info$group, levels = group_labels)
    point_col <- group_cols[f]
    legend_labels <- group_labels
    legend_cols <- group_cols
  }

  if (color_by == "cluster") {
    extra_info <- cluster_raters(
      x,
      k = cluster_k,
      k_max = cluster_k_max,
      scale = cluster_scale,
      probs = probs
    )

    levs <- levels(extra_info$cluster)

    if (is.null(cluster_cols)) {
      cluster_cols <- grDevices::hcl.colors(length(levs), palette = "Dark 3")
    }

    if (!is.character(cluster_cols) || length(cluster_cols) < length(levs)) {
      stop("`cluster_cols` must be a character vector with at least as many elements as the number of clusters.",
           call. = FALSE)
    }

    f <- factor(extra_info$cluster, levels = levs)
    point_col <- cluster_cols[f]
    legend_labels <- paste("Cluster", levs)
    legend_cols <- cluster_cols[seq_along(levs)]
  }

  graphics::plot(a_mean, b_mean,
                 xlim = xlim,
                 ylim = ylim,
                 xlab = xlab,
                 ylab = ylab,
                 main = main,
                 pch = pch,
                 cex = cex,
                 col = point_col,
                 ...)

  graphics::abline(v = 0, lty = 2, col = "gray")
  graphics::abline(h = 0, lty = 2, col = "gray")

  if (isTRUE(show_intervals)) {
    graphics::segments(x0 = a_lo, y0 = b_mean, x1 = a_hi, y1 = b_mean, col = point_col)
    graphics::segments(x0 = a_mean, y0 = b_lo, x1 = a_mean, y1 = b_hi, col = point_col)
  }

  if (isTRUE(show_labels)) {
    graphics::text(a_mean, b_mean,
                   labels = seq_len(n),
                   pos = 3,
                   cex = 0.8,
                   col = point_col)
  }

  if (!is.null(legend_labels)) {
    graphics::legend(legend_pos,
                     legend = legend_labels,
                     col = legend_cols,
                     pch = pch,
                     bty = "n")
  }

  out <- data.frame(
    rater = seq_len(n),
    a_mean = a_mean,
    a_lo = a_lo,
    a_hi = a_hi,
    b_mean = b_mean,
    b_lo = b_lo,
    b_hi = b_hi
  )

  if (color_by == "group") {
    out$group <- extra_info$group
  }

  if (color_by == "cluster") {
    out$cluster <- extra_info$cluster
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
