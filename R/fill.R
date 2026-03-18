# Fill missing components of init and construct chain-specific initial values:
# if init is NULL, generate default initial values
# if n.chains = 1, return one completed init list
# if n.chains > 1 and init is a single named list, treat it as a template and generate perturbed chain-specific init lists
# if n.chains > 1 and init is a list of init lists, fill missing components separately for each chain
# observed entries of z are perturbed, while missing positions remain NA
# A is perturbed and then adjusted to remain strictly increasing
.fill_init <- function(init, w, D, n.chains,
                       sd_z = 0.05, sd_z0 = 0.05, sd_a = 0.05, sd_b = 0.05,
                       sd_u = 0.05, sd_A = 0.02) {
  K <- length(sort(unique(w[!is.na(w)])))
  m <- nrow(w)
  n <- ncol(w)

  valid_names <- c("z", "z0", "a", "b", "u", "A")

  .is_named_init_list <- function(x) {
    is.list(x) && (length(x) == 0L || (!is.null(names(x)) && all(names(x) != "")))
  }

  .is_missing_chain_init <- function(x) {
    is.null(x) || (is.list(x) && length(x) == 0L)
  }

  .check_unknown_names <- function(x) {
    if (is.null(x) || length(x) == 0L) return(invisible(NULL))
    nm <- names(x)
    if (is.null(nm)) return(invisible(NULL))
    bad <- setdiff(nm, valid_names)
    if (length(bad) > 0L) {
      stop(
        sprintf("Unknown init component(s): %s", paste(bad, collapse = ", ")),
        call. = FALSE
      )
    }
    invisible(NULL)
  }

  .validate_A <- function(A, K) {
    if (is.null(A)) return(invisible(NULL))
    if (!is.numeric(A) || anyNA(A)) {
      stop("`A` must be a numeric vector without missing values.", call. = FALSE)
    }
    if (length(A) != (K - 1L)) {
      stop(sprintf("`A` must have length K-1 = %d.", K - 1L), call. = FALSE)
    }
    if (any(diff(A) <= 0)) {
      stop("`A` must be strictly increasing.", call. = FALSE)
    }
    invisible(NULL)
  }

  .validate_one_init <- function(init_one, m, n, K) {
    if (is.null(init_one) || length(init_one) == 0L) return(invisible(NULL))

    .check_unknown_names(init_one)

    if (!is.null(init_one$z)) {
      if (!is.matrix(init_one$z) || !all(dim(init_one$z) == c(m, n))) {
        stop(sprintf("`z` must be a %d x %d matrix.", m, n), call. = FALSE)
      }
    }
    if (!is.null(init_one$z0)) {
      if (!is.numeric(init_one$z0) || length(init_one$z0) != m) {
        stop(sprintf("`z0` must be a numeric vector of length %d.", m), call. = FALSE)
      }
    }
    if (!is.null(init_one$a)) {
      if (!is.numeric(init_one$a) || length(init_one$a) != n) {
        stop(sprintf("`a` must be a numeric vector of length %d.", n), call. = FALSE)
      }
    }
    if (!is.null(init_one$b)) {
      if (!is.numeric(init_one$b) || length(init_one$b) != n) {
        stop(sprintf("`b` must be a numeric vector of length %d.", n), call. = FALSE)
      }
    }
    if (!is.null(init_one$u)) {
      if (!is.numeric(init_one$u) || length(init_one$u) != m) {
        stop(sprintf("`u` must be a numeric vector of length %d.", m), call. = FALSE)
      }
    }
    if (!is.null(init_one$A)) {
      .validate_A(init_one$A, K)
    }

    invisible(NULL)
  }

  .make_default_A <- function(w, D, K, m, n) {
    if (!requireNamespace("ordinal", quietly = TRUE)) {
      stop(
        "Package 'ordinal' is required to construct default cutpoint initial values via `clmm()`. Please install it first.",
        call. = FALSE
      )
    }

    dat <- data.frame(
      y = as.vector(w),
      rater = factor(rep(seq_len(n), each = m)),
      patient = factor(rep(seq_len(m), times = n))
    )

    if (!is.null(D)) {
      dat$D <- rep(D, times = n)
    }

    dat <- dat[!is.na(dat$y), , drop = FALSE]
    dat$y <- ordered(dat$y, levels = seq_len(K))

    fit.cut <- if (!is.null(D)) {
      ordinal::clmm(
        y ~ D + (1 | rater) + (1 | patient),
        data = dat,
        link = "probit",
        Hess = TRUE,
        control = ordinal::clmm.control(
          maxIter = 100,
          maxLineIter = 100,
          gradTol = 1e-4
        ),
        threshold = "flexible"
      )
    } else {
      ordinal::clmm(
        y ~ 1 + (1 | rater) + (1 | patient),
        data = dat,
        link = "probit",
        Hess = TRUE,
        control = ordinal::clmm.control(
          maxIter = 100,
          maxLineIter = 100,
          gradTol = 1e-4
        ),
        threshold = "flexible"
      )
    }

    as.numeric(fit.cut$alpha)
  }

  .make_template_no_A <- function(w, m, n) {
    z <- matrix(1, m, n)
    if (anyNA(w)) z[is.na(w)] <- NA_real_

    list(
      z = z,
      z0 = rep(0, m),
      a = stats::rnorm(n, 0, 1),
      b = stats::rexp(n, rate = 0.1),
      u = stats::rnorm(m, 0, 1)
    )
  }

  .fill_missing_from_template <- function(init_one, template) {
    out <- template

    if (is.null(init_one) || length(init_one) == 0L) {
      return(out)
    }

    nm <- names(init_one)
    if (is.null(nm)) nm <- character(0)

    out[nm] <- init_one[nm]

    if (!is.null(out$z) && anyNA(w)) {
      out$z[is.na(w)] <- NA_real_
    }

    out
  }

  .perturb_init <- function(init_one, w,
                            sd_z, sd_z0, sd_a, sd_b, sd_u, sd_A) {
    out <- init_one

    if (!is.null(out$z)) {
      obs <- !is.na(w)
      out$z[obs] <- out$z[obs] + stats::rnorm(sum(obs), 0, sd_z)
      if (anyNA(w)) out$z[is.na(w)] <- NA_real_
    }

    if (!is.null(out$z0)) {
      out$z0 <- out$z0 + stats::rnorm(length(out$z0), 0, sd_z0)
    }

    if (!is.null(out$a)) {
      out$a <- out$a + stats::rnorm(length(out$a), 0, sd_a)
    }

    if (!is.null(out$b)) {
      out$b <- out$b + stats::rnorm(length(out$b), 0, sd_b)
      out$b <- pmax(out$b, 1e-6)
    }

    if (!is.null(out$u)) {
      out$u <- out$u + stats::rnorm(length(out$u), 0, sd_u)
    }

    if (!is.null(out$A)) {
      out$A <- sort(out$A + stats::rnorm(length(out$A), 0, sd_A))
      if (length(out$A) >= 2L) {
        for (j in 2:length(out$A)) {
          if (out$A[j] <= out$A[j - 1]) {
            out$A[j] <- out$A[j - 1] + 1e-6
          }
        }
      }
    }

    out
  }

  .extract_first_A <- function(init_obj, K) {
    if (is.null(init_obj)) return(NULL)

    if (.is_named_init_list(init_obj)) {
      .validate_one_init(init_obj, m, n, K)
      if (!is.null(init_obj$A)) return(init_obj$A)
      return(NULL)
    }

    for (i in seq_along(init_obj)) {
      init_i <- init_obj[[i]]
      if (is.null(init_i) || length(init_i) == 0L) next
      .validate_one_init(init_i, m, n, K)
      if (!is.null(init_i$A)) return(init_i$A)
    }

    NULL
  }

  .normalize_init_to_chain_list <- function(init, n.chains) {
    out <- vector("list", n.chains)

    if (is.null(init)) {
      return(out)
    }

    if (.is_named_init_list(init)) {
      out[[1]] <- init
      return(out)
    }

    if (!is.list(init)) {
      stop("`init` must be NULL, a named list, or a list of chain-specific init lists.", call. = FALSE)
    }

    if (length(init) > n.chains) {
      stop("Length of chain-specific `init` cannot exceed `n.chains`.", call. = FALSE)
    }

    for (i in seq_along(init)) {
      out[[i]] <- init[[i]]
    }

    out
  }

  chain_init <- .normalize_init_to_chain_list(init, n.chains)

  for (i in seq_len(n.chains)) {
    .validate_one_init(chain_init[[i]], m, n, K)
  }

  template <- .make_template_no_A(w, m, n)

  first_A <- .extract_first_A(init, K)
  if (is.null(first_A)) {
    template$A <- .make_default_A(w, D, K, m, n)
  } else {
    template$A <- first_A
  }

  out <- vector("list", n.chains)
  base_init <- NULL

  if (n.chains == 1L) {
    out[[1]] <- .fill_missing_from_template(chain_init[[1]], template)
    return(out)
  }

  for (i in seq_len(n.chains)) {
    init_i <- chain_init[[i]]

    if (.is_missing_chain_init(init_i)) {
      if (is.null(base_init)) {
        base_init <- template
      }

      out[[i]] <- .perturb_init(
        base_init, w,
        sd_z, sd_z0, sd_a, sd_b, sd_u, sd_A
      )
    } else {
      out[[i]] <- .fill_missing_from_template(init_i, template)

      if (is.null(base_init)) {
        base_init <- out[[i]]
      }
    }
  }

  out
}



# Fill missing prior hyperparmeters with default values:
# if prior is NULL, return the full default prior list
# if prior is partially specified, replace only the corresponding default values
.fill_prior <- function(prior) {
  default_prior <- list(mu_a = 0,
                        tau_a = 1,
                        mu_b = 0,
                        tau_b = 1,
                        delta = 5,
                        gamma = 0.1,
                        tau = 0.1,
                        mu_u = 0,
                        tau_u = 1)
  if (is.null(prior)) return(default_prior)
  user_prior <- names(prior)
  if (is.null(user_prior)) user_prior <- character(0)

  out <- default_prior
  out[user_prior] <- prior[user_prior]
  out
}
