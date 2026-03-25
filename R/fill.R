# Internal helper to fill missing initial values and construct chain-specific MCMC initial values.
# Builds one default template, determines whether a default A must be generated, creates default initial values for all chains, and overlays user-supplied values onto the corresponding chain defaults.
.fill_init <- function(init, w, D, n.chains,
                       sd_z = 0.05, sd_z0 = 0.05, sd_a = 0.05, sd_b = 0.05, sd_u = 0.05, sd_A = 0.02) {
  K <- length(unique(w[!is.na(w)]))
  m <- nrow(w)
  n <- ncol(w)

  .is_named_init_list <- function(x) {
    is.list(x) && (length(x) == 0L || (!is.null(names(x)) && all(names(x) != "")))
  }

  .normalize_init_to_chain_list <- function(init, n.chains) {
    out <- vector("list", n.chains)
    if (is.null(init)) {
      return(out)
    }
    # Single named init list: use only for chain 1, others remain empty
    if (.is_named_init_list(init)) {
      out[[1]] <- init
      return(out)
    }
    # Chain-specific init lists
    for (i in seq_along(init)) {
      out[[i]] <- init[[i]]
    }
    out
  }

  .make_default_A <- function(w, D, K, m, n) {
    if (!requireNamespace("ordinal", quietly = TRUE)) stop("Package 'ordinal' is required to construct default cutpoint initial values via `clmm()`. Please install it first.",
                                                           call. = FALSE)
    dat <- data.frame(y = as.vector(w), rater = factor(rep(seq_len(n), each = m)), patient = factor(rep(seq_len(m), times = n)))

    if (!is.null(D)) {
      dat$D <- rep(D, times = n)
    }

    dat <- dat[!is.na(dat$y), , drop = FALSE]
    dat$y <- ordered(dat$y, levels = seq_len(K))

    fit.cut <- if (!is.null(D)) {ordinal::clmm(y ~ D + (1 | rater) + (1 | patient),
                                               data = dat,
                                               link = "probit",
                                               Hess = TRUE,
                                               control = ordinal::clmm.control(maxIter = 100, maxLineIter = 100, gradTol = 1e-4),
                                               threshold = "flexible")
    } else {ordinal::clmm(y ~ 1 + (1 | rater) + (1 | patient),
                          data = dat,
                          link = "probit",
                          Hess = TRUE,
                          control = ordinal::clmm.control(maxIter = 100, maxLineIter = 100, gradTol = 1e-4),
                          threshold = "flexible")
    }
    as.numeric(fit.cut$alpha)
  }

  .make_template_no_A <- function(w, m, n) {
    z <- matrix(1, m, n)
    if (anyNA(w)) z[is.na(w)] <- NA_real_

    list(z = z,
         z0 = rep(0, m),
         a = stats::rnorm(n, 0, 1),
         b = stats::rexp(n, rate = 0.1),
         u = stats::rnorm(m, 0, 1))
  }

  .count_supplied_A <- function(chain_init) {
    sum(vapply(chain_init, function(x) !is.null(x) && !is.null(x$A), logical(1)))
  }

  .extract_first_A <- function(chain_init) {
    for (i in seq_along(chain_init)) {
      if (!is.null(chain_init[[i]]) && !is.null(chain_init[[i]]$A)) {
        return(chain_init[[i]]$A)
      }
    }
    NULL
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
    }

    if (!is.null(out$u)) {
      out$u <- out$u + stats::rnorm(length(out$u), 0, sd_u)
    }

    if (!is.null(out$A)) {
      out$A <- sort(out$A + stats::rnorm(length(out$A), 0, sd_A))
      if (length(out$A) >= 2L) {
        for (j in 2:length(out$A)) {
          if (out$A[j] <= out$A[j - 1L]) {
            out$A[j] <- out$A[j - 1L] + 1e-6
          }
        }
      }
    }

    out
  }

  .overlay_init <- function(default_init, user_init, w) {
    out <- default_init

    if (is.null(user_init) || length(user_init) == 0L) {
      return(out)
    }

    nm <- names(user_init)
    if (is.null(nm)) nm <- character(0)

    out[nm] <- user_init[nm]

    # Always keep missing positions of z as NA
    if (!is.null(out$z) && anyNA(w)) {
      out$z[is.na(w)] <- NA_real_
    }
    out
  }

  # Step 1: normalize user init to a chain list
  chain_init <- .normalize_init_to_chain_list(init, n.chains)

  # Step 2: build one default template without A
  template <- .make_template_no_A(w, m, n)

  # Step 3: determine the base A
  # If A is supplied for all chains, use the first supplied A as the template A.
  # Otherwise, generate A via clmm().
  n_A_supplied <- .count_supplied_A(chain_init)

  if (n_A_supplied == n.chains) {
    A_base <- .extract_first_A(chain_init)
  } else {
    A_base <- .make_default_A(w, D, K, m, n)
  }

  template$A <- A_base

  # Step 4: create default chain-specific initial values
  default_chain_init <- vector("list", n.chains)
  default_chain_init[[1]] <- template

  if (n.chains >= 2L) {
    for (i in 2:n.chains) {
      default_chain_init[[i]] <- .perturb_init(template, w,
                                               sd_z, sd_z0, sd_a, sd_b, sd_u, sd_A)
    }
  }

  # Step 5: overlay user-supplied values onto chain defaults
  out <- vector("list", n.chains)
  for (i in seq_len(n.chains)) {
    out[[i]] <- .overlay_init(default_chain_init[[i]], chain_init[[i]], w)
  }

  out
}


# Internal helper to fill missing prior hyperparameters with default values.
# Starts from the default prior list and replaces only the components supplied by the user.
.fill_prior <- function(prior) {
  default_prior <- list(mu_a = 0, tau_a = 1, mu_b = 0, tau_b = 1, delta = 5, gamma = 0.1, tau = 0.1, mu_u = 0, tau_u = 1)
  if (is.null(prior)) return(default_prior)
  user_prior <- names(prior)
  if (is.null(user_prior)) user_prior <- character(0)

  out <- default_prior
  out[user_prior] <- prior[user_prior]
  out
}
