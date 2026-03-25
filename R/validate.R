# Internal helper to validate the ordinal rating matrix w.
# Returns TRUE invisibly if all checks pass; otherwise throws an error.
.validate_w <- function(w) {
  if (!is.matrix(w)) stop("w must be a matrix.")
  if (!is.numeric(w)) stop("w must be a numeric (integer-valued, with possible NA).")
  if (any(is.nan(w))) stop("w must not contain NaN.")
  x <- w[!is.na(w)] # Store non-missing entries in x
  if (length(x) == 0L) stop("w has only NA values; cannot determine rating levels.")
  if (any(!is.finite(x))) stop("w must not contain Inf, or -Inf (NA is allowed).")
  if (any(x != floor(x))) stop("All non-missing entries of w must be integers.")
  if (any(x < 1)) stop("All non-missing entries of w must be >= 1.")
  lev <- sort(unique(as.integer(x)))
  if (!identical(lev, seq.int(from = 1L, to = max(lev)))) stop("Non-missing entries of w must contain consecutive integer levels starting at 1 (i.e., 1, 2, ..., K with no gaps).")
  # 8) K(등급 수) 최소 조건
  K <- length(lev)
  if (K < 3) stop(sprintf("w must have at least 3 observed ordinal levels (found K = %d).", K))

  invisible(TRUE)
}



# Internal helper to validate the disease-status vector D.
# Returns TRUE invisibly if all checks pass; otherwise throws an error.
.validate_D <- function(D, w) {
  if (is.null(D)) return(invisible(TRUE))
  if (!is.atomic(D) || !is.null(dim(D))) stop("D must be a vector (no dimensions).")
  if (length(D) != nrow(w)) stop(sprintf("Length of D must equal nrow(w) (got length(D)=%d, nrow(w)=%d).", length(D), nrow(w)))
  if (anyNA(D) || any(is.nan(D)) || any(is.infinite(D))) stop("D must not contain NA/NaN/Inf.")
  if (is.logical(D)) D <- as.integer(D) # Allow logical vectors by converting them to 0/1
  if (!is.numeric(D) || any(D != floor(D))) stop("D must be integer-valued (0/1).")
  if (!all(D %in% c(0L, 1L))) stop("D must contain only 0 and 1.")
  if (!all(c(0L, 1L) %in% D)) stop("D must contain both classes 0 and 1.")

  invisible(TRUE)
}



# Internal helper to validate MCMC tuning parameters.
# Returns TRUE invisibly if all checks pass; otherwise throws an error.
.validate_mcmc <- function(n.iter, n.burnin, n.thin) {
  # Internal helper to validate n.iter, n.burnin, and n.thin
  .check_scalar_integer <- function(x, name) {
    if (length(x) != 1L) stop(sprintf("%s must be a single value.", name))
    if (!is.numeric(x)) stop(sprintf("%s must be numeric.", name))
    if (is.na(x) || is.nan(x) || is.infinite(x)) stop(sprintf("%s must be finite (not NA/NaN/Inf).", name))
    if (x != floor(x)) stop(sprintf("%s must be an integer.", name))
    as.integer(x)
  }
  n.iter <- .check_scalar_integer(n.iter, "n.iter")
  n.burnin <- .check_scalar_integer(n.burnin, "n.burnin")
  n.thin <- .check_scalar_integer(n.thin, "n.thin")

  if (n.iter < 1L) stop("n.iter must be >= 1.")
  if (n.burnin < 0L) stop("n.burnin must be >= 0.")
  if (n.thin < 1L) stop("n.thin must be >= 1.")

  if (n.iter <= n.burnin) stop("n.iter must be > n.burnin.")
  n.samples <- (n.iter - n.burnin) %/% n.thin
  if (n.samples < 1L) stop("MCMC settings yield zero post-burnin samples; increase n.iter or decrease n.burnin/n.thin.")

  invisible(TRUE)
}



# Internal helper to validate initial values for MCMC chains.
# Checks whether init has a valid single-chain or multi-chain structure, and validates the supplied components z, z0, a, b, u, and A.
# Returns TRUE invisibly if valid; otherwise throws an error.
.validate_init_chains <- function(init, w, n.chains) {
  K <- length(unique(w[!is.na(w)])) # w 검증됐으므로 NA를 제외한 unique한 값의 개

  if (is.null(init)) return(invisible(TRUE))
  if (length(n.chains) != 1L || !is.numeric(n.chains) || is.na(n.chains) || is.nan(n.chains) || is.infinite(n.chains) || n.chains != floor(n.chains) || n.chains < 1L) {
    stop("n.chains must be a single integer >= 1.")
  }

  chains <- as.integer(n.chains)

  m <- nrow(w)
  n <- ncol(w)

  allowed <- c("z", "z0", "a", "b", "u", "A")

  .check_finite_numeric <- function(x, name, where = "init") { # where: 다중체인에서 어느 체인에서 잘못됐는지 정확히 에러 메시지 전달
    if (!is.numeric(x)) stop(sprintf("%s$%s must be numeric.", where, name))
    if (any(is.nan(x), na.rm = TRUE)) stop(sprintf("%s$%s must not contain NaN.", where, name))
    if (anyNA(x)) stop(sprintf("%s$%s must not contain NA.", where, name))
    if (any(is.infinite(x))) stop(sprintf("%s$%s must not contain Inf/-Inf.", where, name))

    invisible(TRUE)
  }
  # Validate initial values for a single chain
  .validate_one_init <- function(init, where = "init") {
    if (is.null(init)) return(invisible(TRUE))
    if (!is.list(init)) stop(sprintf("%s must be a named list.", where))
    if (length(init) == 0L) return(invisible(TRUE))
    if (is.null(names(init)) || any(names(init) == "")) stop(sprintf("%s must be a named list with non-empty names.", where))
    bad <- setdiff(names(init), allowed)
    if (length(bad) > 0L) stop(sprintf("%s has unknown element(s): %s", where, paste(bad, collapse = ", ")))

    # Check z
    if (!is.null(init$z)) {
      if (!is.matrix(init$z)) stop(sprintf("%s$z must be a matrix.", where))
      if (!identical(dim(init$z), c(m, n))) stop(sprintf("%s$z must have dim = c(nrow(w), ncol(w)).", where))
      if (!is.numeric(init$z)) stop(sprintf("%s$z must be numeric.", where))

      obs <- !is.na(w)
      miss <- is.na(w)

      if (any(is.nan(init$z[obs]), na.rm = TRUE)) stop(sprintf("%s$z must not contain NaN on observed positions.", where))
      if (any(is.infinite(init$z[obs]))) stop(sprintf("%s$z must not contain Inf/-Inf on observed positions.", where))
      if (anyNA(init$z[obs])) stop(sprintf("%s$z must not contain NA on observed positions.", where))

      if (any(is.nan(init$z[miss]), na.rm = TRUE)) stop(sprintf("%s$z must not contain NaN on missing positions.", where))
      if (any(is.infinite(init$z[miss]), na.rm = TRUE)) stop(sprintf("%s$z must not contain Inf/-Inf on missing positions.", where))
    }


    # Check z0
    if (!is.null(init$z0)) {
      if (!is.atomic(init$z0) || !is.null(dim(init$z0))) stop(sprintf("%s$z0 must be a vector (no dimensions).", where))
      if (length(init$z0) != m) stop(sprintf("%s$z0 must have length nrow(w).", where))
      .check_finite_numeric(init$z0, "z0", where)
    }

    # Check a
    if (!is.null(init$a)) {
      if (!is.atomic(init$a) || !is.null(dim(init$a))) stop(sprintf("%s$a must be a vector (no dimensions).", where))
      if (length(init$a) != n) stop(sprintf("%s$a must have length ncol(w).", where))
      .check_finite_numeric(init$a, "a", where)
    }

    # Check b
    if (!is.null(init$b)) {
      if (!is.atomic(init$b) || !is.null(dim(init$b))) stop(sprintf("%s$b must be a vector (no dimensions).", where))
      if (length(init$b) != n) stop(sprintf("%s$b must have length ncol(w).", where))
      .check_finite_numeric(init$b, "b", where)
    }

    # Check u
    if (!is.null(init$u)) {
      if (!is.atomic(init$u) || !is.null(dim(init$u))) stop(sprintf("%s$u must be a vector (no dimensions).", where))
      if (length(init$u) != m) stop(sprintf("%s$u must have length nrow(w).", where))
      .check_finite_numeric(init$u, "u", where)
    }

    # Check A
    if (!is.null(init$A)) {
      if (!is.atomic(init$A) || !is.null(dim(init$A))) stop(sprintf("%s$A must be a vector (no dimensions).", where))
      if (!is.numeric(init$A)) stop(sprintf("%s$A must be numeric.", where))
      if (length(init$A) != K - 1L) stop(sprintf("%s$A must have length K-1.", where))
      .check_finite_numeric(init$A, "A", where)
      if (length(init$A) >= 2L && any(diff(init$A) <= 0)) stop(sprintf("%s$A must be strictly increasing.", where))
    }

    invisible(TRUE)
  }

  # Check whether init is a single-chain init or a multiple-chain init
  .is_named_init_list <- function(x) {
    is.list(x) && (length(x) == 0L || (!is.null(names(x)) && all(names(x) != "")))
  }

  # Single-chain case
  if (chains == 1L) {
    .validate_one_init(init, where = "init")
    return(invisible(TRUE))
  }

  # Multi-chain case
  if (.is_named_init_list(init)) {
    .validate_one_init(init, where = "init")
    return(invisible(TRUE))
  }

  # Validate the outer list for chain-specific init lists
  if (!is.list(init)) stop(paste("For multiple chains, init must be either", "(1) a single named list, or", "(2) an unnamed outer list of chain-specific init lists."))
  if (!is.null(names(init)) && any(names(init) != "")) stop("For multiple chains, a chain-specific init must be an unnamed outer list.")

  if (length(init) > chains) stop(sprintf("For multiple chains, a chain-specific init must have length less than or equal to chains (got %d, expected at most %d).", length(init), chains))

  for (i in seq_along(init)) {
    .validate_one_init(init[[i]], where = sprintf("init[[%d]]", i))
  }

  invisible(TRUE)
}

# Internal helper to validate the prior hyperparameter list.
# Checks whether prior is NULL or a named list with allowed component names, and validates each supplied hyperparameter as a single finite numeric value.
# Parameters that represent precisions or scale-like quantities must be positive.
# Returns TRUE invisibly if all checks pass; otherwise throws an error.
.validate_prior <- function(prior) {
  if (is.null(prior)) return(invisible(TRUE))
  if (!is.list(prior)) stop("prior must be NULL or a list.")
  if (length(prior) == 0L) return(invisible(TRUE))
  if (is.null(names(prior)) || any(names(prior) == "")) stop("prior must be a named list with non-empty names.")
  allowed <- c("mu_a", "tau_a", "mu_b", "tau_b", "delta", "gamma", "tau", "mu_u", "tau_u")
  prior_name <- names(prior)
  if (is.null(prior_name)) prior_name <- character(0)
  bad <- setdiff(prior_name, allowed)
  if (length(bad) > 0L) stop(sprintf("prior has unknown element(s): %s", paste(bad, collapse = ", ")))

  .check_scalar_numeric <- function(x, name, positive = FALSE) {
    if (length(x) != 1L) stop(sprintf("prior$%s must be a single value.", name))
    if (!is.numeric(x)) stop(sprintf("prior$%s must be numeric.", name))
    if (is.na(x) || is.nan(x) || is.infinite(x)) stop(sprintf("prior$%s must be finite (not NA/NaN/Inf).", name))
    if (positive && x <= 0) stop(sprintf("prior$%s must be > 0.", name))
    invisible(TRUE)
  }
  # 각 piror 인자 검사
  if (!is.null(prior[["mu_a"]]))
    .check_scalar_numeric(prior[["mu_a"]], "mu_a", positive = FALSE)
  if (!is.null(prior[["tau_a"]]))
    .check_scalar_numeric(prior[["tau_a"]], "tau_a", positive = TRUE)
  if (!is.null(prior[["mu_b"]]))
    .check_scalar_numeric(prior[["mu_b"]], "mu_b", positive = FALSE)
  if (!is.null(prior[["tau_b"]]))
    .check_scalar_numeric(prior[["tau_b"]], "tau_b", positive = TRUE)
  if (!is.null(prior[["delta"]]))
    .check_scalar_numeric(prior[["delta"]], "delta", positive = TRUE)
  if (!is.null(prior[["gamma"]]))
    .check_scalar_numeric(prior[["gamma"]], "gamma", positive = TRUE)
  if (!is.null(prior[["tau"]]))
    .check_scalar_numeric(prior[["tau"]], "tau", positive = TRUE)
  if (!is.null(prior[["mu_u"]]))
    .check_scalar_numeric(prior[["mu_u"]], "mu_u", positive = FALSE)
  if (!is.null(prior[["tau_u"]]))
    .check_scalar_numeric(prior[["tau_u"]], "tau_u", positive = TRUE)

  invisible(TRUE)
}
