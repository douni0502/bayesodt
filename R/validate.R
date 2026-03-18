# Validate rating matrix:
# must be numeric matrix
# NA allowed
# observed entries must be finite integers >= 1
# levels must be consecutive 1, 2, ..., K with K >= 3
.validate_w <- function(w) {
  if (!is.matrix(w)) stop("w must be a matrix.")
  if (!is.numeric(w)) stop("w must be a numeric (integer-valued, with possible NA).")
  if (any(is.nan(w))) stop("w must not contain NaN.")
  x <- w[!is.na(w)]
  if (length(x) == 0L) stop("w has only NA values; cannot determine rating levels.")
  # 4) NaN/Inf/-Inf 검사 (NA는 허용)
  if (any(!is.finite(x))) stop("w must not contain Inf, or -Inf (NA is allowed).")
  # 5) 정수인지 검사
  if (any(x != floor(x))) stop("All non-missing entries of w must be integers.")
  # 6) 1 이상인지 검사
  if (any(x < 1)) stop("All non-missing entries of w must be >= 1.")
  # 7) 등급이 연속인지 검사 -> {1, 2, ..., K} 꼴이어야 함
  lev <- sort(unique(as.integer(x)))
  if (!identical(lev, seq.int(from = 1L, to = max(lev)))) stop("Non-missing entries of w must contain consecutive integer levels starting at 1 (i.e., 1, 2, ..., K with no gaps).")
  # 8) K(등급 수) 최소 조건
  K <- length(lev)
  if (K < 3) stop(sprintf("w must have at least 3 observed ordinal levels (found K = %d).", K))

  invisible(TRUE)
}

# Validate binary disease status vector D:
# must be a dimensionless atomic vector
# length must match nrow(w)
# must not contain NA, NaN, Inf, or -Inf
# only 0 and 1 are allwed
# both classes 0 and 1 msut be a present
.validate_D <- function(D, w) {
  if (is.null(D)) return(invisible(TRUE))
  # 1) D는 벡터인지 확인
  if (!is.atomic(D) || !is.null(dim(D))) stop("D must be a vector (no dimensions).")
  # 2) 길이 체크 w 행의 길이와 같은지
  if (length(D) != nrow(w)) stop(sprintf("Length of D must equal nrow(w) (got length(D)=%d, nrow(w)=%d).", length(D), nrow(w)))
  # 3) NA, NaN, Inf 금지
  if (anyNA(D) || any(is.nan(D)) || any(is.infinite(D))) stop("D must not contain NA/NaN/Inf.")
  # 4) 0 또는 1만 허용
  if (is.logical(D)) D <- as.integer(D)
  if (!is.numeric(D) || any(D != floor(D))) stop("D must be integer-valued (0/1).")
  if (!all(D %in% c(0L, 1L))) stop("D must contain only 0 and 1.")
  # 5) 0과 1 둘 다 최소 하나 이상 존재해야 됨.
  if (!all(c(0L, 1L) %in% D)) stop("D must contain both classes 0 and 1.")

  invisible(TRUE)
}


# Validate MCMC tuning parameters:
# n.iter, n.burnin, and n.thin must each be a single finite integer
# n.iter must be >= 1
# n.burnin must be >= 0
# n.thin must be >= 1
# n.iter must be greater than n.burnin
# at least one post-burin sample must be retained
.validate_mcmc <- function(n.iter, n.burnin, n.thin) {
  # n.iter, n.burnin, n.thin 공통적으로 검사하는 함수
  .check_scalar_integer <- function(x, name) {
    # 길이가 1인 스칼라인지
    if (length(x) != 1L) stop(sprintf("%s must be a single value.", name))
    # 숫자인지 체크
    if (!is.numeric(x)) stop(sprintf("%s must be numeric.", name))
    # NA/NaN/Inf 아닌지
    if (is.na(x) || is.nan(x) || is.infinite(x)) stop(sprintf("%s must be finite (not NA/NaN/Inf).", name))
    # 정수인지
    if (x != floor(x)) stop(sprintf("%s must be an integer.", name))
    # 정수로 변환해서 반환
    as.integer(x)
  }
  n.iter <- .check_scalar_integer(n.iter, "n.iter")
  n.burnin <- .check_scalar_integer(n.burnin, "n.burnin")
  n.thin <- .check_scalar_integer(n.thin, "n.thin")

  # 범위 맞는지 확인
  if (n.iter < 1L) stop("n.iter must be >= 1.")
  if (n.burnin < 0L) stop("n.burnin must be >= 0.")
  if (n.thin < 1L) stop("n.thin must be >= 1.")

  # burnin보다 iter가 당연히 더 커야함
  if (n.iter <= n.burnin) stop("n.iter must be > n.burnin.")
  # 사후 sample이 최소 하나가 되어야 함.
  n.keep <- (n.iter - n.burnin) %/% n.thin
  if (n.keep < 1L) stop("MCMC settings yield zero post-burnin samples; increase n.iter or decrease n.burnin/n.thin.")

  invisible(TRUE)
}

# Validate init:
# NULL is allowed
# for a single chain, init must be a named list
# for multiple chains, init can be either:
# (1) a single named list, or
# (2) an unnamed list of chain-specific named lists of length = n.chains (i.e., the outer list must be unnamed)
# allowed element names are: z, z0, a, b, u, A
# each provided element must satisfy dimension and finiteness constraints
.validate_init_chains <- function(init, w, n.chains) {
  K <- length(unique(w[!is.na(w)])) # w 검증됐으므로 NA를 제외한 unique한 값의 개
  # NULL 은 허용
  if (is.null(init)) return(invisible(TRUE))
  # 체인이 1이상의 정수인지 검사
  if (length(n.chains) != 1L || !is.numeric(n.chains) || is.na(n.chains) || is.nan(n.chains) || is.infinite(n.chains) || n.chains != floor(n.chains) || n.chains < 1L) {
    stop("n.chains must be a single integer >= 1.")
  }
  chains <- as.integer(n.chains)

  m <- nrow(w)
  n <- ncol(w)

  # 허용된 parameter
  allowed <- c("z", "z0", "a", "b", "u", "A")

  # 내부함수: numeric인지 유한한지
  .check_finite_numeric <- function(x, name, where = "init") {
    if (!is.numeric(x)) stop(sprintf("%s$%s must be numeric.", where, name))
    if (any(is.nan(x), na.rm = TRUE)) stop(sprintf("%s$%s must not contain NaN.", where, name))
    if (anyNA(x)) stop(sprintf("%s$%s must not contain NA.", where, name))
    if (any(is.infinite(x))) stop(sprintf("%s$%s must not contain Inf/-Inf.", where, name))

    invisible(TRUE)
  }

  # 내부 함수: init 한 세트 검사
  .validate_one_init <- function(init_one, where = "init") {
    # 리스트인지 체크
    if (!is.list(init_one)) stop(sprintf("%s must be a named list.", where))
    # 빈 리리트 허용
    if (length(init_one) == 0L) return(invisible(TRUE))
    # 이름이 무조건 있어야 되고 이름이 있더라도 공백이면 에러
    if (is.null(names(init_one)) || any(names(init_one) == "")) stop(sprintf("%s must be a named list with non-empty names.", where))
    # 허용된 이름만 있는지 체크
    bad <- setdiff(names(init_one), allowed)
    if (length(bad) > 0L) stop(sprintf("%s has unknown element(s): %s", where, paste(bad, collapse = ", ")))

    # z 검사
    if (!is.null(init_one$z)) {
      if (!is.matrix(init_one$z)) stop(sprintf("%s$z must be a matrix.", where))
      if (!identical(dim(init_one$z), c(m, n))) stop(sprintf("%s$z must have dim = c(nrow(w), ncol(w)).", where))
      if (!is.numeric(init_one$z)) stop(sprintf("%s$z must be numeric.", where))

      obs <- !is.na(w)
      miss <- is.na(w)

      if (any(is.nan(init_one$z[obs]), na.rm = TRUE)) stop(sprintf("%s$z must not contain NaN on observed positions.", where))
      if (any(is.infinite(init_one$z[obs]))) stop(sprintf("%s$z must not contain Inf/-Inf on observed positions.", where))
      if (anyNA(init_one$z[obs])) stop(sprintf("%s$z must not contain NA on observed positions.", where))

      # missing 위치는 NA 허용
      if (any(is.nan(init_one$z[miss]), na.rm = TRUE)) stop(sprintf("%s$z must not contain NaN on missing positions.", where))
      if (any(is.infinite(init_one$z[miss]), na.rm = TRUE)) stop(sprintf("%s$z must not contain Inf/-Inf on missing positions.", where))
    }


    # z0 검사
    if (!is.null(init_one$z0)) {
      if (!is.atomic(init_one$z0) || !is.null(dim(init_one$z0))) stop(sprintf("%s$z0 must be a vector (no dimensions).", where))
      if (length(init_one$z0) != m) stop(sprintf("%s$z0 must have length nrow(w).", where))
      .check_finite_numeric(init_one$z0, "z0", where)
    }

    # a 검사
    if (!is.null(init_one$a)) {
      if (!is.atomic(init_one$a) || !is.null(dim(init_one$a))) stop(sprintf("%s$a must be a vector (no dimensions).", where))
      if (length(init_one$a) != n) stop(sprintf("%s$a must have length ncol(w).", where))
      .check_finite_numeric(init_one$a, "a", where)
    }

    # b 검사
    if (!is.null(init_one$b)) {
      if (!is.atomic(init_one$b) || !is.null(dim(init_one$b))) stop(sprintf("%s$b must be a vector (no dimensions).", where))
      if (length(init_one$b) != n) stop(sprintf("%s$b must have length ncol(w).", where))
      .check_finite_numeric(init_one$b, "b", where)
    }

    # u 검사
    if (!is.null(init_one$u)) {
      if (!is.atomic(init_one$u) || !is.null(dim(init_one$u))) stop(sprintf("%s$u must be a vector (no dimensions).", where))
      if (length(init_one$u) != m) stop(sprintf("%s$u must have length nrow(w).", where))
      .check_finite_numeric(init_one$u, "u", where)
    }

    # A 검사
    if (!is.null(init_one$A)) {
      if (!is.atomic(init_one$A) || !is.null(dim(init_one$A))) stop(sprintf("%s$A must be a vector (no dimensions).", where))
      if (!is.numeric(init_one$A)) stop(sprintf("%s$A must be numeric.", where))
      if (length(init_one$A) != K - 1L) stop(sprintf("%s$A must have length K-1.", where))
      .check_finite_numeric(init_one$A, "A", where)
      if (length(init_one$A) >= 2L && any(diff(init_one$A) <= 0)) stop(sprintf("%s$A must be strictly increasing.", where))
    }

    invisible(TRUE)
  }

  # 내부 함수: init 하나의 세트인지 여러 개의 세트가 담긴 리스트인지 검사
  .is_named_init_list <- function(x) {
    is.list(x) && (length(x) == 0L || (!is.null(names(x)) && all(names(x) != "")))
  }

  # 단일 체인인 경우 검사
  if (chains == 1L) {
    .validate_one_init(init, where = "init")
    return(invisible(TRUE))
  }

  # 다중 체인인 경우
  # case 1) 하나의 체인 리스트만 주어진 경우
  if (.is_named_init_list(init)) {
    .validate_one_init(init, where = "init")
    return(invisible(TRUE))
  }

  # case 2) 이름이 없는 바깥 리스트 안에 체인별 init 리스트들이 주어진 경우
  if (!is.list(init)) stop(paste("For multiple chains, init must be either", "(1) a single named list, or", "(2) an unnamed outer list of chain-specific init lists."))
  if (!is.null(names(init)) && any(names(init) != "")) stop("For multiple chains, a chain-specific init must be an unnamed outer list.")

  # 체인의 수와 리스트 안의 체인 리스트 개수가 맞는지 검사
  if (length(init) != chains) stop(sprintf(paste("For multiple chains, a chain-specific init must have length equal to chains", "(got %d, expected %d)."), length(init), chains))

  for (i in seq_len(n.chains)) {
    .validate_one_init(init[[i]], where = sprintf("init[[%d]]", i))
  }

  invisible(TRUE)
}

# Validate prior hyperparameters:
# NULL is allowed
# prior must be a named list with non-empty names
# only known prior names are allowed
# each supplied prior value must be a single finite numeric value
# precision/rate-like parameters must be strictly positive
.validate_prior <- function(prior) {
  # NULL 허용 (기본값 생성은 다른 함수에서 처리)
  if (is.null(prior)) return(invisible(TRUE))
  # 1) list인지
  if (!is.list(prior)) stop("prior must be NULL or a list.")
  # 빈 list 허용
  if (length(prior) == 0L) return(invisible(TRUE))
  # prior에 이름이 있어야 함. 그리고 빈 이름 ""도 금지
  if (is.null(names(prior)) || any(names(prior) == "")) stop("prior must be a named list with non-empty names.")
  # 2) 허용된 이름만 있는지 체크
  allowed <- c("mu_a", "tau_a", "mu_b", "tau_b", "delta", "gamma", "tau", "mu_u", "tau_u")
  prior_name <- names(prior)
  if (is.null(prior_name)) prior_name <- character(0)
  bad <- setdiff(prior_name, allowed)
  if (length(bad) > 0L) stop(sprintf("prior has unknown element(s): %s", paste(bad, collapse = ", ")))

  # 내부 함수: 유한한 실수인지 체크
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
