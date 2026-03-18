# ---- validate_mcmc ----

testthat::test_that(".validate_mcmc accepts valid settings", {
  testthat::expect_true(.validate_mcmc(n.iter = 4000, n.burnin = 2000, n.thin = 5))
})

testthat::test_that(".validate_mcmc rejects non-scalar inputs", {
  testthat::expect_error(
    .validate_mcmc(n.iter = c(10, 20), n.burnin = 0, n.thin = 1),
    "n.iter must be a single value"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = c(0, 1), n.thin = 1),
    "n.burnin must be a single value"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = 0, n.thin = c(1, 2)),
    "n.thin must be a single value"
  )
})

testthat::test_that(".validate_mcmc rejects non-numeric", {
  testthat::expect_error(
    .validate_mcmc(n.iter = "10", n.burnin = 0, n.thin = 1),
    "n.iter must be numeric"
  )
})

testthat::test_that(".validate_mcmc rejects NA/NaN/Inf", {
  testthat::expect_error(
    .validate_mcmc(n.iter = NA, n.burnin = 0, n.thin = 1),
    "must be numeric"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = NaN, n.burnin = 0, n.thin = 1),
    "must be finite"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = Inf, n.burnin = 0, n.thin = 1),
    "must be finite"
  )
})

testthat::test_that(".validate_mcmc rejects non-integer values", {
  testthat::expect_error(
    .validate_mcmc(n.iter = 10.5, n.burnin = 0, n.thin = 1),
    "n.iter must be an integer"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = 0.2, n.thin = 1),
    "n.burnin must be an integer"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = 0, n.thin = 1.1),
    "n.thin must be an integer"
  )
})

testthat::test_that(".validate_mcmc checks ranges", {
  testthat::expect_error(
    .validate_mcmc(n.iter = 0, n.burnin = 0, n.thin = 1),
    "n.iter must be >= 1"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = -1, n.thin = 1),
    "n.burnin must be >= 0"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = 0, n.thin = 0),
    "n.thin must be >= 1"
  )
})

testthat::test_that(".validate_mcmc requires n.iter > n.burnin", {
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = 10, n.thin = 1),
    "n.iter must be > n.burnin"
  )
  testthat::expect_error(
    .validate_mcmc(n.iter = 10, n.burnin = 11, n.thin = 1),
    "n.iter must be > n.burnin"
  )
})

testthat::test_that(".validate_mcmc requires at least one post-burnin sample", {
  # (n.iter - n.burnin) %/% n.thin = 0
  testthat::expect_error(
    .validate_mcmc(n.iter = 11, n.burnin = 10, n.thin = 2),
    "zero post-burnin samples"
  )
})


# ---- validate_w ----

testthat::test_that(".validate_w accepts valid w", {
  w <- matrix(c(1,2,3, 2,3,1, NA,2,3), nrow = 3)
  testthat::expect_true(.validate_w(w))
})

testthat::test_that(".validate_w rejects non-matrix", {
  w <- data.frame(a = 1:3, b = 1:3)
  testthat::expect_error(.validate_w(w), "w must be a matrix")
})

testthat::test_that(".validate_w rejects non-numeric", {
  w <- matrix(c("1","2","3"), nrow=1)
  testthat::expect_error(.validate_w(w), "w must be a numeric")
})

testthat::test_that(".validate_w rejects all-NA", {
  w <- matrix(NA_real_, nrow = 2, ncol = 2)
  testthat::expect_error(.validate_w(w), "w has only NA values")
})



testthat::test_that(".validate_w rejects NaN/Inf", {
  w1 <- matrix(c(1, 2, 3, NaN), nrow=1)
  testthat::expect_error(.validate_w(w1), "must not contain NaN, Inf, or -Inf")

  w2 <- matrix(c(1, 2, 3, Inf), nrow=1)
  testthat::expect_error(.validate_w(w2), "must not contain NaN, Inf, or -Inf")
})

testthat::test_that(".validate_w rejects non-integers", {
  w <- matrix(c(1, 2.5, 3), nrow=1)
  testthat::expect_error(.validate_w(w), "must be integers")
})

testthat::test_that(".validate_w rejects values < 1", {
  w <- matrix(c(0, 1, 2), nrow=1)
  testthat::expect_error(.validate_w(w), "must be >= 1")
})

testthat::test_that(".validate_w rejects non-consecutive levels", {
  # {1, 3} gap 있음
  w <- matrix(c(1, 3, 3, NA), nrow=2)
  testthat::expect_error(.validate_w(w), "consecutive integer levels")
})

testthat::test_that(".validate_w requires at least 3 levels", {
  w <- matrix(c(1, 2, 2, 1, NA), nrow=1)
  testthat::expect_error(.validate_w(w), "at least 3 observed ordinal levels")
})


# ---- validate_D -----
testthat::test_that(".validate_D accepts valid D", {
  w <- matrix(1:12, nrow = 4)         # nrow(w) = 4
  D <- c(0L, 1L, 0L, 1L)
  testthat::expect_true(.validate_D(D, w))
})

testthat::test_that(".validate_D rejects non-vector (has dimensions)", {
  w <- matrix(1:12, nrow = 4)
  D <- matrix(c(0, 1, 0, 1), nrow = 2)  # dim 있음
  testthat::expect_error(.validate_D(D, w), "D must be a vector")
})

testthat::test_that(".validate_D rejects wrong length", {
  w <- matrix(1:12, nrow = 4)
  D <- c(0L, 1L, 0L)  # length 3 != 4
  testthat::expect_error(.validate_D(D, w), "Length of D must equal nrow\\(w\\)")
})

testthat::test_that(".validate_D rejects NA/NaN/Inf", {
  w <- matrix(1:12, nrow = 4)

  D1 <- c(0L, 1L, NA_integer_, 1L)
  testthat::expect_error(.validate_D(D1, w), "must not contain NA/NaN/Inf")

  D2 <- c(0, 1, NaN, 1)
  testthat::expect_error(.validate_D(D2, w), "must not contain NA/NaN/Inf")

  D3 <- c(0, 1, Inf, 1)
  testthat::expect_error(.validate_D(D3, w), "must not contain NA/NaN/Inf")
})

testthat::test_that(".validate_D rejects non-integer numeric", {
  w <- matrix(1:12, nrow = 4)
  D <- c(0, 1, 0.5, 1)
  testthat::expect_error(.validate_D(D, w), "integer-valued")
})

testthat::test_that(".validate_D rejects values not in {0,1}", {
  w <- matrix(1:12, nrow = 4)
  D <- c(0L, 1L, 2L, 1L)
  testthat::expect_error(.validate_D(D, w), "contain only 0 and 1")
})

testthat::test_that(".validate_D requires both classes 0 and 1", {
  w <- matrix(1:12, nrow = 4)

  D0 <- c(0L, 0L, 0L, 0L)
  testthat::expect_error(.validate_D(D0, w), "contain both classes 0 and 1")

  D1 <- c(1L, 1L, 1L, 1L)
  testthat::expect_error(.validate_D(D1, w), "contain both classes 0 and 1")
})

testthat::test_that(".validate_D allows logical and coerces to 0/1", {
  w <- matrix(1:12, nrow = 4)
  D <- c(FALSE, TRUE, FALSE, TRUE)
  testthat::expect_true(.validate_D(D, w))
})


# -------validate_init---------
testthat::test_that(".validate_init accepts NULL", {
  w <- matrix(1:12, nrow = 3)
  K <- 4
  testthat::expect_true(.validate_init(NULL, w, K))
})
# 정상적인 init인 경우
testthat::test_that(".validate_init accepts valid init", {
  w <- matrix(1:12, nrow = 3)
  m <- nrow(w)
  n <- ncol(w)
  K <- 4

  init <- list(z  = matrix(0, m, n),
               z0 = rep(0, m),
               a  = rep(0, n),
               b  = rep(1, n),
               u  = rep(0, m),
               A  = c(-1, 0, 1, Inf))

  testthat::expect_true(.validate_init(init, w, K))
})

# list가 아닌 경우
testthat::test_that(".validate_init rejects non-list", {
  w <- matrix(1:12, nrow = 3)
  testthat::expect_error(.validate_init(123, w, 4),
                         "init must be NULL or a list")
})

# init에 허용되지 않는 이름일 경우
testthat::test_that(".validate_init rejects unknown names", {
  w <- matrix(1:12, nrow = 3)
  init <- list(foo = 1)
  testthat::expect_error(.validate_init(init, w, 4),
                         "unknown element")
})

# z 테스트
testthat::test_that(".validate_init checks z dimension", {
  w <- matrix(1:12, nrow = 3)
  init <- list(z = matrix(0, 2, 2))  # wrong dim
  testthat::expect_error(.validate_init(init, w, 4),
                         "init\\$z must have dim")
})

# z0 길이가 잘못된 경우
testthat::test_that(".validate_init checks z0 length", {
  w <- matrix(1:12, nrow = 3)
  m <- nrow(w)
  n <- ncol(w)

  init <- list(
    z  = matrix(0, m, n),
    z0 = c(0, 0)
  )

  testthat::expect_error(
    .validate_init(init, w, 4),
    "init\\$z0 must have length"
  )
})


# a 길이가 잘못된 경우
testthat::test_that(".validate_init checks a length", {
  w <- matrix(1:12, nrow = 3)
  init <- list(a = c(0,0,0))  # ncol(w)=4
  testthat::expect_error(.validate_init(init, w, 4),
                         "init\\$a must have length")
})

# A 길이가 잘못된 경우
testthat::test_that(".validate_init checks A length", {
  w <- matrix(1:12, nrow = 3)
  init <- list(A = c(-1,0,Inf))  # length 3 but K=4
  testthat::expect_error(.validate_init(init, w, 4),
                         "init\\$A must have length")
})

# A 마지막 원소가 Inf가 아닌 경우
testthat::test_that(".validate_init requires A last element Inf", {
  w <- matrix(1:12, nrow = 3)
  init <- list(A = c(-1,0,1,2))  # last not Inf
  testthat::expect_error(.validate_init(init, w, 4),
                         "last element must be Inf")
})

# A 중간에 Inf가 있는 경우
testthat::test_that(".validate_init forbids Inf except last", {
  w <- matrix(1:12, nrow = 3)
  init <- list(A = c(-1,Inf,1,Inf))
  testthat::expect_error(.validate_init(init, w, 4),
                         "Inf only as the last element")
})

# A strictly increasing이 아닌 경우
testthat::test_that(".validate_init checks A strictly increasing", {
  w <- matrix(1:12, nrow = 3)
  init <- list(A = c(-1, 2, 1, Inf))  # not increasing
  testthat::expect_error(.validate_init(init, w, 4),
                         "strictly increasing")
})


#------- .validate_prior
testthat::test_that(".validate_prior accepts NULL", {
  testthat::expect_true(.validate_prior(NULL))
})

testthat::test_that(".validate_prior rejects non-list", {
  testthat::expect_error(.validate_prior(123), "prior must be NULL or a list")
  testthat::expect_error(.validate_prior("x"), "prior must be NULL or a list")
})

testthat::test_that(".validate_prior rejects unknown names", {
  pr <- list(mu_a = 0, badname = 1)
  testthat::expect_error(.validate_prior(pr), "prior has unknown element")
})

testthat::test_that(".validate_prior accepts valid prior (full)", {
  pr <- list(
    mu_a = 0,
    tau_a = 1,
    mu_b = 0,
    tau_b = 1,
    delta = 5,
    gamma = 0.1,
    tau = 0.1,
    mu_u = 0,
    tau_u = 1
  )
  testthat::expect_true(.validate_prior(pr))
})

testthat::test_that(".validate_prior accepts partial prior (subset allowed)", {
  pr <- list(tau_a = 1, gamma = 0.1)
  testthat::expect_true(.validate_prior(pr))
})

testthat::test_that(".validate_prior enforces scalar length-1", {
  pr <- list(tau_a = c(1, 2))
  testthat::expect_error(.validate_prior(pr), "prior\\$tau_a must be a single value")
})

testthat::test_that(".validate_prior enforces numeric", {
  pr <- list(mu_b = "0")
  testthat::expect_error(.validate_prior(pr), "prior\\$mu_b must be numeric")
})

testthat::test_that(".validate_prior rejects NA/NaN/Inf", {
  pr1 <- list(tau_a = NA_real_)
  testthat::expect_error(.validate_prior(pr1), "must be finite")

  pr2 <- list(tau_a = NaN)
  testthat::expect_error(.validate_prior(pr2), "must be finite")

  pr3 <- list(tau_a = Inf)
  testthat::expect_error(.validate_prior(pr3), "must be finite")
})

testthat::test_that(".validate_prior enforces positivity where required", {
  # tau_* must be > 0
  testthat::expect_error(.validate_prior(list(tau_a = 0)), "prior\\$tau_a must be > 0")
  testthat::expect_error(.validate_prior(list(tau_b = -1)), "prior\\$tau_b must be > 0")
  testthat::expect_error(.validate_prior(list(tau = 0)), "prior\\$tau must be > 0")
  testthat::expect_error(.validate_prior(list(tau_u = -0.1)), "prior\\$tau_u must be > 0")

  # gamma, delta must be > 0
  testthat::expect_error(.validate_prior(list(gamma = 0)), "prior\\$gamma must be > 0")
  testthat::expect_error(.validate_prior(list(delta = -1)), "prior\\$delta must be > 0")
})

testthat::test_that(".validate_prior allows mu_* to be negative", {
  pr <- list(mu_a = -10, mu_b = -0.5, mu_u = -3)
  testthat::expect_true(.validate_prior(pr))
})
