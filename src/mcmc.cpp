#include <Rcpp.h>
// [[Rcpp::depends(RcppTN)]]
#include <RcppTN.h>
// [[Rcpp::depends(fntl)]]
#include <fntl.h>
using namespace Rcpp;

// Update latent ordinal variables Z_ij.
// Observed entries are sampled from truncated N(a_j + b_j * u_i, 1) according to the observed category; missing entries are sampled from the corresponding untruncated normal distribution.
NumericMatrix sample_z_ij(NumericMatrix z, NumericVector a, NumericVector b, NumericVector u, NumericVector k_lower, NumericVector k_upper,
                          List W_list, SEXP W_miss) {
  RNGScope scope;

  const int K = W_list.size();
  double* zptr = REAL(z);

  for (int kk = 0; kk < K; ++kk) {
    List W = W_list[kk];

    IntegerVector ii = W["i"];
    IntegerVector jj = W["j"];
    IntegerVector idx = W["idx"];

    const int nk = idx.size();
    const double low  = k_lower[kk];
    const double high = k_upper[kk];

    for (int t = 0; t < nk; ++t) {
      const int i = ii[t] - 1;
      const int j = jj[t] - 1;

      const double mu = a[j] + b[j] * u[i];

      // RcppTN truncated normal sampler: N(mu, 1) truncated to (low, high)
      const double draw = RcppTN::rtn1(mu, 1.0, low, high);

      const int pos = idx[t] - 1;
      zptr[pos] = draw;
    }
  }

  // Missing entries are updated from the corresponding untruncated N(mu, 1)
  if (W_miss != R_NilValue) {
    List Wm(W_miss);
    IntegerVector ii = Wm["i"];
    IntegerVector jj = Wm["j"];
    IntegerVector idx = Wm["idx"];

    const int nm = idx.size();
    for (int t = 0; t < nm; ++t) {
      const int i = ii[t] - 1;
      const int j = jj[t] - 1;

      const double mu = a[j] + b[j] * u[i];
      const int pos = idx[t] - 1;

      zptr[pos] = R::rnorm(mu, 1.0);
    }
  }

  return z;
}

// Update latent disease severities u_i when D is unknown.
// Samples u_i from its Gaussian full conditional and then centers the sampled vector to have mean 0 for identifiability.
NumericVector sample_u_i_unknown(const NumericMatrix& z, const NumericVector& a, const NumericVector& b, const double mu_u, const double tau_u) {
  RNGScope scope;

  const int m = z.nrow();
  const int n = z.ncol();

  // den = sum(b^2) + tau_u
  double sum_b2 = 0.0;
  for (int j = 0; j < n; ++j) sum_b2 += b[j] * b[j];
  const double den = sum_b2 + tau_u;

  // ab = sum(a*b)
  double ab = 0.0;
  for (int j = 0; j < n; ++j) ab += a[j] * b[j];

  // mean_u[i] = ((z %*% b)[i] - ab + tau_u * mu_u) / den
  // and sample u_i ~ N(mean_u[i], 1/den)
  NumericVector u(m);
  const double sd = std::sqrt(1.0 / den);
  const double add_const = (-ab + tau_u * mu_u) / den;

  for (int i = 0; i < m; ++i) {
    double zb_i = 0.0;
    // zb_i = sum_j z(i,j) * b[j]
    for (int j = 0; j < n; ++j) {
      zb_i += z(i, j) * b[j];
    }
    const double mean_i = zb_i / den + add_const;
    u[i] = R::rnorm(mean_i, sd);
  }

  // center: u <- u - mean(u)
  double ubar = 0.0;
  for (int i = 0; i < m; ++i) ubar += u[i];
  ubar /= (double)m;
  for (int i = 0; i < m; ++i) u[i] -= ubar;

  return u;
}


// Update latent disease severities u_i when D is known.
// Samples u_i from its Gaussian full conditional given z and z0, then applies a common shift so that mean(Phi(u_i)) matches mean(D).
NumericVector sample_u_i_known(const NumericMatrix& z, const NumericVector& a, const NumericVector& b, const NumericVector& z0, const NumericVector& D,
                               const double mu_u, const double tau_u,
                               const double tol = 1e-8, const unsigned int max_iter = 200) {
  RNGScope scope;

  const int m = z.nrow();
  const int n = z.ncol();

  // den = sum(b^2) + 1 + tau_u
  double sum_b2 = 0.0;
  for (int j = 0; j < n; ++j) sum_b2 += b[j] * b[j];
  const double den = sum_b2 + 1.0 + tau_u;

  // ab = sum(a*b)
  double ab = 0.0;
  for (int j = 0; j < n; ++j) ab += a[j] * b[j];

  // target = mean(D)
  double target = 0.0;
  for (int i = 0; i < m; ++i) target += D[i];
  target /= (double)m;

  // sample u_i
  NumericVector u(m);
  const double sd = std::sqrt(1.0 / den);
  const double add_const = (-ab + tau_u * mu_u) / den;

  for (int i = 0; i < m; ++i) {
    double zb_i = 0.0;
    for (int j = 0; j < n; ++j) {
      zb_i += z(i, j) * b[j];
    }
    const double mean_i = (zb_i + z0[i]) / den + add_const;
    u[i] = R::rnorm(mean_i, sd);
  }

  // f(c) = mean(Phi(u + c)) - target  (monotone increasing)
  fntl::dfd f = [&](double c) {
    double s = 0.0;
    for (int i = 0; i < m; ++i) {
      s += R::pnorm5(u[i] + c, 0.0, 1.0, 1, 0);
    }
    return s / (double)m - target;
  };

  // bracket 찾기 (R 코드와 동일한 전략)
  double low = -20.0, high = 20.0;
  double flow = f(low), fhigh = f(high);

  if (flow * fhigh > 0.0) {
    low = -50.0; high = 50.0;
    flow = f(low); fhigh = f(high);
  }

  // 그래도 안 잡히면 조금 더 확장 (안전장치)
  if (flow * fhigh > 0.0) {
    double L = 50.0;
    for (int k = 0; k < 20; ++k) {
      L *= 2.0;
      low = -L; high = L;
      flow = f(low); fhigh = f(high);
      if (flow * fhigh <= 0.0) break;
      if (L > 1e6) break;
    }
    if (flow * fhigh > 0.0) {
      stop("Failed to bracket root for shift. target may be unreachable numerically.");
    }
  }

  // fntl Brent root-finding
  fntl::findroot_args args;
  args.tol = tol;
  args.maxiter = max_iter;

  auto out = fntl::findroot_brent(f, low, high, args);

  // 상태 체크 (OK 아니면 에러)
  if (out.status != fntl::findroot_status::OK) {
    stop(out.message);
  }

  const double ss = out.root;

  // shift u <- u + ss
  for (int i = 0; i < m; ++i) u[i] += ss;

  return u;
}




// Update latent disease-status variables z0_i.
// For D_i = 1, sample z0_i from N(u_i, 1) truncated to (0, Inf); for D_i = 0, sample z0_i from N(u_i, 1) truncated to (-Inf, 0).
// idx1 and idx0 are 1-based indices from R.
NumericVector sample_z0(NumericVector u,
                        IntegerVector idx1, IntegerVector idx0) {
  RNGScope scope;
  int m = u.size();
  NumericVector z0(m);

  // idx1, idx0 는 R의 which() 결과라 1-based
  for (int t = 0; t < idx1.size(); ++t) {
    int i = idx1[t] - 1;
    double mu = u[i];
    z0[i] = RcppTN::rtn1(mu, 1.0, 0.0, R_PosInf);     // (0, Inf)
  }

  for (int t = 0; t < idx0.size(); ++t) {
    int i = idx0[t] - 1;
    double mu = u[i];
    z0[i] = RcppTN::rtn1(mu, 1.0, R_NegInf, 0.0);     // (-Inf, 0)
  }

  return z0;
}

// Update cutpoints A_k for the ordinal model.
// Each cutpoint is sampled uniformly from its full conditional interval, determined by (i) neighboring cutpoints, (ii) the spacing constraint delta, and (iii) the current latent z values in adjacent categories.
// W_list contains 1-based linear indices for entries in each observed category.
NumericVector sample_A(NumericVector A, int K, NumericMatrix z, double delta,
                       List W_list) {
  RNGScope scope;

  // z is column-major; we will use linear indexing
  double* zptr = REAL(z);

  // k별 z.max / z.min
  NumericVector z_max(K, R_NegInf);
  NumericVector z_min(K, R_PosInf);

  for (int kk = 0; kk < K; ++kk) {
    List W = W_list[kk];
    if (!W.containsElementNamed("idx")) stop("Each W_list[[k]] must contain 'idx'.");

    IntegerVector idx = W["idx"]; // 1-based linear indices
    int nk = idx.size();
    if (nk > 0) {
      double mx = R_NegInf;
      double mn = R_PosInf;
      for (int t = 0; t < nk; ++t) {
        int pos = idx[t] - 1;     // to 0-based
        double val = zptr[pos];
        if (val > mx) mx = val;
        if (val < mn) mn = val;
      }
      z_max[kk] = mx;
      z_min[kk] = mn;
    }
  }
  // --- Update cutpoints A[1:(K-1)] (0-based: A[0..K-2]) ---
  // Edge-case handling:
  // K==2: only one cutpoint. No neighbor constraints (A[1] is Inf).
  // K==3: two cutpoints. Has k=1 and k=K-1 updates only.
  // K>=4: full general loop.

  auto draw_unif_checked = [&](double L, double U) -> double {
    if (L > U) stop("Invalid bounds in sample_A_rcpp: L > U (L=%f, U=%f).", L, U);
    if (L == U) return L;
    return R::runif(L, U);
  };

  // k = 1  (0-based index 0)
  {
    double neighbor_low = (K >= 3) ? (A[1] - delta) : (-delta);
    double neighbor_up  = (K >= 3) ? (A[1])         : ( delta);

    double L = std::max({-delta, neighbor_low, (double)z_max[0]});
    double U = std::min({ delta, neighbor_up,  (double)z_min[1]}); // z_min[2] in R -> z_min[1] here
    A[0] = draw_unif_checked(L, U);
  }

  // middle k = 2 .. K-2  (0-based: 1 .. K-3) only if K>=4
  if (K >= 4) {
    for (int k = 1; k <= K - 3; ++k) { // k corresponds to R's k=2..K-2
      double L = std::max({ (double)A[k - 1],
                          (double)(A[k + 1] - delta),
                          (double)z_max[k] });

      double U = std::min({ (double)(A[k - 1] + delta),
                          (double)A[k + 1],
                                   (double)z_min[k + 1] });

      A[k] = draw_unif_checked(L, U);
    }
  }

  // last cutpoint k = K-1  (0-based index K-2), exists when K>=2 always
  if (K >= 2) {
    int last = K - 2;

    double left  = (K >= 3) ? A[last - 1] : (-delta);
    double leftU = (K >= 3) ? (A[last - 1] + delta) : (delta);

    double L = std::max({ (double)left, (double)z_max[last] });
    double U = std::min({ (double)leftU, (double)A[K - 1], (double)z_min[K - 1] }); // A[K-1] is Inf
    A[last] = draw_unif_checked(L, U);
  }

  return A;
}



List MCMC(NumericMatrix z, NumericVector z0, NumericVector u, Nullable<NumericVector> D, NumericVector a, NumericVector b, NumericVector A,
          double mu_a, double tau_a, double mu_b, double tau_b, double mu_u, double tau_u, double gamma, double delta, double tau,
          int K, List W_list, SEXP W_miss, Nullable<IntegerVector> idx1, Nullable<IntegerVector> idx0,
          double tol_shift = 1e-8, unsigned int max_iter_shift = 200) {
  RNGScope scope;

  const int m = z.nrow();
  const int n = z.ncol();

  // cutpoint bounds
  NumericVector k_lower(K), k_upper(K);
  k_lower[0] = R_NegInf;
  for (int kk = 1; kk < K; ++kk) k_lower[kk] = A[kk - 1];
  for (int kk = 0; kk < K; ++kk) k_upper[kk] = A[kk];

  const bool known = D.isNotNull();
  NumericVector Dv;
  IntegerVector idx1v;
  IntegerVector idx0v;

  if (known) {
    Dv = as<NumericVector>(D);
    if (Dv.size() != m) stop("D must have length m");

    idx1v = IntegerVector(idx1);
    idx0v = IntegerVector(idx0);
  }

  // ---- z0, z, u 업데이트 ----
  if (!known) {
    z = sample_z_ij(z, a, b, u, k_lower, k_upper, W_list, W_miss);
    u = sample_u_i_unknown(z, a, b, mu_u, tau_u);
  } else {
    z0 = sample_z0(u, idx1v, idx0v);
    z = sample_z_ij(z, a, b, u, k_lower, k_upper, W_list, W_miss);
    u = sample_u_i_known(z, a, b, z0, Dv, mu_u, tau_u, tol_shift, max_iter_shift);
  }

  // ---- a 업데이트 ----
  double sum_u = 0.0;
  for (int i = 0; i < m; ++i) sum_u += u[i];

  NumericVector colsum_z(n);
  for (int j = 0; j < n; ++j) {
    double s = 0.0;
    for (int i = 0; i < m; ++i) s += z(i, j);
    colsum_z[j] = s;
  }

  const double den_a = (double)m + tau_a;
  const double sd_a = std::sqrt(1.0 / den_a);
  for (int j = 0; j < n; ++j) {
    const double mean_aj = (colsum_z[j] - b[j] * sum_u) / den_a;
    a[j] = R::rnorm(mean_aj, sd_a);
  }

  // ---- b 업데이트 ----
  double sum_u2 = 0.0;
  for (int i = 0; i < m; ++i) sum_u2 += u[i] * u[i];

  const double den_b = sum_u2 + tau_b;
  const double sd_b = std::sqrt(1.0 / den_b);

  for (int j = 0; j < n; ++j) {
    double cross = 0.0;
    for (int i = 0; i < m; ++i) cross += u[i] * z(i, j);
    const double mean_bj = (cross - a[j] * sum_u + mu_b * tau_b) / den_b;
    b[j] = R::rnorm(mean_bj, sd_b);
  }

  // ---- tau_a 업데이트 (Gamma: shape/rate -> R::rgamma는 shape/scale) ----
  double ss_a = 0.0;
  for (int j = 0; j < n; ++j) {
    const double d = a[j] - mu_a;
    ss_a += d * d;
  }
  {
    const double shape = gamma + 0.5 * (double)n;
    const double rate  = gamma + 0.5 * ss_a;
    tau_a = R::rgamma(shape, 1.0 / rate);
  }

  // ---- mu_b 업데이트 ----
  double sum_b = 0.0;
  for (int j = 0; j < n; ++j) sum_b += b[j];

  {
    const double den_mub = (double)n * tau_b + tau;
    const double mean_mub = (tau_b * sum_b) / den_mub;
    const double sd_mub = std::sqrt(1.0 / den_mub);
    mu_b = R::rnorm(mean_mub, sd_mub);
  }

  // ---- tau_b 업데이트 ----
  double ss_b = 0.0;
  for (int j = 0; j < n; ++j) {
    const double d = b[j] - mu_b;
    ss_b += d * d;
  }
  {
    const double shape = gamma + 0.5 * (double)n;
    const double rate  = gamma + 0.5 * ss_b;
    tau_b = R::rgamma(shape, 1.0 / rate);
  }

  // ---- A 업데이트 ----
  A = sample_A(A, K, z, delta, W_list);

  return List::create(_["z0"] = z0, _["z"] = z, _["u"] = u, _["a"] = a, _["b"] = b, _["A"] = A,
                      _["tau_a"] = tau_a, _["mu_b"] = mu_b, _["tau_b"] = tau_b);
}


// [[Rcpp::export]]
List MCMC_loop(List state, Nullable<NumericVector> D,
               double mu_a, double tau_a, double mu_b, double tau_b, double mu_u, double tau_u, double gamma, double delta, double tau,
               int K, List W_list, SEXP W_miss, Nullable<IntegerVector> idx1, Nullable<IntegerVector> idx0,
               int n_iter, int n_burnin, int n_thin,
               double tol_shift = 1e-8, unsigned int max_iter_shift = 200) {

  RNGScope scope;

  NumericMatrix z = state["z"];
  NumericVector z0 = state["z0"];
  NumericVector u  = state["u"];
  NumericVector a  = state["a"];
  NumericVector b  = state["b"];
  NumericVector A  = state["A"];

  double tau_a_cur = state["tau_a"];
  double mu_b_cur  = state["mu_b"];
  double tau_b_cur = state["tau_b"];

  const int m = z.nrow();
  const int n = z.ncol();
  const int num = (n_iter - n_burnin) / n_thin;

  NumericMatrix A_save(num, K - 1);   // cutpoints
  NumericMatrix a_save(num, n);       // rater bias
  NumericMatrix b_save(num, n);       // rater magnifier/scale
  NumericMatrix u_save(num, m);       // latent patient effects
  NumericVector tau_a_save(num);
  NumericVector mu_b_save(num);
  NumericVector tau_b_save(num);

  int save_idx = 0;

  for (int iter = 1; iter <= n_iter; ++iter) {
    List up = MCMC(z, z0, u, D, a, b, A,
      mu_a, tau_a_cur, mu_b_cur, tau_b_cur, mu_u, tau_u, gamma, delta, tau, K, W_list, W_miss, idx1, idx0,
      tol_shift, max_iter_shift
    );

    // overwrite current state
    z  = as<NumericMatrix>(up["z"]);
    z0 = as<NumericVector>(up["z0"]);
    u  = as<NumericVector>(up["u"]);
    a  = as<NumericVector>(up["a"]);
    b  = as<NumericVector>(up["b"]);
    A  = as<NumericVector>(up["A"]);
    tau_a_cur = as<double>(up["tau_a"]);
    mu_b_cur  = as<double>(up["mu_b"]);
    tau_b_cur = as<double>(up["tau_b"]);

    if (iter > n_burnin && ((iter - n_burnin) % n_thin == 0)) {
      for (int k = 0; k < K - 1; ++k) A_save(save_idx, k) = A[k];
      for (int j = 0; j < n; ++j) {
        a_save(save_idx, j) = a[j];
        b_save(save_idx, j) = b[j];
      }
      for (int i = 0; i < m; ++i) u_save(save_idx, i) = u[i];

      tau_a_save[save_idx] = tau_a_cur;
      mu_b_save[save_idx]  = mu_b_cur;
      tau_b_save[save_idx] = tau_b_cur;

      save_idx++;
    }
  }

  return List::create(_["A"] = A_save, _["a"] = a_save, _["tau_a"] = tau_a_save, _["b"] = b_save, _["mu_b"] = mu_b_save, _["tau_b"] = tau_b_save, _["u"] = u_save, _["num"] = num,
                      _["state_last"] = List::create(_["z"] = z, _["z0"] = z0, _["u"] = u, _["a"] = a, _["b"] = b, _["A"] = A, _["tau_a"] = tau_a_cur, _["mu_b"] = mu_b_cur, _["tau_b"] = tau_b_cur
    )
  );
}
