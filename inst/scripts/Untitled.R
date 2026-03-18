pkgload::unload("bayesodt")        # ← 이게 핵심(기존 DLL/namespace 정리)
Rcpp::compileAttributes()          # RcppExports 갱신(새 export/시그니처 변경 시 필수)
devtools::document()               # roxygen/NAMESPACE 갱신(주석 바꿨으면)
devtools::load_all()               # 다시 컴파일 + 로드


path <- system.file("extdata", "simdata.rds", package="bayesodt")
path


simdata <- readRDS(path)
dim(simdata)

w <- simdata$w
D <- simdata$D



m <- dim(w)[1] # Num. of patients
n <- dim(w)[2] # Num. of raters
init <- list(z = matrix(1, nrow=m, ncol=n), z0 = rep(0, m), a = rnorm(n, 0, 1), b = rexp(n, 0.1), u = rnorm(m, 0, 1), A = c(-1, 0, 1))
init2 <- list(z = matrix(1, nrow=m, ncol=n), z0 = rep(0, m), a = rnorm(n, 0, 1), b = rexp(n, 0.1), u = rnorm(m, 0, 1), A = c(-1, 0, 1, Inf))
init3 <- list(z = matrix(1, nrow=m, ncol=n), z0 = rep(0, m), a = rnorm(n, 0, 1), b = rexp(n, 0.1), u = rnorm(m, 0, 1))


prior <- list(mu_a = 0, tau_a = 1, mu_b = 0, tau_b = 1, delta = 5, gamma = 0.1, tau = 0.1, mu_u = 0, tau_u = 1)


# 결측 만들 때
make_MCAR <- function(w, missing_rate = 0.2) {
  w_miss <- w
  miss_idx <- sample(seq_len(length(w)),
                     size = round(length(w) * missing_rate),
                     replace = FALSE)
  w_miss[miss_idx] <- NA
  w_miss
}
w <- make_MCAR(w, 0.2)






########### bayesodt 초안 실험 (교수님) ######################

set.seed(1)
system.time({
  fit <- BayesODT_draft(w=w, D=D, init=init3, prior=prior)})
str(fit)
head(fit$alpha)
library(coda)
m_alpha <- mcmc(fit$alpha)
m_A <- mcmc(fit$A)
m_B <- mcmc(fit$B)
traceplot(m_alpha)
traceplot(m_A)
traceplot(m_B)

################
set.seed(1)
system.time({
  fit <- BayesODT(w=w, D=D, n.chains=1, prior=prior, n.thin=5, n.burnin=2000, n.iter=4000)})

print(fit)
summary(fit)
summary(fit)$a

roc1 <- bayes_roc(fit, type = "ROC1")
str(roc1)

roc2 <- bayes_roc(fit, type = "ROC2"); str(roc2)
roc_ind <- bayes_roc(fit, type = "individual", rater = 2); str(roc_ind)
roc_sm <- bayes_roc(fit, type = "smoothed"); str(roc_sm)
roc_emp <- bayes_roc(fit, type = "pROC"); str(roc_emp)
roc_emp.sm <- bayes_roc(fit, type = "smoothed.pROC"); str(roc_emp.sm)

# plot
plot(fit, which="ABU")
plot(fit, which="ROC", type="individual", rater=10)
plot(fit, which="ROC", type="individual", rater= c(1,2,10,14,17))
plot(fit, which="Udist")

# mcmc 또는 mcmc.list 형식으로 바꿔주는
m <- as.mcmc(fit, param = "A", index = 1)
coda::traceplot(m)
# traceplot
plot_trace(fit, param = "A")
plot_trace(fit, param = "a", index = c(1,2,3))
plot_trace(fit, param = "b", index = c(1, 2, 3))
plot_trace(fit, param = "u", index = c(5, 10))
plot_trace(fit, param = "tau_a")
plot_trace(fit, param = "mu_b")
plot_trace(fit, param = "tau_b")


# effective sample size
effective_size(fit, param = "A", index = 1:3)
effective_size(fit, param = "a", index = 1:5)


# autocorr plot
plot_autocorr(fit, param = "A", index = 1:3)
plot_autocorr(fit, param = "b", index = c(1, 2, 3))


# Geweke_diag
geweke_diag(fit, param = "A", index = c(1, 2, 3), frac1 = 0.1, frac2 = 0.1)


# heidel_diag
heidel_diag(fit, param = "A", index = c(1,2,3))

# density plot
plot_density(fit, param = "a", index = c(1, 2, 3, 4, 5, 6, 7, 8))


# rater anaylsis
plot_rater_map(fit)
plot_rater_map(fit, color_by = "group", show_intervals = FALSE)
plot_rater_map(fit, color_by = "cluster", show_intervals = FALSE)
head(group_raters(fit))
head(cluster_raters(fit))

cl <- cluster_raters(fit)
attr(cl, "k_selected")
attr(cl, "silhouette_scores")
