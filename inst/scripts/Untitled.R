pkgload::unload("bayesodt")        # ← 이게 핵심(기존 DLL/namespace 정리)
Rcpp::compileAttributes()          # RcppExports 갱신(새 export/시그니처 변경 시 필수)
devtools::document()               # roxygen/NAMESPACE 갱신(주석 바꿨으면)
devtools::load_all()               # 다시 컴파일 + 로드


# path <- system.file("extdata", "simdata.rds", package="bayesodt")
# path
#
#
# simdata <- readRDS(path)
# dim(simdata)
#
# w <- simdata$w
# D <- simdata$D



# m <- dim(w)[1] # Num. of patients
# n <- dim(w)[2] # Num. of raters
# init <- list(z = matrix(1, nrow=m, ncol=n), z0 = rep(0, m), a = rnorm(n, 0, 1), b = rexp(n, 0.1), u = rnorm(m, 0, 1), A = c(-1, 0, 1))
# init2 <- list(z = matrix(1, nrow=m, ncol=n), z0 = rep(0, m), a = rnorm(n, 0, 1), b = rexp(n, 0.1), u = rnorm(m, 0, 1), A = c(-1, 0, 1, Inf))
# init3 <- list(z = matrix(1, nrow=m, ncol=n), z0 = rep(0, m), a = rnorm(n, 0, 1), b = rexp(n, 0.1), u = rnorm(m, 0, 1))


prior <- list(mu_a = 0, tau_a = 1, mu_b = 0, tau_b = 1, delta = 5, gamma = 0.1, tau = 0.1, mu_u = 0, tau_u = 1)

set.seed(33)
data <- simulate_data(m=90, n=50)
W <- data$W
D <- data$D

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
#
# set.seed(1)
# system.time({
#   fit <- BayesODT_draft(w=w, D=D, init=init2, prior=prior)})
#
#
# str(fit)
# head(fit$alpha)
# library(coda)
# m_alpha <- mcmc(fit$alpha)
# m_A <- mcmc(fit$A)
# m_B <- mcmc(fit$B)
# traceplot(m_alpha)
# traceplot(m_A)
# traceplot(m_B)

################
# set.seed(1)
system.time({
  fit <- BayesODT(w=W, D=D, n.chains = 2, n.iter=45000, n.burnin=20000, n.thin = 5)})#, n.chains=1, prior=prior, n.thin=5, n.burnin=2000, n.iter=4000)})

print(fit)
summary(fit)
summary(fit)$a

roc1 <- bayes_roc(fit, type = "ROC1"); str(roc1)

roc2 <- bayes_roc(fit, type = "ROC2"); str(roc2)
roc_ind <- bayes_roc(fit, type = "individual", rater = 2); str(roc_ind)
roc_sm <- bayes_roc(fit, type = "smoothed"); str(roc_sm)
roc_emp <- bayes_roc(fit, type = "pROC"); str(roc_emp)
roc_emp.sm <- bayes_roc(fit, type = "smoothed.pROC"); str(roc_emp.sm)

# plot
plot(fit, which="ABU")
par(mfrow = c(2, 2))
plot(fit, which="ROC", type="ROC1")
plot(fit, which="ROC", type="ROC2")
plot(fit, which="ROC", type= "smoothed")
plot(fit, which="ROC", type="individual", rater= c(1,2,10,14,17))
dev.off()
plot(fit, which="Udist")

# mcmc 또는 mcmc.list 형식으로 바꿔주는
m <- as.mcmc(fit, param="a", index=10)
m_all <- as.mcmc_all(fit, param=c("A", "tau_a", "mu_b", "tau_b"))

library(coda)

gelman.diag(m_all)
plot(m_all)

coda::traceplot(m)
# traceplot
par(mfrow = c(2, 3))
plot_trace(fit, param = "A")
plot_trace(fit, param = "a", index = 17)
plot_trace(fit, param = "b", index = 21)
plot_trace(fit, param = "u", index = 35)
plot_trace(fit, param = "tau_a")
plot_trace(fit, param = "mu_b")
plot_trace(fit, param = "tau_b")
dev.off()



# rater anaylsis
group_raters(fit)
print(m)
plot_rater_map(fit,
               show_intervals = FALSE,
               group = TRUE)


gof <- GoF(fit)
print(gof)

