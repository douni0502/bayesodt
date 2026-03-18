#' Fitting the Bayesain hierarchical model for ordinal classification processes
#'
#' @param w The classification data (m: num. of patients by n: num. of raters) matrix
#' @param D The true disease status data
#' @param n.iter Number of MCMC iterations
#' @param n.burnin Number of burnin iterations
#' @param n.thin Number of the thinning interval
#' @param ROC Optional. ROC=TRUE for plotting ROC curves. Default is FALSE.
#' @param init Setting initial values for z, z0, a, b, u, A. If z0 is empty, the model will not use Model (1) in the manuscript.
#' @param prior Setting hyper-parameters (mu_a, tau_a, mu_b, tau_b, delta, gamma, tau, mu_u, tau_u)
#' @return List of posterior samples of the parameters.
#' @export


BayesODT_draft <- function(w=w, D=D, init=init, prior=prior, n.iter=4000, n.burnin=2000, n.thin=5, ROC=FALSE){

  # REQUIRED LIBRARIES
  #  library(truncnorm)
  #  library(data.table)
  #  library(msm)
  #  library(ordinal)
  #  library(reshape)
  #  library(pROC)
  #  library(rootSolve)
  #  library(utils)


  #pre-processing
  a <- init$a
  b <- init$b
  u <- init$u
  z <- init$z
  z0 <- init$z0




  mu_a <- 0
  #mu_a <- prior$mu_a
  tau_a <- prior$tau_a
  mu_b <- prior$mu_b
  tau_b <- prior$tau_b
  delta <- prior$delta
  gamma <- prior$gamma
  tau <- prior$tau
  mu_u <- 0
  tau_u <- 1
  #  mu_u <- prior$mu_u
  #  tau_u <- prior$tau_u

  m <- length(D)
  n <- length(a)
  if(!is.null(init$A)){
    A <- init$A
  }else{
    fit.cut<-ordinal::clmm(as.factor(matrix(w, ncol=1)) ~ rep(D, n) +(1|c(sapply(1:n, function(x) rep(x, m))))+(1|c(sapply(1:m, function(x) rep(x, n)))), link = "probit", Hess=TRUE,
                  control = ordinal::clmm.control(maxIter = 100,maxLineIter = 100),threshold = "flexible", gradTol=1e-4)
    A <- c(fit.cut$alpha,Inf)
  }
  K <- length(A)

  ##################################
  ## Main MCMC Function ############
  ##################################

  MCMC <- function(z=z, z0=z0, u=u, w=w, D=D, a=a, b=b, A=A, mu_a=mu_a, tau_a=tau_a, mu_b=mu_b, tau_b=tau_b, mu_u = mu_u, tau_u = tau_u){

    if(is.null(z0)){
      # sample z_ij
      for(k in 1:K){
        W <- which(w==k, arr.ind=T)
        if(k==1){
          lower <- -Inf
        }else{
          lower <- A[k-1]
        }
        z[W] <- truncnorm::rtruncnorm(dim(W)[1], a=lower, b= A[k], mean=a[W[,2]]+b[W[,2]]*u[W[,1]], sd=1)
      }

      # sample u_i
      u <- rnorm(m,((z-matrix(a, nrow=m, ncol=n, byrow=TRUE))%*%b+tau_u*mu_u)/(sum(b^2)+tau_u), sqrt(1/(sum(b^2)+tau_u)))
      model <- function(c){
        F1 <- mean(pnorm(u+c)) -  mean(D)
        c(F1 = F1)
      }
      ss <- rootSolve::multiroot(f = model, start = c(1))
      u <- u + ss$root
    }else{
      # sample z_i0
      z0 <- (D)* msm::rtnorm(m, lower=0, upper=Inf, mean=u, sd=1)+(1-D)* msm::rtnorm(m, lower=-Inf, upper=0, mean=u, sd=1)

      # sample z_ij
      for(k in 1:K){
        W <- which(w==k, arr.ind=T)
        if(k==1){
          lower <- -Inf
        }else{
          lower <- A[k-1]
        }
        z[W] <- truncnorm::rtruncnorm(dim(W)[1], a=lower, b= A[k], mean=a[W[,2]]+b[W[,2]]*u[W[,1]], sd=1)
      }

      # sample u_i
      u <- rnorm(m,((z-matrix(a, nrow=m, ncol=n, byrow=TRUE))%*%b+z0+tau_u*mu_u)/(sum(b^2)+1+tau_u), sqrt(1/(sum(b^2)+1+tau_u)))
      model <- function(c){
        F1 <- mean(pnorm(u+c)) -  mean(D)
        c(F1 = F1)
      }
      ss <- rootSolve::multiroot(f = model, start = c(1))
      u <- u + ss$root
    }

    # sample a_j
    a <- rnorm(n,(colSums(z-u%*%t(b))+mu_a*tau_a)/(m+tau_a), sqrt(1/(m+tau_a)))

    # sample b_j
    b <- rnorm(n, (u%*%(z-matrix(a, nrow=m, ncol=n, byrow=T))+mu_b*tau_b)/(sum(u^2)+tau_b), sqrt(1/(sum(u^2)+tau_b)))

    # sample mu_a
    # mu_a <- rnorm(1, tau_a*sum(a)/(n*tau_a+tau) , sqrt(1/(n*tau_a+tau)))

    # sample tau_a
    tau_a <- rgamma(1, gamma+n/2, gamma+0.5*sum((a-mu_a)^2))

    # sample mu_b
    mu_b <- rnorm(1, tau_b*sum(b)/(n*tau_b+tau) , sqrt(1/(n*tau_b+tau)))

    # sample tau_b
    tau_b <- rgamma(1, gamma+n/2, gamma+0.5*sum((b-mu_b)^2))

    # sample A (alpha)
    U <- NULL
    L <- NULL
    # sample A1
    L[1] <- max(-delta,A[2]-delta,max(z[which(w==1, arr.ind=T)]))
    U[1] <- min(delta,A[2],min(z[which(w==2, arr.ind=T)]))
    A[1] <- runif(1, L[1], U[1])
    #A[1] <- A[1]
    # sample A2 to Ak-2
    for(k in 2:(K-2)){
      U[k] <- min(A[k-1]+delta,A[k+1],min(z[which(w==k+1, arr.ind=T)]))
      L[k] <- max(A[k-1],A[k+1]-delta,max(z[which(w==k, arr.ind=T)]))
      A[k] <- runif(1, L[k], U[k])
    }
    U[K-1] <- min(A[K-2]+delta,A[K],min(z[which(w==K, arr.ind=T)]))
    L[K-1] <- max(A[K-2],max(z[which(w==K-1, arr.ind=T)]))
    A[K-1] <- runif(1, L[K-1], U[K-1])

    if(ROC){
      sp <- as.matrix(sapply(1:K, function(x) (1/m)*(pnorm(A[x]-(matrix(a, ncol=m, nrow=n, byrow=FALSE)+outer(b,u))))%*%(1-pnorm(u))/mean(1-pnorm(u)) ) )
      se <- as.matrix(sapply(1:K, function(x) (1/m)*(pnorm(A[x]-(matrix(a, ncol=m, nrow=n, byrow=FALSE)+outer(b,u))))%*%(pnorm(u))/mean(pnorm(u)) ) )

      sp.mean <- se.mean <- NULL
      for(k in 1:K){
        sp.mean[k] <- (1/m)*(pnorm(A[k]-(mean(a)+outer(mean(b),u))))%*%(1-pnorm(u))/mean(1-pnorm(u))
        se.mean[k] <- (1/m)*(pnorm(A[k]-(mean(a)+outer(mean(b),u))))%*%(pnorm(u))/mean(pnorm(u))
      }
    }else{
      sp <- se <- sp.mean <- se.mean <- NULL
    }

    return(list(z0=z0,z=z,
                u=u,
                w=w,
                D=D,
                a=a,
                b=b,
                mu_a=mu_a,
                tau_a=tau_a,
                mu_b=mu_b,
                tau_b=tau_b,
                mu_u=mu_u,
                tau_u=tau_u,
                A=A,
                sp=sp,
                se=se,
                sp.mean=sp.mean,
                se.mean=se.mean))
  }

  update <- list()
  update[[1]] <- list(z0=z0,z=z,
                      u=u,
                      w=w,
                      D=D,
                      a=a,
                      b=b,
                      mu_a=mu_a,
                      tau_a=tau_a,
                      mu_b=mu_b,
                      tau_b=tau_b,
                      mu_u=mu_u,
                      tau_u=tau_u,
                      A=A)


  for(t in 2:n.iter){
    update[[t]] <- MCMC(z0=update[[t-1]]$z0,
                        z=update[[t-1]]$z,
                        u=update[[t-1]]$u,
                        w=update[[t-1]]$w,
                        D=update[[t-1]]$D,
                        a=update[[t-1]]$a,
                        b=update[[t-1]]$b,
                        mu_a=update[[t-1]]$mu_a,
                        tau_a=update[[t-1]]$tau_a,
                        mu_b=update[[t-1]]$mu_b,
                        tau_b=update[[t-1]]$tau_b,
                        A=update[[t-1]]$A,
                        mu_u=update[[t-1]]$mu_u,
                        tau_u=update[[t-1]]$tau_u)
  }

  num <- (n.iter-n.burnin)/n.thin
  out <- list()
  for(ind in 1:num){
    IND <- ind*n.thin+n.burnin
    out[[ind]] <- list(z0=update[[IND]]$z0,
                       z=update[[IND]]$z,
                       a=update[[IND]]$a,
                       b=update[[IND]]$b,
                       u=update[[IND]]$u,
                       A=update[[IND]]$A,
                       mu_a=update[[IND]]$mu_a,
                       mu_b=update[[IND]]$mu_b,
                       tau_a=update[[IND]]$tau_a,
                       tau_b=update[[IND]]$tau_b,
                       mu_a=update[[IND]]$mu_a,
                       mu_b=update[[IND]]$mu_b,
                       tau_a=update[[IND]]$tau_a,
                       tau_b=update[[IND]]$tau_b,
                       sp=update[[IND]]$sp,
                       se=update[[IND]]$se,
                       sp.mean=update[[IND]]$sp.mean,
                       se.mean=update[[IND]]$se.mean)
  }


  #  MU.A <- NULL
  MU.B <- NULL
  TAU.A <- NULL
  TAU.B <- NULL
  alpha <- matrix(nrow=num, ncol=K-1)
  U <- matrix(nrow=num, ncol=m)
  B <- matrix(nrow=num, ncol=n)
  A <- matrix(nrow=num, ncol=n)
  for(i in 1:num){
    alpha[i,] <- out[[i]]$A[1:(K-1)]
    #    MU.A[i] <- out[[i]]$mu_a
    TAU.A[i] <- out[[i]]$tau_a
    MU.B[i] <- out[[i]]$mu_b
    TAU.B[i] <- out[[i]]$tau_b
    U[i,] <- out[[i]]$u
    B[i,] <- out[[i]]$b
    A[i,] <- out[[i]]$a
  }

  if(ROC){
    sp <- se <- array(dim=c(n, K, num))
    sp.mean <- se.mean <- matrix(nrow=num, ncol=K)

    for(i in 1:num){
      sp[,,i] <- (1-out[[i]]$sp)
      se[,,i] <- (1-out[[i]]$se)
    }
    for(i in 1:num){
      sp.mean[i,] <- (1-out[[i]]$sp.mean)
      se.mean[i,] <- (1-out[[i]]$se.mean)
    }
  }else{
    sp <- se <- sp.mean <- se.mean <- NULL
  }

  OUT <- list(alpha=alpha,
              A=A,
              TAU.A=TAU.A,
              B=B,
              MU.B=MU.B,
              TAU.B=TAU.B,
              U=U,
              num=num,
              D=D,
              w=w,
              sp=sp,
              se=se,
              sp.mean=sp.mean,
              se.mean=se.mean)

  OUT$call <- match.call()
  class(OUT) <- "BayesODT"
  return(OUT)
}


