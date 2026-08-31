rm(list = ls())

library(lubridate)
library(TSA)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(readr)
library(ggplot2)
library(xtable)
library(car)
library(ggpubr)
library(astsa)

library(TSAOP)
library(pbivnorm)

source("~/Desktop/Changepoint/aop/TSAOP/CLSEinnov/tools.R")
options("scipen"=10) # show all digits


fit_CLSInnov <- function(X_hour, X_hour_wide, DesignXEst){
  
  K <- ncol(X_hour_wide)
  q <- ncol(DesignXEst)
  
  obj <- make_aop1_cls_optim(Xw = X_hour_wide, DesignX = DesignXEst,
                             lag_max = 50L, lag_tol = 1e-08,
                             jitter = 1e-10,
                             do_grad = TRUE)
  
  cout <- summary(as.factor(X_hour))
  ci_initial <- as.vector(qnorm(cumsum(cout / sum(cout))))[1:(K - 1)]
  phi_initial <- stats::pacf(X_hour, plot = FALSE)$acf[1]
  par_initial <- c(ci_initial[2:(K-1)] - ci_initial[1], -ci_initial[1],
                   rep(0, NCOL(DesignXEst) - 1), phi_initial)
  
  con <- make_constraints_aop1(K = K, q = q, delta = 1e-6,
                               delta0 = 1e-6, eps_rho = 1e-6)
  fit <- constrOptim(theta = par_initial,
                     f = obj$fn, grad = obj$gr,
                     ui = con$ui, ci = con$ci, method = "BFGS", hessian = TRUE,
                     control = list(maxit = 200))
  paramEst <- fit$par
  
  seEst <- sandwichSE(obj=obj, fit=fit, Xw=X_hour_wide, DesignX=DesignXEst)
  # seEst$grad_check
  
  RES_mat <- as.data.frame(cbind(paramEst, seEst$se))
  colnames(RES_mat) <- c("Estimate", "Std. Error")
  rownames(RES_mat) <- c(paste0("c", 2:(K-1)), colnames(DesignXEst), "rho")
  
  RES_mat$z.value <- RES_mat$Estimate / RES_mat$`Std. Error`
  RES_mat$p.value <- 2 * pnorm(-abs(RES_mat$z.value))
  
  RES_mat$ci.lb <- RES_mat$Estimate - 1.96*RES_mat$`Std. Error`
  RES_mat$ci.ub <- RES_mat$Estimate + 1.96*RES_mat$`Std. Error`
  
  return(RES_mat)
  
}

#############################
##### HVAC quality data #####
#############################

data.full <- c(
  6,1,2,6,5,1,5,6,1, 6,2,2,4,4,3,5,3,2, 1,3,4,5,2,5,4,4,6,
  2,2,1,2,6,2,2,4,6, 1,5,1,1,1,6,6,5,1, 6,6,2,2,3,4,6,6,5,
  4,6,4,2,2,2,6,4,6, 5,6,5,2,2,3,2,6,4, 6,4,6,5,2,2,3,2,3,
  1,2,2,6,5,6,2,1,6, 6,5,3,1,1,2,2,3,3, 3,4,6,6,3,4,5,4,1,
  5,4,2,5,5,6,1,1,1, 3,4,6,3,1,6,5,3,6, 3,2,2,1,2,1,6,2,3,
  1,1,2,3,1,1,3,5,4, 3,1,2,4,5,5,6,2,1, 5,2,5,6,5,6,2,1,6,
  6,5,3,4,4,6,5,4,3, 2,6,6,6,3,5,1,1,1, 1,3,4,3,3,3,1,6,2,
  3,2,1,1,1,6,5,5,2, 1,3,1,6,6,2,6,3,4, 3,2,5,2,6,6,6,3,5,
  4,4,3,6,3,3,5,4,3, 5,1,1,1,6,6,6,6,3, 4,6,3,2,5,6,6,5,5,
  5,3,4,2,5,6,5,3,3, 4,4,4,2,6,3,5,4,2, 5,4,2,2,1,6,2,1,4,
  6,4,6,6,6,3,2,4,4, 1,1,2,3,2,4,6,4,3, 1,1,2,4,5,4,6,2,2,
  5,6,5,1,2,3,6,6,6, 6,6,6,6,6,4,5,4,6, 4,6,5,4,4,6,4,3,5,
  5,6,6,6,6,6,6,6,6, 6,6,1,1,1,6,6,5,2, 4,5,1,5,3,3,4,3,4,
  4,1,6,6,3,5,2,3,6, 4,5,6,6,6,3,6,5,2, 6,5,6,6,6,4,6,6,1,
  1,6,6,6,6,5,1,6,6, 5,6,6,6,4,2,5,6,6, 6,5,6,5,5,4,6,3,6,
  3,6,6,5,6,6,5,6,2, 4,2,1,6,1,1
)


tauT <- 262

data.use <- data.full[1:(tauT - 1)]


## Training / forecasting split


n_train <- 240

X_hour_obs  <- data.use[1:n_train]
X_hour_pred <- data.use[(n_train + 1):length(data.use)]

K <- max(data.use)

make_wide <- function(x, K) {
  Xw <- matrix(0, nrow = length(x), ncol = K)
  for (t in seq_along(x)) {
    Xw[t, x[t]] <- 1
  }
  Xw
}

X_hour_wide_obs  <- make_wide(X_hour_obs, K)
X_hour_wide_pred <- make_wide(X_hour_pred, K)

DesignXEst_obs <- matrix(1, nrow = length(X_hour_obs), ncol = 1)
colnames(DesignXEst_obs) <- "Intercept"

DesignXEst_pred <- matrix(1, nrow = length(X_hour_pred), ncol = 1)
colnames(DesignXEst_pred) <- "Intercept"

fit1 <- fit_CLSInnov(
  X_hour = X_hour_obs,
  X_hour_wide = X_hour_wide_obs,
  DesignXEst = DesignXEst_obs
)

round(fit1, 4)
#           Estimate Std. Error z.value p.value  ci.lb  ci.ub
# c2          0.5658     0.0783  7.2307  0.0000 0.4124 0.7192
# c3          0.9753     0.0930 10.4889  0.0000 0.7931 1.1576
# c4          1.2824     0.1024 12.5243  0.0000 1.0817 1.4831
# c5          1.6974     0.1175 14.4481  0.0000 1.4672 1.9277
# Intercept   0.9125     0.1104  8.2631  0.0000 0.6960 1.1289
# rho         0.2122     0.0916  2.3153  0.0206 0.0326 0.3918

save.image("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/HVAC_InnovFit_nocpt_forecast.RData")





load("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/HVAC_InnovFit_nocpt_forecast.RData")

forecast_probs_innov <- function(par, Xw, DesignX, DesignX_future,
                                 lag_max = 50L, lag_tol = 1e-8,
                                 jitter = 1e-10, symmetrize = TRUE) {
  
  Xw <- as.matrix(Xw)
  DesignX <- as.matrix(DesignX)
  DesignX_future <- as.matrix(DesignX_future)
  
  Tt <- nrow(Xw)
  H  <- nrow(DesignX_future)
  K  <- ncol(Xw)
  q  <- ncol(DesignX)
  m  <- K - 1L
  
  cuts <- as.numeric(par[1:(K - 2L)])
  beta <- as.numeric(par[(K - 1L):(K - 2L + q)])
  rho  <- as.numeric(par[K - 2L + q + 1L])
  
  seg <- c(-Inf, 0, cuts, Inf)
  
  DesignX_all <- rbind(DesignX, DesignX_future)
  mst_all <- as.numeric(DesignX_all %*% beta)
  
  Fmat <- pnorm(outer(mst_all, seg, function(mu, s) s - mu))
  Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] -
    Fmat[, 1:K, drop = FALSE]
  
  EXW_all <- Pmat[, 1:m, drop = FALSE]
  EXW_obs <- EXW_all[1:Tt, , drop = FALSE]
  EXW_fut <- EXW_all[(Tt + 1):(Tt + H), , drop = FALSE]
  
  u_obs <- Xw[, 1:m, drop = FALSE] - EXW_obs
  
  kappa_all <- function(i, j) {
    multikappa_pbivnorm(
      i, j,
      msti = mst_all[i], mstj = mst_all[j],
      EXwi = EXW_all[i, ], EXwj = EXW_all[j, ],
      seg  = seg,
      rho  = rho
    )
  }
  
  u_fut <- matrix(0, nrow = H, ncol = m)
  u_pad <- rbind(u_obs, u_fut)
  
  fit_all <- innovations_mv(
    kappa = kappa_all,
    u = u_pad,
    jitter = jitter,
    symmetrize = symmetrize,
    compute_innov = TRUE,
    lag_max = lag_max,
    lag_tol = lag_tol,
    ahead = NULL
  )
  
  Theta_all <- fit_all$Theta
  e_obs <- fit_all$innov[1:Tt, , drop = FALSE]
  
  phat_fut <- matrix(NA_real_, nrow = H, ncol = K)
  
  for (h in seq_len(H)) {
    
    t_idx <- Tt + h
    coefs <- Theta_all[[t_idx]]
    Lt <- length(coefs)
    
    mu_u <- rep(0, m)
    jmax <- min(Lt, t_idx - 1L)
    
    if (jmax >= h) {
      for (j in h:jmax) {
        mu_u <- mu_u + coefs[[j]] %*% e_obs[t_idx - j, ]
      }
    }
    
    p_m_raw <- as.numeric(EXW_fut[h, ]) + as.numeric(mu_u)
    p_full_raw <- c(p_m_raw, 1 - sum(p_m_raw))
    
    p_full <- pmin(1, pmax(0, p_full_raw))
    sprob <- sum(p_full)
    
    if (sprob > 0) {
      p_full <- p_full / sprob
    }
    
    phat_fut[h, ] <- p_full
  }
  
  phat_fut
}

## Rolling one-step-ahead forecast

paramEst <- fit1$Estimate
H <- length(X_hour_pred)
prob_roll <- matrix(NA_real_, nrow = H, ncol = K)

X_hist <- X_hour_obs
Design_hist <- DesignXEst_obs

for (h in seq_len(H)) {
  
  Xw_hist <- make_wide(X_hist, K)
  
  prob_h <- forecast_probs_innov(
    par = paramEst,
    Xw = Xw_hist,
    DesignX = Design_hist,
    DesignX_future = DesignXEst_pred[h, , drop = FALSE]
  )
  
  prob_roll[h, ] <- prob_h[1, ]
  
  X_hist <- c(X_hist, X_hour_pred[h])
  Design_hist <- rbind(Design_hist, DesignXEst_pred[h, , drop = FALSE])
}

score_roll <- as.numeric(prob_roll %*% seq_len(K))

## Plot observed categories + forecast score

x_obs <- seq_along(X_hour_obs)
x_fut <- length(X_hour_obs) + seq_len(H)
x_all <- c(x_obs, x_fut)

y_all <- c(X_hour_obs, X_hour_pred)

score_fit <- as.numeric(Predw %*% seq_len(K))



setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/HVACFigFitline.eps", 
           width = 14, height = 7)

par(cex.lab = 1.45,
    cex.axis = 1.45,
    mfrow = c(1, 1))

plot(x = x_all, y = y_all,
     type = "l",
     col = "grey65",
     lwd = 1,
     xlab = "Time",
     ylab = "HVAC Quality category",
     xaxt = "n",
     yaxt = "n",
     xlim = c(-5, length(y_all) + 5),
     ylim = c(1, K + 1.5),
     xaxs = "i")

axis(1, at = pretty(x_all), labels = pretty(x_all), las = 1)
axis(2, at = 1:K)

lines(x_obs, score_fit, type = "l", lty = "solid", col = "blue", lwd = 2.5)
lines(x_fut, score_roll, type = "l", lty = "solid", col = "red4", lwd = 2.5)
abline(v = length(X_hour_obs), lty = 2, lwd = 2.5, col = "red")

legend(
  "topleft",
  legend = c("Observed HVAC quality category",
             "Average fitted categorical response",
             "Rolling one-step-ahead forecast categorical response"),
  col = c("grey65", "blue", "red4"),
  lty = c("solid", "solid", "solid"),
  lwd = 2.5,
  bty = "n",
  horiz = FALSE,
  seg.len = 3,
  cex=1.45
)

dev.off()
