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

tool_file     <- "~/Desktop/Changepoint/aop/TSAOP/CLSEinnov/tools.R"
tool_file_ar2 <- "~/Desktop/Changepoint/aop/TSAOP/CLSEinnov/tools_ar2.R"

source(tool_file)
source(tool_file_ar2)

options("scipen"=10) # show all digits

fit_CLSInnov_AR2 <- function(X_hour, X_hour_wide, DesignXEst,
                             lag_max = 50L,
                             lag_tol = 1e-8) {
  
  K <- ncol(X_hour_wide)
  q <- ncol(DesignXEst)
  
  obj <- make_aop2_cls_optim(
    Xw = X_hour_wide,
    DesignX = DesignXEst,
    lag_max = lag_max,
    lag_tol = lag_tol,
    jitter = 1e-10,
    do_grad = TRUE
  )
  
  cout <- tabulate(X_hour, nbins = K)
  cumprob <- cumsum(cout / sum(cout))
  cumprob <- pmin(pmax(cumprob, 1e-6), 1 - 1e-6)
  
  ci_initial <- qnorm(cumprob)[1:(K - 1)]
  
  cuts_initial <- ci_initial[2:(K - 1)] - ci_initial[1]
  cuts_initial <- pmax(cuts_initial, 1e-5)
  cuts_initial <- sort(cuts_initial)
  
  theta_initial <- c(-ci_initial[1], rep(0, q - 1L))
  
  yy <- as.numeric(scale(X_hour, center = TRUE, scale = FALSE))
  
  phi_initial <- tryCatch({
    arfit <- stats::arima(yy, order = c(2, 0, 0), include.mean = FALSE)
    tmp <- as.numeric(arfit$coef[1:2])
    if (!is_stationary_ar2(tmp[1], tmp[2])) c(0.3, 0.2) else tmp
  }, error = function(e) {
    c(0.3, 0.2)
  })
  
  par_initial <- c(cuts_initial, theta_initial, phi_initial)
  
  con <- make_constraints_aop2(
    K = K,
    q = q,
    delta = 1e-6,
    delta0 = 1e-6,
    eps_ar = 1e-6
  )
  
  fit <- constrOptim(
    theta = par_initial,
    f = obj$fn,
    grad = obj$gr,
    ui = con$ui,
    ci = con$ci,
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 300, reltol = 1e-7)
  )
  
  
  seEst <- sandwichSE_AOP2(
    obj = obj,
    fit = fit,
    Xw = X_hour_wide,
    DesignX = DesignXEst,
    lag_max = lag_max,
    lag_tol = lag_tol
  )
  
  se <- seEst$se
  
  RES_mat <- as.data.frame(cbind(fit$par, se))
  colnames(RES_mat) <- c("Estimate", "Std. Error")
  rownames(RES_mat) <- c(
    paste0("c", 2:(K - 1)),
    colnames(DesignXEst),
    "phi1",
    "phi2"
  )
  
  RES_mat$z.value <- RES_mat$Estimate / RES_mat$`Std. Error`
  RES_mat$p.value <- 2 * pnorm(-abs(RES_mat$z.value))
  RES_mat$ci.lb <- RES_mat$Estimate - 1.96 * RES_mat$`Std. Error`
  RES_mat$ci.ub <- RES_mat$Estimate + 1.96 * RES_mat$`Std. Error`
  
  list(
    RES_mat = RES_mat,
    paramEst = fit$par,
    se = se,
    seEst = seEst,
    convergence = fit$convergence,
    value = fit$value,
    par_initial = par_initial,
    fit = fit,
    obj = obj
  )
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

####################################################
#####      Training / forecasting split        #####
####################################################

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

####################################################
#####      Fit innovations-based CLS           #####
####################################################

library(TSAOP)
library(pbivnorm)

lag_max_fit <- 50L
lag_tol_fit <- 1e-8


t1 <- Sys.time()

fit_ar2 <- fit_CLSInnov_AR2(
  X_hour = X_hour_obs,
  X_hour_wide = X_hour_wide_obs,
  DesignXEst = DesignXEst_obs,
  lag_max = lag_max_fit,
  lag_tol = lag_tol_fit
)

t2 <- Sys.time()


t2 - t1
# Time difference of 4.332215 mins

RES_ar2 <- fit_ar2$RES_mat
round(RES_ar2, 4)

#           Estimate Std. Error z.value p.value   ci.lb  ci.ub
# c2          0.5674     0.0789  7.1880  0.0000  0.4127 0.7221
# c3          0.9781     0.0938 10.4266  0.0000  0.7943 1.1620
# c4          1.2867     0.1037 12.4084  0.0000  1.0835 1.4900
# c5          1.6987     0.1174 14.4642  0.0000  1.4685 1.9289
# Intercept   0.9210     0.1091  8.4434  0.0000  0.7072 1.1348
# phi1        0.2250     0.0920  2.4456  0.0145  0.0447 0.4054
# phi2       -0.0610     0.0774 -0.7880  0.4307 -0.2127 0.0907

