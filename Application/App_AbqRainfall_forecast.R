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

# library(TSAOP)
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

# 0. Read and subset data
Abq <- read_csv("~/Downloads/Abq19310301_20201211.csv", col_types = cols(DATE = col_date(format = "%Y-%m-%d")))
Abq <- as.data.frame(Abq[,c("DATE", "PRCP")])
# select data range by changing starting and ending values
pstart <- which(as.character(Abq$DATE)== "1990-01-01")
pend   <- which(as.character(Abq$DATE)== "2000-12-31")
mydat <- Abq[pstart:pend,]

# 1. Missing values Detection and Imputation

# 1.1 Date missing
FullSeq <- data.frame(DateFull=seq.Date(from = min(mydat$DATE), to = max(mydat$DATE), by = 1))
Missing <- FullSeq$DateFull[!FullSeq$DateFull %in% mydat$DATE]
Missing

# Imputation (& delete leap days)
if(length(Missing)>0){
  # Leave the missing values slots
  mydat <- merge(mydat, FullSeq, by.x= "DATE", by.y="DateFull", all.y = TRUE)
  # To obtain an exact period of 365 days,
  # delete leap days observations (February 29).
  p2 <- which(format(mydat$DATE, format = "%m%d") == "0229")
  mydat <- mydat[-p2,]
  # calculation 365 day's average among 11 years
  DaysMean <- data.frame(date=format(seq.Date(from = as.Date("2019-01-01", "%Y-%m-%d"),
                                              to   = as.Date("2019-12-31", "%Y-%m-%d"),
                                              by = 1),
                                     format="%m%d"),
                         mean=colMeans(t(matrix(mydat$PRCP, nrow=365)), na.rm = T))
  # Impute day's average for missing date
  p3 <- which(is.na(mydat$PRCP))
  # by days mean
  for(i in 1:length(p2)){
    pp <- which(DaysMean$date == format(mydat$DATE[p2[i]], format="%m%d"))
    mydat$PRCP[p2[i]] <- DaysMean$mean[pp]
  }
  # # by Overall mean
  # p4 <- which(is.na(mydat$PRCP))
  # mydat$PRCP[p3] <- mean(mydat$PRCP, na.rm = T)
}else{
  # To obtain an exact period of 365 days,
  # delete leap days observations (February 29).
  p2 <- which(format(mydat$DATE, format = "%m%d") == "0229")
  mydat <- mydat[-p2,]
}

# 1.2 Value Missing
p4 <- which(is.na(mydat$PRCP))
p4

# 1.3 Discretized to 3 categories
mydat$Catg <- 1
mydat$Catg[which(mydat$PRCP > 0 & mydat$PRCP < 0.2)] <- 2
mydat$Catg[which(mydat$PRCP >=0.2)]                  <- 3



####################################################
##### Getting 1990-1999 data for model fitting #####
####################################################

flag <- which(mydat$DATE > "1999-12-31")[1]

mydat_obs  <- mydat[1:(flag - 1), ]
mydat_pred <- mydat[flag:NROW(mydat), ]

X_hour <- mydat$Catg

Ts <- length(X_hour)
K  <- max(X_hour)

make_wide <- function(x, K) {
  Xw <- matrix(0, nrow = length(x), ncol = K)
  for (t in seq_along(x)) {
    Xw[t, x[t]] <- 1
  }
  Xw
}

X_hour_wide <- make_wide(X_hour, K)

ss1 <- 365
ss2 <- 183

TrendValue <- (1:Ts) / Ts
BValue <- cos(2 * pi * (1:Ts) / ss1)
DValue <- sin(2 * pi * (1:Ts) / ss1)
EValue <- cos(2 * pi * (1:Ts) / ss2)
FValue <- sin(2 * pi * (1:Ts) / ss2)

CpLocY <- 2366
CpValueY <- c(rep(0, CpLocY - 1), rep(1, Ts - CpLocY + 1))

DesignXEst <- cbind(1, TrendValue, BValue, DValue, EValue, FValue, CpValueY)
colnames(DesignXEst) <- c("Intercept", "alpha", "B", "D", "E", "F", "Delta")

## Training / forecasting split

X_hour_obs       <- X_hour[1:(flag - 1)]
X_hour_pred      <- X_hour[flag:Ts]

X_hour_wide_obs  <- X_hour_wide[1:(flag - 1), , drop = FALSE]
X_hour_wide_pred <- X_hour_wide[flag:Ts, , drop = FALSE]

DesignXEst_obs   <- DesignXEst[1:(flag - 1), , drop = FALSE]
DesignXEst_pred  <- DesignXEst[flag:Ts, , drop = FALSE]

fit <- fit_CLSInnov(X_hour = X_hour_obs, X_hour_wide = X_hour_wide_obs, DesignXEst = DesignXEst_obs)

round(fit, 4)
#           Estimate Std. Error z.value p.value   ci.lb   ci.ub
# c2          0.8170     0.0452 18.0888  0.0000  0.7285  0.9056
# Intercept  -0.7272     0.0763 -9.5338  0.0000 -0.8767 -0.5777
# alpha      -0.7392     0.2254 -3.2798  0.0010 -1.1810 -0.2975
# B          -0.1875     0.0488 -3.8406  0.0001 -0.2832 -0.0918
# D          -0.1368     0.0497 -2.7521  0.0059 -0.2343 -0.0394
# E           0.0864     0.0492  1.7548  0.0793 -0.0101  0.1829
# F           0.1654     0.0499  3.3152  0.0009  0.0676  0.2632
# Delta       0.3182     0.1286  2.4740  0.0134  0.0661  0.5702
# rho         0.4145     0.0295 14.0373  0.0000  0.3566  0.4724

save.image("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/AbqRainfall_InnovFit_forecast.RData")




load("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/AbqRainfall_InnovFit_forecast.RData")

## Fitted probabilities on the training period

par     <- fit$Estimate
Xw      <- X_hour_wide_obs
DesignX <- DesignXEst_obs

Ts_obs <- nrow(Xw)
K      <- ncol(Xw)
m      <- K - 1L
num_par <- length(par)

cuts  <- as.numeric(par[1:(K - 2L)])
theta <- as.numeric(par[(K - 1L):(num_par - 1L)])
rho   <- as.numeric(par[num_par])

seg <- c(-Inf, 0, cuts, Inf)

mst <- as.numeric(DesignX %*% theta)

Fmat <- pnorm(outer(mst, seg, function(mu, s) s - mu))
Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] -
  Fmat[, 1:K, drop = FALSE]

EXW <- Pmat[, 1:m, drop = FALSE]
u <- Xw[, 1:m, drop = FALSE] - EXW

kappa <- function(i, j) {
  multikappa_pbivnorm(
    i, j,
    msti = mst[i],
    mstj = mst[j],
    EXwi = EXW[i, ],
    EXwj = EXW[j, ],
    seg  = seg,
    rho  = rho
  )
}

InnovRes <- innovations_mv(
  kappa = kappa,
  u     = u,
  jitter = 1e-10,
  symmetrize = TRUE,
  compute_innov = TRUE,
  lag_max = 50L,
  lag_tol = 1e-08,
  ahead   = NULL
)

PredX <- EXW + InnovRes$uhat
Predw <- cbind(PredX, 1 - rowSums(PredX))

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


## Rolling one-step-ahead forecasting

paramEst <- fit$Estimate

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

## Plot observed rainfall + forecast score

df_obs <- data.frame(
  date = mydat_obs$DATE,
  cat  = as.numeric(X_hour_obs),
  Ehat = as.numeric(Predw %*% seq_len(K))
)

df_fut <- data.frame(
  date = mydat_pred$DATE,
  cat  = as.numeric(X_hour_pred),
  Ehat = score_roll
)

df_all <- data.frame(
  date = mydat$DATE,
  cat  = as.numeric(mydat$Catg)
)

yrs <- seq(year(min(df_all$date)), year(max(df_all$date)), by = 1)
at_dates <- as.Date(paste0(yrs, "-01-01"))

setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/RainfallFigFitline_innov.eps",
           width = 16, height = 8)

par(cex.lab = 2,
    cex.axis = 2,
    mar = c(5, 5, 4, 0.5),
    mfrow = c(1, 1))

plot(range(df_all$date), c(1, K + 0.7),
     type = "n",
     xlab = "Year",
     ylab = "Categorized rainfall",
     xaxt = "n",
     yaxt = "n",
     xaxs = "i")

axis(1, at = at_dates, labels = yrs, las = 1)
axis(2, at = 1:K)

segments(df_all$date, 1, df_all$date, df_all$cat, col = "grey65", lwd = 1.5)

lines(df_obs$date, df_obs$Ehat, type = "l", lty = "solid", col = "blue", lwd = 1.5)
lines(df_fut$date, df_fut$Ehat, type = "l", lty = "solid", col = "red4", lwd = 2)
abline(v = min(df_fut$date), lty = 2, lwd = 2, col = "red")

legend(
  "topleft",
  legend = c("Categorized rainfall category",
             "Average fitted categorical response",
             "Rolling one-step-ahead forecast rainfall categories"),
  col = c("grey65", "blue", "red4"),
  lty = c("solid", "solid", "solid"),
  lwd = c(1.5, 1.5),
  bty = "n",
  horiz = FALSE,
  seg.len = 3,
  cex = 2
)

dev.off()
