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

pit_diag <- function(pred, y, H = 20,
                     method = c("expected", "randomized"),
                     seed = NULL,
                     main = "PIT histogram",
                     xlab = "u (PIT)",
                     ylab = "Relative frequency",
                     x_breaks = c(0, .25, .5, .75, 1),
                     x_digits = 2,
                     x_cex = 1.5, x_las = 1,
                     left_margin = 6.2,
                     ylab_line = 4) {
  
  method <- match.arg(method)
  pred <- as.matrix(pred)
  Ti <- nrow(pred)
  
  if (is.factor(y)) y <- as.integer(y)
  
  mar <- par("mar")
  mar[2] <- max(mar[2], left_margin)
  par(mar = mar, xaxs = "i")
  
  Pdist <- t(apply(cbind(0, pred), 1L, cumsum))
  idx <- seq_len(Ti)
  
  lower <- Pdist[cbind(idx, y)]
  upper <- Pdist[cbind(idx, y + 1L)]
  
  if (method == "expected") {
    
    u_grid <- (1:H) / H
    
    Ftubar <- vapply(u_grid, function(u) {
      below  <- (u < lower)
      above  <- (u >= upper)
      inside <- !(below | above)
      w <- pmax(upper - lower, .Machine$double.eps)
      mean(0 + 1 * above + inside * (u - lower) / w)
    }, numeric(1))
    
    counts <- diff(c(0, Ftubar))
    
  } else {
    
    if (!is.null(seed)) set.seed(seed)
    U <- lower + runif(Ti) * pmax(upper - lower, .Machine$double.eps)
    brks <- seq(0, 1, length.out = H + 1)
    counts <- as.numeric(table(cut(U, breaks = brks,
                                   include.lowest = TRUE))) / Ti
  }
  
  p0 <- 1 / H
  se <- sqrt(p0 * (1 - p0) / Ti)
  ytop <- max(c(counts, p0 + 3 * se))
  
  bp <- barplot(counts,
                col = "lightblue",
                border = "grey30",
                ylim = c(0, ytop),
                main = main,
                cex.main = 1.5,
                yaxt = "n",
                ylab = "",
                xaxt = "n",
                xlab = "")
  
  axis(2, at = pretty(c(0, ytop)))
  mtext(ylab, side = 2, line = ylab_line, cex = 1.5)
  
  abline(h = p0, col = "red", lwd = 2, lty = 2)
  abline(h = p0 + 1.96 * se, col = "blue", lwd = 2, lty = 2)
  abline(h = p0 - 1.96 * se, col = "blue", lwd = 2, lty = 2)
  
  w <- if (length(bp) >= 2) (bp[2] - bp[1]) else 1
  edge0 <- bp[1] - 0.5 * w
  edge1 <- bp[length(bp)] + 0.5 * w
  at <- edge0 + x_breaks * (edge1 - edge0)
  
  axis(1, at = at,
       labels = format(x_breaks, nsmall = x_digits),
       cex.axis = x_cex,
       las = x_las,
       tck = -0.015)
  
  mtext(xlab, side = 1, line = 3.2, cex = 1.5)
  
  invisible(list(relfreq = counts,
                 H = H,
                 n = Ti,
                 p0 = p0,
                 se = se,
                 breaks_pos = at))
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


X_hour <- mydat$Catg

Ts <- length(X_hour)
K  <- max(X_hour)

make_wide <- function(x, K) {
  Xw <- matrix(0, nrow = length(x), ncol = K)
  for (t in seq_along(x)) {Xw[t, x[t]] <- 1}
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



##### Innovation-based CLS #####

fit_innov <- fit_CLSInnov(X_hour = X_hour, X_hour_wide = X_hour_wide, DesignXEst = DesignXEst)
paramEst_innov <- fit_innov$Estimate
RES_innov <- fit_innov
RES_innov
#           Estimate Std. Error z.value p.value   ci.lb   ci.ub
# c2          0.8134     0.0436 18.6521  0.0000  0.7280  0.8989
# Intercept  -0.7282     0.0730 -9.9784  0.0000 -0.8712 -0.5851
# alpha      -0.7260     0.2124 -3.4179  0.0006 -1.1423 -0.3097
# B          -0.1804     0.0463 -3.8961  0.0001 -0.2711 -0.0896
# D          -0.1570     0.0469 -3.3516  0.0008 -0.2489 -0.0652
# E           0.0722     0.0469  1.5373  0.1242 -0.0198  0.1642
# F           0.1566     0.0470  3.3321  0.0009  0.0645  0.2487
# Delta       0.3121     0.1268  2.4617  0.0138  0.0636  0.5605
# rho         0.4076     0.0284 14.3675  0.0000  0.3520  0.4632

## Diagnostics

par     <- paramEst_innov
X_hour  <- X_hour
Xw      <- X_hour_wide
DesignX <- DesignXEst

Ts <- nrow(Xw)
K  <- ncol(Xw)
m  <- K - 1L
num_par <- length(par)

cuts  <- as.numeric(par[1:(K - 2L)])
theta <- as.numeric(par[(K - 1L):(num_par - 1L)])
rho   <- as.numeric(par[num_par])

seg <- c(-Inf, 0, cuts, Inf)

mst <- as.numeric(DesignX %*% theta)

Fmat <- pnorm(outer(mst, seg, function(mu, s) s - mu))
Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] - Fmat[, 1:K, drop = FALSE]
EXW <- Pmat[, 1:m, drop = FALSE]
u <- Xw[, 1:m, drop = FALSE] - EXW

kappa <- function(i, j) {
  multikappa_pbivnorm(
    i, j,
    msti = mst[i], mstj = mst[j],
    EXwi = EXW[i, ], EXwj = EXW[j, ],
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
Predw_innov <- cbind(PredX, 1 - rowSums(PredX))

resid  <- Xw - Predw_innov
V_list <- InnovRes$V
ErrErr <- InnovRes$innov
mylag  <- InnovRes$mylag

p_dim <- m
std_innov <- matrix(NA_real_, Ts, p_dim)

for (t in seq_len(Ts)) {
  e_t <- as.numeric(ErrErr[t, 1:p_dim])
  V_t <- matrix(V_list[[t]], nrow = p_dim, ncol = p_dim)
  R_t <- chol(V_t)
  z_t <- backsolve(R_t, e_t, transpose = TRUE)
  std_innov[t, ] <- z_t
}

std_innov_diag <- std_innov[-1, , drop = FALSE]
Predw_innov_diag <- Predw_innov[-1, , drop = FALSE]
Y_innov_diag <- X_hour[-1]


##### CLS-lag-1 fit #####

fit_lag1 <- aopts(y = X_hour, X = DesignXEst, method = "clse")

paramEst_lag1 <- fit_lag1$par
paramSE_lag1  <- fit_lag1$se

RES_lag1 <- cbind(
  Estimate   = paramEst_lag1,
  `CI lower` = paramEst_lag1 - 1.96 * paramSE_lag1,
  `CI upper` = paramEst_lag1 + 1.96 * paramSE_lag1,
  `paramSE`  = paramSE_lag1
)
rownames(RES_lag1) <- c(paste0("c", 2:(K - 1)), colnames(DesignXEst), "rho")
round(RES_lag1, 4)
#           Estimate CI lower CI upper paramSE
# c2          0.8115   0.7252   0.8978  0.0440
# Intercept  -0.7277  -0.8670  -0.5884  0.0711
# alpha      -0.7266  -1.1332  -0.3201  0.2074
# B          -0.1791  -0.2665  -0.0917  0.0446
# D          -0.1569  -0.2461  -0.0677  0.0455
# E           0.0713  -0.0173   0.1600  0.0452
# F           0.1554   0.0662   0.2446  0.0455
# Delta       0.3115   0.0696   0.5534  0.1234
# rho         0.4134   0.3586   0.4683  0.0280

## Diagnostics

aop_marginal_probs <- function(par, x_cur, K) {
  
  x_cur <- as.numeric(x_cur)
  q <- length(x_cur)
  
  cuts <- as.numeric(par[1:(K - 2L)])
  beta <- as.numeric(par[(K - 1L):(K - 2L + q)])
  
  seg <- c(-Inf, 0, cuts, Inf)
  mu_cur <- sum(x_cur * beta)
  
  probs <- pnorm(seg[-1] - mu_cur) -
    pnorm(seg[-length(seg)] - mu_cur)
  
  probs <- pmax(probs, 0)
  probs / sum(probs)
}

rect_bivnorm <- function(lower, upper, mean, rho) {
  
  as.numeric(
    mvtnorm::pmvnorm(
      lower = lower,
      upper = upper,
      mean = mean,
      sigma = matrix(c(1, rho, rho, 1), 2, 2)
    )
  )
}

lag1_cond_probs <- function(par, y_prev, x_prev, x_cur, K) {
  
  x_prev <- as.numeric(x_prev)
  x_cur  <- as.numeric(x_cur)
  q <- length(x_cur)
  
  cuts <- as.numeric(par[1:(K - 2L)])
  beta <- as.numeric(par[(K - 1L):(K - 2L + q)])
  rho  <- as.numeric(par[K - 2L + q + 1L])
  
  seg <- c(-Inf, 0, cuts, Inf)
  
  mu_prev <- sum(x_prev * beta)
  mu_cur  <- sum(x_cur * beta)
  
  denom <- pnorm(seg[y_prev + 1L] - mu_prev) -
    pnorm(seg[y_prev] - mu_prev)
  
  probs <- rep(NA_real_, K)
  
  for (k in seq_len(K)) {
    probs[k] <- rect_bivnorm(
      lower = c(seg[y_prev], seg[k]),
      upper = c(seg[y_prev + 1L], seg[k + 1L]),
      mean = c(mu_prev, mu_cur),
      rho = rho
    ) / denom
  }
  
  probs <- pmax(probs, 0)
  probs / sum(probs)
}


par     <- paramEst_lag1
X_hour  <- X_hour
Xw      <- X_hour_wide
DesignX <- DesignXEst

Ts_obs <- nrow(Xw)
K      <- ncol(Xw)
m      <- K - 1L

Predw_lag1 <- matrix(NA_real_, nrow = Ts, ncol = K)

Predw_lag1[1, ] <- aop_marginal_probs(
  par = par, 
  x_cur = DesignX[1, , drop = TRUE],
  K = K
)

for (t in 2:Ts) {
  Predw_lag1[t, ] <- lag1_cond_probs(
    par    = par,
    y_prev = X_hour[t - 1L],
    x_prev = DesignX[t - 1L, , drop = TRUE],
    x_cur  = DesignX[t, , drop = TRUE],
    K      = K
  )
}

Predw_diag <- Predw_lag1[-1, , drop = FALSE]
Xw_diag    <- Xw[-1, , drop = FALSE]
Y_diag     <- X_hour[-1]


std_lag1 <- matrix(NA_real_, nrow = nrow(Predw_diag), ncol = m)

for (i in seq_len(nrow(Predw_diag))) {
  
  p_i <- as.numeric(Predw_diag[i, 1:m])
  e_i <- as.numeric(Xw_diag[i, 1:m] - p_i)
  
  V_i <- diag(p_i, nrow = m) - tcrossprod(p_i)
  R_i <- chol(V_i + diag(1e-10, m))
  
  z_i <- backsolve(R_i, e_i, transpose = TRUE)
  std_lag1[i, ] <- z_i
}

Predw_diag_lag1 = Predw_diag







# setEPS()
# postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/RainfallFigResid_ACF.eps",
#            width = 20, height = 4)
# 
# par(mfrow = c(1, 4),
#     mar = c(5.5, 6, 3.5, 6),
#     cex.main = 1.7,
#     cex.lab = 1.7,
#     cex.axis = 1.7)
# 
# acf(std_innov_diag[, 1], main = "Category 1 (innovations-based CLS)", ylab = "Sample ACF")
# acf(std_innov_diag[, 2], main = "Category 2 (innovations-based CLS)", ylab = "Sample ACF")
# 
# acf(std_lag1[, 1], main = "Category 1 (CLS-lag-1)", ylab = "Sample ACF")
# acf(std_lag1[, 2], main = "Category 2 (CLS-lag-1)", ylab = "Sample ACF")
# 
# 
# dev.off()
# 
# 
# 
# 
# 
# setEPS()
# postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/RainfallFigResid_PIT.eps",
#            width = 25, height = 10)
# 
# par(mfrow = c(1, 2), mar = c(5.5, 6, 3.5, 6),
#     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
# 
# pit_diag(Predw_innov_diag, Y_innov_diag,
#          H = 20,
#          left_margin = 6.5,
#          method = "expected",
#          main = "Innovations-based CLS",
#          xlab = "Probability Integral Transform (U)",
#          x_breaks = seq(from = 0, to = 1, by = 0.1))
# 
# 
# pit_diag(Predw_diag_lag1, Y_diag,
#          H = 20,
#          left_margin = 6.5,
#          method = "expected",
#          main = "CLS-lag-1",
#          xlab = "Probability Integral Transform (U)",
#          x_breaks = seq(from = 0, to = 1, by = 0.1))
# 
# dev.off()



setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/RainfallFigResid.eps",
           width = 20, height = 8)

par(mfrow = c(2, 3), mar = c(5.5, 6, 3.5, 6),
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

acf(std_innov_diag[, 1], main = "Category 1 (innovations-based CLS)", ylab = "Sample ACF")
acf(std_innov_diag[, 2], main = "Category 2 (innovations-based CLS)", ylab = "Sample ACF")
pit_diag(Predw_innov_diag, Y_innov_diag,
         H = 20,
         left_margin = 6.5,
         method = "expected",
         main = "Innovations-based CLS",
         xlab = "Probability Integral Transform (U)",
         x_breaks = seq(from = 0, to = 1, by = 0.1))

acf(std_lag1[, 1], main = "Category 1 (CLS-lag-1)", ylab = "Sample ACF")
acf(std_lag1[, 2], main = "Category 2 (CLS-lag-1)", ylab = "Sample ACF")
pit_diag(Predw_diag_lag1, Y_diag,
         H = 20,
         left_margin = 6.5,
         method = "expected",
         main = "CLS-lag-1",
         xlab = "Probability Integral Transform (U)",
         x_breaks = seq(from = 0, to = 1, by = 0.1))

dev.off()
