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


fit1 <- aopts(y = X_hour_obs, X = DesignXEst_obs, method = "clse")
fit1$convergence

paramEst <- fit1$par
paramSE  <- fit1$se

RES_lag1 <- cbind(
  Estimate = paramEst,
  `CI lower` = paramEst - 1.96 * paramSE,
  `CI upper` = paramEst + 1.96 * paramSE,
  `paramSE` = paramSE
)

rownames(RES_lag1) <- c(paste0("c", 2:(K - 1)),
                        colnames(DesignXEst_obs),
                        "rho")

round(RES_lag1, 4)


aop_marginal_probs <- function(par, x_cur, K) {
  
  x_cur <- as.numeric(x_cur)
  q <- length(x_cur)
  
  cuts <- as.numeric(par[1:(K - 2L)])
  beta <- as.numeric(par[(K - 1L):(K - 2L + q)])
  
  seg <- c(-Inf, 0, cuts, Inf)
  mu_cur <- sum(x_cur * beta)
  
  p <- pnorm(seg[-1] - mu_cur) -
    pnorm(seg[-length(seg)] - mu_cur)
  
  p <- pmax(p, 0)
  p / sum(p)
}


rect_bivnorm <- function(lower, upper, mean, rho) {
  as.numeric(
    mvtnorm::pmvnorm(
      lower = lower,
      upper = upper,
      mean  = mean,
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
  mu_cur  <- sum(x_cur  * beta)
  
  denom <- pnorm(seg[y_prev + 1L] - mu_prev) -
    pnorm(seg[y_prev] - mu_prev)
  
  probs <- rep(NA_real_, K)
  
  for (k in seq_len(K)) {
    probs[k] <- rect_bivnorm(
      lower = c(seg[y_prev],      seg[k]),
      upper = c(seg[y_prev + 1L], seg[k + 1L]),
      mean  = c(mu_prev, mu_cur),
      rho   = rho
    ) / denom
  }
  
  probs <- pmax(probs, 0)
  probs / sum(probs)
}


## Diagnostics on training

par     <- paramEst
X_hour  <- X_hour_obs
Xw      <- X_hour_wide_obs
DesignX <- DesignXEst_obs

Ts <- nrow(Xw)
K  <- ncol(Xw)
m  <- K - 1L

Predw_lag1 <- matrix(NA_real_, nrow = Ts, ncol = K)

Predw_lag1[1, ] <- aop_marginal_probs(
  par = par,
  x_cur = DesignX[1, , drop = TRUE],
  K = K
)

for (t in 2:Ts) {
  Predw_lag1[t, ] <- lag1_cond_probs(
    par = par,
    y_prev = X_hour[t - 1],
    x_prev = DesignX[t - 1, , drop = TRUE],
    x_cur  = DesignX[t, , drop = TRUE],
    K = K
  )
}

Predw_diag <- Predw_lag1[-1, , drop = FALSE]
Xw_diag    <- Xw[-1, , drop = FALSE]
Y_diag     <- X_hour[-1]


std_lag1 <- matrix(NA_real_,
                   nrow = nrow(Predw_diag),
                   ncol = m)

for (i in seq_len(nrow(Predw_diag))) {
  
  p_i <- as.numeric(Predw_diag[i, 1:m])
  e_i <- as.numeric(Xw_diag[i, 1:m] - p_i)
  
  V_i <- diag(p_i, nrow = m) - tcrossprod(p_i)
  R_i <- chol(V_i + diag(1e-10, m))
  
  z_i <- backsolve(R_i, e_i, transpose = TRUE)
  std_lag1[i, ] <- z_i
}


pit_diag <- function(pred, y, H = 20,
                     method = c("expected", "randomized"),
                     seed = NULL,
                     main = "PIT histogram",
                     xlab = "u (PIT)",
                     ylab = "Relative frequency",
                     x_breaks = c(0, .25, .5, .75, 1),
                     x_digits = 2,
                     x_cex = 1.7, x_las = 1,
                     left_margin = 6.2,
                     ylab_line = 4) {
  
  method <- match.arg(method)
  pred <- as.matrix(pred)
  Ti <- nrow(pred)
  
  if (is.factor(y)) y <- as.integer(y)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  mar <- op$mar
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
                cex.main = 2.7,
                yaxt = "n",
                ylab = "",
                xaxt = "n",
                xlab = "")
  
  axis(2, at = pretty(c(0, ytop)))
  mtext(ylab, side = 2, line = ylab_line, cex = 1.7)
  
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
  
  mtext(xlab, side = 1, line = 3.2, cex = 1.7)
  
  invisible(list(relfreq = counts,
                 H = H,
                 n = Ti,
                 p0 = p0,
                 se = se,
                 breaks_pos = at))
}



## Residual ACF + PIT diagnostics

setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/applications/HVACFigResidDiag_lag1.eps",
           width = 20, height = 10)

par(mfrow = c(2, 3),
    mar = c(5.5, 6, 3.5, 6),
    cex.main = 2.7,
    cex.lab = 2.7,
    cex.axis = 2.7)

acf(std_lag1[, 1],
    main = paste0("Level ", 1),
    ylab = "Sample ACF")

acf(std_lag1[, 2],
    main = paste0("Level ", 2),
    ylab = "Sample ACF")

acf(std_lag1[, 3],
    main = paste0("Level ", 3),
    ylab = "Sample ACF")

acf(std_lag1[, 4],
    main = paste0("Level ", 4),
    ylab = "Sample ACF")

acf(std_lag1[, 5],
    main = paste0("Level ", 5),
    ylab = "Sample ACF")


pit_diag(Predw_diag, Y_diag,
         H = 20,
         left_margin = 6.5,
         method = "expected",
         main = "PIT histograms",
         xlab = "Probability Integral Transform (U)",
         x_breaks = seq(from = 0, to = 1, by = 0.1))

dev.off()





## Plot observed categories + lag1 scores

H <- length(X_hour_pred)

prob_roll_lag1 <- matrix(NA_real_, nrow = H, ncol = K)

X_hist <- X_hour_obs
Design_hist <- DesignXEst_obs

for (h in seq_len(H)) {
  
  prob_roll_lag1[h, ] <- lag1_cond_probs(
    par = paramEst,
    y_prev = X_hist[length(X_hist)],
    x_prev = Design_hist[nrow(Design_hist), , drop = TRUE],
    x_cur  = DesignXEst_pred[h, , drop = TRUE],
    K = K
  )
  
  X_hist <- c(X_hist, X_hour_pred[h])
  Design_hist <- rbind(
    Design_hist,
    DesignXEst_pred[h, , drop = FALSE]
  )
}

score_fit_lag1 <- as.numeric(Predw_lag1 %*% seq_len(K))
score_roll_lag1 <- as.numeric(prob_roll_lag1 %*% seq_len(K))

x_obs <- seq_along(X_hour_obs)
x_fut <- length(X_hour_obs) + seq_len(H)
x_all <- c(x_obs, x_fut)

y_all <- c(X_hour_obs, X_hour_pred)

setEPS()
postscript(
  "~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/applications/HVACFigFitline_lag1.eps",
  width = 16,
  height = 8
)

par(
  cex.lab = 1.45,
  cex.axis = 1.45,
  mfrow = c(1, 1)
)

plot(
  x = x_all,
  y = y_all,
  type = "l",
  col = "grey25",
  lwd = 1,
  xlab = "Time",
  ylab = "Quality levels",
  xaxt = "n",
  yaxt = "n",
  xlim = c(-5, length(y_all) + 5),
  ylim = c(1, K + 1),
  xaxs = "i"
)

axis(1, at = pretty(x_all), labels = pretty(x_all), las = 1)
axis(2, at = 1:K)

lines(
  x_obs,
  score_fit_lag1,
  type = "l",
  lty = "solid",
  col = "blue",
  lwd = 2.5
)
lines(
  x_fut,
  score_roll_lag1,
  type = "l",
  lty = "solid",
  col = "red4",
  lwd = 2.5
)
abline(
  v = length(X_hour_obs),
  lty = 2,
  lwd = 2.5,
  col = "red"
)

legend(
  "topleft",
  legend = c("Average fitted quality level",
             "Rolling one-step-ahead forecast quality level"),
  col = c("blue", "red4"),
  lty = c("solid", "solid"),
  lwd = 2.5,
  bty = "n",
  horiz = FALSE,
  seg.len = 3,
  cex = 1.45
)

dev.off()