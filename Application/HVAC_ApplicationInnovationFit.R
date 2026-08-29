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

save.image("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/applications/HVAC_InnovFit_nocpt_forecast.RData")





load("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/applications/HVAC_InnovFit_nocpt_forecast.RData")

ParEst <- fit1$Estimate

## Diagnostics on training period

par     <- fit1$Estimate
X_hour  <- X_hour_obs
Xw      <- X_hour_wide_obs
DesignX <- DesignXEst_obs

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
Predw <- cbind(PredX, 1 - rowSums(PredX))

resid  <- Xw - Predw
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


## Residual ACF + PIT diagnostics

setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/applications/HVACFigResidDiag_innov.eps", 
           width = 20, height = 10)

par(mfrow = c(2, 3),
    mar = c(5.5, 6, 3.5, 6),
    cex.main = 2.7, cex.lab = 2.7, cex.axis = 2.7)

acf(std_innov[, 1],
    main = paste0("Level ", 1),
    ylab = "Sample ACF")

acf(std_innov[, 2],
    main = paste0("Level ", 2),
    ylab = "Sample ACF")

acf(std_innov[, 3],
    main = paste0("Level ", 3),
    ylab = "Sample ACF")

acf(std_innov[, 4],
    main = paste0("Level ", 4),
    ylab = "Sample ACF")

acf(std_innov[, 5],
    main = paste0("Level ", 5),
    ylab = "Sample ACF")


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
    counts <- as.numeric(table(cut(U, breaks = brks, include.lowest = TRUE))) / Ti
  }
  
  p0 <- 1 / H
  se <- sqrt(p0 * (1 - p0) / Ti)
  ytop <- max(c(counts, p0 + 3 * se))
  
  bp <- barplot(counts, col = "lightblue", border = "grey30",
                ylim = c(0, ytop), main = main, cex.main = 2.7,
                yaxt = "n", ylab = "", xaxt = "n", xlab = "")
  
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
       cex.axis = x_cex, las = x_las, tck = -0.015)
  
  mtext(xlab, side = 1, line = 3.2, cex = 1.7)
  
  invisible(list(relfreq = counts, H = H, n = Ti, p0 = p0, se = se,
                 breaks_pos = at))
}


pit_diag(Predw, X_hour_obs, H = 20, left_margin = 6.5,
         method = "expected",
         main = "PIT histograms",
         xlab = "Probability Integral Transform (U)",
         x_breaks = seq(from = 0, to = 1, by = 0.1))

dev.off()


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


## Forecast evaluation

y <- as.numeric(X_hour_pred)
eps <- 1e-12

Ymat <- matrix(0, nrow = H, ncol = K)
Ymat[cbind(seq_len(H), y)] <- 1

Pcdf <- t(apply(prob_roll, 1, cumsum))[, 1:(K - 1), drop = FALSE]
Ycdf <- sapply(1:(K - 1), function(k) as.numeric(y <= k))

forecast_eval_roll <- data.frame(
  LogScore = -mean(log(pmax(prob_roll[cbind(seq_len(H), y)], eps))),
  BrierScore = mean(rowSums((prob_roll - Ymat)^2)),
  RPS = mean(rowSums((Pcdf - Ycdf)^2))
)

round(forecast_eval_roll, 4)



## Plot observed categories + forecast score

x_obs <- seq_along(X_hour_obs)
x_fut <- length(X_hour_obs) + seq_len(H)
x_all <- c(x_obs, x_fut)

y_all <- c(X_hour_obs, X_hour_pred)

score_fit <- as.numeric(Predw %*% seq_len(K))



setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/applications/HVACFigFitline.eps", 
           width = 14, height = 7)

par(cex.lab = 1.45,
    cex.axis = 1.45,
    mfrow = c(1, 1))

plot(x = x_all, y = y_all,
     type = "l",
     col = "grey65",
     lwd = 1,
     xlab = "Time",
     ylab = "Quality levels",
     xaxt = "n",
     yaxt = "n",
     xlim = c(-5, length(y_all) + 5),
     ylim = c(1, K + 1),
     xaxs = "i")

axis(1, at = pretty(x_all), labels = pretty(x_all), las = 1)
axis(2, at = 1:K)

lines(x_obs, score_fit,
      type = "l",
      lty = "solid",
      col = "blue",
      lwd = 2.5)
lines(x_fut, score_roll,
      type = "l",
      lty = "solid",
      col = "red4",
      lwd = 2.5)
abline(v = length(X_hour_obs), lty = 2, lwd = 2.5, col = "red")

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
  cex=1.45
)

dev.off()