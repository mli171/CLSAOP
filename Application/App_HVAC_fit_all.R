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

X_hour <- data.full[1:(tauT - 1)]

K <- max(X_hour)

make_wide <- function(x, K) {
  Xw <- matrix(0, nrow = length(x), ncol = K)
  for (t in seq_along(x)) {Xw[t, x[t]] <- 1}
  Xw
}

X_hour_wide  <- make_wide(X_hour, K)
DesignXEst <- matrix(1, nrow = length(X_hour), ncol = 1)
colnames(DesignXEst) <- "Intercept"




##### Innovation-based CLS #####

fit_innov <- fit_CLSInnov(X_hour=X_hour, X_hour_wide=X_hour_wide, DesignXEst=DesignXEst)
paramEst_innov <- fit_innov$Estimate
RES_innov <- fit_innov
RES_innov
#            Estimate Std. Error   z.value      p.value      ci.lb     ci.ub
# c2        0.5738583 0.07648479  7.502907 6.241762e-14 0.42394812 0.7237685
# c3        0.9908480 0.09077948 10.914889 9.784881e-28 0.81292019 1.1687758
# c4        1.3218429 0.09994602 13.225568 6.246346e-40 1.12594873 1.5177371
# c5        1.7761323 0.11634346 15.266284 1.282770e-52 1.54809911 2.0041655
# Intercept 0.9569722 0.11028407  8.677338 4.051381e-18 0.74081540 1.1731289
# rho       0.2219107 0.09331081  2.378188 1.739794e-02 0.03902148 0.4047999

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
#            Estimate   CI lower  CI upper    paramSE
# c2        0.5738623 0.42291286 0.7248116 0.07701500
# c3        0.9908601 0.81181152 1.1699087 0.09135131
# c4        1.3218865 1.12518303 1.5185901 0.10035894
# c5        1.7757495 1.54950696 2.0019920 0.11542986
# Intercept 0.9574057 0.74555808 1.1692534 0.10808553
# rho       0.2232137 0.05776497 0.3886624 0.08441261

## Diagnostics


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


par     <- paramEst_lag1
X_hour  <- X_hour
Xw      <- X_hour_wide
DesignX <- DesignXEst

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







setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/HVACFigResid_ACF.eps",
           width = 20, height = 5)


par(mfrow = c(2, 5),
    mar = c(4.2, 3.4, 2.8, 0.6),
    mgp = c(2.1, 0.7, 0),
    cex.main = 1.25,
    cex.lab = 1.25,
    cex.axis = 1.15)

acf(std_innov_diag[, 1], main = "Category 1 (innovations-based CLS)", ylab = "Sample ACF", lag.max = 15)
acf(std_innov_diag[, 2], main = "Category 2 (innovations-based CLS)", ylab = "Sample ACF", lag.max = 15)
acf(std_innov_diag[, 3], main = "Category 3 (innovations-based CLS)", ylab = "Sample ACF", lag.max = 15)
acf(std_innov_diag[, 4], main = "Category 4 (innovations-based CLS)", ylab = "Sample ACF", lag.max = 15)
acf(std_innov_diag[, 5], main = "Category 5 (innovations-based CLS)", ylab = "Sample ACF", lag.max = 15)

acf(std_lag1[, 1], main = "Category 1 (CLS-lag-1)", ylab = "Sample ACF", lag.max = 15)
acf(std_lag1[, 2], main = "Category 2 (CLS-lag-1)", ylab = "Sample ACF", lag.max = 15)
acf(std_lag1[, 3], main = "Category 3 (CLS-lag-1)", ylab = "Sample ACF", lag.max = 15)
acf(std_lag1[, 4], main = "Category 4 (CLS-lag-1)", ylab = "Sample ACF", lag.max = 15)
acf(std_lag1[, 5], main = "Category 5 (CLS-lag-1)", ylab = "Sample ACF", lag.max = 15)

dev.off()





setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/HVACFigResid_PIT.eps",
           width = 25, height = 10)

par(mfrow = c(1, 2), mar = c(5.5, 6, 3.5, 6),
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

pit_diag(Predw_innov_diag, Y_innov_diag,
         H = 20,
         left_margin = 6.5,
         method = "expected",
         main = "Innovations-based CLS",
         xlab = "Probability Integral Transform (U)",
         x_breaks = seq(from = 0, to = 1, by = 0.1))

pit_diag(Predw_diag_lag1, Y_diag,
         H = 20, left_margin = 6.5,
         method = "expected",
         main = "CLS-lag-1",
         xlab = "Probability Integral Transform (U)",
         x_breaks = seq(from = 0, to = 1, by = 0.1))

dev.off()


save.image("~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/application_final/HVAC_FullFit.RData")

