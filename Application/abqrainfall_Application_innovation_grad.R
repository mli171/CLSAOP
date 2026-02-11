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



X_hour = mydat$Catg

Ts = length(X_hour)
K  = length(unique(X_hour))

X_hour_wide = matrix(0, nrow=length(X_hour), ncol=K)
for (t in 1:Ts) {X_hour_wide[t,X_hour[t]] <- 1}


ss1 <- 365
ss2 <- 183
TrendValue <- 1:Ts/Ts
BValue <- cos(2*pi*(1:Ts)/ss1)
DValue <- sin(2*pi*(1:Ts)/ss1)
EValue <- cos(2*pi*(1:Ts)/ss2)
FValue <- sin(2*pi*(1:Ts)/ss2)
CpLocY <- 2366
CpValueY <- c(rep(0, CpLocY-1), rep(1, Ts-CpLocY+1))


DesignXEst1 <- cbind(1, TrendValue, BValue, DValue, EValue, FValue)
colnames(DesignXEst1) <- c("Intercept", "alpha", "B", "D", "E", "F")
DesignXEst2 <- cbind(1, TrendValue, BValue, DValue, EValue, FValue, CpValueY)
colnames(DesignXEst2) <- c("Intercept", "alpha", "B", "D", "E", "F", "Delta")

fit1 <- fit_CLSInnov(X_hour=X_hour, X_hour_wide=X_hour_wide, DesignXEst=DesignXEst1)
round(fit1, 4)
#           Estimate Std. Error  z.value p.value   ci.lb   ci.ub
# c2          0.8137     0.0439  18.5453  0.0000  0.7277  0.8997
# Intercept  -0.8204     0.0650 -12.6268  0.0000 -0.9478 -0.6931
# alpha      -0.2819     0.1194  -2.3611  0.0182 -0.5159 -0.0479
# B          -0.1796     0.0463  -3.8816  0.0001 -0.2702 -0.0889
# D          -0.1607     0.0469  -3.4287  0.0006 -0.2525 -0.0688
# E           0.0769     0.0469   1.6392  0.1012 -0.0151  0.1689
# F           0.1632     0.0468   3.4855  0.0005  0.0714  0.2549
# rho         0.4117     0.0283  14.5551  0.0000  0.3562  0.4671

fit2 <- fit_CLSInnov(X_hour=X_hour, X_hour_wide=X_hour_wide, DesignXEst=DesignXEst2)
round(fit2, 4)
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








par     <- fit2$Estimate
Xl      <- X_hour
Xw      <- X_hour_wide
DesignX <- DesignXEst2


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
EXW  <- Pmat[, 1:m, drop = FALSE]

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

# Standardized one-step-ahead prediction residuals

# dimension of innovations
p_dim <- m
std_innov <- matrix(NA_real_, Ts, p_dim)

for (t in seq_len(Ts)) {
  e_t <- as.numeric(ErrErr[t, 1:p_dim])
  V_t <- matrix(V_list[[t]], nrow = p_dim, ncol = p_dim)
  R_t <- chol(V_t)  # upper triangular, V_t = t(R_t) R_t
  ## Solve t(R_t) z_t = e_t  => z_t = (t(R_t))^{-1} e_t
  z_t <- backsolve(R_t, e_t, transpose = FALSE)
  std_innov[t, ] <- z_t
}

setEPS()
postscript("~/Desktop/Changepoint/aop/AMOCAOP/CLS/grad_app/application/abqrainfall/ApplyFigure3_new.eps", width = 20, height = 5)


par(mfrow=c(1,3),
    mar = c(5.5, 6, 3.5, 6),
    cex.main=1.7, cex.lab=1.7, cex.axis=1.7)

acf(std_innov[,1],
    main = paste0("One-step-ahead probability residuals for category ", 1),
    ylab = 'Sample ACF')
acf(std_innov[,2],
    main = paste0("One-step-ahead probability residuals for category ", 2),
    ylab = 'Sample ACF')


pit_diag <- function(pred, y, H = 20,
                     method = c("expected","randomized"),
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
  pred <- as.matrix(pred); Ti <- nrow(pred)
  if (is.factor(y)) y <- as.integer(y)

  # widen left margin for y label
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  mar <- op$mar; mar[2] <- max(mar[2], left_margin)
  par(mar = mar, xaxs = "i")

  # PIT components
  Pdist <- t(apply(cbind(0, pred), 1L, cumsum))
  idx <- seq_len(Ti)
  lower <- Pdist[cbind(idx, y)]
  upper <- Pdist[cbind(idx, y + 1L)]

  if (method == "expected") {
    u_grid <- (1:H) / H
    Ftubar <- vapply(u_grid, function(u) {
      below  <- (u <  lower)
      above  <- (u >= upper)
      inside <- !(below | above)
      w <- pmax(upper - lower, .Machine$double.eps)
      mean(0 + 1*above + inside * (u - lower)/w)
    }, numeric(1))
    counts <- diff(c(0, Ftubar))
  } else {
    if (!is.null(seed)) set.seed(seed)
    U <- lower + runif(Ti) * pmax(upper - lower, .Machine$double.eps)
    brks <- seq(0, 1, length.out = H + 1)
    counts <- as.numeric(table(cut(U, breaks = brks, include.lowest = TRUE))) / Ti
  }

  p0 <- 1/H; se <- sqrt(p0*(1-p0)/Ti)
  ytop <- max(c(counts, p0 + 3*se))

  # draw bars (suppress axes/labels; we add custom ones)
  bp <- barplot(counts, col = "lightblue", border = "grey30",
                ylim = c(0, ytop), main = main, cex.main=1.7,
                yaxt = "n", ylab = "", xaxt = "n", xlab = "")

  # y axis + label
  axis(2, at = pretty(c(0, ytop)))
  mtext(ylab, side = 2, line = ylab_line, cex=1.7)

  # reference lines
  abline(h = p0, col = "red", lwd = 2, lty = 2)
  abline(h = p0 + 1.96*se, col = "blue", lwd = 2, lty = 2)
  abline(h = p0 - 1.96*se, col = "blue", lwd = 2, lty = 2)

  # map 0–1 to the barplot coordinate system for proper tick placement
  w <- if (length(bp) >= 2) (bp[2] - bp[1]) else 1
  edge0 <- bp[1] - 0.5 * w
  edge1 <- bp[length(bp)] + 0.5 * w
  at <- edge0 + x_breaks * (edge1 - edge0)

  axis(1, at = at,
       labels = format(x_breaks, nsmall = x_digits),
       cex.axis = x_cex, las = x_las, tck = -0.015)
  mtext(xlab, side = 1, line = 3.2, cex = 1)

  invisible(list(relfreq = counts, H = H, n = Ti, p0 = p0, se = se,
                 breaks_pos = at))
}


pit_diag(Predw, X_hour, H = 50, left_margin = 6.5,
         method = "expected",
         main = "PIT histograms",
         xlab = "Probability Integral Transform (U)",
         x_breaks = seq(from=0, to=1, by=0.1))

dev.off()
