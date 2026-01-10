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
# library(qqplotr)
# library(latex2exp)
library(astsa)

library(TSAOP)

options("scipen"=10) # show all digits


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
DesignXEst2 <- cbind(1, TrendValue, BValue, DValue, EValue, FValue, CpValueY)


FitRes1 <- aopts(y=X_hour, X=DesignXEst1, method="clse")
FitRes2 <- aopts(y=X_hour, X=DesignXEst2, method="clse")

FitRes1$convergence
FitRes2$convergence


ParEst1 <- FitRes1$par
ParEst2 <- FitRes2$par

ParEst1.se <- FitRes1$se
ParEst2.se <- FitRes2$se

round(ParEst1.se, 4)
round(ParEst2.se, 4)

round(cbind(ParEst1, ParEst1-1.96*ParEst1.se, ParEst1+1.96*ParEst1.se), 4)
round(cbind(ParEst2, ParEst2-1.96*ParEst2.se, ParEst2+1.96*ParEst2.se), 4)





##### Univariate Expectation

seg_est1 <- c(-Inf, 0, ParEst1[1:(K-2)], Inf)
# Univariate Expectation
if(length(ParEst1)==K){
  mst1 <- rep(0, Ts)
}else{
  mst1 <- DesignXEst1%*%ParEst1[(K-1):(length(ParEst1)-1)]
}

EXL1 <- rep(0, Ts)
for(t in 1:Ts){EXL1[t] <- K - sum(pnorm(seg_est1-mst1[t])[2:K])}


seg_est2 <- c(-Inf, 0, ParEst2[1:(K-2)], Inf)
# Univariate Expectation
if(length(ParEst2)==K){
  mst2 <- rep(0, Ts)
}else{
  mst2 <- DesignXEst2%*%ParEst2[(K-1):(length(ParEst2)-1)]
}

EXL2 <- rep(0, Ts)
for(t in 1:Ts){EXL2[t] <- K - sum(pnorm(seg_est2-mst2[t])[2:K])}

##################################################
#-------------- Application Figure 1  -----------# Time Series with Expected value as mean structure
##################################################

setEPS()
postscript("grad_app/application/abqrainfall/ApplyFigure1.eps", width = 10, height = 6)

par(cex.lab=1.2, cex.axis=1.2, mfrow=c(1,1))
plot(x=1:Ts, y=X_hour, type="l", col="darkgrey",
     xlab="Year", ylab="Categorized Rainfall",
     xaxt="n", yaxt="n", 
     xlim=c(-100, Ts+100), ylim=c(1,K+0.2), xaxs="i")
tmp <- 1990:2001
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(as.numeric(format(mydat$DATE, "%Y")) == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1)
axis(2, at = c(1:K))

lines(1:Ts, EXL1, type="l", lty="solid", col="orange", lwd=2)
lines(1:(CpLocY-1), EXL2[1:(CpLocY-1)], type="l", lty="solid", col="blue", lwd=2)
lines(CpLocY:Ts, EXL2[CpLocY:Ts], type="l", lty="dashed", col="red", lwd=2)
legend(
  x = 500, y=3.25, 
  legend = c("Without changepoint","Before changepoint","After changepoint"),
  col = c("orange","blue","red"),
  lty = c("solid","solid","dashed"),
  lwd = 2,
  bty = "n",
  horiz = TRUE,
  # x.intersp = 0.3,   # tighter gap between line and text (default ~1)
  seg.len   = 2.5,   # shorter legend lines
  # inset     = 0.01
)
dev.off()






















# ---------- helpers ----------

# PSD-safe inverse (uses Cholesky; tiny jitter fallback)
solve_psd <- function(M) {
  M <- (M + t(M)) / 2
  p <- nrow(M)
  tryCatch(chol2inv(chol(M)),
           error = function(e) solve(M + 1e-10 * diag(p)))
}

# Project each row onto the probability simplex {x>=0, sum x = 1}
.project_simplex_row <- function(v) {
  u <- sort(v, decreasing = TRUE)
  sv <- cumsum(u)
  rho <- max(which(u > (sv - 1) / seq_along(u)))
  theta <- (sv[rho] - 1) / rho
  pmax(v - theta, 0)
}
project_rows_to_simplex <- function(M) t(apply(M, 1L, .project_simplex_row))

# ---------- multikappa: cov(w_i, w_j) block for K-1 working cats ----------

multikappa <- function(i, j, msti, mstj, EXwi, EXwj, seg, rho) {
  # seg: cutpoints (-Inf, c1, c2, ..., c_{K-1}, Inf)
  # EXwi/EXwj: length K-1 vectors of category probs at times i,j (excluding last cat)
  K  <- length(seg) - 1L
  out <- matrix(NA_real_, nrow = K-1L, ncol = K-1L)
  
  if (i == j) {
    # same time: Cov(1{Y=k},1{Y=k'}) = diag(p) - p p^T, take (K-1)x(K-1) block
    out <- - EXwi %*% t(EXwi)
    diag(out) <- diag(out) + EXwi
    return(out)
  }
  
  # different times: latent bivar normal with corr = rho^|i-j|
  lag  <- abs(i - j)
  r    <- rho^lag
  Sigma <- matrix(c(1, r, r, 1), 2, 2)   # correlation matrix
  meanvec <- c(msti, mstj)
  
  for (k in 1:(K-1L)) for (kp in 1:(K-1L)) {
    lower <- c(seg[k],    seg[kp])
    upper <- c(seg[k+1L], seg[kp+1L])
    pij <- mvtnorm::pmvnorm(lower = lower, upper = upper,
                            mean = meanvec, corr = Sigma,
                            algorithm = mvtnorm::GenzBretz(abseps = 1e-8, maxpts = 1e5))
    out[k, kp] <- as.numeric(pij)
  }
  out - EXwi %*% t(EXwj)
}

# ---------- Innovations (truncated) ----------

multiInnov <- function(mst, EXW, seg, rho, num_coef = NULL, mytol = 1e-6, mylag = 0) {
  # mst: mean of latent Z_t (length Ts)
  # EXW: Ts x (K-1) probs (excluding last category)
  Ts <- length(mst)
  K  <- length(seg) - 1L
  if (is.null(num_coef)) num_coef <- Ts
  
  V  <- vector(mode = "list", length = Ts)
  HQ <- array(NA_real_, c(Ts, num_coef, K-1L, K-1L))  # Θ coefficients
  
  # t=1 (index 1)
  V[[1]] <- multikappa(1, 1, mst[1], mst[1], EXW[1, ], EXW[1, ], seg, rho)
  
  # t=2 (index 2)
  HQ[1, 1, , ] <- multikappa(2, 1, mst[2], mst[1], EXW[2, ], EXW[1, ], seg, rho) %*% solve_psd(V[[1]])
  tmp <- HQ[1, 1, , ] %*% V[[1]] %*% t(HQ[1, 1, , ])
  V[[2]] <- multikappa(2, 2, mst[2], mst[2], EXW[2, ], EXW[2, ], seg, rho) - tmp
  
  # t >= 3
  tmpcheck <- matrix(100, nrow = K-1L, ncol = K-1L)
  for (t in 2:(Ts-1)) {
    if (any(tmpcheck > mytol) && mylag == 0) {
      # k = 0
      HQ[t, t, , ] <- multikappa(t+1, 1, mst[t+1], mst[1], EXW[t+1, ], EXW[1, ], seg, rho) %*% solve_psd(V[[1]])
      # k > 0
      for (k in 1:(t-1)) {
        tmp <- matrix(0, K-1L, K-1L)
        for (j in 0:(k-1)) {
          tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[k, k-j, , ])
        }
        HQ[t, t-k, , ] <- (multikappa(t+1, k+1, mst[t+1], mst[k+1], EXW[t+1, ], EXW[k+1, ], seg, rho) - tmp) %*% solve_psd(V[[k+1]])
      }
      # Prediction MSE
      tmp <- matrix(0, K-1L, K-1L)
      for (j in 0:(t-1)) {
        tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[t, t-j, , ])
      }
      V[[t+1]] <- multikappa(t+1, t+1, mst[t+1], mst[t+1], EXW[t+1, ], EXW[t+1, ], seg, rho) - tmp
      
      tmpcheck <- abs(HQ[t-1, t-1, , ])
      if (all(abs(HQ[t, t, , ]) < mytol)) {
        mylag <- t
        HQ <- HQ[, 1:mylag, , , drop = FALSE]
      }
    } else {
      if (t <= mylag) {
        # k = 0
        HQ[t, t, , ] <- multikappa(t+1, 1, mst[t+1], mst[1], EXW[t+1, ], EXW[1, ], seg, rho) %*% solve_psd(V[[1]])
        # k > 0
        for (k in 1:(t-1)) {
          tmp <- matrix(0, K-1L, K-1L)
          for (j in 0:(k-1)) {
            tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[k, k-j, , ])
          }
          HQ[t, t-k, , ] <- (multikappa(t+1, k+1, mst[t+1], mst[k+1], EXW[t+1, ], EXW[k+1, ], seg, rho) - tmp) %*% solve_psd(V[[k+1]])
        }
        tmp <- matrix(0, K-1L, K-1L)
        for (j in 0:(t-1)) {
          tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[t, t-j, , ])
        }
        V[[t+1]] <- multikappa(t+1, t+1, mst[t+1], mst[t+1], EXW[t+1, ], EXW[t+1, ], seg, rho) - tmp
        tmpcheck <- abs(HQ[t-1, t-1, , ])
      } else {
        # only effective lags
        for (k in (t-mylag):(t-1)) {
          tmp <- matrix(0, K-1L, K-1L)
          for (j in max(k-mylag, 0):(k-1)) {
            if ((t-j) < mylag) {
              tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[k, k-j, , ])
            }
          }
          HQ[t, t-k, , ] <- (multikappa(t+1, k+1, mst[t+1], mst[k+1], EXW[t+1, ], EXW[k+1, ], seg, rho) - tmp) %*% solve_psd(V[[k+1]])
        }
        tmp <- matrix(0, K-1L, K-1L)
        for (j in (t-mylag):(t-1)) {
          if ((t-j) < mylag) {
            tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[t, t-j, , ])
          }
        }
        V[[t+1]] <- multikappa(t+1, t+1, mst[t+1], mst[t+1], EXW[t+1, ], EXW[t+1, ], seg, rho) - tmp
        tmpcheck <- abs(HQ[t-1, mylag, , ])
      }
    }
  }
  list(HQ = HQ, V = V, mylag = mylag)
}

# ---------- Innovations (exact, all lags) ----------

multiInnovALL <- function(mst, EXW, seg, rho) {
  Ts <- length(mst)
  K  <- length(seg) - 1L
  
  V  <- vector(mode = "list", length = Ts)
  HQ <- array(NA_real_, c(Ts, Ts, K-1L, K-1L))
  
  # t=1
  V[[1]] <- multikappa(1, 1, mst[1], mst[1], EXW[1, ], EXW[1, ], seg, rho)
  
  # t=2
  HQ[1, 1, , ] <- multikappa(2, 1, mst[2], mst[1], EXW[2, ], EXW[1, ], seg, rho) %*% solve_psd(V[[1]])
  tmp <- HQ[1, 1, , ] %*% V[[1]] %*% t(HQ[1, 1, , ])
  V[[2]] <- multikappa(2, 2, mst[2], mst[2], EXW[2, ], EXW[2, ], seg, rho) - tmp
  
  # t>=3
  for (t in 2:(Ts-1)) {
    HQ[t, t, , ] <- multikappa(t+1, 1, mst[t+1], mst[1], EXW[t+1, ], EXW[1, ], seg, rho) %*% solve_psd(V[[1]])
    for (k in 1:(t-1)) {
      tmp <- matrix(0, K-1L, K-1L)
      for (j in 0:(k-1)) {
        tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[k, k-j, , ])
      }
      HQ[t, t-k, , ] <- (multikappa(t+1, k+1, mst[t+1], mst[k+1], EXW[t+1, ], EXW[k+1, ], seg, rho) - tmp) %*% solve_psd(V[[k+1]])
    }
    tmp <- matrix(0, K-1L, K-1L)
    for (j in 0:(t-1)) {
      tmp <- tmp + HQ[t, t-j, , ] %*% V[[j+1]] %*% t(HQ[t, t-j, , ])
    }
    V[[t+1]] <- multikappa(t+1, t+1, mst[t+1], mst[t+1], EXW[t+1, ], EXW[t+1, ], seg, rho) - tmp
  }
  list(HQ = HQ, V = V, mylag = 0L)
}

# ---------- One-step prediction & innovations (K-1 working cats) ----------

MultiOneStepPred <- function(EXW, XwKminus1, mylag, HQ) {
  Ts  <- nrow(EXW)
  Err <- XwKminus1 - EXW   # error in probs (K-1 columns)
  
  PredErr <- matrix(0, Ts, ncol(Err))
  for (t in 1:(Ts-1)) {
    predlag <- if (mylag > 0) min(mylag, t) else t
    for (j in 1:predlag) {
      PredErr[t+1, ] <- PredErr[t+1, ] + (HQ[t, j, , ] %*% (Err[t+1-j, ] - PredErr[t+1-j, ]))
    }
  }
  
  ErrErr <- Err - PredErr            # innovations (K-1)
  PredX  <- PredErr + EXW            # one-step-ahead predicted probs (K-1)
  list(ErrErr = ErrErr, PredX = PredX)
}

# ---------- Top-level: run innovations and return predictions & residuals ----------

ResWInnov <- function(par, Xl, Xw, DesignX, UseAllLag = FALSE, num_coef = NULL, mylag = 0) {
  Ts <- length(Xl)
  K  <- nlevels(factor(Xl))
  num_par <- length(par)
  
  # thresholds
  ci  <- c(0, par[1:(K-2)])
  seg <- c(-Inf, ci, Inf)
  
  # regression + ar
  theta <- par[(K-1):(num_par-1)]
  rho   <- par[num_par]
  
  # latent means
  mst <- as.vector(DesignX %*% theta)
  
  # conditional category probs (K-1 working cats)
  tmp <- t(matrix(rep(seg, Ts), nrow = length(seg), ncol = Ts)) -
    matrix(rep(mst, length(seg)), nrow = Ts)
  EXW <- pnorm(tmp[, 2:(K+1)]) - pnorm(tmp[, 1:K])
  EXW <- EXW[, 1:(K-1), drop = FALSE]
  
  # innovations recursion
  if (UseAllLag) {
    inv <- multiInnovALL(mst = mst, EXW = EXW, seg = seg, rho = rho)
  } else {
    inv <- multiInnov(mst = mst, EXW = EXW, seg = seg, rho = rho,
                      num_coef = num_coef, mytol = 1e-6, mylag = mylag)
  }
  HQ <- inv$HQ; mylag <- inv$mylag; V <- inv$V
  
  # one-step prediction for probs of first K-1 cats
  predK1 <- MultiOneStepPred(EXW, Xw[, 1:(K-1), drop = FALSE], mylag, HQ)
  
  # append last category, then project to valid simplex
  Predw <- cbind(predK1$PredX, 1 - rowSums(predK1$PredX))
  Predw <- project_rows_to_simplex(Predw)   # keep everything in [0,1] & rowsum=1
  
  list(
    Predw      = Predw,             # Ts x K, one-step-ahead predicted probs
    residw     = Xw - Predw,        # Ts x K, residuals (indicators - predicted probs)
    innovation = predK1$ErrErr,     # Ts x (K-1), innovations (working block)
    HQ         = HQ,
    V          = V,
    mylag      = mylag
  )
}

RR <- ResWInnov(par = ParEst2, Xl = X_hour, Xw = X_hour_wide, DesignX = DesignXEst2)
dim(RR$Predw)       # Ts x K
dim(RR$residw)      # Ts x K
dim(RR$innovation)  # Ts x (K-1)

range(RR$Predw)           # now in [0,1]
summary(rowSums(RR$Predw))# exactly 1 for each row
colMeans(RR$innovation)   # ≈ 0 (ignore first few transients)


setEPS()
postscript("grad_app/application/abqrainfall/ApplyFigure2.eps", width = 17, height = 10)


par(mfrow=c(2,2), 
    mar = c(5.5, 6, 3.5, 6),
    cex.main=1.7, cex.lab=1.7, cex.axis=1.7)
for(k in 1:K){
  acf(RR$residw[,k], 
      main = paste0("One-step-ahead probability residuals for category ", k),
      ylab = 'Sample ACF')
}

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
  mtext(xlab, side = 1, line = 2.6, cex = 1.7)
  
  invisible(list(relfreq = counts, H = H, n = Ti, p0 = p0, se = se,
                 breaks_pos = at))
}



pit_diag(RR$Predw, X_hour, H = 50, left_margin = 6.5,
         method = "expected",
         main = "PIT histograms",
         xlab = "Probability Integral Transform (U)", 
         x_breaks = seq(from=0, to=1, by=0.1))

dev.off()
