rm(list = ls())

options("scipen"=10) # show all digits

library(Rcpp)
library(RcppArmadillo)
library(VGAM)

library(mvtnorm)
library(foreach)
library(doMC)

Rcpp::sourceCpp("tools/myfuncClipPaper.cpp")
Rcpp::sourceCpp("tools/CLSE_Tools_1_test.cpp")


#####  Parameter Estimation by LSE (Li and lu, 2022)
aopLSE = function(X_hour, DesignXEst, start=NULL){
  
  Ts = length(X_hour)
  K = length(levels(factor(X_hour)))
  
  # adapt design matrix for LSE code
  if(sum(DesignXEst)==1){
    DesignXLSE = matrix(0, nrow=Ts, ncol=1)
  }else{
    DesignXLSE = as.matrix(DesignXEst)
  }
  
  #------------------------ Initial Value Calculation
  
  if(is.null(start)){
    cout = summary(factor(X_hour))
    ciInit = as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
    p = which(!is.finite(ciInit))
    if(length(p)!=0){
      if(p==1){
        ciInit[p] = -6
      }else{
        ciInit[p] =  6
      }
    }
    if(sum(DesignXLSE)==0){
      # stationary
      parInit = ciInit
    }else{
      # non-stationary
      friInit = rep(0, NCOL(DesignXLSE))
      parInit <- c(ciInit, friInit)
    }
  }
  
  # create multivariate form of categorical data
  X_hour_wide = matrix(0, nrow=Ts, ncol=K)
  for (h in 1:Ts) {X_hour_wide[h,X_hour[h]] = 1}
  
  #------------------------ Parameter Estimation
  resNW = NW_cpp(par=parInit, X=X_hour_wide, DesignX=DesignXLSE, stepsize=0.5, conv=1e-05)
  MargEst = resNW$par
  MySampleCov = sum(diag(t(X_hour_wide[1:(Ts-1),])%*%X_hour_wide[2:Ts,]))/Ts
  OptimSecond2 = optim(par=acf(X_hour, plot=F)$acf[2,1,1], fn=diffMOM,
                       method = "Brent", lower=-0.99, upper=0.99,
                       control=list(reltol=1e-05), Marg_est=MargEst,
                       K=K, DesignX=DesignXEst, MySampleCov=MySampleCov)
  phiEst = OptimSecond2$par
  paramEst = c(MargEst, phiEst)
  paramEst = c(paramEst[2:(K-1)] - paramEst[1], - paramEst[1], paramEst[K:length(paramEst)])
  
  paramCvg = resNW$convergence + OptimSecond2$convergence
  
  return(list(paramEst=paramEst, paramCvg=paramCvg))
}

##### multivariate innovation algorithm (Exact)
multiInnovALL = function(mst, EXW, seg, rho){
  
  Ts = length(mst)
  K = length(seg)-1
  
  V = vector(mode="list", length=Ts)
  HQ = array(NA, c(Ts, Ts, K-1, K-1))
  
  # t=0 with index 1
  # t=1 with index 2
  
  #--------------------- t=0 ---------------------#
  V[[1]] = multikappa(i=1, j=1, msti=mst[1], mstj=mst[1], EXwi=EXW[1,], EXwj=EXW[1,], seg, rho)
  
  #--------------------- t=1 ---------------------#
  HQ[1,1,,] = multikappa(i=2, j=1, msti=mst[2], mstj=mst[1], EXwi=EXW[2,], EXwj=EXW[1,], seg, rho) %*% solve(V[[1]]) # k=0
  tmp = HQ[1,1,,] %*% V[[1]] %*% t(HQ[1,1,,])
  V[[2]] = multikappa(i=2, j=2, msti=mst[2], mstj=mst[2], EXwi=EXW[2,], EXwj=EXW[2,], seg, rho) - tmp
  
  #--------------------- t>=2 --------------------#
  mylag = 0
  tmpcheck = matrix(100, nrow=K-1, ncol=K-1)
  for(t in 2:(Ts-1)){
    # k = 0
    HQ[t,t,,] = multikappa(i=t+1, j=1, msti=mst[t+1], mstj=mst[1], EXwi=EXW[t+1,], EXwj=EXW[1,], seg, rho) %*% solve(V[[1]])
    # k > 0
    for(k in 1:(t-1)){
      tmp = matrix(0, nrow=K-1, ncol=K-1)
      for(j in 0:(k-1)){
        tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[k,k-j,,])
      }
      HQ[t,t-k,,] = (multikappa(i=t+1, j=k+1, msti=mst[t+1], mstj=mst[k+1], EXwi=EXW[t+1,], EXwj=EXW[k+1,], seg, rho) - tmp) %*% solve(V[[k+1]])
    }
    # Prediction MSE
    tmp = matrix(0, nrow=K-1, ncol=K-1)
    for(j in 0:(t-1)){
      tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[t,t-j,,])
    }
    V[[t+1]] = multikappa(i=t+1, j=t+1, msti=mst[t+1], mstj=mst[t+1], EXwi=EXW[t+1,], EXwj=EXW[t+1,], seg, rho) - tmp
  }
  
  RES = list(HQ=HQ, V=V, mylag=mylag)
  return(RES)
}

##### Vector form conditional least squares estimation by multivariate innovation algorithm
CLSWInnov = function(par, Xl, Xw, DesignX, UseAllLag=FALSE, num_coef=NULL, mylag=0){
  
  Ts = length(Xl)
  K = nlevels(factor(Xl))
  num_par = length(par)
  
  ci = c(0, par[1:(K-2)])
  seg = c(-Inf, ci, Inf)
  theta = par[(K-1):(num_par-1)]
  rho = par[num_par]
  
  ## mean vector
  mst = DesignX%*%theta
  
  ## Conditional expectation
  tmp = t(matrix(rep(seg, Ts), nrow=length(seg), ncol=Ts)) - matrix(rep(mst, length(seg)), nrow=Ts)
  EXW = pnorm(tmp[,2:(K+1)]) - pnorm(tmp[,1:K])
  EXW = EXW[,1:(K-1)]
  
  # Innovation function starts from here
  if(UseAllLag){
    ## all lags version
    InnovRes = multiInnovALL(mst=mst, EXW=EXW, seg=seg, rho=rho)
    Predw = MultiOneStepPred(EXW, Xw[,1:(K-1)], InnovRes$mylag, InnovRes$HQ)
    Predw = cbind(Predw$PredX, 1-rowSums(Predw$PredX))
  }else{
    ## necessary lags version
    InnovRes = multiInnov_cpp(mst = mst, EXW = EXW, seg = seg, rho = rho, 
                              numcoef=num_coef, mytol=1e-05, mylag=mylag)
    Predw = MultiOneStepPred_cpp(EXW = EXW, Xw = Xw[,1:(K-1)], HQ = InnovRes$HQ, 
                                 numcoef=num_coef, mylag=InnovRes$mylag)
    Predw = cbind(Predw$PredX, 1-rowSums(Predw$PredX))
  }
  
  res = sum((Xw - Predw)^2)
  
  # cat("\n=====================")
  # cat("\n par = ", par)
  # cat("\n res = ", res)
  
  return(res)
}


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

options("scipen"=10)


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


#------------------- Parameter Estimation   ------------------------#
cout <- summary(as.factor(X_hour))
ci_initial <- as.vector(qnorm(cumsum(cout / sum(cout))))[1:(K - 1)]
phi_initial <- stats::pacf(X_hour, plot = FALSE)$acf[1]
par_initial <- c(ci_initial[2:(K-1)] - ci_initial[1], -ci_initial[1],
                 rep(0, NCOL(DesignXEst1) - 1), phi_initial)

if (K == 3) {
  constrLSE = rbind(c(rep(0, length(par_initial)-1),  1),
                    c(rep(0, length(par_initial)-1), -1)) # rho
} else {
  constrLSE <- matrix(0, nrow = K - 2, ncol = length(par_initial))
  constrLSE[1,1] <- 1
  for (ii in 2:(K-2)) constrLSE[ii, c(ii-1, ii)] <- c(-1, 1) # ci monotone
  constrLSE = rbind(constrLSE,
                    c(rep(0, ncol(constrLSE)-1),  1),
                    c(rep(0, ncol(constrLSE)-1), -1)) # rho bounds
}
constrLSE_ci <- c(rep(0, nrow(constrLSE) - 2), -1, -1)
## optimization
OptimWInnov1 = constrOptim(par_initial, f=CLSWInnov,
                          method = "Nelder-Mead", ui=constrLSE,
                          ci=constrLSE_ci, hessian=F,
                          control= list(reltol=1e-05),
                          Xl=X_hour, Xw=X_hour_wide, DesignX=DesignXEst1,
                          UseAllLag=FALSE, num_coef=Ts, mylag=0)
OptimWInnov1$par
# [1]  0.81394319 -0.82091886 -0.28142605 -0.18235601 -0.15789591  0.07620144  0.16528796  0.41188530
OptimWInnov1$counts[1] # Number of function evaluation
# function 
# 295 
OptimWInnov1$convergence
# [1] 0





#------------------- Parameter Estimation   ------------------------#
cout <- summary(as.factor(X_hour))
ci_initial <- as.vector(qnorm(cumsum(cout / sum(cout))))[1:(K - 1)]
phi_initial <- stats::pacf(X_hour, plot = FALSE)$acf[1]
par_initial <- c(ci_initial[2:(K-1)] - ci_initial[1], -ci_initial[1],
                 rep(0, NCOL(DesignXEst2) - 1), phi_initial)

if (K == 3) {
  constrLSE = rbind(c(rep(0, length(par_initial)-1),  1),
                    c(rep(0, length(par_initial)-1), -1)) # rho
} else {
  constrLSE <- matrix(0, nrow = K - 2, ncol = length(par_initial))
  constrLSE[1,1] <- 1
  for (ii in 2:(K-2)) constrLSE[ii, c(ii-1, ii)] <- c(-1, 1) # ci monotone
  constrLSE = rbind(constrLSE,
                    c(rep(0, ncol(constrLSE)-1),  1),
                    c(rep(0, ncol(constrLSE)-1), -1)) # rho bounds
}
constrLSE_ci <- c(rep(0, nrow(constrLSE) - 2), -1, -1)
## optimization
OptimWInnov2 = constrOptim(par_initial, f=CLSWInnov,
                          method = "Nelder-Mead", ui=constrLSE,
                          ci=constrLSE_ci, hessian=F,
                          control= list(reltol=1e-05),
                          Xl=X_hour, Xw=X_hour_wide, DesignX=DesignXEst2,
                          UseAllLag=FALSE, num_coef=Ts, mylag=0)
OptimWInnov2$par
# [1]  0.81149549 -0.73383461 -0.69748410 -0.18328755 -0.15562774  0.07269736  0.15913022  0.29352053  0.41096039
OptimWInnov2$counts[1] # Number of function evaluation
# function 
# 410 
OptimWInnov2$convergence
# [1] 0


# long time running, save results for easy use.
# save(OptimWInnov1, OptimWInnov2, 
#      file="grad_app/application/abqrainfall/CLS_innovation.RData")

load("grad_app/application/abqrainfall/CLS_innovation.RData")

par = OptimWInnov2$par
Xl = X_hour
Xw = X_hour_wide
DesignX = DesignXEst2
num_coef = Ts

Ts = length(Xl)
K = nlevels(factor(Xl))
num_par = length(par)

ci = c(0, par[1:(K-2)])
seg = c(-Inf, ci, Inf)
theta = par[(K-1):(num_par-1)]
rho = par[num_par]

## mean vector
mst = DesignX%*%theta

## Conditional expectation
tmp = t(matrix(rep(seg, Ts), nrow=length(seg), ncol=Ts)) - matrix(rep(mst, length(seg)), nrow=Ts)
EXW = pnorm(tmp[,2:(K+1)]) - pnorm(tmp[,1:K])
EXW = EXW[,1:(K-1)]

## necessary lags version
InnovRes = multiInnov_cpp(mst = mst, EXW = EXW, seg = seg, rho = rho, 
                          numcoef=num_coef, mytol=1e-05, mylag=0)
out_pred = MultiOneStepPred_cpp(EXW = EXW, Xw = Xw[,1:(K-1)], HQ = InnovRes$HQ, 
                                numcoef=num_coef, mylag=InnovRes$mylag)
Predw = cbind(out_pred$PredX, 1 - rowSums(out_pred$PredX))

resid  = Xw - Predw          # prediction errors on probability scale
V_list = InnovRes$V          # list of (K-1)x(K-1) innovation cov matrices
ErrErr = out_pred$ErrErr     # Ts x (K-1): innovations u_t - uhat_t

# Standardized one-step-ahead prediction residuals

# dimension of innovations
p_dim <- K - 1                       
std_innov <- matrix(NA_real_, Ts, p_dim)

for (t in seq_len(Ts)) {
  e_t <- ErrErr[t, 1:p_dim]
  V_raw <- V_list[[t]] # (K-1) x (K-1)
  
  ## force to (K-1)x(K-1) matrix [2 by 2]
  V_t <- matrix(V_raw, nrow = p_dim, ncol = p_dim)
  
  # Cholesky: V_t = L_t %*% t(L_t)
  L_t <- chol(V_t)
  
  # solve L_t z_t = e_t
  # z_t = L_t^{-1} e_t
  z_t <- backsolve(L_t, e_t, transpose = FALSE)
  
  std_innov[t, ] <- z_t
}


setEPS()
postscript("grad_app/application/abqrainfall/ApplyFigure3.eps", width = 20, height = 5)


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
