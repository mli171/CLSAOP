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







par = ParEst2

Xl=X_hour
Xw=X_hour_wide
DesignX=DesignXEst2


Ts = length(Xl)
K = NCOL(Xw)

num_par = length(par)

seg = c(-Inf, 0, par[1:(K-2)], Inf)
theta = par[(K-1):(num_par-1)]
rho = par[num_par]

# mean vector
mst = DesignX%*%theta
CondEX = matrix(0, nrow=Ts, ncol=K)

# t = 1 (using marginal expectation)
tmpprob = pnorm(seg[2:(K+1)]-mst[1]) - pnorm(seg[1:K]-mst[1])
CondEX[1,] = tmpprob

# t > 2 (using conditonal expectation with lag=1)
corr_nom = rho^abs(matrix(1:2 - 1, nrow = 2, ncol = 2, byrow = TRUE) - (1:2 - 1))
diag(corr_nom) = 1
for(t in 2:Ts){
  prevX = Xl[t-1]
  denom = pnorm(seg[prevX+1]-mst[t-1]) - pnorm(seg[prevX]-mst[t-1])
  mu_nom = mst[t:(t-1)]
  for(k in 1:K){
    a_nom = seg[c(k, prevX)]
    b_nom = seg[c(k, prevX)+1]
    nom = pmvnorm(lower=a_nom, upper=b_nom, mean=mu_nom, corr=corr_nom)
    CondEX[t,k] = nom/denom
  }
}

resid = Xw - CondEX



setEPS()
postscript("grad_app/application/abqrainfall/ApplyFigure2.eps", width = 17, height = 10)


par(mfrow=c(2,2), 
    mar = c(5.5, 6, 3.5, 6),
    cex.main=1.7, cex.lab=1.7, cex.axis=1.7)
for(k in 1:K){
  acf(resid[,k], 
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
  mtext(xlab, side = 1, line = 2.8, cex = 1.7)
  
  invisible(list(relfreq = counts, H = H, n = Ti, p0 = p0, se = se,
                 breaks_pos = at))
}



pit_diag(CondEX, X_hour, H = 50, left_margin = 6.5,
         method = "expected",
         main = "PIT histograms",
         xlab = "Probability Integral Transform (U)", 
         x_breaks = seq(from=0, to=1, by=0.1))

dev.off()

