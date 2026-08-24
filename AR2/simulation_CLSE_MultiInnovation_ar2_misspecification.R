rm(list = ls())

options(scipen = 10)

library(TSAOP)
library(pbivnorm)
library(foreach)
library(doParallel)

tool_file_ar1 <- "~/Desktop/Changepoint/aop/TSAOP/CLSEinnov/tools.R"
tool_file_ar2 <- "~/Desktop/Changepoint/aop/TSAOP/CLSEinnov/tools_ar2.R"

source(tool_file_ar1)
source(tool_file_ar2)

extract_numeric_param <- function(fit, target_len = 3L) {
  
  candidate_names <- c(
    "paramEst", "par", "coef", "coefficients",
    "estimate", "est", "Est", "Estimate"
  )
  
  for (nm in candidate_names) {
    if (!is.null(fit[[nm]])) {
      tmp <- fit[[nm]]
      
      if (is.data.frame(tmp) || is.matrix(tmp)) {
        tmp <- as.numeric(tmp)
      } else {
        tmp <- as.numeric(tmp)
      }
      
      tmp <- tmp[is.finite(tmp)]
      
      if (length(tmp) >= target_len) {
        return(tmp[1:target_len])
      }
    }
  }
  
}

extract_value <- function(fit) {
  
  candidate_names <- c("value", "obj.val", "objective", "obj", "SSE", "sse")
  
  for (nm in candidate_names) {
    if (!is.null(fit[[nm]])) {
      tmp <- as.numeric(fit[[nm]])
      tmp <- tmp[is.finite(tmp)]
      if (length(tmp) >= 1L) return(tmp[1])
    }
  }
  
  NA_real_
}

extract_convergence <- function(fit) {
  
  candidate_names <- c("convergence", "convg", "conv", "code")
  
  for (nm in candidate_names) {
    if (!is.null(fit[[nm]])) {
      tmp <- as.numeric(fit[[nm]])
      tmp <- tmp[is.finite(tmp)]
      if (length(tmp) >= 1L) return(tmp[1])
    }
  }
  
  0
}

fit_CLSInnov_AR1_raw <- function(X_hour, X_hour_wide, DesignXEst,
                                 lag_max = 50L,
                                 lag_tol = 1e-8) {
  
  K <- ncol(X_hour_wide)
  q <- ncol(DesignXEst)
  
  obj <- make_aop1_cls_optim(
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
  
  theta_initial <- c(
    -ci_initial[1],
    rep(0, q - 1L)
  )
  
  rho_initial <- tryCatch({
    as.numeric(stats::pacf(X_hour, plot = FALSE)$acf[1])
  }, error = function(e) {
    0.3
  })
  
  if (!is.finite(rho_initial)) rho_initial <- 0.3
  rho_initial <- max(min(rho_initial, 0.8), -0.8)
  
  par_initial <- c(cuts_initial, theta_initial, rho_initial)
  
  con <- make_constraints_aop1(
    K = K,
    q = q,
    delta = 1e-6,
    delta0 = 1e-6,
    eps_rho = 1e-6
  )
  
  fit <- constrOptim(
    theta = par_initial,
    f = obj$fn,
    grad = obj$gr,
    ui = con$ui,
    ci = con$ci,
    method = "BFGS",
    hessian = FALSE,
    control = list(maxit = 300, reltol = 1e-7)
  )
  
  list(
    paramEst = fit$par,
    convergence = fit$convergence,
    value = fit$value,
    par_initial = par_initial,
    fit = fit
  )
}


fit_CLSLag1_AR1_raw <- function(X_hour, DesignXEst) {
  
  fit <- aopts(
    y = X_hour,
    X = DesignXEst,
    method = "clse"
  )
  
  paramEst <- extract_numeric_param(fit, target_len = 3L)
  value <- extract_value(fit)
  convergence <- extract_convergence(fit)
  
  list(
    paramEst = paramEst,
    convergence = convergence,
    value = value,
    fit = fit
  )
}







## Simulation start ##

Ts <- 500
K  <- 3

ciT    <- 0.8615
thetaT <- 0.4308

phi1T <- 0.5
phi2T <- 0.2
rhoT  <- c(phi1T, phi2T)

DesignXT <- matrix(1, nrow = Ts)

## Correct AR(2) parameter vector
parT_ar2 <- c(ciT, thetaT, phi1T, phi2T)

## Reference values for misspecified AR(1) fits.
rho1_ref <- phi1T / (1 - phi2T)
parT_ar1_ref <- c(ciT, thetaT, rho1_ref)

param_name_ar2 <- c("c2.AR2", "theta0.AR2", "phi1.AR2", "phi2.AR2")
param_name_ar1_innov <- c("c2.AR1.innov", "theta0.AR1.innov", "rho.AR1.innov")
param_name_ar1_lag1  <- c("c2.AR1.lag1",  "theta0.AR1.lag1",  "rho.AR1.lag1")

nSim <- 1000
myseed <- 20251103

lag_max_fit <- 50L
lag_tol_fit <- 1e-8

ncores <- 10

out_dir <- "~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/ar2/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

tim_start <- Sys.time()

cl <- parallel::makeCluster(ncores)
on.exit(parallel::stopCluster(cl), add = TRUE)
doParallel::registerDoParallel(cl)

parallel::clusterExport(
  cl,
  c(
    "tool_file_ar1", "tool_file_ar2",
    "extract_numeric_param", "extract_value", "extract_convergence",
    "fit_CLSInnov_AR1_raw",
    "fit_CLSLag1_AR1_raw",
    "Ts", "K", "DesignXT",
    "ciT", "thetaT", "rhoT",
    "myseed",
    "parT_ar2", "parT_ar1_ref",
    "lag_max_fit", "lag_tol_fit"
  ),
  envir = environment()
)

parallel::clusterEvalQ(cl, {
  library(TSAOP)
  library(pbivnorm)
  source(tool_file_ar1)
  source(tool_file_ar2)
  NULL
})


RES <- foreach(
  iii = 1:nSim,
  .combine = "rbind",
  .inorder = TRUE
) %dopar% {
  
  sim <- ClipSimulation_AR2(
    ci = ciT,
    theta = thetaT,
    rho = rhoT,
    K = K,
    Ts = Ts,
    DesignX = DesignXT,
    seed = myseed + iii
  )
  
  fit_ar2_out <- tryCatch({
    
    t1 <- Sys.time()
    
    fit_ar2 <- fit_CLSInnov_AR2(
      X_hour = sim$X_hour,
      X_hour_wide = sim$X_hour_wide,
      DesignXEst = DesignXT,
      lag_max = lag_max_fit,
      lag_tol = lag_tol_fit
    )
    
    t2 <- Sys.time()
    
    n_fun <- if (!is.null(fit_ar2$fit$counts)) {
      as.numeric(fit_ar2$fit$counts["function"])
    } else {
      NA
    }
    
    n_gr <- if (!is.null(fit_ar2$fit$counts)) {
      as.numeric(fit_ar2$fit$counts["gradient"])
    } else {
      NA
    }
    
    c(
      convg.AR2 = fit_ar2$convergence,
      time.AR2 = as.numeric(difftime(t2, t1, units = "secs")),
      obj.AR2 = fit_ar2$value,
      fit_ar2$paramEst,
      n.fun.AR2 = n_fun,
      n.grad.AR2 = n_gr,
      err.AR2 = ""
    )
    
  }, error = function(e) {
    
    c(
      convg.AR2 = 999,
      time.AR2 = NA,
      obj.AR2 = NA,
      rep(NA, length(parT_ar2)),
      n.fun.AR2 = NA,
      n.grad.AR2 = NA,
      err.AR2 = conditionMessage(e)
    )
  })
  
  
  fit_ar1_innov_out <- tryCatch({
    
    t1 <- Sys.time()
    
    fit_ar1_innov <- fit_CLSInnov_AR1_raw(
      X_hour = sim$X_hour,
      X_hour_wide = sim$X_hour_wide,
      DesignXEst = DesignXT,
      lag_max = lag_max_fit,
      lag_tol = lag_tol_fit
    )
    
    t2 <- Sys.time()
    
    n_fun <- if (!is.null(fit_ar1_innov$fit$counts)) {
      as.numeric(fit_ar1_innov$fit$counts["function"])
    } else {
      NA
    }
    
    n_gr <- if (!is.null(fit_ar1_innov$fit$counts)) {
      as.numeric(fit_ar1_innov$fit$counts["gradient"])
    } else {
      NA
    }
    
    c(
      convg.AR1.innov = fit_ar1_innov$convergence,
      time.AR1.innov = as.numeric(difftime(t2, t1, units = "secs")),
      obj.AR1.innov = fit_ar1_innov$value,
      fit_ar1_innov$paramEst,
      n.fun.AR1.innov = n_fun,
      n.grad.AR1.innov = n_gr,
      err.AR1.innov = ""
    )
    
  }, error = function(e) {
    
    c(
      convg.AR1.innov = 999,
      time.AR1.innov = NA,
      obj.AR1.innov = NA,
      rep(NA, length(parT_ar1_ref)),
      n.fun.AR1.innov = NA,
      n.grad.AR1.innov = NA,
      err.AR1.innov = conditionMessage(e)
    )
  })
  
  fit_ar1_lag1_out <- tryCatch({
    
    t1 <- Sys.time()
    
    fit_ar1_lag1 <- fit_CLSLag1_AR1_raw(
      X_hour = sim$X_hour,
      DesignXEst = DesignXT
    )
    
    t2 <- Sys.time()
    
    c(
      convg.AR1.lag1 = fit_ar1_lag1$convergence,
      time.AR1.lag1 = as.numeric(difftime(t2, t1, units = "secs")),
      obj.AR1.lag1 = fit_ar1_lag1$value,
      fit_ar1_lag1$paramEst,
      err.AR1.lag1 = ""
    )
    
  }, error = function(e) {
    
    c(
      convg.AR1.lag1 = 999,
      time.AR1.lag1 = NA,
      obj.AR1.lag1 = NA,
      rep(NA, length(parT_ar1_ref)),
      err.AR1.lag1 = conditionMessage(e)
    )
  })
  
  c(
    seed = myseed + iii,
    fit_ar2_out,
    fit_ar1_innov_out,
    fit_ar1_lag1_out
  )
}

tim_end <- Sys.time()


RES <- as.data.frame(RES)

colnames(RES) <- c(
  "seed",
  "convg.AR2", "time.AR2", "obj.AR2",
  param_name_ar2,
  "n.fun.AR2", "n.grad.AR2", "err.AR2",
  "convg.AR1.innov", "time.AR1.innov", "obj.AR1.innov",
  param_name_ar1_innov,
  "n.fun.AR1.innov", "n.grad.AR1.innov", "err.AR1.innov",
  "convg.AR1.lag1", "time.AR1.lag1", "obj.AR1.lag1",
  param_name_ar1_lag1,
  "err.AR1.lag1"
)

char_cols <- c("err.AR2", "err.AR1.innov", "err.AR1.lag1")
num_cols <- setdiff(colnames(RES), char_cols)

RES[num_cols] <- lapply(RES[num_cols], function(x) as.numeric(as.character(x)))

print(table(RES$convg.AR2))
print(table(RES$convg.AR1.innov))
print(table(RES$convg.AR1.lag1))

print(table(RES$err.AR2))
print(table(RES$err.AR1.innov))
print(table(RES$err.AR1.lag1))

print(mean(RES$time.AR2, na.rm = TRUE))
print(mean(RES$time.AR1.innov, na.rm = TRUE))
print(mean(RES$time.AR1.lag1, na.rm = TRUE))

print(tim_end - tim_start)

print(mean(RES$obj.AR2, na.rm = TRUE))
print(mean(RES$obj.AR1.innov, na.rm = TRUE))
print(mean(RES$obj.AR1.lag1, na.rm = TRUE))

print(head(RES))


RES_AR2 <- RES[RES$convg.AR2 == 0, ]
RES_AR1_INNOV <- RES[RES$convg.AR1.innov == 0, ]
RES_AR1_LAG1 <- RES[RES$convg.AR1.lag1 == 0, ]

RES_ALL <- RES[
  RES$convg.AR2 == 0 &
    RES$convg.AR1.innov == 0 &
    RES$convg.AR1.lag1 == 0,
]


## Correct AR(2) innovation fit summary


EstTable_AR2 <- matrix(NA_real_, nrow = 4, ncol = length(parT_ar2))
colnames(EstTable_AR2) <- c("c2", "theta0", "phi1", "phi2")
rownames(EstTable_AR2) <- c("True", "Mean", "Rel.bias", "SD")

EstTable_AR2[1, ] <- parT_ar2

if (nrow(RES_AR2) > 0) {
  EstTable_AR2[2, ] <- colMeans(RES_AR2[, param_name_ar2], na.rm = TRUE)
  EstTable_AR2[3, ] <- (EstTable_AR2[2, ] - EstTable_AR2[1, ]) / EstTable_AR2[1, ]
  EstTable_AR2[4, ] <- apply(RES_AR2[, param_name_ar2], 2, sd, na.rm = TRUE)
}


EstTable_AR1_INNOV <- matrix(NA_real_, nrow = 4, ncol = length(parT_ar1_ref))
colnames(EstTable_AR1_INNOV) <- c("c2", "theta0", "rho")
rownames(EstTable_AR1_INNOV) <- c("Reference", "Mean", "Rel.diff", "SD")

EstTable_AR1_INNOV[1, ] <- parT_ar1_ref

if (nrow(RES_AR1_INNOV) > 0) {
  EstTable_AR1_INNOV[2, ] <- colMeans(
    RES_AR1_INNOV[, param_name_ar1_innov],
    na.rm = TRUE
  )
  EstTable_AR1_INNOV[3, ] <- (
    EstTable_AR1_INNOV[2, ] - EstTable_AR1_INNOV[1, ]
  ) / EstTable_AR1_INNOV[1, ]
  EstTable_AR1_INNOV[4, ] <- apply(
    RES_AR1_INNOV[, param_name_ar1_innov],
    2,
    sd,
    na.rm = TRUE
  )
}


EstTable_AR1_LAG1 <- matrix(NA_real_, nrow = 4, ncol = length(parT_ar1_ref))
colnames(EstTable_AR1_LAG1) <- c("c2", "theta0", "rho")
rownames(EstTable_AR1_LAG1) <- c("Reference", "Mean", "Rel.diff", "SD")

EstTable_AR1_LAG1[1, ] <- parT_ar1_ref

if (nrow(RES_AR1_LAG1) > 0) {
  EstTable_AR1_LAG1[2, ] <- colMeans(
    RES_AR1_LAG1[, param_name_ar1_lag1],
    na.rm = TRUE
  )
  EstTable_AR1_LAG1[3, ] <- (
    EstTable_AR1_LAG1[2, ] - EstTable_AR1_LAG1[1, ]
  ) / EstTable_AR1_LAG1[1, ]
  EstTable_AR1_LAG1[4, ] <- apply(
    RES_AR1_LAG1[, param_name_ar1_lag1],
    2,
    sd,
    na.rm = TRUE
  )
}

print(round(EstTable_AR2, 4))
print(round(EstTable_AR1_INNOV, 4))
print(round(EstTable_AR1_LAG1, 4))




ObjTable <- c(
  Mean.obj.AR2 = mean(RES_AR2$obj.AR2, na.rm = TRUE),
  Mean.obj.AR1.innov = mean(RES_AR1_INNOV$obj.AR1.innov, na.rm = TRUE),
  Mean.obj.AR1.lag1 = mean(RES_AR1_LAG1$obj.AR1.lag1, na.rm = TRUE),
  Mean.diff.AR1.innov.minus.AR2 = mean(
    RES_ALL$obj.AR1.innov - RES_ALL$obj.AR2,
    na.rm = TRUE
  ),
  Median.diff.AR1.innov.minus.AR2 = median(
    RES_ALL$obj.AR1.innov - RES_ALL$obj.AR2,
    na.rm = TRUE
  ),
  Prop.AR2.obj.smaller.than.AR1.innov = mean(
    RES_ALL$obj.AR2 < RES_ALL$obj.AR1.innov,
    na.rm = TRUE
  )
)

print(round(ObjTable, 4))
