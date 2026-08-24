rm(list = ls())

library(TSAOP)
library(pbivnorm)
library(foreach)
library(doParallel)

tool_file_ar1 <- "~/Desktop/Changepoint/aop/TSAOP/CLSEinnov/tools.R"
tool_file_ar2 <- "~/Desktop/Changepoint/aop/TSAOP/CLSEinnov/tools_ar2.R"

source(tool_file_ar1)
source(tool_file_ar2)

options(scipen = 10)

############################################################
## Simulation setting
############################################################

Ts <- 500
K  <- 3

ciT    <- 0.8615
thetaT <- 0.4308

phi1T <- 0.5
phi2T <- 0.2
rhoT  <- c(phi1T, phi2T)

DesignXT <- matrix(1, nrow = Ts)
parT <- c(ciT, thetaT, phi1T, phi2T)

param_name <- c(paste0("c", 2:(K - 1)), "theta0", "phi1", "phi2")

nSim <- 1000
myseed <- 20251103

lag_max_fit <- 50L
lag_tol_fit <- 1e-8

ncores <- 10

############################################################
## Parallel setup
############################################################

tim_start <- Sys.time()

cl <- parallel::makeCluster(ncores)
on.exit(parallel::stopCluster(cl), add = TRUE)
doParallel::registerDoParallel(cl)

parallel::clusterExport(
  cl,
  c(
    "tool_file_ar1", "tool_file_ar2",
    "Ts", "K", "DesignXT",
    "ciT", "thetaT", "rhoT",
    "myseed", "parT",
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

############################################################
## Run simulation
############################################################

RES <- foreach(
  iii = 1:nSim,
  .combine = "rbind",
  .inorder = TRUE
) %dopar% {
  
  out <- tryCatch({
    
    sim <- ClipSimulation_AR2(
      ci = ciT,
      theta = thetaT,
      rho = rhoT,
      K = K,
      Ts = Ts,
      DesignX = DesignXT,
      seed = myseed + iii
    )
    
    t1 <- Sys.time()
    
    fit <- fit_CLSInnov_AR2(
      X_hour = sim$X_hour,
      X_hour_wide = sim$X_hour_wide,
      DesignXEst = DesignXT,
      lag_max = lag_max_fit,
      lag_tol = lag_tol_fit
    )
    
    t2 <- Sys.time()
    
    n_fun <- if (!is.null(fit$fit$counts)) as.numeric(fit$fit$counts["function"]) else NA
    n_gr  <- if (!is.null(fit$fit$counts)) as.numeric(fit$fit$counts["gradient"]) else NA
    
    c(
      seed = myseed + iii,
      convg = fit$convergence,
      time = as.numeric(difftime(t2, t1, units = "secs")),
      obj.val = fit$value,
      fit$paramEst,
      n.fun = n_fun,
      n.grad = n_gr,
      err = ""
    )
    
  }, error = function(e) {
    
    c(
      seed = myseed + iii,
      convg = 999,
      time = NA,
      obj.val = NA,
      rep(NA, length(parT)),
      n.fun = NA,
      n.grad = NA,
      err = conditionMessage(e)
    )
  })
  
  out
}

tim_end <- Sys.time()




RES <- as.data.frame(RES)
colnames(RES) <- c(
  "seed", "convg", "time", "obj.val",
  param_name,
  "n.fun", "n.grad",
  "err"
)

num_cols <- setdiff(colnames(RES), "err")
RES[num_cols] <- lapply(RES[num_cols], function(x) as.numeric(as.character(x)))

print(table(RES$convg))
print(table(RES$err))
print(mean(RES$time, na.rm = TRUE))
print(tim_end - tim_start)
print(head(RES))




RES0 <- RES[RES$convg == 0, ]

EstTable <- matrix(NA_real_, nrow = 4, ncol = length(parT))
colnames(EstTable) <- param_name
rownames(EstTable) <- c("True", "Mean", "Rel.bias", "SD")

EstTable[1, ] <- parT

if (nrow(RES0) > 0) {
  EstTable[2, ] <- colMeans(RES0[, param_name], na.rm = TRUE)
  EstTable[3, ] <- (EstTable[2, ] - EstTable[1, ]) / EstTable[1, ]
  EstTable[4, ] <- apply(RES0[, param_name], 2, sd, na.rm = TRUE)
}

print(round(EstTable, 4))


save(
  RES,
  RES0,
  EstTable,
  file = paste0(
    "~/Desktop/Changepoint/aop/AMOCAOP/CLS/revision/ar2/",
    "CLSE_AR2_Innovation_T", Ts,
    "_nSim", nSim,
    "_lag", lag_max_fit,
    "_phi1_", phi1T,
    "_phi2_", phi2T,
    ".RData"
  )
)