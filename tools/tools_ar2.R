BIG <- 1e30

is_stationary_ar2 <- function(phi1, phi2, eps = 1e-8) {
  is.finite(phi1) && is.finite(phi2) &&
    abs(phi2) < 1 - eps &&
    (phi1 + phi2) < 1 - eps &&
    (phi2 - phi1) < 1 - eps
}

ar2_innov_var_unit <- function(phi1, phi2) {
  1 - phi2^2 - phi1^2 * (1 + phi2) / (1 - phi2)
}

ClipSimulation_AR2 <- function(ci, theta, rho, K, Ts, DesignX, seed = NULL) {
  
  if (length(ci) != K - 2) stop("Number of cut points and categories NOT match.")
  if (length(theta) != NCOL(DesignX)) stop("Design matrix and theta do not match.")
  if (length(rho) != 2L) stop("rho must be c(phi1, phi2) for AR(2).")
  
  phi1 <- rho[1]
  phi2 <- rho[2]
  
  if (!is_stationary_ar2(phi1, phi2)) stop("AR(2) is not stationary.")
  
  sig2e <- ar2_innov_var_unit(phi1, phi2)
  if (!is.finite(sig2e) || sig2e <= 0) stop("Invalid AR(2) innovation variance.")
  
  if (!is.null(seed)) set.seed(seed)
  
  mst <- as.vector(DesignX %*% theta)
  
  et <- as.numeric(arima.sim(
    n = Ts,
    list(ar = c(phi1, phi2)),
    sd = sqrt(sig2e)
  ))
  
  Z <- et + mst
  
  ci2 <- c(0, ci)  # c1 fixed at 0
  X_hour <- sapply(as.vector(Z), function(x) length(which(x > ci2)) + 1L)
  
  X_hour_wide <- matrix(0L, nrow = Ts, ncol = K)
  for (t in 1:Ts) X_hour_wide[t, X_hour[t]] <- 1L
  
  list(Z = Z, X_hour = X_hour, X_hour_wide = X_hour_wide)
}

multikappa_pbivnorm_ar2 <- function(i, j, msti, mstj, EXwi, EXwj, seg, acf_vec) {
  
  K <- length(seg) - 1L
  m <- K - 1L
  
  if (i == j) {
    p <- as.numeric(EXwi)
    return(diag(p, m, m) - tcrossprod(p))
  }
  
  h <- abs(i - j)
  r <- acf_vec[h + 1L]
  
  if (!is.finite(r) || abs(r) >= 1) {
    return(matrix(NA_real_, nrow = m, ncol = m))
  }
  
  lo1 <- seg[1:m]       - msti
  up1 <- seg[2:(m + 1)] - msti
  lo2 <- seg[1:m]       - mstj
  up2 <- seg[2:(m + 1)] - mstj
  
  x_up <- rep(up1, each = m)
  x_lo <- rep(lo1, each = m)
  y_up <- rep(up2, times = m)
  y_lo <- rep(lo2, times = m)
  
  Fuu <- pbiv_safe(x_up, y_up, r)
  Flu <- pbiv_safe(x_lo, y_up, r)
  Ful <- pbiv_safe(x_up, y_lo, r)
  Fll <- pbiv_safe(x_lo, y_lo, r)
  
  pij <- Fuu - Flu - Ful + Fll
  EWW <- matrix(pij, nrow = m, byrow = TRUE)
  
  EWW - tcrossprod(EXwi, EXwj)
}

CLSWInnov_AR2 <- function(par, Xw, DesignX,
                          lag_max = 50L,
                          lag_tol = 1e-8,
                          jitter = 1e-10) {
  
  Xw <- as.matrix(Xw)
  DesignX <- as.matrix(DesignX)
  
  Ts <- nrow(Xw)
  K  <- ncol(Xw)
  q  <- ncol(DesignX)
  m  <- K - 1L
  
  p_cuts <- K - 2L
  
  cuts  <- as.numeric(par[1:p_cuts])
  theta <- as.numeric(par[p_cuts + (1:q)])
  phi1  <- as.numeric(par[p_cuts + q + 1L])
  phi2  <- as.numeric(par[p_cuts + q + 2L])
  
  if (any(!is.finite(cuts)) || is.unsorted(cuts, strictly = TRUE)) return(BIG)
  if (cuts[1] <= 0) return(BIG)
  if (!is_stationary_ar2(phi1, phi2)) return(BIG)
  
  sig2e <- ar2_innov_var_unit(phi1, phi2)
  if (!is.finite(sig2e) || sig2e <= 0) return(BIG)
  
  acf_vec <- tryCatch(
    stats::ARMAacf(ar = c(phi1, phi2), lag.max = Ts - 1L),
    error = function(e) NULL
  )
  if (is.null(acf_vec) || any(!is.finite(acf_vec))) return(BIG)
  
  cuts_all <- c(0, cuts)
  seg <- c(-Inf, cuts_all, Inf)
  
  mst <- as.numeric(DesignX %*% theta)
  
  Fmat <- stats::pnorm(outer(mst, seg, function(mu, c) c - mu))
  Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] - Fmat[, 1:K, drop = FALSE]
  EXW  <- Pmat[, 1:m, drop = FALSE]
  u    <- Xw[, 1:m, drop = FALSE] - EXW
  
  kappa <- function(i, j) {
    multikappa_pbivnorm_ar2(
      i, j,
      msti = mst[i],
      mstj = mst[j],
      EXwi = EXW[i, ],
      EXwj = EXW[j, ],
      seg = seg,
      acf_vec = acf_vec
    )
  }
  
  fit <- tryCatch(
    innovations_mv(
      kappa = kappa,
      u = u,
      jitter = jitter,
      symmetrize = TRUE,
      compute_innov = TRUE,
      lag_max = lag_max,
      lag_tol = lag_tol,
      ahead = NULL
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(BIG)
  
  PredX <- EXW + fit$uhat
  Predw <- cbind(PredX, 1 - rowSums(PredX))
  
  val <- sum((Xw - Predw)^2)
  if (!is.finite(val)) val <- BIG
  
  val
}

make_constraints_aop2 <- function(K, q,
                                  delta = 1e-6,
                                  delta0 = 1e-6,
                                  eps_ar = 1e-6) {
  
  p_cuts <- K - 2L
  p <- p_cuts + q + 2L
  
  phi1_idx <- p - 1L
  phi2_idx <- p
  
  ui <- NULL
  ci <- NULL
  
  ## c2 >= delta0
  if (p_cuts >= 1L) {
    row <- rep(0, p)
    row[1] <- 1
    ui <- rbind(ui, row)
    ci <- c(ci, delta0)
  }
  
  ## ordered cuts
  if (p_cuts >= 2L) {
    for (j in 1:(p_cuts - 1L)) {
      row <- rep(0, p)
      row[j + 1L] <-  1
      row[j]      <- -1
      ui <- rbind(ui, row)
      ci <- c(ci, delta)
    }
  }
  
  ## phi2 <= 1 - eps
  row <- rep(0, p)
  row[phi2_idx] <- -1
  ui <- rbind(ui, row)
  ci <- c(ci, -(1 - eps_ar))
  
  ## phi2 >= -1 + eps
  row <- rep(0, p)
  row[phi2_idx] <- 1
  ui <- rbind(ui, row)
  ci <- c(ci, -1 + eps_ar)
  
  ## phi1 + phi2 <= 1 - eps
  row <- rep(0, p)
  row[phi1_idx] <- -1
  row[phi2_idx] <- -1
  ui <- rbind(ui, row)
  ci <- c(ci, -(1 - eps_ar))
  
  ## phi2 - phi1 <= 1 - eps
  row <- rep(0, p)
  row[phi1_idx] <-  1
  row[phi2_idx] <- -1
  ui <- rbind(ui, row)
  ci <- c(ci, -1 + eps_ar)
  
  list(ui = ui, ci = ci)
}

# fit_CLSInnov_AR2 <- function(X_hour, X_hour_wide, DesignXEst,
#                              lag_max = 50L,
#                              lag_tol = 1e-8) {
#   
#   K <- ncol(X_hour_wide)
#   q <- ncol(DesignXEst)
#   
#   cout <- tabulate(X_hour, nbins = K)
#   cumprob <- cumsum(cout / sum(cout))
#   cumprob <- pmin(pmax(cumprob, 1e-6), 1 - 1e-6)
#   
#   ci_initial <- qnorm(cumprob)[1:(K - 1)]
#   
#   cuts_initial <- ci_initial[2:(K - 1)] - ci_initial[1]
#   cuts_initial <- pmax(cuts_initial, 1e-5)
#   cuts_initial <- sort(cuts_initial)
#   
#   theta_initial <- c(
#     -ci_initial[1],
#     rep(0, q - 1L)
#   )
#   
#   yy <- as.numeric(scale(X_hour, center = TRUE, scale = FALSE))
#   
#   phi_initial <- tryCatch({
#     arfit <- stats::arima(yy, order = c(2, 0, 0), include.mean = FALSE)
#     tmp <- as.numeric(arfit$coef[1:2])
#     if (!is_stationary_ar2(tmp[1], tmp[2])) c(0.3, 0.2) else tmp
#   }, error = function(e) {
#     c(0.3, 0.2)
#   })
#   
#   par_initial <- c(cuts_initial, theta_initial, phi_initial)
#   
#   con <- make_constraints_aop2(K = K, q = q)
#   
#   fit <- constrOptim(
#     theta = par_initial,
#     f = CLSWInnov_AR2,
#     # grad = NULL,
#     grad = num_grad_CLSWInnov_AR2,
#     ui = con$ui,
#     ci = con$ci,
#     method = "BFGS",
#     hessian = FALSE,
#     control = list(maxit = 300, reltol = 1e-7),
#     Xw = X_hour_wide,
#     DesignX = DesignXEst,
#     lag_max = lag_max,
#     lag_tol = lag_tol,
#     jitter = 1e-10
#   )
#   
#   list(
#     paramEst = fit$par,
#     convergence = fit$convergence,
#     value = fit$value,
#     par_initial = par_initial
#   )
# }


num_grad_CLSWInnov_AR2 <- function(theta, Xw, DesignX,
                                   lag_max = 50L,
                                   lag_tol = 1e-8,
                                   jitter = 1e-10,
                                   eps = 1e-5) {
  
  theta <- as.numeric(theta)
  p <- length(theta)
  g <- numeric(p)
  
  f0 <- CLSWInnov_AR2(
    par = theta,
    Xw = Xw,
    DesignX = DesignX,
    lag_max = lag_max,
    lag_tol = lag_tol,
    jitter = jitter
  )
  
  if (!is.finite(f0) || f0 >= BIG / 10) {
    return(rep(0, p))
  }
  
  for (j in seq_len(p)) {
    
    step <- eps * max(1, abs(theta[j]))
    
    theta_p <- theta
    theta_m <- theta
    
    theta_p[j] <- theta_p[j] + step
    theta_m[j] <- theta_m[j] - step
    
    fp <- CLSWInnov_AR2(
      par = theta_p,
      Xw = Xw,
      DesignX = DesignX,
      lag_max = lag_max,
      lag_tol = lag_tol,
      jitter = jitter
    )
    
    fm <- CLSWInnov_AR2(
      par = theta_m,
      Xw = Xw,
      DesignX = DesignX,
      lag_max = lag_max,
      lag_tol = lag_tol,
      jitter = jitter
    )
    
    if (is.finite(fp) && is.finite(fm) && fp < BIG / 10 && fm < BIG / 10) {
      g[j] <- (fp - fm) / (2 * step)
    } else if (is.finite(fp) && fp < BIG / 10) {
      g[j] <- (fp - f0) / step
    } else if (is.finite(fm) && fm < BIG / 10) {
      g[j] <- (f0 - fm) / step
    } else {
      g[j] <- 0
    }
  }
  
  g[!is.finite(g)] <- 0
  g
}








########## analytic gradient

ar2_acf_deriv <- function(phi1, phi2, maxlag) {
  
  rho <- numeric(maxlag + 1L)
  d1  <- numeric(maxlag + 1L)  # d rho_h / d phi1
  d2  <- numeric(maxlag + 1L)  # d rho_h / d phi2
  
  rho[1] <- 1
  d1[1]  <- 0
  d2[1]  <- 0
  
  if (maxlag >= 1L) {
    rho[2] <- phi1 / (1 - phi2)
    d1[2]  <- 1 / (1 - phi2)
    d2[2]  <- phi1 / (1 - phi2)^2
  }
  
  if (maxlag >= 2L) {
    for (h in 2:maxlag) {
      idx  <- h + 1L
      idx1 <- h
      idx2 <- h - 1L
      
      rho[idx] <- phi1 * rho[idx1] + phi2 * rho[idx2]
      d1[idx]  <- rho[idx1] + phi1 * d1[idx1] + phi2 * d1[idx2]
      d2[idx]  <- rho[idx2] + phi1 * d2[idx1] + phi2 * d2[idx2]
    }
  }
  
  list(rho = rho, dphi1 = d1, dphi2 = d2)
}

make_du_aop2 <- function(cuts, theta, DesignX) {
  
  DesignX <- as.matrix(DesignX)
  Ts <- nrow(DesignX)
  q  <- ncol(DesignX)
  
  cuts  <- as.numeric(cuts)
  theta <- as.numeric(theta)
  
  p_cuts <- length(cuts)
  K <- p_cuts + 2L
  m <- K - 1L
  
  p <- p_cuts + q + 2L
  
  cuts_all <- c(0, cuts)
  seg <- c(-Inf, cuts_all, Inf)
  
  mu <- as.numeric(DesignX %*% theta)
  
  pdf_seg <- dnorm(outer(mu, seg,      function(a, b) b - a))
  pdf_cut <- dnorm(outer(mu, cuts_all, function(a, b) b - a))
  
  ## marginal mean derivs
  dm_dmu <- pdf_seg[, 1:m, drop = FALSE] -
    pdf_seg[, 2:(m + 1L), drop = FALSE]
  
  du_arr <- array(0, dim = c(Ts, m, p))
  
  ## cutpoint parameters
  if (p_cuts > 0L) {
    for (t in 1:Ts) {
      dm_dc <- matrix(0, nrow = m, ncol = p_cuts)
      
      for (r in 1:p_cuts) {
        j <- r + 1L # free cut points from c2
        val <- pdf_cut[t, j]
        
        if (j <= m) {
          dm_dc[j, r] <- dm_dc[j, r] + val # lower category
        }
        
        if ((j + 1L) <= m) {
          dm_dc[j + 1L, r] <- dm_dc[j + 1L, r] - val # upper category
        }
      }
      
      ## u_t = W_t - E(W_t)
      du_arr[t, , 1:p_cuts] <- -dm_dc
    }
  }
  
  ## theta parameters
  theta_idx <- p_cuts + (1:q)
  
  for (t in 1:Ts) {
    dm_dtheta <- outer(dm_dmu[t, ], DesignX[t, ])
    du_arr[t, , theta_idx] <- -dm_dtheta
  }
  
  du_arr
}

make_dkappa_aop2 <- function(cuts, theta, phi1, phi2, DesignX, eps_s = 1e-12) {
  
  DesignX <- as.matrix(DesignX)
  Ts <- nrow(DesignX)
  q  <- ncol(DesignX)
  
  cuts  <- as.numeric(cuts)
  theta <- as.numeric(theta)
  phi1  <- as.numeric(phi1)
  phi2  <- as.numeric(phi2)
  
  p_cuts <- length(cuts)
  K <- p_cuts + 2L
  m <- K - 1L
  
  ## cuts + theta + phi1 + phi2
  p <- p_cuts + q + 2L
  phi1_idx <- p - 1L
  phi2_idx <- p
  
  cuts_all <- c(0, cuts)
  seg <- c(-Inf, cuts_all, Inf)
  
  mu <- as.numeric(DesignX %*% theta)
  
  ## marginal probabilities
  Fmat <- pnorm(outer(mu, seg, function(a,b) b - a))
  Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] - Fmat[, 1:K, drop = FALSE]
  EXW  <- Pmat[, 1:m, drop = FALSE]
  
  ## derivative of marginal mean
  dm <- array(0, dim = c(Ts, m, p))
  
  pdf_seg <- dnorm(outer(mu, seg,      function(a,b) b - a))
  pdf_cut <- dnorm(outer(mu, cuts_all, function(a,b) b - a))
  
  dm_dmu <- pdf_seg[, 1:m, drop = FALSE] - pdf_seg[, 2:(m + 1L), drop = FALSE]
  
  theta_idx <- p_cuts + (1:q)
  for (t in 1:Ts) {
    dm[t, , theta_idx] <- outer(dm_dmu[t, ], DesignX[t, ])
  }
  
  if (p_cuts > 0L) {
    for (t in 1:Ts) {
      dm_dc <- matrix(0, nrow = m, ncol = p_cuts)
      for (r in 1:p_cuts) {
        j <- r + 1L
        val <- pdf_cut[t, j]
        if (j <= m) dm_dc[j, r] <- dm_dc[j, r] + val
        if ((j + 1L) <= m) dm_dc[j + 1L, r] <- dm_dc[j + 1L, r] - val
      }
      dm[t, , 1:p_cuts] <- dm_dc
    }
  }
  
  ## AR(2) autocorrelation and derivatives
  acf_deriv <- ar2_acf_deriv(phi1 = phi1, phi2 = phi2, maxlag = Ts - 1L)
  rho_vec   <- acf_deriv$rho
  dr_phi1   <- acf_deriv$dphi1
  dr_phi2   <- acf_deriv$dphi2
  
  ## BVN derivative helpers
  Phi2_safe <- function(a, b, r) {
    if (is.infinite(a) && a < 0) return(0)
    if (is.infinite(b) && b < 0) return(0)
    if (is.infinite(a) && a > 0 && is.infinite(b) && b > 0) return(1)
    if (is.infinite(a) && a > 0) return(pnorm(b))
    if (is.infinite(b) && b > 0) return(pnorm(a))
    pbivnorm::pbivnorm(a, b, r)
  }
  
  s_safe <- function(r) {
    s2 <- 1 - r^2
    if (!is.finite(s2) || s2 <= 0) s2 <- eps_s
    sqrt(s2)
  }
  
  dPhi2_da <- function(a, b, r) {
    if (!is.finite(a)) return(0)
    if (is.infinite(b) && b < 0) return(0)
    s <- s_safe(r)
    dnorm(a) * pnorm((b - r * a) / s)
  }
  
  dPhi2_db <- function(a, b, r) {
    if (!is.finite(b)) return(0)
    if (is.infinite(a) && a < 0) return(0)
    s <- s_safe(r)
    dnorm(b) * pnorm((a - r * b) / s)
  }
  
  phi2_pdf <- function(a, b, r) {
    if (!is.finite(a) || !is.finite(b)) return(0)
    s2 <- 1 - r^2
    if (!is.finite(s2) || s2 <= 0) s2 <- eps_s
    z <- (a^2 - 2 * r * a * b + b^2) / s2
    (1 / (2 * pi * sqrt(s2))) * exp(-0.5 * z)
  }
  
  dPhi2_dr <- function(a, b, r) phi2_pdf(a, b, r)
  
  rect_partials <- function(al, au, bl, bu, r) {
    
    Fuu <- Phi2_safe(au, bu, r)
    Flu <- Phi2_safe(al, bu, r)
    Ful <- Phi2_safe(au, bl, r)
    Fll <- Phi2_safe(al, bl, r)
    P   <- Fuu - Flu - Ful + Fll
    
    dau <-  dPhi2_da(au, bu, r) - dPhi2_da(au, bl, r)
    dal <- -dPhi2_da(al, bu, r) + dPhi2_da(al, bl, r)
    dbu <-  dPhi2_db(au, bu, r) - dPhi2_db(al, bu, r)
    dbl <- -dPhi2_db(au, bl, r) + dPhi2_db(al, bl, r)
    dr  <-  dPhi2_dr(au, bu, r) - dPhi2_dr(al, bu, r) -
      dPhi2_dr(au, bl, r) + dPhi2_dr(al, bl, r)
    
    dmu_i <- -(dau + dal)
    dmu_j <- -(dbu + dbl)
    
    list(
      P = P,
      dmu_i = dmu_i,
      dmu_j = dmu_j,
      dr = dr,
      dau = dau,
      dal = dal,
      dbu = dbu,
      dbl = dbl
    )
  }
  
  dkappa <- function(i, j) {
    
    i <- as.integer(i)
    j <- as.integer(j)
    
    if (i < j) {
      A <- dkappa(j, i)
      return(aperm(A, c(2, 1, 3)))
    }
    
    mi  <- as.numeric(EXW[i, ])
    mj  <- as.numeric(EXW[j, ])
    dmi <- matrix(dm[i, , ], nrow = m, ncol = p)
    dmj <- matrix(dm[j, , ], nrow = m, ncol = p)
    
    out <- array(0, dim = c(m, m, p))
    
    ## same-time covariance
    if (i == j) {
      for (r in 1:p) {
        dmr <- dmi[, r]
        out[, , r] <- diag(dmr, m, m) -
          tcrossprod(dmr, mi) -
          tcrossprod(mi, dmr)
      }
      return(out)
    }
    
    h <- abs(i - j)
    r_ij <- rho_vec[h + 1L]
    
    if (!is.finite(r_ij) || abs(r_ij) >= 1) {
      return(array(0, dim = c(m, m, p)))
    }
    
    drij_dphi1 <- dr_phi1[h + 1L]
    drij_dphi2 <- dr_phi2[h + 1L]
    
    dE_dmu_i <- matrix(0, m, m)
    dE_dmu_j <- matrix(0, m, m)
    dE_dr    <- matrix(0, m, m)
    dE_dcuts <- array(0, dim = c(m, m, p_cuts))
    
    for (k in 1:m) {
      al <- seg[k]     - mu[i]
      au <- seg[k + 1] - mu[i]
      
      for (l in 1:m) {
        bl <- seg[l]     - mu[j]
        bu <- seg[l + 1] - mu[j]
        
        part <- rect_partials(al, au, bl, bu, r_ij)
        
        dE_dmu_i[k, l] <- part$dmu_i
        dE_dmu_j[k, l] <- part$dmu_j
        dE_dr[k, l]    <- part$dr
        
        if (p_cuts > 0L) {
          for (rcut in 1:p_cuts) {
            ku <- rcut + 1L
            kl <- rcut + 2L
            
            val <- 0
            if (k == ku) val <- val + part$dau
            if (kl <= m && k == kl) val <- val + part$dal
            if (l == ku) val <- val + part$dbu
            if (kl <= m && l == kl) val <- val + part$dbl
            
            dE_dcuts[k, l, rcut] <- val
          }
        }
      }
    }
    
    ## cut parameters
    if (p_cuts > 0L) {
      for (r in 1:p_cuts) {
        out[, , r] <- dE_dcuts[, , r] -
          tcrossprod(dmi[, r], mj) -
          tcrossprod(mi, dmj[, r])
      }
    }
    
    ## theta parameters
    Xi <- DesignX[i, ]
    Xj <- DesignX[j, ]
    
    for (ell in 1:q) {
      idx <- p_cuts + ell
      
      dE_ell <- dE_dmu_i * Xi[ell] + dE_dmu_j * Xj[ell]
      
      out[, , idx] <- dE_ell -
        tcrossprod(dmi[, idx], mj) -
        tcrossprod(mi, dmj[, idx])
    }
    
    ## phi1 and phi2 parameters
    out[, , phi1_idx] <- dE_dr * drij_dphi1
    out[, , phi2_idx] <- dE_dr * drij_dphi2
    
    out
  }
  
  dkappa
}


make_aop2_cls_optim <- function(Xw, DesignX,
                                lag_max = 50L,
                                lag_tol = 1e-8,
                                jitter = 1e-10,
                                symmetrize = TRUE,
                                do_grad = TRUE) {
  
  Xw <- as.matrix(Xw)
  DesignX <- as.matrix(DesignX)
  
  Ts <- nrow(Xw)
  K  <- ncol(Xw)
  q  <- ncol(DesignX)
  m  <- K - 1L
  
  p_cuts <- K - 2L
  
  ## cuts + theta + phi1 + phi2
  p <- p_cuts + q + 2L
  
  stopifnot(nrow(DesignX) == Ts, K >= 3L)
  
  unpack <- function(par) {
    par <- as.numeric(par)
    stopifnot(length(par) == p)
    
    cuts  <- par[1:p_cuts]
    theta <- par[p_cuts + (1:q)]
    phi1  <- par[p_cuts + q + 1L]
    phi2  <- par[p_cuts + q + 2L]
    
    list(cuts = cuts, theta = theta, phi1 = phi1, phi2 = phi2)
  }
  
  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$val <- NULL
  cache$PredX <- NULL
  cache$EXW <- NULL
  cache$u <- NULL
  cache$kappa <- NULL
  cache$upar <- NULL
  cache$grad <- NULL
  cache$lag_max_used <- NULL
  
  eval_value <- function(par) {
    
    par <- as.numeric(par)
    
    if (!is.null(cache$par) &&
        length(par) == length(cache$par) &&
        all(par == cache$par) &&
        !is.null(cache$val)) {
      return(cache$val)
    }
    
    upar  <- unpack(par)
    cuts  <- upar$cuts
    theta <- upar$theta
    phi1  <- upar$phi1
    phi2  <- upar$phi2
    
    if (any(!is.finite(cuts)) ||
        is.unsorted(cuts, strictly = TRUE) ||
        cuts[1] <= 0) {
      cache$par <- par
      cache$val <- BIG
      cache$grad <- NULL
      return(BIG)
    }
    
    if (!is_stationary_ar2(phi1, phi2)) {
      cache$par <- par
      cache$val <- BIG
      cache$grad <- NULL
      return(BIG)
    }
    
    sig2e <- ar2_innov_var_unit(phi1, phi2)
    if (!is.finite(sig2e) || sig2e <= 0) {
      cache$par <- par
      cache$val <- BIG
      cache$grad <- NULL
      return(BIG)
    }
    
    cuts_all <- c(0, cuts)
    seg <- c(-Inf, cuts_all, Inf)
    
    mst <- as.numeric(DesignX %*% theta)
    
    Fmat <- stats::pnorm(outer(mst, seg, function(mu, c) c - mu))
    Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] - Fmat[, 1:K, drop = FALSE]
    EXW  <- Pmat[, 1:m, drop = FALSE]
    u    <- Xw[, 1:m, drop = FALSE] - EXW
    
    acf_deriv <- ar2_acf_deriv(phi1 = phi1, phi2 = phi2, maxlag = Ts - 1L)
    rho_vec <- acf_deriv$rho
    
    kappa <- function(i, j) {
      multikappa_pbivnorm_ar2(
        i, j,
        msti = mst[i],
        mstj = mst[j],
        EXwi = EXW[i, ],
        EXwj = EXW[j, ],
        seg = seg,
        acf_vec = rho_vec
      )
    }
    
    fit <- tryCatch(
      innovations_mv(
        kappa = kappa,
        u = u,
        jitter = jitter,
        symmetrize = symmetrize,
        compute_innov = TRUE,
        lag_max = lag_max,
        lag_tol = lag_tol,
        ahead = NULL
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      cache$par <- par
      cache$val <- BIG
      cache$grad <- NULL
      return(BIG)
    }
    
    PredX <- EXW + fit$uhat
    Predw <- cbind(PredX, 1 - rowSums(PredX))
    
    val <- sum((Xw - Predw)^2)
    if (!is.finite(val)) val <- BIG
    
    cache$par <- par
    cache$val <- val
    cache$PredX <- PredX
    cache$EXW <- EXW
    cache$u <- u
    cache$kappa <- kappa
    cache$upar <- upar
    cache$grad <- NULL
    cache$lag_max_used <- fit$lag_max_used
    
    val
  }
  
  fn <- function(par) eval_value(par)
  
  gr <- if (do_grad) {
    
    function(par) {
      
      par <- as.numeric(par)
      val <- eval_value(par)
      
      if (!is.finite(val) || val >= BIG / 10) {
        return(rep(0, p))
      }
      
      if (!is.null(cache$grad) &&
          !is.null(cache$par) &&
          length(par) == length(cache$par) &&
          all(par == cache$par)) {
        return(cache$grad)
      }
      
      upar <- cache$upar
      
      cuts  <- upar$cuts
      theta <- upar$theta
      phi1  <- upar$phi1
      phi2  <- upar$phi2
      
      du_arr <- make_du_aop2(
        cuts = cuts,
        theta = theta,
        DesignX = DesignX
      )
      
      dkappa <- make_dkappa_aop2(
        cuts = cuts,
        theta = theta,
        phi1 = phi1,
        phi2 = phi2,
        DesignX = DesignX
      )
      
      fitg <- tryCatch(
        innovations_mv_grad(
          kappa = cache$kappa,
          u = cache$u,
          dkappa = dkappa,
          du = du_arr,
          jitter = jitter,
          symmetrize = symmetrize,
          lag_max = cache$lag_max_used,
          lag_tol = NULL
        ),
        error = function(e) NULL
      )
      
      if (is.null(fitg)) {
        cache$grad <- rep(0, p)
        return(cache$grad)
      }
      
      g <- grad_cls_sse(
        Xw_full = Xw,
        PredX = cache$PredX,
        duhat = fitg$duhat,
        du_arr = du_arr
      )
      
      g[!is.finite(g)] <- 0
      cache$grad <- g
      g
    }
    
  } else NULL
  
  list(fn = fn, gr = gr, unpack = unpack)
}


fit_CLSInnov_AR2 <- function(X_hour, X_hour_wide, DesignXEst,
                             lag_max = 50L,
                             lag_tol = 1e-8) {
  
  K <- ncol(X_hour_wide)
  q <- ncol(DesignXEst)
  
  obj <- make_aop2_cls_optim(
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
  
  yy <- as.numeric(scale(X_hour, center = TRUE, scale = FALSE))
  
  phi_initial <- tryCatch({
    arfit <- stats::arima(yy, order = c(2, 0, 0), include.mean = FALSE)
    tmp <- as.numeric(arfit$coef[1:2])
    if (!is_stationary_ar2(tmp[1], tmp[2])) c(0.3, 0.2) else tmp
  }, error = function(e) {
    c(0.3, 0.2)
  })
  
  par_initial <- c(cuts_initial, theta_initial, phi_initial)
  
  con <- make_constraints_aop2(
    K = K,
    q = q,
    delta = 1e-6,
    delta0 = 1e-6,
    eps_ar = 1e-6
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