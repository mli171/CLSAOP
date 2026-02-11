
BIG <- 1e30

innovations_mv <- function(
    kappa,
    u,
    jitter = 1e-10,
    symmetrize = TRUE,
    compute_innov = !is.null(u),
    lag_max = NULL,
    lag_tol = 1e-8,
    ahead = 1L) {
  
  if (!is.function(kappa)) stop("kappa must be a function of (i, j) returning an m x m matrix.")
  
  compute_ahead <- !is.null(ahead)
  if (compute_ahead) {
    ahead <- as.integer(ahead)
    if (ahead < 1L) stop("ahead must be >= 1 or NULL.")
  }
  
  u <- as.matrix(u)
  Tt <- nrow(u)
  m  <- ncol(u)
  max_n <- Tt - 1L
  
  # initial cap before we possibly shrink to tol-based lag
  lag_cap0 <- if (is.null(lag_max)) max_n else {
    lag_max <- as.integer(lag_max)
    if (lag_max < 1L) stop("lag_max must be >= 1 or NULL.")
    min(lag_max, max_n)
  }
  
  do_tol_stop <- !is.null(lag_tol)
  if (do_tol_stop) {
    lag_tol <- as.numeric(lag_tol)
    if (!is.finite(lag_tol) || lag_tol < 0) stop("lag_tol must be finite and >= 0, or NULL.")
  }
  
  getK <- function(i, j) {
    A <- as.matrix(kappa(i, j))
    if (!all(dim(A) == c(m, m))) stop("kappa(i,j) must return an m x m matrix.")
    if (symmetrize && i == j) A <- (A + t(A)) / 2
    A
  }
  
  make_spd <- function(A) {
    if (symmetrize) A <- (A + t(A))/2
    A + jitter * diag(m)
  }
  
  Theta <- vector("list", Tt)
  V     <- vector("list", Tt)
  Vinv  <- vector("list", Tt)
  
  Theta[[1]] <- list()
  V[[1]]     <- make_spd(getK(1, 1))
  Vinv[[1]]  <- chol2inv(chol(V[[1]]))
  
  # Lcurr is the max lag length we WILL compute/store at each t
  Lcurr <- lag_cap0
  lag_fixed <- FALSE
  lag_used_path <- integer(Tt)  # optional diagnostic: lag used at each t
  lag_used_path[1] <- 0L
  
  if (max_n >= 1L && Lcurr >= 1L) {
    for (n in 1:max_n) {
      # time index t = n+1
      Lt <- min(n, Lcurr)           # number of lags actually computed at time t
      Theta_n <- vector("list", Lt) # store only lags 1..Lt
      
      k_start <- n - Lt             # 0-based k start so that lag = n-k <= Lt
      
      for (k in k_start:(n - 1)) {
        ell <- n - k                # lag index in 1..Lt
        term <- getK(n + 1, k + 1)
        
        if (k >= 1) {
          # j must keep both:
          #   (1) ell_tj = n - j <= Lt  => j >= n - Lt
          #   (2) ell_kj = k - j <= L_k (<= Lcurr once fixed/capped) => j >= k - Lcurr
          # We can safely use Lt here because:
          # - before lag_fixed, Theta[[k+1]] length is min(k, Lcurr)
          # - after lag_fixed, Lcurr is the fixed lag
          j0 <- max(0L, n - Lt, k - Lt)
          if (j0 <= k - 1) {
            for (j in j0:(k - 1)) {
              ell_tj <- n - j        # in 1..Lt
              ell_kj <- k - j        # must exist in Theta[[k+1]]
              term <- term - Theta_n[[ell_tj]] %*% V[[j + 1]] %*% t(Theta[[k + 1]][[ell_kj]])
            }
          }
        }
        
        Theta_n[[ell]] <- term %*% Vinv[[k + 1]]
      }
      
      # V at time t=n+1 uses only the retained lags
      Vn <- getK(n + 1, n + 1)
      j0v <- max(0L, n - Lt)
      if (j0v <= n - 1) {
        for (j in j0v:(n - 1)) {
          ell <- n - j
          A  <- Theta_n[[ell]]
          Vn <- Vn - A %*% V[[j + 1]] %*% t(A)
        }
      }
      Vn <- make_spd(Vn)
      
      Theta[[n + 1]] <- Theta_n
      V[[n + 1]]     <- Vn
      Vinv[[n + 1]]  <- chol2inv(chol(Vn))
      lag_used_path[n + 1] <- Lt
      
      # tol-based "freeze lag length" check:
      # only meaningful when we computed the FULL n lags (i.e., Lt == n)
      if (!lag_fixed && do_tol_stop && Lt == n) {
        if (max(abs(Theta_n[[n]])) < lag_tol) {
          Lcurr <- n
          lag_fixed <- TRUE
        }
      }
    }
  }
  
  lag_max_used <- Lcurr
  
  # innovations / prediction
  uhat <- innov <- NULL
  if (compute_innov) {
    uhat  <- matrix(0, nrow = Tt, ncol = m)
    innov <- matrix(0, nrow = Tt, ncol = m)
    
    for (t in 1:Tt) {
      if (t == 1) {
        uhat[t, ] <- 0
      } else {
        coefs <- Theta[[t]]
        jmax  <- min(length(coefs), t - 1L)
        pred  <- rep(0, m)
        if (jmax >= 1L) {
          for (j in 1:jmax) pred <- pred + coefs[[j]] %*% innov[t - j, ]
        }
        uhat[t, ] <- pred
      }
      innov[t, ] <- u[t, ] - uhat[t, ]
    }
  }
  
  uhat_ahead <- NULL
  if (compute_ahead && compute_innov) {
    if (ahead == 1L) {
      uhat_ahead <- uhat
    } else {
      uhat_ahead <- matrix(0, nrow = Tt, ncol = m)
      for (t in 1:Tt) {
        if (t <= ahead) {
          uhat_ahead[t, ] <- 0
        } else {
          origin <- t - ahead
          coefs  <- Theta[[origin + 1L]]
          jmax   <- min(length(coefs), t - 1L)
          pred   <- rep(0, m)
          if (jmax >= ahead) {
            for (j in ahead:jmax) pred <- pred + coefs[[j]] %*% innov[t - j, ]
          }
          uhat_ahead[t, ] <- pred
        }
      }
    }
  }
  
  list(
    Theta = Theta,
    V = V,
    uhat = uhat,
    innov = innov,
    uhat_ahead = uhat_ahead,
    lag_max_used = lag_max_used,
    lag_used_path = lag_used_path
  )
}


pbiv_safe <- function(x, y, r) {
  x <- as.numeric(x); y <- as.numeric(y)
  out <- numeric(length(x))
  out[x == -Inf | y == -Inf] <- 0
  out[x == Inf & y == Inf] <- 1
  idx <- (x == Inf & is.finite(y))
  if (any(idx)) out[idx] <- stats::pnorm(y[idx])
  idx <- (y == Inf & is.finite(x))
  if (any(idx)) out[idx] <- stats::pnorm(x[idx])
  fin <- is.finite(x) & is.finite(y)
  if (any(fin)) out[fin] <- pbivnorm::pbivnorm(x[fin], y[fin], r)
  out
}

multikappa_pbivnorm <- function(i, j, msti, mstj, EXwi, EXwj, seg, rho) {
  
  msti <- as.numeric(msti); mstj <- as.numeric(mstj)
  rho  <- as.numeric(rho)
  seg  <- as.numeric(seg)
  EXwi <- as.numeric(EXwi)
  EXwj <- as.numeric(EXwj)
  
  K <- length(seg) - 1L
  m <- K - 1L
  
  if (i == j) {
    p <- EXwi
    return(diag(p, m, m) - tcrossprod(p))
  }
  
  r <- rho^abs(i - j)
  
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

CLSWInnov <- function(par, Xw, DesignX,
                      lag_max = 50L, lag_tol = 1e-6, jitter = 1e-10) {
  
  Xw <- as.matrix(Xw)
  DesignX <- as.matrix(DesignX)
  
  Ts <- nrow(Xw)
  K  <- ncol(Xw)
  q  <- ncol(DesignX)
  m  <- K - 1L
  
  cuts  <- as.numeric(par[1:(K - 2L)])                 # c2..c_{K-1}
  theta <- as.numeric(par[(K - 1L):(K - 2L + q)])
  rho   <- as.numeric(par[K - 2L + q + 1L])
  
  ## feasibility checks
  if (!is.finite(rho) || abs(rho) >= 1) return(BIG)
  if (any(!is.finite(cuts)) || is.unsorted(cuts, strictly = TRUE)) return(BIG)
  if (cuts[1] <= 0) return(BIG)   # replicate exp/cumsum implication: c2>0
  
  cuts_all <- c(0, cuts)
  seg <- c(-Inf, cuts_all, Inf)
  
  mst <- as.numeric(DesignX %*% theta)
  Fmat <- stats::pnorm(outer(mst, seg, function(mu, c) c - mu))
  Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] - Fmat[, 1:K, drop = FALSE]
  EXW  <- Pmat[, 1:m, drop = FALSE]
  u    <- Xw[, 1:m, drop = FALSE] - EXW
  
  kappa <- function(i, j) {
    multikappa_pbivnorm(i, j,
                        msti = mst[i], mstj = mst[j],
                        EXwi = EXW[i, ], EXwj = EXW[j, ],
                        seg  = seg, rho = rho)
  }
  
  fit <- tryCatch(
    innovations_mv(kappa = kappa, u = u,
                   jitter = jitter, symmetrize = TRUE,
                   compute_innov = TRUE,
                   lag_max = lag_max, lag_tol = lag_tol,
                   ahead = NULL),
    error = function(e) NULL
  )
  if (is.null(fit)) return(BIG)
  
  PredX <- EXW + fit$uhat
  Predw <- cbind(PredX, 1 - rowSums(PredX))
  
  val <- sum((Xw - Predw)^2)
  if (!is.finite(val)) val <- BIG
  val
}

unpack_CLSWInnov_par <- function(par, K, q) {
  stopifnot(length(par) == (K - 2L + q + 1L))
  cuts  <- as.numeric(par[1:(K - 2L)])
  theta <- as.numeric(par[(K - 1L):(K - 2L + q)])
  rho   <- as.numeric(par[K - 2L + q + 1L])
  list(cuts = cuts, theta = theta, rho = rho)
}

make_du_aop1 <- function(cuts, theta, DesignX) {
  
  DesignX <- as.matrix(DesignX)
  Ts <- nrow(DesignX)
  q  <- ncol(DesignX)
  
  cuts  <- as.numeric(cuts)
  theta <- as.numeric(theta)
  
  p_cuts <- length(cuts)
  K <- p_cuts + 2L
  m <- K - 1L
  p <- p_cuts + q + 1L
  
  cuts_all <- c(0, cuts)
  seg <- c(-Inf, cuts_all, Inf)
  
  mu <- as.numeric(DesignX %*% theta)
  
  pdf_seg <- dnorm(outer(mu, seg,      function(a,b) b - a))      # Ts x (K+1) for all boundary points
  pdf_cut <- dnorm(outer(mu, cuts_all, function(a,b) b - a))      # Ts x (K-1) for all cutpoints parameters
  
  # dEXW_dmu=phi(seg_k-mu)-phi(seg_{k+1}-mu) (Ts x m)
  dm_dmu <- pdf_seg[, 1:m, drop=FALSE] - pdf_seg[, 2:(m+1), drop=FALSE]
  
  du_arr <- array(0, dim = c(Ts, m, p))
  
  # cutpoint parameters (affecting adjacent categories up to category m=K-1)
  if (p_cuts > 0L) {
    for (t in 1:Ts) {
      # dm_dc[k,r] = dEXW_dcuts
      dm_dc <- matrix(0, nrow = m, ncol = p_cuts)
      for (r in 1:p_cuts) {
        j <- r + 1L
        val <- pdf_cut[t, j]
        if (j <= m) dm_dc[j, r] <- dm_dc[j, r] + val
        if ((j + 1L) <= m) dm_dc[j + 1L, r] <- dm_dc[j + 1L, r] - val
      }
      du_arr[t, , 1:p_cuts] <- -dm_dc
    }
  }
  
  # theta parameters: dEXW_dtheta = dEXW_dmu*Xt
  theta_idx <- p_cuts + (1:q)
  for (t in 1:Ts) {
    dm_dtheta <- outer(dm_dmu[t, ], DesignX[t, ])
    du_arr[t, , theta_idx] <- -dm_dtheta
  }
  
  # rho block is zero (ut only depends on marginal probabilities, not on rho)
  du_arr
}


make_dkappa_aop1 <- function(cuts, theta, rho, DesignX, eps_s = 1e-12) {
  
  DesignX <- as.matrix(DesignX)
  Ts <- nrow(DesignX)
  q  <- ncol(DesignX)
  
  cuts  <- as.numeric(cuts)
  theta <- as.numeric(theta)
  rho   <- as.numeric(rho)
  
  p_cuts <- length(cuts)
  K <- p_cuts + 2L
  m <- K - 1L
  p <- p_cuts + q + 1L # cuts + theta + rho
  rho_idx <- p
  
  cuts_all <- c(0, cuts)
  seg <- c(-Inf, cuts_all, Inf)
  
  mu <- as.numeric(DesignX %*% theta)
  
  # Marginal probabilities EXW
  Fmat <- pnorm(outer(mu, seg, function(a,b) b - a))
  Pmat <- Fmat[, 2:(K+1), drop=FALSE] - Fmat[, 1:K, drop=FALSE]
  EXW  <- Pmat[, 1:m, drop=FALSE]
  
  # dm = dEXW_dpar
  dm <- array(0, dim = c(Ts, m, p))
  
  pdf_seg <- dnorm(outer(mu, seg,      function(a,b) b - a))
  pdf_cut <- dnorm(outer(mu, cuts_all, function(a,b) b - a))
  # dEXW_dmu=phi(seg_k-mu)-phi(seg_{k+1}-mu) (Ts x m)
  dm_dmu  <- pdf_seg[, 1:m, drop=FALSE] - pdf_seg[, 2:(m+1), drop=FALSE]
  
  # dmu_dtheta = Xt
  theta_idx <- p_cuts + (1:q)
  for (t in 1:Ts) dm[t, , theta_idx] <- outer(dm_dmu[t, ], DesignX[t, ])
  
  # cutpoint parameters (affecting adjacent categories up to category m=K-1)
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
  
  # BVN derivative helpers for rectangle probability derivatives
  
  # numerical safer
  Phi2_safe <- function(a, b, r) {
    if (is.infinite(a) && a < 0) return(0)
    if (is.infinite(b) && b < 0) return(0)
    if (is.infinite(a) && a > 0 && is.infinite(b) && b > 0) return(1)
    if (is.infinite(a) && a > 0) return(pnorm(b))
    if (is.infinite(b) && b > 0) return(pnorm(a))
    pbivnorm::pbivnorm(a, b, r)
  }
  # numerical safer
  s_safe <- function(r) {
    s2 <- 1 - r^2
    if (!is.finite(s2) || s2 <= 0) s2 <- eps_s
    sqrt(s2)
  }
  
  # dPhi2_da = phi(a) * Phi((b - r*a)/sqrt(1-r^2))
  dPhi2_da <- function(a, b, r) {
    if (!is.finite(a)) return(0)
    if (is.infinite(b) && b < 0) return(0)
    s <- s_safe(r)
    dnorm(a) * pnorm((b - r*a)/s)
  }
  # dPhi2_db = phi(b) * Phi((a - r*b)/sqrt(1-r^2))
  dPhi2_db <- function(a, b, r) {
    if (!is.finite(b)) return(0)
    if (is.infinite(a) && a < 0) return(0)
    s <- s_safe(r)
    dnorm(b) * pnorm((a - r*b)/s)
  }
  # dPhi2_dr = phi2(a,b,r)
  phi2_pdf <- function(a, b, r) {
    if (!is.finite(a) || !is.finite(b)) return(0)
    s2 <- 1 - r^2
    if (!is.finite(s2) || s2 <= 0) s2 <- eps_s
    z <- (a^2 - 2*r*a*b + b^2) / s2
    (1/(2*pi*sqrt(s2))) * exp(-0.5*z)
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
    
    list(P=P, dmu_i=dmu_i, dmu_j=dmu_j, dr=dr,
         dau=dau, dal=dal, dbu=dbu, dbl=dbl)
  }
  
  # dkappa_dpar: m x m x p derivative array
  dkappa <- function(i, j) {
    i <- as.integer(i); j <- as.integer(j)
    
    # symmetric
    if (i < j) {
      A <- dkappa(j, i)
      return(aperm(A, c(2, 1, 3)))
    }
    
    mi  <- as.numeric(EXW[i, ])
    mj  <- as.numeric(EXW[j, ])
    dmi <- matrix(dm[i, , ], nrow = m, ncol = p)
    dmj <- matrix(dm[j, , ], nrow = m, ncol = p)
    
    out <- array(0, dim = c(m, m, p))
    
    # kappa(i,i)
    if (i == j) {
      for (r in 1:p) {
        dmr <- dmi[, r]
        out[,,r] <- diag(dmr, m, m) - tcrossprod(dmr, mi) - tcrossprod(mi, dmr)
      }
      return(out)
    }
    
    # kappa(i,j) where i != j
    h <- abs(i - j)
    r_ij <- rho^h
    drij_drho <- if (h == 0L) 0 else (h * rho^(h - 1L))
    
    dE_dmu_i <- matrix(0, m, m)
    dE_dmu_j <- matrix(0, m, m)
    dE_dr    <- matrix(0, m, m)
    dE_dcuts <- array(0, dim = c(m, m, p_cuts))
    
    for (k in 1:m) {
      al <- seg[k]   - mu[i]
      au <- seg[k+1] - mu[i]
      for (l in 1:m) {
        bl <- seg[l]   - mu[j]
        bu <- seg[l+1] - mu[j]
        
        part <- rect_partials(al, au, bl, bu, r_ij)
        
        dE_dmu_i[k, l] <- part$dmu_i
        dE_dmu_j[k, l] <- part$dmu_j
        dE_dr[k, l]    <- part$dr
        
        ## cut cuts[r] = c_{r+1} is boundary between categories ku=r+1 and kl=r+2.
        ##   Zi upper bound (au) when k==ku   -> contributes to dau
        ##   Zi lower bound (al) when k==kl   -> contributes to dal
        ##   Zj upper bound (bu) when l==ku   -> contributes to dbu
        ##   Zj lower bound (bl) when l==kl   -> contributes to dbl
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
    
    # cut parameters
    if (p_cuts > 0L) {
      for (r in 1:p_cuts) {
        out[, , r] <- dE_dcuts[, , r] -
          tcrossprod(dmi[, r], mj) - tcrossprod(mi, dmj[, r])
      }
    }
    
    # theta parameters
    Xi <- DesignX[i, ]
    Xj <- DesignX[j, ]
    for (ell in 1:q) {
      idx <- p_cuts + ell
      dE_ell <- dE_dmu_i * Xi[ell] + dE_dmu_j * Xj[ell]
      out[, , idx] <- dE_ell -
        tcrossprod(dmi[, idx], mj) - tcrossprod(mi, dmj[, idx])
    }
    
    # rho parameters
    out[, , rho_idx] <- dE_dr * drij_drho
    
    out
  }
  
  dkappa
}


innovations_mv_grad <- function(
    kappa,
    u,
    dkappa,
    du,  # du[t, , r] = d u_t / d par[r]
    jitter = 1e-10,
    symmetrize = TRUE,
    lag_max = NULL,
    lag_tol = 1e-8) {
  
  if (!is.function(kappa))  stop("kappa must be a function(i,j)->m x m.")
  if (!is.function(dkappa)) stop("dkappa must be a function(i,j)->m x m x p.")
  
  u <- as.matrix(u)
  Tt <- nrow(u)
  m  <- ncol(u)
  max_n <- Tt - 1L
  
  # initial cap (memory/compute cap)
  lag_cap0 <- if (is.null(lag_max)) max_n else {
    lag_max <- as.integer(lag_max)
    if (lag_max < 1L) stop("lag_max must be >= 1 or NULL.")
    min(lag_max, max_n)
  }
  
  do_tol_stop <- !is.null(lag_tol)
  if (do_tol_stop) {
    lag_tol <- as.numeric(lag_tol)
    if (!is.finite(lag_tol) || lag_tol < 0) stop("lag_tol must be finite and >= 0, or NULL.")
  }
  
  tmp <- as.array(dkappa(1L, 1L))
  if (length(dim(tmp)) != 3L || any(dim(tmp)[1:2] != c(m, m)))
    stop("dkappa(1,1) must return array dim (m,m,p).")
  p <- dim(tmp)[3L]
  
  du_arr <- as.array(du)
  if (length(dim(du_arr)) != 3L || any(dim(du_arr) != c(Tt, m, p)))
    stop("du must be array dim (T,m,p).")
  
  getK <- function(i, j) {
    A <- as.matrix(kappa(i, j))
    if (symmetrize && i == j) A <- (A + t(A)) / 2
    A
  }
  getdK <- function(i, j) {
    A <- as.array(dkappa(i, j))
    if (symmetrize && i == j) A <- 0.5 * (A + aperm(A, c(2,1,3)))
    A
  }
  getdu <- function(t) matrix(du_arr[t, , ], nrow = m, ncol = p)
  
  make_spd <- function(A) {
    if (symmetrize) A <- (A + t(A))/2
    A + jitter * diag(m)
  }
  
  # storage
  Theta  <- vector("list", Tt)
  V      <- vector("list", Tt)
  Vinv   <- vector("list", Tt)
  dTheta <- vector("list", Tt)
  dV     <- vector("list", Tt)
  
  lag_used_path <- integer(Tt)
  lag_used_path[1] <- 0L
  
  Theta[[1]]  <- list()
  dTheta[[1]] <- list()
  
  V[[1]]    <- make_spd(getK(1, 1))
  Vinv[[1]] <- chol2inv(chol(V[[1]]))
  dV[[1]]   <- getdK(1, 1)  # no jitter derivative
  
  # adaptive lag length
  Lcurr <- lag_cap0
  lag_fixed <- FALSE
  
  if (max_n >= 1L && Lcurr >= 1L) {
    for (n in 1:max_n) {
      # time t = n+1
      Lt <- min(n, Lcurr)
      lag_used_path[n + 1L] <- Lt
      
      Theta_n  <- vector("list", Lt)
      dTheta_n <- vector("list", Lt)
      
      k_start <- n - Lt  # keep k in {n-Lt, ..., n-1}
      
      for (k in k_start:(n - 1L)) {
        ell <- n - k  # lag index in 1..Lt
        
        term  <- getK(n + 1L, k + 1L)
        dterm <- getdK(n + 1L, k + 1L)
        
        if (k >= 1L) {
          Lk <- length(Theta[[k + 1L]])          # lags stored at time k+1
          j0 <- max(0L, n - Lt, k - Lk)          # ensure both lag indices exist
          if (j0 <= (k - 1L)) {
            for (j in j0:(k - 1L)) {
              ell_tj <- n - j                    # in 1..Lt
              ell_kj <- k - j                    # in 1..Lk
              
              A_n <- Theta_n[[ell_tj]]
              A_k <- Theta[[k + 1L]][[ell_kj]]
              term <- term - A_n %*% V[[j + 1L]] %*% t(A_k)
              
              dA_n <- dTheta_n[[ell_tj]]                 # m x m x p
              dA_k <- dTheta[[k + 1L]][[ell_kj]]         # m x m x p
              dV_j <- dV[[j + 1L]]                       # m x m x p
              
              for (r in 1:p) {
                dterm[,,r] <- dterm[,,r] -
                  (dA_n[,,r] %*% V[[j + 1L]] %*% t(A_k)) -
                  (A_n %*% dV_j[,,r] %*% t(A_k)) -
                  (A_n %*% V[[j + 1L]] %*% t(dA_k[,,r]))
              }
            }
          }
        }
        
        # Theta = term * V_k^{-1}
        Theta_mat <- term %*% Vinv[[k + 1L]]
        Theta_n[[ell]] <- Theta_mat
        
        dTheta_mat <- array(0, dim = c(m, m, p))
        dV_k <- dV[[k + 1L]]
        for (r in 1:p) {
          # d( term V^{-1}) = (dterm - Theta dV) V^{-1}
          dTheta_mat[,,r] <- (dterm[,,r] - Theta_mat %*% dV_k[,,r]) %*% Vinv[[k + 1L]]
        }
        dTheta_n[[ell]] <- dTheta_mat
      }
      
      # Vn uses only retained lags j in {n-Lt,...,n-1}
      Vn  <- getK(n + 1L, n + 1L)
      dVn <- getdK(n + 1L, n + 1L)
      
      j0v <- max(0L, n - Lt)
      if (j0v <= (n - 1L)) {
        for (j in j0v:(n - 1L)) {
          ell <- n - j
          A  <- Theta_n[[ell]]
          Vn <- Vn - A %*% V[[j + 1L]] %*% t(A)
          
          dA   <- dTheta_n[[ell]]
          dV_j <- dV[[j + 1L]]
          for (r in 1:p) {
            dVn[,,r] <- dVn[,,r] -
              (dA[,,r] %*% V[[j + 1L]] %*% t(A)) -
              (A %*% dV_j[,,r] %*% t(A)) -
              (A %*% V[[j + 1L]] %*% t(dA[,,r]))
          }
        }
      }
      
      Vn <- make_spd(Vn)
      
      Theta[[n + 1L]]  <- Theta_n
      dTheta[[n + 1L]] <- dTheta_n
      
      V[[n + 1L]]    <- Vn
      Vinv[[n + 1L]] <- chol2inv(chol(Vn))
      dV[[n + 1L]]   <- dVn
      
      # tol rule: once the boundary lag (oldest retained) is small, fix Lcurr forever
      if (!lag_fixed && do_tol_stop && Lt >= 1L) {
        if (max(abs(Theta_n[[Lt]])) < lag_tol) {
          Lcurr <- Lt
          lag_fixed <- TRUE
        }
      }
    }
  }
  
  # now forward recursion for uhat, innov and their derivatives (A.18)-(A.19)
  uhat   <- matrix(0, nrow = Tt, ncol = m)
  innov  <- matrix(0, nrow = Tt, ncol = m)
  duhat  <- array(0, dim = c(Tt, m, p))
  dinnov <- array(0, dim = c(Tt, m, p))
  
  uhat[1, ] <- 0
  innov[1, ] <- u[1, ]
  dinnov[1, , ] <- getdu(1L)
  
  if (Tt >= 2L) {
    for (t in 2:Tt) {
      coefs <- Theta[[t]]
      jmax  <- min(length(coefs), t - 1L)
      
      pred <- rep(0, m)
      if (jmax >= 1L) {
        for (j in 1:jmax) pred <- pred + coefs[[j]] %*% innov[t - j, ]
      }
      uhat[t, ]  <- pred
      innov[t, ] <- u[t, ] - pred
      
      dpred <- matrix(0, nrow = m, ncol = p)
      du_t  <- getdu(t)
      
      if (jmax >= 1L) {
        for (j in 1:jmax) {
          H   <- coefs[[j]]
          dH  <- dTheta[[t]][[j]]      # m x m x p
          e_l <- innov[t - j, ]
          de  <- dinnov[t - j, , ]     # m x p
          for (r in 1:p) {
            dpred[, r] <- dpred[, r] + (dH[,,r] %*% e_l) + (H %*% de[, r])
          }
        }
      }
      
      duhat[t, , ]  <- dpred
      dinnov[t, , ] <- du_t - dpred
    }
  }
  
  list(
    uhat = uhat,
    innov = innov,
    duhat = duhat,
    dinnov = dinnov,
    lag_max_used = Lcurr,
    lag_used_path = lag_used_path
  )
}




grad_cls_sse <- function(Xw_full, PredX, duhat, du_arr) {
  
  # Dimensions:
  #   Xw_full : T x K one-hot matrix W_t
  #   PredX   : T x m predicted probabilities for categories 1..m, where m=K-1
  #   duhat   : T x m x p, duhat[t,k,r] = duhat_dpar_r{t,k}
  #   du_arr  : T x m x p, du_arr[t,k,r] = du_par_r{t,k}
  #
  # Recall:
  #   u_t = X_t - EXW_t
  #   du_t = - d(EXW_t)
  #   PredX_t = EXW_t + uhat_t
  #   dPredX_t = d(EXW_t) + d(uhat_t) = (-du_t) + duhat_t
  
  Ts <- nrow(PredX)
  m  <- ncol(PredX)
  K  <- ncol(Xw_full)
  stopifnot(K == m + 1L)
  
  p <- dim(duhat)[3]
  g <- numeric(p)
  
  residX <- Xw_full[, 1:m, drop = FALSE] - PredX
  residL <- Xw_full[, K] - (1 - rowSums(PredX))
  
  for (t in 1:Ts) {
    dPredX_t <- -du_arr[t, , ] + duhat[t, , ]
    w <- residX[t, ] - residL[t]
    # Per-time gradient contribution: s_t = -2 * (dPredX_t)' * w_t
    g <- g - 2 * as.numeric(crossprod(dPredX_t, w))
  }
  
  g
}

# create constraints for constrOptim usage (no need if use optim)
make_constraints_aop1 <- function(K, q, delta = 1e-6, delta0 = 1e-6, eps_rho = 1e-6) {
  p_cuts <- K - 2L
  p <- p_cuts + q + 1L
  rho_idx <- p
  
  ui <- NULL
  ci <- NULL
  
  ## c2 >= delta0
  if (p_cuts >= 1L) {
    row <- rep(0, p); row[1] <- 1
    ui <- rbind(ui, row); ci <- c(ci, delta0)
  }
  
  ## c_{j+1} - c_j >= delta
  if (p_cuts >= 2L) {
    for (j in 1:(p_cuts - 1L)) {
      row <- rep(0, p)
      row[j + 1L] <-  1
      row[j]      <- -1
      ui <- rbind(ui, row)
      ci <- c(ci, delta)
    }
  }
  
  ## rho <= 1-eps  -> -rho >= -(1-eps)
  row <- rep(0, p); row[rho_idx] <- -1
  ui <- rbind(ui, row); ci <- c(ci, -(1 - eps_rho))
  
  ## rho >= -1+eps
  row <- rep(0, p); row[rho_idx] <-  1
  ui <- rbind(ui, row); ci <- c(ci, -1 + eps_rho)
  
  list(ui = ui, ci = ci)
}

make_aop1_cls_optim <- function(Xw, DesignX,
                                lag_max = 50L,
                                lag_tol = 1e-8,
                                jitter  = 1e-10,
                                symmetrize = TRUE,
                                do_grad = TRUE) {
  
  Xw <- as.matrix(Xw)
  DesignX <- as.matrix(DesignX)
  
  Ts <- nrow(Xw)
  K  <- ncol(Xw)
  q  <- ncol(DesignX)
  m  <- K - 1L
  p_cuts <- K - 2L
  p  <- p_cuts + q + 1L
  
  stopifnot(nrow(DesignX) == Ts, K >= 3L)
  
  # find cut, theta, rho given pre-defined order
  unpack <- function(par) {
    par <- as.numeric(par)
    stopifnot(length(par) == p)
    cuts  <- par[1:p_cuts]
    theta <- par[p_cuts + (1:q)]
    rho   <- par[p]
    list(cuts = cuts, theta = theta, rho = rho)
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
  cache$seg <- NULL
  cache$mst <- NULL
  
  eval_value <- function(par) {
    par <- as.numeric(par)
    
    if (!is.null(cache$par) && length(par) == length(cache$par) &&
        all(par == cache$par) && !is.null(cache$val)) {
      return(cache$val)
    }
    
    upar  <- unpack(par)
    cuts  <- upar$cuts
    theta <- upar$theta
    rho   <- upar$rho
    
    if (!is.finite(rho) || abs(rho) >= 1) {
      cache$par <- par; cache$val <- BIG; cache$grad <- NULL
      return(BIG)
    }
    if (any(!is.finite(cuts)) || is.unsorted(cuts, strictly = TRUE) || cuts[1] <= 0) {
      cache$par <- par; cache$val <- BIG; cache$grad <- NULL
      return(BIG)
    }
    
    cuts_all <- c(0, cuts)
    seg <- c(-Inf, cuts_all, Inf)
    
    mst <- as.numeric(DesignX %*% theta)
    Fmat <- stats::pnorm(outer(mst, seg, function(mu, c) c - mu))
    Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] - Fmat[, 1:K, drop = FALSE]
    EXW  <- Pmat[, 1:m, drop = FALSE]
    u    <- Xw[, 1:m, drop = FALSE] - EXW
    
    kappa <- function(i, j) {
      multikappa_pbivnorm(i, j,
                          msti = mst[i], mstj = mst[j],
                          EXwi = EXW[i, ], EXwj = EXW[j, ],
                          seg  = seg, rho = rho)
    }
    
    fit <- tryCatch(
      innovations_mv(kappa = kappa, u = u,
                     jitter = jitter, symmetrize = symmetrize,
                     compute_innov = TRUE,
                     lag_max = lag_max, lag_tol = lag_tol,
                     ahead = NULL),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      cache$par <- par; cache$val <- BIG; cache$grad <- NULL
      return(BIG)
    }
    cache$lag_max_used <- fit$lag_max_used
    
    PredX <- EXW + fit$uhat
    Predw <- cbind(PredX, 1 - rowSums(PredX))
    
    val <- sum((Xw - Predw)^2)
    if (!is.finite(val)) val <- BIG
    
    cache$par   <- par
    cache$val   <- val
    cache$PredX <- PredX
    cache$EXW   <- EXW
    cache$u     <- u
    cache$kappa <- kappa
    cache$upar  <- upar
    cache$grad  <- NULL
    cache$seg   <- seg
    cache$mst   <- mst
    
    val
  }
  
  fn <- function(par) eval_value(par)
  
  gr <- if (do_grad) {
    function(par) {
      par <- as.numeric(par)
      val <- eval_value(par)
      if (!is.finite(val) || val >= BIG/10) return(rep(0, p))
      
      if (!is.null(cache$grad) && !is.null(cache$par) &&
          length(par) == length(cache$par) && all(par == cache$par)) {
        return(cache$grad)
      }
      
      upar <- cache$upar
      cuts  <- upar$cuts
      theta <- upar$theta
      rho   <- upar$rho
      
      du_arr <- make_du_aop1(cuts = cuts, theta = theta, DesignX = DesignX)
      dkappa <- make_dkappa_aop1(cuts = cuts, theta = theta, rho = rho, DesignX = DesignX)
      
      fitg <- tryCatch(
        innovations_mv_grad(kappa = cache$kappa,
                            u = cache$u,
                            dkappa = dkappa,
                            du = du_arr,
                            jitter = jitter,
                            symmetrize = symmetrize,
                            lag_max = cache$lag_max_used,
                            lag_tol = NULL),
        error = function(e) NULL
      )
      if (is.null(fitg)) {
        cache$grad <- rep(0, p)
        return(cache$grad)
      }
      
      g <- grad_cls_sse(Xw_full = Xw,
                        PredX   = cache$PredX,
                        duhat   = fitg$duhat,
                        du_arr  = du_arr)
      
      g[!is.finite(g)] <- 0
      cache$grad <- g
      g
    }
  } else NULL
  
  list(fn = fn, gr = gr, unpack = unpack)
}

# per-time gradient contributions for sandwitch variance estimator
score_cls_se <- function(Xw_full, PredX, duhat, du_arr) {
  Ts <- nrow(PredX)
  m  <- ncol(PredX)
  K  <- ncol(Xw_full)
  stopifnot(K == m + 1L)
  
  p <- dim(duhat)[3]
  S <- matrix(0, nrow = Ts, ncol = p)
  
  residX <- Xw_full[, 1:m, drop = FALSE] - PredX
  residL <- Xw_full[, K] - (1 - rowSums(PredX))
  
  for (t in 1:Ts) {
    dPredX_t <- -du_arr[t, , ] + duhat[t, , ]
    w <- residX[t, ] - residL[t]
    S[t, ] <- -2 * as.numeric(crossprod(dPredX_t, w))
  }
  S
}


# re-run everthing to calculate Score at one parameter set
aop1_score_matrix <- function(par, Xw, DesignX,
                              lag_max = 50L, lag_tol = 1e-8,
                              jitter = 1e-10, symmetrize = TRUE) {
  
  Xw <- as.matrix(Xw)
  DesignX <- as.matrix(DesignX)
  
  Ts <- nrow(Xw)
  K  <- ncol(Xw)
  q  <- ncol(DesignX)
  m  <- K - 1L
  p_cuts <- K - 2L
  p <- p_cuts + q + 1L
  
  cuts  <- as.numeric(par[1:p_cuts])
  theta <- as.numeric(par[p_cuts + (1:q)])
  rho   <- as.numeric(par[p])
  
  cuts_all <- c(0, cuts)
  seg <- c(-Inf, cuts_all, Inf)
  
  mst <- as.numeric(DesignX %*% theta)
  Fmat <- stats::pnorm(outer(mst, seg, function(mu, c) c - mu))
  Pmat <- Fmat[, 2:(K + 1L), drop = FALSE] - Fmat[, 1:K, drop = FALSE]
  EXW  <- Pmat[, 1:m, drop = FALSE]
  u    <- Xw[, 1:m, drop = FALSE] - EXW
  
  kappa_fun <- function(i, j) {
    multikappa_pbivnorm(i, j,
                        msti = mst[i], mstj = mst[j],
                        EXwi = EXW[i, ], EXwj = EXW[j, ],
                        seg  = seg, rho = rho)
  }
  
  # 1) run objective innovations recursion to get uhat and the chosen lag
  fit <- innovations_mv(kappa = kappa_fun, u = u,
                        jitter = jitter, symmetrize = symmetrize,
                        compute_innov = TRUE,
                        lag_max = lag_max, lag_tol = lag_tol,
                        ahead = NULL)
  
  PredX <- EXW + fit$uhat
  
  # 2) build derivatives at this parameter
  du_arr     <- make_du_aop1(cuts = cuts, theta = theta, DesignX = DesignX)
  dkappa_fun <- make_dkappa_aop1(cuts = cuts, theta = theta, rho = rho, DesignX = DesignX)
  force(dkappa_fun)
  
  # 3) run differentiated recursion with EXACT same lag length as objective
  fitg <- innovations_mv_grad(kappa = kappa_fun,
                              u = u,
                              dkappa = dkappa_fun,
                              du = du_arr,
                              jitter = jitter,
                              symmetrize = symmetrize,
                              lag_max = fit$lag_max_used,
                              lag_tol = NULL)   # IMPORTANT: do NOT re-trigger tol logic
  
  score_cls_se(Xw_full = Xw, PredX = PredX, duhat = fitg$duhat, du_arr = du_arr)
}


# sandwich SE from Hessian (from fit) and score matrix (from calculation)
#      1) compute score matrix S (T x p), rows s_t
#      2) M = sum s_t s_t' = S'%*% S
#      3) H^{-1}
#      4) vcov = H^{-1} M H^{-1}
#      5) se = sqrt(diag(vcov))
sandwichSE <- function(obj, fit, Xw, DesignX,
                       lag_max=50L, lag_tol=1e-8,
                       jitter=1e-10, symmetrize=TRUE) {
  
  H <- fit$hessian
  stopifnot(!is.null(H))
  
  par_hat <- fit$par
  # Score contributions
  S <- aop1_score_matrix(par_hat, Xw, DesignX,
                         lag_max=lag_max, lag_tol=lag_tol,
                         jitter=jitter, symmetrize=symmetrize)
  M <- crossprod(S) # sum s_t s_t'
  # maybe possible to symmetrize H numerically if any issues exist, then take inverse
  # H <- 0.5*(H+t(H)))
  Hinv <- solve(H)
  
  # Sandwich variance
  V <- Hinv %*% M %*% Hinv
  se <- sqrt(pmax(0, diag(V)))
  
  # check to see if sum on per-time score equals to the gradient calculated one.
  grad_check <- max(abs(colSums(S) - obj$gr(par_hat)))
  
  list(se=se, vcov=V, scores=S, meat=M, hessian=H, grad_check=grad_check)
}
