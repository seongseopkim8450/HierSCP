
library(ggplot2)
library(ppmSuite)
library(bcp)
library(mvtnorm)
setwd("C://Users/User/Desktop/HierSCP_code/corefunctions_code") ### YOUR PATH

#install.packages("RcppArmadillo")
library(RcppArmadillo)

for (f in sort(list.files(pattern = "^[0-9].*\\.R$"))) source(f)


`%||%` <- function(a, b) if (!is.null(a)) a else b

# Single canonical definition used by every bench_* file.
.mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) Inf else mean(x)
}

# SETTLED upper-level changepoints for the shape-driven family (S1-S3).
#   cluster 1 -> tau = 40            (K = 2)
#   cluster 2 -> tau = 55            (K = 2)
#   cluster 3 -> tau = 40, 85        (K = 3)

SHAPE_CP_UPPER <- list(c(50L), c(42L), c(40L, 85L))

# Paths used by the parallel replication runner 
.HIERSCP_CPP    <- "hiercpd_core.cpp"    # C++ kernel
.HIERSCP_BRIDGE <- "15_rcpp_bridge.R"    # set_precomp_datx -> G_DATX cube bridge


###  SECTION 1.DATA GENERATORS  (S1-S3)
simulate_shape_driven <- function(
    scenario = c("S1_in_model", "S2_robust_contaminated", "S3_external_shape"),
    J = 9,
    T_len = 120,
    C = 3,
    cp_upper = list(c(50), c(42), c(40, 85)),
    M_basis = 8,
    m_min = 10,
    seed = 2026,
    return_scale = c("standardized", "raw"),
    attach_model_objects = TRUE
) {
  scenario <- match.arg(scenario)
  return_scale <- match.arg(return_scale)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (C != 3L) {
    stop("This simulation design is written for C = 3 groups.")
  }
  if (J < C) {
    stop("J must be at least C.")
  }
  if (length(cp_upper) != C) {
    stop("cp_upper must be a list of length C.")
  }
  if (T_len < 3L * m_min) {
    stop("T_len is too short relative to m_min.")
  }
  
  x <- seq(0, 1, length.out = T_len)
  
  scenario_cfg <- switch(
    scenario,
    
    S1_in_model = list(
      dgp = "hierscp",
      target_snr = 3.0,
      cp_jitter = 3L,
      hetero_noise = FALSE,
      hetero_log_sd = 0,
      outlier_prob = 0,
      outlier_mag = 0,
      outlier_df = 3,
      description = "Exact HierSCP shape-switching DGP with lower-level changepoint jitter."
    ),
    
    S2_robust_contaminated = list(
      dgp = "hierscp",
      target_snr = 1.5,
      cp_jitter = 5L,
      hetero_noise = TRUE,
      hetero_log_sd = 0.35,
      outlier_prob = 0.05,
      outlier_mag = 6,
      outlier_df = 3,
      description = "HierSCP mean DGP under low SNR, heteroscedastic noise, and heavy-tailed outliers."
    ),
    
    S3_external_shape = list(
      dgp = "external",
      target_snr = 2.5,
      cp_jitter = 3L,
      hetero_noise = TRUE,
      hetero_log_sd = 0.05,
      outlier_prob = 0.05,
      outlier_mag = 6,
      outlier_df = 3,
      amp_min = 0.95,
      amp_max = 1.05,
      phase_max = 0.005,
      description = "misspecified C1 piecewise-polynomial shape Data Generating Process"
    )
  )
  
  validate_cps <- function(cps, T_len, m_min) {
    cps <- as.integer(cps)
    
    if (length(cps) == 0L) {
      return(integer(0))
    }
    
    if (any(cps <= 1L) || any(cps > T_len)) {
      stop("All changepoints must lie in {2, ..., T_len}.")
    }
    
    cps <- sort(unique(cps))
    
    seg_start <- c(1L, cps)
    seg_end <- c(cps - 1L, T_len)
    seg_len <- seg_end - seg_start + 1L
    
    if (any(seg_len < m_min)) {
      stop("At least one segment violates m_min.")
    }
    
    cps
  }
  
  state_from_cps <- function(cps, T_len) {
    s <- rep(1L, T_len)
    
    if (length(cps) > 0L) {
      for (k in seq_along(cps)) {
        s[cps[k]:T_len] <- as.integer(k + 1L)
      }
    }
    
    s
  }
  
  jitter_cps <- function(cps, T_len, m_min, jitter) {
    cps <- as.integer(cps)
    
    if (length(cps) == 0L || jitter <= 0L) {
      return(validate_cps(cps, T_len, m_min))
    }
    
    out <- cps
    
    for (k in seq_along(cps)) {
      proposal <- cps[k] + sample(seq.int(-jitter, jitter), size = 1L)
      
      left_bound <- if (k == 1L) {
        1L + m_min
      } else {
        out[k - 1L] + m_min
      }
      
      right_bound <- if (k < length(cps)) {
        cps[k + 1L] - m_min
      } else {
        T_len - m_min + 1L
      }
      
      out[k] <- max(left_bound, min(right_bound, proposal))
    }
    
    validate_cps(out, T_len, m_min)
  }
  
  cumtrapz <- function(x, y) {
    if (length(x) != length(y)) {
      stop("x and y must have the same length.")
    }
    
    c(0, cumsum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2))
  }
  
  shape_profile_H <- function(theta, x) {
    M <- length(theta)
    
    Phi <- outer(
      x,
      seq_len(M),
      FUN = function(xx, mm) sqrt(2) * cos(mm * pi * xx)
    )
    
    Z <- as.vector(Phi %*% theta)
    Z2 <- Z^2
    
    Q <- tail(cumtrapz(x, Z2), 1)
    
    if (!is.finite(Q) || Q < 1e-10) {
      stop("Degenerate GP profile: integral of Z^2 is too close to zero.")
    }
    
    int_0_to_x <- cumtrapz(x, Z2)
    tail_int <- Q - int_0_to_x
    I <- cumtrapz(x, tail_int)
    
    H <- I / Q
    H <- pmin(pmax(H, 0), 1)
    
    list(H = H, Z = Z, Q = Q)
  }
  
  hierscp_shape_value <- function(x_value, H_value, delta1, delta2, b, g) {
    linear_coef <- delta1 * g + b * (delta1 + delta2) / 2
    linear_coef * x_value - delta2 * b * H_value
  }
  
  external_shape_value <- function(cc, k, x_value, x_cp, atom = NULL) {
    if (cc == 1L) {
      x1 <- x_cp[1]
      
      # G1 state 1: decreasing / convex (-,+)
      # Keep the reversal visible, but avoid an excessively sharp V-shape.
      q1 <- 0.45
      s1 <- 2 * q1 * x1 + 0.10
      
      r11 <- function(x) {
        0.80 - s1 * x + q1 * x^2
      }
      
      if (k == 1L) {
        return(r11(x_value))
      }
      
      # G1 state 2: increasing / concave (+,-)
      # Moderate positive slope at the left boundary, still positive at x = 1.
      dx <- pmax(x_value - x1, 0)
      
      return(
        r11(x1) + 0.35 * dx - 0.25 * dx^2
      )
    }
    
    if (cc == 2L) {
      x1 <- x_cp[1]
      
      # G2 state 1: increasing / convex (+,+)
      r21 <- function(x) {
        0.20 + 0.30 * x + 0.45 * x^2
      }
      
      if (k == 1L) {
        return(r21(x_value))
      }
      
      # G2 state 2: increasing / concave (+,-)
      # C1 connection: derivative at the start of state 2 equals derivative of state 1 at x1.
      d1 <- 0.30 + 0.90 * x1
      dx <- pmax(x_value - x1, 0)
      
      return(
        r21(x1) + d1 * dx - 0.55 * dx^2
      )
    }
    
    if (cc == 3L) {
      if (length(x_cp) < 2L) {
        stop("G3 requires two changepoints.")
      }
      
      x1 <- x_cp[1]
      x2 <- x_cp[2]
      
      if (x2 <= x1) {
        stop("For G3, x_cp[2] must be larger than x_cp[1].")
      }
      
      L1 <- x1
      L2 <- x2 - x1
      L3 <- 1 - x2
      
      if (L1 <= 0 || L2 <= 0 || L3 <= 0) {
        stop("Invalid G3 segment lengths.")
      }
      
      # State 1: increasing / concave (+,-)
      # Mild growth. Avoid explosive early increase.
      y0  <- -0.10
      s10 <-  0.28
      s11 <-  0.10
      
      a1 <- (s10 - s11) / L1
      
      r31 <- function(x) {
        y0 + s10 * x - 0.5 * a1 * x^2
      }
      
      y1 <- r31(x1)
      
      # State 2: decreasing / concave (-,-)
      # This is the important fix:
      # the slope starts mildly negative and becomes clearly more negative.
      s20 <- -0.04
      s21 <- -0.55
      
      a2 <- (s20 - s21) / L2
      
      if (!is.finite(a2) || a2 <= 0) {
        stop("Invalid G3 state-2 curvature.")
      }
      
      r32 <- function(x) {
        dx <- pmax(x - x1, 0)
        y1 + s20 * dx - 0.5 * a2 * dx^2
      }
      
      y2 <- r32(x2)
      
      if (k == 1L) {
        return(r31(x_value))
      }
      
      if (k == 2L) {
        return(r32(x_value))
      }
      
      # State 3: decreasing / convex (-,+)
      # Starts with a mild negative slope and flattens toward almost zero.
      # This prevents the third segment from crashing downward.
      s30 <- -0.12
      s31 <- -0.01
      
      q3 <- (s31 - s30) / L3
      
      if (!is.finite(q3) || q3 <= 0) {
        stop("Invalid G3 state-3 curvature.")
      }
      
      dx <- pmax(x_value - x2, 0)
      
      return(
        y2 + s30 * dx + 0.5 * q3 * dx^2
      )
    }
    
    stop("Invalid group index.")
  }
  
  make_cluster_assign <- function(J, C) {
    base_size <- J %/% C
    remainder <- J %% C
    
    sizes <- rep(base_size, C)
    if (remainder > 0L) {
      sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
    }
    
    rep(seq_len(C), times = sizes)
  }
  
  make_atoms <- function() {
    list(
      data.frame(
        delta1 = c(-1L, +1L),
        delta2 = c(+1L, -1L),
        b = c(1.25, 1.05),
        g = c(0.35, 0.45)
      ),
      data.frame(
        delta1 = c(+1L, +1L),
        delta2 = c(+1L, -1L),
        b = c(1.10, 0.95),
        g = c(0.30, 0.50)
      ),
      data.frame(
        delta1 = c(+1L, -1L, -1L),
        delta2 = c(-1L, -1L, +1L),
        b = c(1.20, 2.65, 1.80),
        g = c(0.55, 0.75, 1.25)
      )
    )
  }
  
  cp_upper <- lapply(cp_upper, validate_cps, T_len = T_len, m_min = m_min)
  K_by_cluster <- vapply(cp_upper, length, integer(1)) + 1L
  
  cluster_assign <- make_cluster_assign(J, C)
  atoms <- make_atoms()
  
  for (cc in seq_len(C)) {
    if (nrow(atoms[[cc]]) != K_by_cluster[cc]) {
      stop("The number of shape atoms does not match K_by_cluster.")
    }
  }
  
  S_upper <- matrix(1L, nrow = C, ncol = T_len)
  for (cc in seq_len(C)) {
    S_upper[cc, ] <- state_from_cps(cp_upper[[cc]], T_len)
  }
  
  cp_lower <- vector("list", J)
  S_lower <- matrix(1L, nrow = J, ncol = T_len)
  
  for (j in seq_len(J)) {
    cc <- cluster_assign[j]
    
    cp_lower[[j]] <- jitter_cps(
      cps = cp_upper[[cc]],
      T_len = T_len,
      m_min = m_min,
      jitter = scenario_cfg$cp_jitter
    )
    
    S_lower[j, ] <- state_from_cps(cp_lower[[j]], T_len)
  }
  
  theta0 <- matrix(NA_real_, nrow = C, ncol = M_basis)
  theta <- matrix(NA_real_, nrow = J, ncol = M_basis)
  
  for (cc in seq_len(C)) {
    for (m in seq_len(M_basis)) {
      theta0[cc, m] <- rnorm(1, mean = 0, sd = 1 / m)
    }
  }
  
  for (j in seq_len(J)) {
    cc <- cluster_assign[j]
    
    for (m in seq_len(M_basis)) {
      theta[j, m] <- rnorm(1, mean = theta0[cc, m], sd = 0.35 / m)
    }
  }
  
  H <- matrix(NA_real_, nrow = J, ncol = T_len)
  Z <- matrix(NA_real_, nrow = J, ncol = T_len)
  Q <- rep(NA_real_, J)
  
  for (j in seq_len(J)) {
    profile <- shape_profile_H(theta[j, ], x)
    H[j, ] <- profile$H
    Z[j, ] <- profile$Z
    Q[j] <- profile$Q
  }
  
  alpha <- 0
  beta_c <- c(0, 0.25, -0.25)
  gamma_j <- rnorm(J, mean = 0, sd = 0.10)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  amp_j <- runif(
    J,
    min = scenario_cfg$amp_min %||% 0.85,
    max = scenario_cfg$amp_max %||% 1.15
  )
  
  phase_j <- runif(
    J,
    min = -(scenario_cfg$phase_max %||% 0.03),
    max =  (scenario_cfg$phase_max %||% 0.03)
  )
  
  mu <- matrix(NA_real_, nrow = J, ncol = T_len)
  f_true <- matrix(NA_real_, nrow = J, ncol = T_len)
  
  for (j in seq_len(J)) {
    cc <- cluster_assign[j]
    
    x_tilde <- pmin(pmax(x + phase_j[j], 0), 1)
    
    x_cp <- if (length(cp_lower[[j]]) > 0L) {
      x_tilde[cp_lower[[j]]]
    } else {
      numeric(0)
    }
    
    for (t in seq_len(T_len)) {
      k <- S_lower[j, t]
      atom <- atoms[[cc]][k, ]
      
      if (scenario_cfg$dgp == "hierscp") {
        f_val <- hierscp_shape_value(
          x_value = x[t],
          H_value = H[j, t],
          delta1 = atom$delta1,
          delta2 = atom$delta2,
          b = atom$b,
          g = atom$g
        )
      } else if (scenario_cfg$dgp == "external") {
        f_val <- amp_j[j] * external_shape_value(
          cc = cc,
          k = k,
          x_value = x_tilde[t],
          x_cp = x_cp,
          atom = atom
        )
      } else {
        stop("Unknown DGP type.")
      }
      
      f_true[j, t] <- f_val
      mu[j, t] <- alpha + beta_c[cc] + gamma_j[j] + f_val
    }
  }
  
  hetero_mult <- if (scenario_cfg$hetero_noise) {
    rlnorm(J, meanlog = 0, sdlog = scenario_cfg$hetero_log_sd)
  } else {
    rep(1, J)
  }
  
  sigma_j <- rep(NA_real_, J)
  xi_true <- matrix(0L, nrow = J, ncol = T_len)
  eps <- matrix(NA_real_, nrow = J, ncol = T_len)
  Y_raw <- matrix(NA_real_, nrow = J, ncol = T_len)
  
  for (j in seq_len(J)) {
    signal_sd <- sd(mu[j, ])
    
    if (!is.finite(signal_sd) || signal_sd < 1e-8) {
      signal_sd <- 1
    }
    
    sigma_j[j] <- hetero_mult[j] * signal_sd / scenario_cfg$target_snr
    
    xi <- rbinom(
      n = T_len,
      size = 1L,
      prob = scenario_cfg$outlier_prob
    )
    
    xi_true[j, ] <- as.integer(xi)
    
    eps_j <- rnorm(T_len, mean = 0, sd = sigma_j[j])
    
    outlier_idx <- which(xi_true[j, ] == 1L)
    
    if (length(outlier_idx) > 0L) {
      shock <- rt(length(outlier_idx), df = scenario_cfg$outlier_df)
      
      if (scenario_cfg$outlier_df > 2) {
        shock <- shock / sqrt(scenario_cfg$outlier_df / (scenario_cfg$outlier_df - 2))
      }
      
      eps_j[outlier_idx] <- scenario_cfg$outlier_mag * sigma_j[j] * shock
    }
    
    eps[j, ] <- eps_j
    Y_raw[j, ] <- mu[j, ] + eps_j
  }
  
  row_mean <- rowMeans(Y_raw)
  row_sd <- apply(Y_raw, 1, sd)
  row_sd[row_sd < 1e-10] <- 1
  
  Y_std <- sweep(sweep(Y_raw, 1, row_mean, "-"), 1, row_sd, "/")
  mu_std <- sweep(sweep(mu, 1, row_mean, "-"), 1, row_sd, "/")
  f_std <- sweep(sweep(f_true, 1, row_mean, "-"), 1, row_sd, "/")
  
  shape_truth <- lapply(seq_len(C), function(cc) {
    out <- atoms[[cc]][, c("delta1", "delta2"), drop = FALSE]
    rownames(out) <- paste0("state_", seq_len(nrow(out)))
    out
  })
  names(shape_truth) <- paste0("G", seq_len(C))
  
  tau_upper <- lapply(cp_upper, function(cps) c(1L, cps))
  tau_lower <- lapply(cp_lower, function(cps) c(1L, cps))
  
  truth <- list(
    scenario = scenario,
    scenario_description = scenario_cfg$description,
    scenario_config = scenario_cfg,
    
    J = J,
    C = C,
    T_len = T_len,
    x = x,
    M_basis = M_basis,
    m_min = m_min,
    
    cluster_assign = cluster_assign,
    K_by_cluster = K_by_cluster,
    
    cp_upper = cp_upper,
    cp_lower = cp_lower,
    tau_upper = tau_upper,
    tau_lower = tau_lower,
    S_upper = S_upper,
    S_lower = S_lower,
    
    atoms = atoms,
    shape_truth = shape_truth,
    
    theta0 = theta0,
    theta = theta,
    H = H,
    Z = Z,
    Q = Q,
    
    alpha = alpha,
    beta_c = beta_c,
    gamma_j = gamma_j,
    amp_j = amp_j,
    phase_j = phase_j,
    
    f_true = f_true,
    f_std = f_std,
    mu = mu,
    mu_std = mu_std,
    
    Y_raw = Y_raw,
    Y_std = Y_std,
    eps = eps,
    sigma_j = sigma_j,
    sigma2_j = sigma_j^2,
    xi_true = xi_true,
    
    standardization = list(
      row_mean = row_mean,
      row_sd = row_sd
    ),
    
    seed = seed
  )
  
  if (
    attach_model_objects &&
    exists("create_model_spec", mode = "function") &&
    exists("precompute_basis", mode = "function") &&
    exists("precompute_global_x", mode = "function")
  ) {
    model <- create_model_spec(
      matrix(0, J, T_len),
      C = C,
      K_init = 1L,
      K_min = 1L,
      K_max = max(K_by_cluster) + 2L,
      M = M_basis,
      m_min = m_min,
      type = "continuous"
    )
    
    precomp <- precompute_basis(M_basis, n_grid = 600)
    precomp <- precompute_global_x(T_len, M_basis, precomp)
    
    precomp$std_info <- list(
      means = row_mean,
      sds = row_sd,
      Y_original = Y_raw
    )
    
    truth$model <- model
    truth$precomp <- precomp
  }
  
  Y_out <- if (return_scale == "standardized") Y_std else Y_raw
  
  list(
    Y = Y_out,
    truth = truth
  )
}


#####  SECTION 2.METRICS  

## metric functions compatible with CCP-PPM, ICP-PPM, and BCP 

ppm_normalize_cps <- function(cps, T_len) {
  if (is.null(cps) || length(cps) == 0L) return(integer(0))
  
  cps <- suppressWarnings(as.integer(round(unlist(cps))))
  cps <- cps[is.finite(cps)]
  cps <- cps[cps >= 2L & cps <= T_len]
  
  sort(unique(cps))
}

ppm_cps_to_segvec <- function(cps, T_len) {
  cps <- ppm_normalize_cps(cps, T_len)
  
  seg <- rep(1L, T_len)
  
  if (length(cps) == 0L) {
    return(seg)
  }
  
  for (i in seq_along(cps)) {
    seg[cps[i]:T_len] <- as.integer(i + 1L)
  }
  
  seg
}

ppm_metric_ari <- function(seg1, seg2) {
  n <- length(seg1)
  stopifnot(n == length(seg2))
  
  l1 <- sort(unique(seg1))
  l2 <- sort(unique(seg2))
  
  nij <- matrix(0, length(l1), length(l2))
  
  for (i in seq_along(l1)) {
    for (j in seq_along(l2)) {
      nij[i, j] <- sum(seg1 == l1[i] & seg2 == l2[j])
    }
  }
  
  ai <- rowSums(nij)
  bj <- colSums(nij)
  
  s_nij <- sum(choose(nij, 2))
  s_ai  <- sum(choose(ai, 2))
  s_bj  <- sum(choose(bj, 2))
  
  n2 <- choose(n, 2)
  
  if (n2 == 0) return(1)
  
  ex <- s_ai * s_bj / n2
  mx <- 0.5 * (s_ai + s_bj)
  
  if (mx == ex) {
    1
  } else {
    (s_nij - ex) / (mx - ex)
  }
}

ppm_metric_hausdorff <- function(set1, set2) {
  set1 <- as.integer(set1)
  set2 <- as.integer(set2)
  
  if (length(set1) == 0L && length(set2) == 0L) return(0)
  if (length(set1) == 0L || length(set2) == 0L) return(Inf)
  
  d12 <- max(vapply(set1, function(x) min(abs(x - set2)), numeric(1)))
  d21 <- max(vapply(set2, function(x) min(abs(x - set1)), numeric(1)))
  
  max(d12, d21)
}

ppm_metric_covering <- function(seg_true, seg_est) {
  T_len <- length(seg_true)
  total <- 0
  
  for (lt in sort(unique(seg_true))) {
    S_t <- which(seg_true == lt)
    best <- 0
    
    for (le in sort(unique(seg_est))) {
      S_e <- which(seg_est == le)
      
      inter <- length(intersect(S_t, S_e))
      uni   <- length(union(S_t, S_e))
      
      if (uni > 0) {
        best <- max(best, inter / uni)
      }
    }
    
    total <- total + length(S_t) * best
  }
  
  total / T_len
}

ppm_metric_f1 <- function(true_cps, est_cps, margin = 12) {
  true_cps <- as.integer(true_cps)
  est_cps  <- as.integer(est_cps)
  
  if (length(true_cps) == 0L && length(est_cps) == 0L) {
    return(list(f1 = 1, precision = 1, recall = 1, TP = 0L, FP = 0L, FN = 0L))
  }
  
  if (length(true_cps) == 0L) {
    return(list(f1 = 0, precision = 0, recall = 0, TP = 0L, FP = length(est_cps), FN = 0L))
  }
  
  if (length(est_cps) == 0L) {
    return(list(f1 = 0, precision = 0, recall = 0, TP = 0L, FP = 0L, FN = length(true_cps)))
  }
  
  matched_true <- logical(length(true_cps))
  matched_est  <- logical(length(est_cps))
  
  pairs <- expand.grid(
    ti = seq_along(true_cps),
    ei = seq_along(est_cps)
  )
  
  pairs$dist <- abs(true_cps[pairs$ti] - est_cps[pairs$ei])
  pairs <- pairs[order(pairs$dist), , drop = FALSE]
  
  for (r in seq_len(nrow(pairs))) {
    ti <- pairs$ti[r]
    ei <- pairs$ei[r]
    d  <- pairs$dist[r]
    
    if (d > margin) break
    
    if (!matched_true[ti] && !matched_est[ei]) {
      matched_true[ti] <- TRUE
      matched_est[ei]  <- TRUE
    }
  }
  
  TP <- sum(matched_true)
  FP <- sum(!matched_est)
  FN <- sum(!matched_true)
  
  precision <- if (TP + FP > 0) TP / (TP + FP) else 0
  recall    <- if (TP + FN > 0) TP / (TP + FN) else 0
  
  f1 <- if (precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }
  
  list(
    f1 = f1,
    precision = precision,
    recall = recall,
    TP = TP,
    FP = FP,
    FN = FN
  )
}

ppm_mean_finite <- function(x) {
  if (all(!is.finite(x))) return(Inf)
  mean(x[is.finite(x)])
}

### metric functions compatible with HierSCP  

.normalize_cps <- function(cps, T_len) {
  if (is.null(cps) || length(cps) == 0L) return(integer(0))
  
  cps <- suppressWarnings(as.integer(round(unlist(cps))))
  cps <- cps[is.finite(cps)]
  
  # In this simulation, changepoints are "start indices of the new segment".
  ### (Important) Therefore time point 1 is not a changepoint and should never be evaluated.
  cps <- cps[cps >= 2L & cps <= T_len]
  
  sort(unique(cps))
}

metric_ari <- function(seg1, seg2) {
  n <- length(seg1)
  stopifnot(n == length(seg2))
  
  l1 <- sort(unique(seg1))
  l2 <- sort(unique(seg2))
  
  nij <- matrix(0, length(l1), length(l2))
  
  for (i in seq_along(l1)) {
    for (j2 in seq_along(l2)) {
      nij[i, j2] <- sum(seg1 == l1[i] & seg2 == l2[j2])
    }
  }
  
  ai <- rowSums(nij)
  bj <- colSums(nij)
  
  s_nij <- sum(choose(nij, 2))
  s_ai  <- sum(choose(ai, 2))
  s_bj  <- sum(choose(bj, 2))
  
  n2 <- choose(n, 2)
  
  if (n2 == 0) return(1)
  
  ex <- s_ai * s_bj / n2
  mx <- 0.5 * (s_ai + s_bj)
  
  if (mx == ex) 1 else (s_nij - ex) / (mx - ex)
}

metric_hausdorff <- function(set1, set2) {
  set1 <- as.integer(set1)
  set2 <- as.integer(set2)
  
  if (length(set1) == 0L && length(set2) == 0L) return(0)
  if (length(set1) == 0L || length(set2) == 0L) return(Inf)
  
  d12 <- max(sapply(set1, function(x) min(abs(x - set2))))
  d21 <- max(sapply(set2, function(x) min(abs(x - set1))))
  
  max(d12, d21)
}

metric_covering <- function(seg_true, seg_est) {
  T_len <- length(seg_true)
  total <- 0
  
  for (lt in sort(unique(seg_true))) {
    S_t <- which(seg_true == lt)
    best <- 0
    
    for (le in sort(unique(seg_est))) {
      S_e <- which(seg_est == le)
      
      inter <- length(intersect(S_t, S_e))
      uni   <- length(union(S_t, S_e))
      
      if (uni > 0) best <- max(best, inter / uni)
    }
    
    total <- total + length(S_t) * best
  }
  
  total / T_len
}

metric_f1 <- function(true_cps, est_cps, margin = 10) {
  true_cps <- as.integer(true_cps)
  est_cps  <- as.integer(est_cps)
  
  if (length(true_cps) == 0L && length(est_cps) == 0L) {
    return(list(
      f1 = 1,
      precision = 1,
      recall = 1,
      TP = 0L,
      FP = 0L,
      FN = 0L
    ))
  }
  
  if (length(true_cps) == 0L) {
    return(list(
      f1 = 0,
      precision = 0,
      recall = 0,
      TP = 0L,
      FP = length(est_cps),
      FN = 0L
    ))
  }
  
  if (length(est_cps) == 0L) {
    return(list(
      f1 = 0,
      precision = 0,
      recall = 0,
      TP = 0L,
      FP = 0L,
      FN = length(true_cps)
    ))
  }
  
  matched_true <- logical(length(true_cps))
  matched_est  <- logical(length(est_cps))
  
  pairs <- expand.grid(
    ti = seq_along(true_cps),
    ei = seq_along(est_cps)
  )
  
  pairs$dist <- abs(true_cps[pairs$ti] - est_cps[pairs$ei])
  pairs <- pairs[order(pairs$dist), , drop = FALSE]
  
  for (r in seq_len(nrow(pairs))) {
    ti <- pairs$ti[r]
    ei <- pairs$ei[r]
    d  <- pairs$dist[r]
    
    if (d > margin) break
    
    if (!matched_true[ti] && !matched_est[ei]) {
      matched_true[ti] <- TRUE
      matched_est[ei]  <- TRUE
    }
  }
  
  TP <- sum(matched_true)
  FP <- sum(!matched_est)
  FN <- sum(!matched_true)
  
  precision <- if (TP + FP > 0) TP / (TP + FP) else 0
  recall    <- if (TP + FN > 0) TP / (TP + FN) else 0
  
  f1 <- if (precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }
  
  list(
    f1 = f1,
    precision = precision,
    recall = recall,
    TP = TP,
    FP = FP,
    FN = FN
  )
}

cps_to_segvec <- function(cps, T_len) {
  cps <- .normalize_cps(cps, T_len)
  
  seg <- rep(1L, T_len)
  
  if (length(cps) == 0L) return(seg)
  
  # IMPORTANT:
  # In simulate_hierscp_continuous_3scenarios(),
  # cp = t means state k + 1 starts at time t.
  for (i in seq_along(cps)) {
    seg[cps[i]:T_len] <- as.integer(i + 1L)
  }
  
  seg
}

.nvi <- function(seg_true, seg_est) {
  n <- length(seg_true)
  if (n != length(seg_est)) stop(".nvi: segmentations must have equal length.")
  if (n <= 1L) return(0)
  lt <- sort(unique(seg_true)); le <- sort(unique(seg_est))
  pi_ <- tabulate(match(seg_true, lt), nbins = length(lt)) / n
  qj_ <- tabulate(match(seg_est,  le), nbins = length(le)) / n
  vi <- 0
  for (i in seq_along(lt)) {
    ti <- which(seg_true == lt[i])
    for (j in seq_along(le)) {
      nij <- sum(seg_est[ti] == le[j])
      if (nij > 0L) {
        pij <- nij / n
        vi <- vi - pij * (log(pij / pi_[i]) + log(pij / qj_[j]))
      }
    }
  }
  denom <- log(n)
  if (!is.finite(denom) || denom <= 0) return(0)
  val <- vi / denom
  if (!is.finite(val)) return(NA_real_)
  min(max(val, 0), 1)
}

# Thin API-symmetric aliases (engine is shared; competitors score RAW segments,
# HierSCP scores SIGN-COLLAPSED segments — same metric, different input).
metric_nvi     <- function(seg_true, seg_est) .nvi(seg_true, seg_est)   # HierSCP path
ppm_metric_nvi <- function(seg_true, seg_est) .nvi(seg_true, seg_est)   # competitor path


#### SECTION 3.HierSCP POSTERIOR SUMMARY

compute_conditional_posterior_summary <- function(mcmc_result, Y,
                                                  prob_lower = 0.025,
                                                  prob_upper = 0.975,
                                                  cp_window = NULL,
                                                  verbose = TRUE) {
  
  samples <- mcmc_result$samples
  model   <- mcmc_result$model
  precomp <- mcmc_result$precomp
  
  J     <- model$J
  T_len <- model$T_len
  C     <- model$C
  M     <- model$M
  B     <- samples$n_saved
  
  if (B < 1) stop("No saved samples found.")
  
  # CP windowed mode bandwidth: default = floor(m_min / 4)
  if (is.null(cp_window)) {
    cp_window <- as.integer(floor(model$m_min / 4))
  }
  cp_window <- max(cp_window, 0L)
  if (verbose) cat(sprintf("CP location mode estimation: windowed mode (w=%d)\n", cp_window))
  
  # === 데이터 준비 ===
  std_info <- precomp$std_info
  if (!is.null(std_info)) {
    Y_work <- sweep(sweep(Y, 1, std_info$means, "-"), 1, std_info$sds, "/")
  } else {
    Y_work <- Y
  }
  
  # === Group 멤버십 ===
  group_members <- model$group_members
  if (is.null(group_members)) {
    n_per <- J %/% C
    group_members <- lapply(1:C, function(g) {
      start <- (g - 1L) * n_per + 1L
      end <- if (g < C) g * n_per else J
      as.integer(start:end)
    })
  }
  series_to_group <- integer(J)
  for (g in seq_along(group_members)) {
    for (jj in group_members[[g]]) series_to_group[jj] <- g
  }
  
  if (verbose) cat(sprintf("\n=== Shape-Driven Conditional MAP Posterior Averaging ===\n"))
  if (verbose) cat(sprintf("Total post-burnin samples: %d, Groups: %d\n", B, C))
  
  has_atoms <- !is.null(samples$atoms_store) && length(samples$atoms_store) > 0
  
  # =================================================================
  # Phase 1: Shape Sequence Joint Mode
  #
  # 핵심: (K, 형상 시퀀스)를 joint로 추정한다.
  # 각 iteration의 shape sequence key = "K:(d1,d2)_1,...,(d1,d2)_K"
  # → 사후 최빈 shape sequence → K_hat, (δ̂1,δ̂2)_k 동시 결정
  #
  # 장점:
  #   - K와 형상이 mutually consistent
  #   - 2단계 필터링(K→shape)의 정보 손실 없음
  #   - HierSCP의 "shape-driven" 정체성과 일치
  # =================================================================
  group_info <- vector("list", C)
  K_hat <- integer(C)
  shape_class <- vector("list", C)
  delta_mean  <- vector("list", C)
  shape_filter_map <- vector("list", C)
  
  for (g in 1:C) {
    j_in_g <- group_members[[g]]
    
    K_g <- integer(B)
    cluster_of_g <- integer(B)
    coherent <- logical(B)
    
    for (b in 1:B) {
      clusters_b <- samples$cluster[b, j_in_g]
      unique_cl <- unique(clusters_b)
      
      if (length(unique_cl) == 1L) {
        coherent[b] <- TRUE
        cc <- unique_cl[1]
      } else {
        coherent[b] <- FALSE
        cc <- as.integer(names(which.max(table(clusters_b))))
      }
      
      cluster_of_g[b] <- cc
      K_g[b] <- samples$K[b, cc]
    }
    
    K_g_coherent <- K_g[coherent]
    if (length(K_g_coherent) == 0) {
      warning(sprintf("Group %d: no coherent iterations. Using all.", g))
      coherent[] <- TRUE
    }
    
    # --- Shape sequence 수집 ---
    seq_keys <- character(B)  # 전체 iteration의 shape key
    for (b in 1:B) {
      if (!coherent[b]) { seq_keys[b] <- NA; next }
      cc <- cluster_of_g[b]
      K_b <- K_g[b]
      
      if (has_atoms && !is.null(samples$atoms_store[[b]]) &&
          !is.null(samples$atoms_store[[b]][[cc]])) {
        shapes <- character(K_b)
        valid <- TRUE
        for (k in 1:K_b) {
          atom_k <- if (k <= length(samples$atoms_store[[b]][[cc]])) {
            samples$atoms_store[[b]][[cc]][[k]]
          } else NULL
          if (is.null(atom_k)) { valid <- FALSE; break }
          d1 <- sign(as.numeric(atom_k$gamma1))
          d2 <- sign(as.numeric(atom_k$gamma2))
          if (d1 == 0) d1 <- 1L; if (d2 == 0) d2 <- -1L
          shapes[k] <- sprintf("(%d,%d)", d1, d2)
        }
        if (valid) {
          seq_keys[b] <- sprintf("%d:%s", K_b, paste(shapes, collapse=","))
        } else {
          seq_keys[b] <- sprintf("%d:NA", K_b)
        }
      } else {
        # atoms 미저장 → K만으로 key 생성 (fallback)
        seq_keys[b] <- sprintf("%d:NA", K_b)
      }
    }
    
    # --- Joint mode: shape sequence ---
    valid_keys <- seq_keys[coherent & !is.na(seq_keys)]
    seq_tab <- sort(table(valid_keys), decreasing = TRUE)
    
    modal_seq_key <- names(seq_tab)[1]
    modal_freq <- as.numeric(seq_tab[1]) / length(valid_keys)
    
    # Parse modal sequence
    parts <- strsplit(modal_seq_key, ":")[[1]]
    K_hat_g <- as.integer(parts[1])
    K_hat[g] <- K_hat_g
    
    # Extract (δ1, δ2) from modal key
    shape_class[[g]] <- vector("list", K_hat_g)
    delta_mean[[g]]  <- vector("list", K_hat_g)
    shape_filter_map[[g]] <- character(K_hat_g)
    
    if (parts[2] != "NA" && K_hat_g > 0) {
      shape_strs <- strsplit(parts[2], ",")[[1]]
      # shape_strs might be like "(1" "1)" "(-1" "-1)" for K=2
      # Need to re-parse: "(d1,d2),(d1,d2)"
      shape_pairs <- regmatches(parts[2], gregexpr("\\([^)]+\\)", parts[2]))[[1]]
      
      for (k in seq_along(shape_pairs)) {
        pair <- gsub("[()]", "", shape_pairs[k])
        dd <- as.integer(strsplit(pair, ",")[[1]])
        d1_modal <- dd[1]; d2_modal <- dd[2]
        
        shape_key <- get_shape_key(d1_modal, d2_modal)
        shape_filter_map[[g]][k] <- shape_key
        
        # Atom 통계량 수집 (B_g 전체에서)
        d1_all <- integer(0); d2_all <- integer(0)
        b_all <- numeric(0); g_all <- numeric(0)
        B_g_coh <- which(coherent & K_g >= k)
        for (bi in B_g_coh) {
          cc_bi <- cluster_of_g[bi]
          if (has_atoms && !is.null(samples$atoms_store[[bi]][[cc_bi]]) &&
              k <= length(samples$atoms_store[[bi]][[cc_bi]])) {
            atom_bi <- samples$atoms_store[[bi]][[cc_bi]][[k]]
            if (!is.null(atom_bi)) {
              d1_all <- c(d1_all, sign(as.numeric(atom_bi$gamma1)))
              d2_all <- c(d2_all, sign(as.numeric(atom_bi$gamma2)))
              b_all <- c(b_all, as.numeric(atom_bi$shape_beta))
              g_all <- c(g_all, as.numeric(atom_bi$shape_gamma))
            }
          }
        }
        
        d1_prob_pos <- if (length(d1_all) > 0) mean(d1_all > 0) else 0.5
        d2_prob_pos <- if (length(d2_all) > 0) mean(d2_all > 0) else 0.5
        
        shape_class[[g]][[k]] <- list(
          shape_key     = shape_key,
          d1_sign       = d1_modal,
          d2_sign       = d2_modal,
          d1_confidence = max(d1_prob_pos, 1 - d1_prob_pos),
          d2_confidence = max(d2_prob_pos, 1 - d2_prob_pos),
          label         = shape_labels_short()[shape_key]
        )
        
        delta_mean[[g]][[k]] <- list(
          d1_prob_pos = d1_prob_pos,
          d2_prob_pos = d2_prob_pos,
          shape_beta_mean  = if (length(b_all) > 0) mean(b_all) else NA,
          shape_gamma_mean = if (length(g_all) > 0) mean(g_all) else NA,
          n_atoms = length(d1_all),
          method = "atom_direct"
        )
      }
    } else {
      # Fallback: atoms 없을 때 K만으로 구성
      for (k in 1:K_hat_g) {
        shape_filter_map[[g]][k] <- "unknown"
        shape_class[[g]][[k]] <- list(
          shape_key="unknown", d1_sign=1, d2_sign=-1,
          d1_confidence=0.5, d2_confidence=0.5,
          label="?"
        )
        delta_mean[[g]][[k]] <- list(
          d1_prob_pos=0.5, d2_prob_pos=0.5,
          shape_beta_mean=NA, shape_gamma_mean=NA,
          n_atoms=0L, method="fallback"
        )
      }
    }
    
    # B_g = shape sequence 일치 iteration
    B_g <- which(coherent & seq_keys == modal_seq_key)
    
    # K 사후분포 (보고용)
    k_tab_g <- sort(table(K_g[coherent]), decreasing = TRUE)
    k_freq_g <- 100 * length(B_g) / sum(coherent)
    
    group_info[[g]] <- list(
      K_hat        = K_hat_g,
      K_freq       = k_freq_g,
      B_g          = B_g,
      cluster_of_g = cluster_of_g,
      K_g          = K_g,
      coherent     = coherent,
      k_posterior  = k_tab_g,
      n_coherent   = sum(coherent),
      n_filtered   = length(B_g),
      modal_seq_key = modal_seq_key,
      modal_seq_freq = modal_freq,
      seq_table    = seq_tab
    )
    
    if (verbose) {
      cat(sprintf("  Group %d: modal_seq='%s' (%.1f%%), K_hat=%d, B_g=%d\n",
                  g, modal_seq_key, 100*modal_freq, K_hat_g, length(B_g)))
      # Top 3 sequences
      for (i in 1:min(3, length(seq_tab))) {
        cat(sprintf("    seq '%s': %d (%.1f%%)\n",
                    names(seq_tab)[i], as.integer(seq_tab[i]),
                    100 * as.numeric(seq_tab[i]) / length(valid_keys)))
      }
      # K marginal
      cat(sprintf("    K marginal: "))
      for (i in 1:min(4, length(k_tab_g))) {
        cat(sprintf("K=%s: %d (%.1f%%)  ", names(k_tab_g)[i],
                    as.integer(k_tab_g[i]),
                    100 * as.numeric(k_tab_g[i]) / sum(coherent)))
      }
      cat("\n")
    }
  }
  
  # =================================================================
  # Phase 2: Per-series B_j 결정 + μ 재구성
  #
  # B_j = {b ∈ B_g : #cp(S_lower[b,j,]) = K_hat_g}
  # 하위 CP 수가 K_hat_g인 iteration만 사용하여 상위-하위 일치성 보장
  # =================================================================
  B_j_list <- vector("list", J)
  for (j in 1:J) {
    g <- series_to_group[j]
    B_g <- group_info[[g]]$B_g
    K_hat_g <- group_info[[g]]$K_hat
    
    B_j <- integer(0)
    for (b in B_g) {
      n_cp_j <- length(extract_changepoints(samples$S_lower[b, j, ]))
      if (n_cp_j == K_hat_g) B_j <- c(B_j, b)
    }
    
    if (length(B_j) == 0) B_j <- B_g   # fallback
    B_j_list[[j]] <- B_j
  }
  
  # μ 재구성: 필요한 iteration 합집합
  #  - B_j (조건부 μ용): 형상 일치 + 하위 CP 수 일치 표본
  #  - coherent 전체 (K-BMA μ용): 클러스터 일관된 모든 표본 (D로 거르지 않음)
  #    → K-BMA = (1/B) Σ_b μ^(b)(t), 전체 표본 평균이 곧 형상 시퀀스 marginalize
  coherent_all <- integer(0)
  for (g in 1:C) {
    coh_g <- which(group_info[[g]]$coherent)
    coherent_all <- union(coherent_all, coh_g)
  }
  coherent_all <- sort(coherent_all)
  all_iters <- sort(unique(c(unlist(B_j_list), coherent_all)))
  if (verbose) cat(sprintf("Reconstructing mu_j(t) for %d unique iterations (incl. %d coherent for K-BMA)...\n",
                           length(all_iters), length(coherent_all)))
  
  mu_store <- list()
  for (idx in seq_along(all_iters)) {
    b <- all_iters[idx]
    params_b <- reconstruct_params_from_sample(b, samples, model)
    state_b  <- reconstruct_state_from_sample(b, samples, model)
    
    mu_mat <- matrix(0, J, T_len)
    for (j in 1:J) {
      mu_mat[j, ] <- tryCatch(
        compute_mu_all(j, state_b, params_b, precomp, model),
        error = function(e) {
          rep(params_b$alpha + params_b$beta[state_b$cluster[j]], T_len)
        }
      )
    }
    mu_store[[as.character(b)]] <- mu_mat
    
    if (verbose && idx %% 100 == 0) {
      cat(sprintf("  ... %d / %d iterations processed\n", idx, length(all_iters)))
    }
  }
  if (verbose) cat(sprintf("  Done. %d iterations reconstructed.\n", length(all_iters)))
  
  # =================================================================
  # Phase 3: Per-series 사후 통계량 (B_j 기반)
  # =================================================================
  mu_mean <- matrix(0, J, T_len)
  mu_lo   <- matrix(0, J, T_len)
  mu_hi   <- matrix(0, J, T_len)
  mu_sd   <- matrix(0, J, T_len)
  
  # (1) 조건부 μ: B_j (형상 시퀀스 일치 D̂ 조건부) — 기존 방식 보존
  for (j in 1:J) {
    B_j <- B_j_list[[j]]
    n_j <- length(B_j)
    if (n_j == 0) next
    
    mu_j_mat <- matrix(0, n_j, T_len)
    for (bi in seq_along(B_j)) {
      mu_j_mat[bi, ] <- mu_store[[as.character(B_j[bi])]][j, ]
    }
    
    mu_mean[j, ] <- colMeans(mu_j_mat, na.rm = TRUE)
    mu_lo[j, ]   <- apply(mu_j_mat, 2, quantile, probs = prob_lower, na.rm = TRUE)
    mu_hi[j, ]   <- apply(mu_j_mat, 2, quantile, probs = prob_upper, na.rm = TRUE)
    mu_sd[j, ]   <- apply(mu_j_mat, 2, sd, na.rm = TRUE)
  }
  
  # =================================================================
  # Phase 3b: K-BMA μ (전체 coherent 표본 평균 — 형상 시퀀스 marginalize)
  #
  #   μ_BMA(t) = (1/|B_coh,g|) Σ_{b ∈ B_coh,g} μ^(b)(t)
  #   D(형상 시퀀스)로 거르지 않음 → 모형 가중치 P(D|Y)가 표본 빈도에 녹아
  #   약분되어 전체 평균이 곧 K-BMA. CI는 형상 시퀀스 불확실성까지 반영.
  #
  #   series j 는 자기 그룹 g 의 coherent 표본을 사용 (그룹 내 클러스터 일관).
  # =================================================================
  mu_mean_bma <- matrix(0, J, T_len)
  mu_lo_bma   <- matrix(0, J, T_len)
  mu_hi_bma   <- matrix(0, J, T_len)
  mu_sd_bma   <- matrix(0, J, T_len)
  
  for (j in 1:J) {
    g <- series_to_group[j]
    B_coh <- which(group_info[[g]]$coherent)
    B_coh <- B_coh[as.character(B_coh) %in% names(mu_store)]
    n_coh <- length(B_coh)
    if (n_coh == 0) {
      # fallback: 조건부 μ 복사
      mu_mean_bma[j, ] <- mu_mean[j, ]
      mu_lo_bma[j, ]   <- mu_lo[j, ]
      mu_hi_bma[j, ]   <- mu_hi[j, ]
      mu_sd_bma[j, ]   <- mu_sd[j, ]
      next
    }
    
    mu_j_bma <- matrix(0, n_coh, T_len)
    for (bi in seq_along(B_coh)) {
      mu_j_bma[bi, ] <- mu_store[[as.character(B_coh[bi])]][j, ]
    }
    
    mu_mean_bma[j, ] <- colMeans(mu_j_bma, na.rm = TRUE)
    mu_lo_bma[j, ]   <- apply(mu_j_bma, 2, quantile, probs = prob_lower, na.rm = TRUE)
    mu_hi_bma[j, ]   <- apply(mu_j_bma, 2, quantile, probs = prob_upper, na.rm = TRUE)
    mu_sd_bma[j, ]   <- apply(mu_j_bma, 2, sd, na.rm = TRUE)
  }
  
  # =================================================================
  # Phase 4: CP 사후 확률 + 위치 MAP
  # =================================================================
  
  # --- 하위 CP 사후확률 (per-series, B_j 기반) ---
  cp_prob <- matrix(0, J, T_len)
  for (j in 1:J) {
    B_j <- B_j_list[[j]]
    n_j <- length(B_j)
    if (n_j == 0) next
    for (b in B_j) {
      S_j <- samples$S_lower[b, j, ]
      for (t in 2:T_len) {
        if (!is.na(S_j[t]) && !is.na(S_j[t-1]) && S_j[t] != S_j[t-1]) {
          cp_prob[j, t] <- cp_prob[j, t] + 1
        }
      }
    }
    cp_prob[j, ] <- cp_prob[j, ] / n_j
  }
  
  # --- 상위 CP 위치 MAP (per-group, B_g 기반) ---
  cp_summary <- vector("list", C)
  for (g in 1:C) {
    K_hat_g <- group_info[[g]]$K_hat
    B_g <- group_info[[g]]$B_g
    cluster_of_g <- group_info[[g]]$cluster_of_g
    
    if (K_hat_g <= 1 || length(B_g) == 0) {
      cp_summary[[g]] <- NULL
      next
    }
    
    tau_mat <- matrix(NA, length(B_g), K_hat_g)
    for (bi in seq_along(B_g)) {
      b <- B_g[bi]
      cc <- cluster_of_g[b]
      tau_b <- extract_changepoints(samples$S_upper[b, cc, ])
      if (length(tau_b) == K_hat_g) tau_mat[bi, ] <- tau_b
    }
    
    cp_locs <- list()
    for (k in 2:K_hat_g) {
      tau_k <- tau_mat[, k]
      tau_k <- tau_k[!is.na(tau_k)]
      if (length(tau_k) == 0) next
      cp_locs[[k-1]] <- list(
        k      = k,
        mean   = mean(tau_k),
        median = median(tau_k),
        sd     = sd(tau_k),
        lo     = quantile(tau_k, prob_lower),
        hi     = quantile(tau_k, prob_upper),
        mode   = windowed_mode(tau_k, w = cp_window, t_range = c(1L, T_len)),
        raw_mode = as.integer(names(which.max(table(tau_k)))),
        point_est = discrete_median(tau_k)
      )
    }
    cp_summary[[g]] <- cp_locs
  }
  
  # --- 하위 CP 위치 MAP (per-series, B_j 기반) ---
  lower_cp_summary <- vector("list", J)
  for (j in 1:J) {
    g <- series_to_group[j]
    K_hat_g <- group_info[[g]]$K_hat
    B_j <- B_j_list[[j]]
    
    if (K_hat_g <= 1 || length(B_j) == 0) {
      lower_cp_summary[[j]] <- NULL
      next
    }
    
    tau_mat_j <- matrix(NA, length(B_j), K_hat_g)
    for (bi in seq_along(B_j)) {
      b <- B_j[bi]
      tau_b <- extract_changepoints(samples$S_lower[b, j, ])
      if (length(tau_b) == K_hat_g) tau_mat_j[bi, ] <- tau_b
    }
    
    cp_locs_j <- list()
    for (k in 2:K_hat_g) {
      tau_k <- tau_mat_j[, k]
      tau_k <- tau_k[!is.na(tau_k)]
      if (length(tau_k) == 0) next
      cp_locs_j[[k-1]] <- list(
        k      = k,
        mean   = mean(tau_k),
        median = median(tau_k),
        sd     = sd(tau_k),
        lo     = quantile(tau_k, prob_lower),
        hi     = quantile(tau_k, prob_upper),
        mode   = windowed_mode(tau_k, w = cp_window, t_range = c(1L, T_len)),
        raw_mode = as.integer(names(which.max(table(tau_k)))),
        point_est = discrete_median(tau_k)
      )
    }
    lower_cp_summary[[j]] <- cp_locs_j
  }
  
  # =================================================================
  # Phase 4b: Partial Model Averaging (PMA)
  #
  # 단조 상태 구조에서 τ_k = "state k가 처음 나타나는 시점"은
  # K에 무관하게 정의됨. 따라서 τ_k 추정에 K≥k인 모든 표본 사용 가능.
  #
  # K_hat은 빈도(Phase 1)로 선택하되,
  # τ_k 위치는 K≥k인 전체 표본에서 windowed mode로 추정.
  # B_g 필터 대비 사용 표본 대폭 증가.
  #
  # Ref: "Monotone state ordering enables partial model averaging
  #       across different model orders, where τ_k is identifiable
  #       for any K ≥ k."
  # =================================================================
  
  if (verbose) cat("\n--- Partial Model Averaging (PMA) ---\n")
  
  # --- Upper CP PMA ---
  pma_upper <- vector("list", C)
  for (g in 1:C) {
    K_hat_g      <- group_info[[g]]$K_hat
    K_g          <- group_info[[g]]$K_g
    coherent     <- group_info[[g]]$coherent
    cluster_of_g <- group_info[[g]]$cluster_of_g
    
    if (K_hat_g <= 1) { pma_upper[[g]] <- NULL; next }
    
    cp_locs_pma <- list()
    for (k in 2:K_hat_g) {
      # Shape-prefix PMA: K≥k이고 첫 k-1개 state 형상이 modal과 일치하는 표본
      B_pma_k <- integer(0)
      for (b in which(coherent & K_g >= k)) {
        cc <- cluster_of_g[b]
        prefix_match <- TRUE
        if (has_atoms && k > 1) {
          for (kk in 1:(k-1)) {
            atom_kk <- NULL
            if (!is.null(samples$atoms_store[[b]][[cc]]) &&
                kk <= length(samples$atoms_store[[b]][[cc]])) {
              atom_kk <- samples$atoms_store[[b]][[cc]][[kk]]
            }
            if (is.null(atom_kk)) { prefix_match <- FALSE; break }
            d1_b <- sign(as.numeric(atom_kk$gamma1))
            d2_b <- sign(as.numeric(atom_kk$gamma2))
            if (d1_b == 0) d1_b <- 1L; if (d2_b == 0) d2_b <- -1L
            b_key <- get_shape_key(d1_b, d2_b)
            if (b_key != shape_filter_map[[g]][kk]) { prefix_match <- FALSE; break }
          }
        }
        if (prefix_match) B_pma_k <- c(B_pma_k, b)
      }
      
      tau_k_vals <- integer(0)
      for (b in B_pma_k) {
        cc <- cluster_of_g[b]
        tau_b <- extract_changepoints(samples$S_upper[b, cc, ])
        if (length(tau_b) >= k) tau_k_vals <- c(tau_k_vals, tau_b[k])
      }
      
      if (length(tau_k_vals) == 0) next
      cp_locs_pma[[k-1]] <- list(
        k       = k,
        n_used  = length(tau_k_vals),
        n_bg    = length(group_info[[g]]$B_g),
        mean    = mean(tau_k_vals),
        median  = median(tau_k_vals),
        sd      = sd(tau_k_vals),
        lo      = quantile(tau_k_vals, prob_lower),
        hi      = quantile(tau_k_vals, prob_upper),
        mode    = windowed_mode(tau_k_vals, w = cp_window, t_range = c(1L, T_len)),
        raw_mode = as.integer(names(which.max(table(tau_k_vals)))),
        point_est = discrete_median(tau_k_vals)
      )
    }
    pma_upper[[g]] <- cp_locs_pma
    
    if (verbose && length(cp_locs_pma) > 0) {
      for (item in cp_locs_pma) {
        cat(sprintf("  G%d τ_%d: PMA n=%d (B_g=%d, +%.0f%%), mode=%d\n",
                    g, item$k, item$n_used, item$n_bg,
                    100*(item$n_used - item$n_bg)/max(1, item$n_bg),
                    item$mode))
      }
    }
  }
  
  # --- Lower CP PMA (per-series) ---
  pma_lower <- vector("list", J)
  for (j in 1:J) {
    g <- series_to_group[j]
    K_hat_g      <- group_info[[g]]$K_hat
    K_g          <- group_info[[g]]$K_g
    coherent     <- group_info[[g]]$coherent
    
    if (K_hat_g <= 1) { pma_lower[[j]] <- NULL; next }
    
    cp_locs_pma_j <- list()
    for (k in 2:K_hat_g) {
      # Shape-prefix PMA for lower level
      cluster_of_g_j <- group_info[[g]]$cluster_of_g
      B_pma_k <- integer(0)
      for (b in which(coherent & K_g >= k)) {
        cc <- cluster_of_g_j[b]
        prefix_match <- TRUE
        if (has_atoms && k > 1) {
          for (kk in 1:(k-1)) {
            atom_kk <- NULL
            if (!is.null(samples$atoms_store[[b]][[cc]]) &&
                kk <= length(samples$atoms_store[[b]][[cc]])) {
              atom_kk <- samples$atoms_store[[b]][[cc]][[kk]]
            }
            if (is.null(atom_kk)) { prefix_match <- FALSE; break }
            d1_b <- sign(as.numeric(atom_kk$gamma1))
            d2_b <- sign(as.numeric(atom_kk$gamma2))
            if (d1_b == 0) d1_b <- 1L; if (d2_b == 0) d2_b <- -1L
            b_key <- get_shape_key(d1_b, d2_b)
            if (b_key != shape_filter_map[[g]][kk]) { prefix_match <- FALSE; break }
          }
        }
        if (prefix_match) B_pma_k <- c(B_pma_k, b)
      }
      
      tau_k_vals <- integer(0)
      for (b in B_pma_k) {
        tau_j_b <- extract_changepoints(samples$S_lower[b, j, ])
        if (length(tau_j_b) >= k) tau_k_vals <- c(tau_k_vals, tau_j_b[k])
      }
      
      if (length(tau_k_vals) == 0) next
      cp_locs_pma_j[[k-1]] <- list(
        k       = k,
        n_used  = length(tau_k_vals),
        mean    = mean(tau_k_vals),
        median  = median(tau_k_vals),
        sd      = sd(tau_k_vals),
        lo      = quantile(tau_k_vals, prob_lower),
        hi      = quantile(tau_k_vals, prob_upper),
        mode    = windowed_mode(tau_k_vals, w = cp_window, t_range = c(1L, T_len)),
        raw_mode = as.integer(names(which.max(table(tau_k_vals)))),
        point_est = discrete_median(tau_k_vals)
      )
    }
    pma_lower[[j]] <- cp_locs_pma_j
  }
  
  # --- BMA cp_prob: 전체 표본에서 CP 확률 (K 필터 없이) ---
  cp_prob_bma <- matrix(0, J, T_len)
  for (j in 1:J) {
    for (b in 1:B) {
      S_j <- samples$S_lower[b, j, ]
      for (t in 2:T_len) {
        if (!is.na(S_j[t]) && !is.na(S_j[t-1]) && S_j[t] != S_j[t-1]) {
          cp_prob_bma[j, t] <- cp_prob_bma[j, t] + 1
        }
      }
    }
    cp_prob_bma[j, ] <- cp_prob_bma[j, ] / B
  }
  
  # =================================================================
  # Phase 5: 스칼라 파라미터 사후 요약
  # =================================================================
  alpha_filtered <- samples$alpha[all_iters]
  beta_filtered  <- samples$beta[all_iters, , drop = FALSE]
  
  alpha_summary <- list(
    mean = mean(alpha_filtered, na.rm = TRUE),
    sd   = sd(alpha_filtered, na.rm = TRUE),
    lo   = quantile(alpha_filtered, prob_lower, na.rm = TRUE),
    hi   = quantile(alpha_filtered, prob_upper, na.rm = TRUE)
  )
  
  beta_summary <- lapply(1:C, function(cc) {
    bv <- beta_filtered[, cc]
    list(mean = mean(bv, na.rm = TRUE), sd = sd(bv, na.rm = TRUE),
         lo = quantile(bv, prob_lower, na.rm = TRUE),
         hi = quantile(bv, prob_upper, na.rm = TRUE))
  })
  
  # =================================================================
  # Phase 6: Group 수준 평균 곡선 + 신용구간
  # =================================================================
  upper_curves <- vector("list", C)
  for (g in 1:C) {
    j_in_g <- group_members[[g]]
    B_g <- group_info[[g]]$B_g
    
    if (length(j_in_g) == 0 || length(B_g) == 0) {
      upper_curves[[g]] <- list(mean = rep(0, T_len),
                                lo = rep(0, T_len), hi = rep(0, T_len))
      next
    }
    
    group_mu_mat <- matrix(NA, length(B_g), T_len)
    for (bi in seq_along(B_g)) {
      b <- B_g[bi]
      bkey <- as.character(b)
      if (!is.null(mu_store[[bkey]])) {
        sub <- mu_store[[bkey]][j_in_g, , drop = TRUE]
        if (is.matrix(sub)) {
          group_mu_mat[bi, ] <- colMeans(sub)
        } else {
          group_mu_mat[bi, ] <- sub
        }
      }
    }
    
    valid_rows <- complete.cases(group_mu_mat)
    if (sum(valid_rows) == 0) {
      upper_curves[[g]] <- list(mean = rep(0, T_len),
                                lo = rep(0, T_len), hi = rep(0, T_len))
      next
    }
    group_mu_mat <- group_mu_mat[valid_rows, , drop = FALSE]
    
    upper_curves[[g]] <- list(
      mean      = colMeans(group_mu_mat, na.rm = TRUE),
      lo        = apply(group_mu_mat, 2, quantile, probs = prob_lower, na.rm = TRUE),
      hi        = apply(group_mu_mat, 2, quantile, probs = prob_upper, na.rm = TRUE),
      sd        = apply(group_mu_mat, 2, sd, na.rm = TRUE),
      member_ids= j_in_g
    )
  }
  
  # =================================================================
  # Phase 6b: Group 수준 K-BMA 곡선 (coherent 전체 표본)
  #   그룹 평균 곡선의 K-BMA 버전. 형상 시퀀스로 거르지 않고 coherent 전체
  #   표본을 평균 → K·형상 불확실성이 그룹 CI에 반영된다.
  # =================================================================
  upper_curves_bma <- vector("list", C)
  for (g in 1:C) {
    j_in_g <- group_members[[g]]
    B_coh  <- which(group_info[[g]]$coherent)
    B_coh  <- B_coh[as.character(B_coh) %in% names(mu_store)]
    
    if (length(j_in_g) == 0 || length(B_coh) == 0) {
      upper_curves_bma[[g]] <- upper_curves[[g]]   # fallback: 조건부 곡선
      next
    }
    
    gmat <- matrix(NA, length(B_coh), T_len)
    for (bi in seq_along(B_coh)) {
      sub <- mu_store[[as.character(B_coh[bi])]][j_in_g, , drop = TRUE]
      gmat[bi, ] <- if (is.matrix(sub)) colMeans(sub) else sub
    }
    valid_rows <- complete.cases(gmat)
    if (sum(valid_rows) == 0) { upper_curves_bma[[g]] <- upper_curves[[g]]; next }
    gmat <- gmat[valid_rows, , drop = FALSE]
    
    upper_curves_bma[[g]] <- list(
      mean      = colMeans(gmat, na.rm = TRUE),
      lo        = apply(gmat, 2, quantile, probs = prob_lower, na.rm = TRUE),
      hi        = apply(gmat, 2, quantile, probs = prob_upper, na.rm = TRUE),
      sd        = apply(gmat, 2, sd, na.rm = TRUE),
      member_ids= j_in_g,
      n_bma     = nrow(gmat)
    )
  }
  
  # =================================================================
  # 결과 반환
  # =================================================================
  K_hat_key  <- paste(K_hat, collapse = ",")
  K_hat_freq <- mean(sapply(group_info, function(gi) gi$K_freq))
  n_filtered <- sum(sapply(group_info, function(gi) gi$n_filtered))
  
  result <- list(
    K_hat            = K_hat,
    K_hat_key        = K_hat_key,
    K_hat_freq       = K_hat_freq,
    n_filtered       = n_filtered,
    group_info       = group_info,
    
    mu_mean          = mu_mean,
    mu_lo            = mu_lo,
    mu_hi            = mu_hi,
    mu_sd            = mu_sd,
    
    # --- K-BMA μ (전체 coherent 표본 평균; 형상 시퀀스 marginalize) ---
    mu_mean_bma      = mu_mean_bma,
    mu_lo_bma        = mu_lo_bma,
    mu_hi_bma        = mu_hi_bma,
    mu_sd_bma        = mu_sd_bma,
    
    cp_prob          = cp_prob,
    cp_prob_bma      = cp_prob_bma,       # [PMA] BMA CP probability (전체 표본)
    cp_summary       = cp_summary,
    lower_cp_summary = lower_cp_summary,
    pma_upper        = pma_upper,          # [PMA] Upper CP partial model averaging
    pma_lower        = pma_lower,          # [PMA] Lower CP partial model averaging
    
    shape_class      = shape_class,
    delta_mean       = delta_mean,
    
    # --- 형상 시퀀스 사후분포 P(D|Y) (그룹별 상위 시퀀스 + 빈도) ---
    shape_seq_posterior = lapply(1:C, function(g) {
      st <- group_info[[g]]$seq_table
      if (is.null(st) || length(st) == 0) return(NULL)
      tot <- sum(st)
      data.frame(
        sequence = names(st),
        count    = as.integer(st),
        prob     = as.numeric(st) / tot,
        stringsAsFactors = FALSE
      )
    }),
    
    alpha_summary    = alpha_summary,
    beta_summary     = beta_summary,
    
    upper_curves     = upper_curves,
    upper_curves_bma = upper_curves_bma,   # [K-BMA] Group-level BMA curves
    
    B_j_list         = B_j_list,
    group_members    = group_members
  )
  
  # === 진행 요약 출력 ===
  if (verbose) {
    cat("\n--- Shape Classification (Shape Sequence Joint Mode) ---\n")
    for (g in 1:C) {
      K_hat_g <- K_hat[g]
      modal_key <- group_info[[g]]$modal_seq_key
      modal_freq <- group_info[[g]]$modal_seq_freq
      cat(sprintf("  Group %d: '%s' (%.0f%%)\n    ", g, modal_key, 100*modal_freq))
      labels <- sapply(1:K_hat_g, function(k) {
        sc <- shape_class[[g]][[k]]
        sprintf("k%d=%s(d1:%.0f%%,d2:%.0f%%)", k, sc$shape_key,
                100*sc$d1_confidence, 100*sc$d2_confidence)
      })
      cat(paste(labels, collapse = " -> "), "\n")
    }
    
    # --- 형상 시퀀스 사후분포 P(D|Y): 상위 후보 ---
    cat("\n--- Shape Sequence Posterior P(D|Y) ---\n")
    for (g in 1:C) {
      st <- group_info[[g]]$seq_table
      if (is.null(st) || length(st) == 0) { cat(sprintf("  Group %d: (none)\n", g)); next }
      tot <- sum(st)
      cat(sprintf("  Group %d (D_hat = '%s'):\n", g, group_info[[g]]$modal_seq_key))
      n_show <- min(4, length(st))
      for (i in 1:n_show) {
        marker <- if (names(st)[i] == group_info[[g]]$modal_seq_key) "  <- D_hat" else ""
        cat(sprintf("    '%s': %d/%d (%.1f%%)%s\n",
                    names(st)[i], as.integer(st[i]), tot,
                    100 * as.numeric(st[i]) / tot, marker))
      }
    }
    
    cat(sprintf("\n--- Upper CP Locations (D_hat-conditional discrete median, w=%d) ---\n", cp_window))
    for (g in 1:C) {
      if (is.null(cp_summary[[g]]) || length(cp_summary[[g]]) == 0) {
        cat(sprintf("  Group %d: no changepoints (K=1)\n", g))
        next
      }
      for (cp_info in cp_summary[[g]]) {
        cat(sprintf("  Group %d, CP %d: tau_hat=%d (discrete median), sd=%.1f, 95%%CI=[%.0f, %.0f]\n",
                    g, cp_info$k, cp_info$point_est, cp_info$sd,
                    cp_info$lo, cp_info$hi))
      }
    }
    
    cat("\n--- Per-series B_j sizes ---\n")
    for (j in 1:J) {
      g <- series_to_group[j]
      n_lcp <- if (!is.null(lower_cp_summary[[j]])) length(lower_cp_summary[[j]]) else 0
      lcp_modes <- ""
      if (n_lcp > 0) {
        lcp_modes <- paste(sapply(lower_cp_summary[[j]], function(x) {
          raw_str <- if (!is.null(x$raw_mode) && x$raw_mode != x$mode) {
            sprintf("%d[raw=%d]", x$mode, x$raw_mode)
          } else {
            sprintf("%d", x$mode)
          }
          raw_str
        }), collapse=",")
      }
      cat(sprintf("  S%d (G%d): B_j=%d, lower_CP_modes=[%s]\n",
                  j, g, length(B_j_list[[j]]), lcp_modes))
    }
    
    # --- K-BMA μ vs 조건부 μ 진단 ---
    # 두 μ 추정의 차이가 클수록 형상 시퀀스 불확실성이 큼. 변화점 근처에 차이가
    # 집중되면 K·형상이 표본마다 흔들린다는 신호. CI 폭 비교도 함께.
    cat("\n--- K-BMA mu vs Conditional mu (diagnostic) ---\n")
    for (j in 1:J) {
      diff_abs <- abs(mu_mean_bma[j, ] - mu_mean[j, ])
      ci_cond  <- mean(mu_hi[j, ] - mu_lo[j, ], na.rm = TRUE)
      ci_bma   <- mean(mu_hi_bma[j, ] - mu_lo_bma[j, ], na.rm = TRUE)
      t_max    <- which.max(diff_abs)
      cat(sprintf("  S%d: max|BMA-cond|=%.3f at t=%d, mean CI width cond=%.3f / BMA=%.3f (x%.2f)\n",
                  j, max(diff_abs, na.rm = TRUE), t_max, ci_cond, ci_bma,
                  if (ci_cond > 0) ci_bma / ci_cond else 1))
    }
  }
  
  result
}


####  SECTION 4.  HierSCP EVALUATION  (truth adapter, shape accuracy, all metrics)

extract_point_estimates <- function(cps_list, T_len = NULL) {
  if (is.null(cps_list) || length(cps_list) == 0L) return(integer(0))
  
  cps <- integer(0)
  
  for (info in cps_list) {
    if (is.null(info)) next
    
    if (!is.null(info$point_est) && is.finite(info$point_est)) {
      cps <- c(cps, as.integer(round(info$point_est)))
    } else if (!is.null(info$mean) && is.finite(info$mean)) {
      cps <- c(cps, as.integer(round(info$mean)))
    } else if (!is.null(info$median) && is.finite(info$median)) {
      cps <- c(cps, as.integer(round(info$median)))
    }
  }
  
  cps <- sort(unique(cps))
  
  if (!is.null(T_len)) {
    cps <- .normalize_cps(cps, T_len)
  }
  
  cps
}

.same_matrix <- function(A, B, tol = 1e-8) {
  if (is.null(A) || is.null(B)) return(FALSE)
  if (!identical(dim(A), dim(B))) return(FALSE)
  
  max(abs(A - B), na.rm = TRUE) < tol
}

adapt_sim_for_metrics <- function(sim_raw) {
  # Already adapted object
  if (is.null(sim_raw$truth)) {
    required <- c(
      "Y", "J", "T_len", "C", "group_members",
      "mu_true", "tau_group", "tau_series",
      "S_upper_true", "S_lower_true", "K_true"
    )
    
    missing <- setdiff(required, names(sim_raw))
    
    if (length(missing) > 0L) {
      stop(
        "sim_raw has no $truth and is missing required adapted fields: ",
        paste(missing, collapse = ", ")
      )
    }
    
    return(sim_raw)
  }
  
  tr <- sim_raw$truth
  
  Y <- sim_raw$Y
  
  J <- if (!is.null(tr$J)) tr$J else nrow(Y)
  T_len <- if (!is.null(tr$T_len)) tr$T_len else ncol(Y)
  C <- if (!is.null(tr$C)) tr$C else length(tr$cp_upper)
  
  cluster_assign <- tr$cluster_assign
  
  if (is.null(cluster_assign)) {
    if (!is.null(tr$state$cluster)) {
      cluster_assign <- tr$state$cluster
    } else {
      n_per <- J %/% C
      cluster_assign <- rep(seq_len(C), each = n_per)
      if (length(cluster_assign) < J) {
        cluster_assign <- c(cluster_assign, rep(C, J - length(cluster_assign)))
      }
    }
  }
  
  group_members <- lapply(seq_len(C), function(g) which(cluster_assign == g))
  
  cp_upper <- if (!is.null(tr$cp_upper)) {
    tr$cp_upper
  } else if (!is.null(tr$tau_upper)) {
    lapply(tr$tau_upper, function(x) .normalize_cps(x, T_len))
  } else {
    stop("Cannot find upper-level true changepoints in sim$truth.")
  }
  
  cp_lower <- if (!is.null(tr$cp_lower)) {
    tr$cp_lower
  } else if (!is.null(tr$tau_lower)) {
    lapply(tr$tau_lower, function(x) .normalize_cps(x, T_len))
  } else {
    stop("Cannot find lower-level true changepoints in sim$truth.")
  }
  
  cp_upper <- lapply(cp_upper, .normalize_cps, T_len = T_len)
  cp_lower <- lapply(cp_lower, .normalize_cps, T_len = T_len)
  
  S_upper_true <- if (!is.null(tr$S_upper)) {
    tr$S_upper
  } else if (!is.null(tr$state$S_upper)) {
    tr$state$S_upper
  } else {
    do.call(rbind, lapply(cp_upper, cps_to_segvec, T_len = T_len))
  }
  
  S_lower_true <- if (!is.null(tr$S_lower)) {
    tr$S_lower
  } else if (!is.null(tr$state$S_lower)) {
    tr$state$S_lower
  } else {
    do.call(rbind, lapply(cp_lower, cps_to_segvec, T_len = T_len))
  }
  
  K_true <- if (!is.null(tr$K_by_cluster)) {
    tr$K_by_cluster
  } else if (!is.null(tr$state$K)) {
    tr$state$K
  } else {
    vapply(cp_upper, length, integer(1)) + 1L
  }
  
  # Choose the truth scale matching sim_raw$Y.
  # Default generator returns standardized Y, but this also supports return_scale = "raw".
  mu_true <- NULL
  
  if (!is.null(tr$Y_std) && !is.null(tr$mu_std) && .same_matrix(Y, tr$Y_std)) {
    mu_true <- tr$mu_std
  } else if (!is.null(tr$Y_raw) && !is.null(tr$mu) && .same_matrix(Y, tr$Y_raw)) {
    mu_true <- tr$mu
  } else if (!is.null(tr$Y_original) && !is.null(tr$mu) && .same_matrix(Y, tr$Y_original)) {
    mu_true <- tr$mu
  } else if (!is.null(tr$mu_std)) {
    mu_true <- tr$mu_std
  } else if (!is.null(tr$mu)) {
    mu_true <- tr$mu
  } else {
    stop("Cannot find true mean curve: expected tr$mu_std or tr$mu.")
  }
  
  scenario_name <- if (!is.null(tr$scenario)) {
    tr$scenario
  } else if (!is.null(tr$snr_level)) {
    tr$snr_level
  } else {
    NA_character_
  }
  
  scenario <- list(
    name = scenario_name,
    description = tr$scenario_description,
    config = tr$scenario_config,
    sigma = if (!is.null(tr$sigma_j)) mean(tr$sigma_j) else if (!is.null(tr$sigma2_j)) sqrt(mean(tr$sigma2_j)) else NA_real_,
    p_out = if (!is.null(tr$xi_true)) mean(tr$xi_true) else if (!is.null(tr$outlier_prob)) tr$outlier_prob else NA_real_,
    jitter_delta = if (!is.null(tr$scenario_config$cp_jitter)) tr$scenario_config$cp_jitter else if (!is.null(tr$cp_jitter)) tr$cp_jitter else NA_integer_
  )
  
  group_specs <- lapply(seq_len(C), function(g) {
    list(name = paste0("Group ", g))
  })
  
  out <- list(
    Y = Y,
    J = J,
    T_len = T_len,
    C = C,
    n_per_group = J %/% C,
    cluster_assign = cluster_assign,
    group_members = group_members,
    
    mu_true = mu_true,
    z_true = tr$xi_true,
    
    tau_group = cp_upper,
    tau_series = cp_lower,
    
    S_upper_true = S_upper_true,
    S_lower_true = S_lower_true,
    
    K_true = as.integer(K_true),
    
    shape_truth = tr$shape_truth,
    atoms = tr$atoms,
    
    scenario = scenario,
    group_specs = group_specs,
    
    raw_truth = tr
  )
  
  out
}

infer_true_shapes <- function(sim) {
  sim_m <- adapt_sim_for_metrics(sim)
  
  C <- sim_m$C
  
  # Preferred route: use exact truth from the simulator.
  if (!is.null(sim_m$shape_truth)) {
    true_shapes <- vector("list", C)
    
    for (g in seq_len(C)) {
      st <- sim_m$shape_truth[[g]]
      
      if (is.data.frame(st) || is.matrix(st)) {
        true_shapes[[g]] <- vector("list", nrow(st))
        
        for (k in seq_len(nrow(st))) {
          true_shapes[[g]][[k]] <- list(
            d1 = as.integer(st[k, "delta1"]),
            d2 = as.integer(st[k, "delta2"])
          )
        }
      } else if (is.list(st)) {
        true_shapes[[g]] <- vector("list", length(st))
        
        for (k in seq_along(st)) {
          x <- st[[k]]
          
          true_shapes[[g]][[k]] <- list(
            d1 = as.integer(if (!is.null(x$d1)) x$d1 else x$delta1),
            d2 = as.integer(if (!is.null(x$d2)) x$d2 else x$delta2)
          )
        }
      } else {
        stop("Unsupported shape_truth format for group ", g)
      }
    }
    
    return(true_shapes)
  }
  
  # Fallback route: infer signs from the true mean curve.
  true_shapes <- vector("list", C)
  
  for (g in seq_len(C)) {
    K_g <- sim_m$K_true[g]
    true_shapes[[g]] <- vector("list", K_g)
    
    j_rep <- sim_m$group_members[[g]][1]
    mu_j <- sim_m$mu_true[j_rep, ]
    seg_g <- sim_m$S_upper_true[g, ]
    
    for (k in seq_len(K_g)) {
      tik <- which(seg_g == k)
      
      if (length(tik) < 4L) {
        true_shapes[[g]][[k]] <- list(d1 = NA_integer_, d2 = NA_integer_)
        next
      }
      
      mu_s <- mu_j[tik]
      x <- seq_along(mu_s)
      
      fit1 <- lm(mu_s ~ x)
      d1 <- if (coef(fit1)[2] > 0) 1L else -1L
      
      fit2 <- lm(mu_s ~ x + I(x^2))
      d2 <- if (coef(fit2)[3] > 0) 1L else -1L
      
      true_shapes[[g]][[k]] <- list(d1 = d1, d2 = d2)
    }
  }
  
  true_shapes
}

.get_post_shape_sign <- function(shape_obj) {
  if (is.null(shape_obj)) {
    return(list(d1 = NA_integer_, d2 = NA_integer_))
  }
  
  d1 <- NA_integer_
  d2 <- NA_integer_
  
  if (!is.null(shape_obj$d1_sign)) d1 <- as.integer(shape_obj$d1_sign)
  if (!is.null(shape_obj$d2_sign)) d2 <- as.integer(shape_obj$d2_sign)
  
  if (is.na(d1) && !is.null(shape_obj$delta1)) d1 <- as.integer(shape_obj$delta1)
  if (is.na(d2) && !is.null(shape_obj$delta2)) d2 <- as.integer(shape_obj$delta2)
  
  if ((is.na(d1) || is.na(d2)) && !is.null(shape_obj$shape_key)) {
    key <- shape_obj$shape_key
    
    map <- list(
      inc_convex  = c(+1L, +1L),
      inc_concave = c(+1L, -1L),
      dec_concave = c(-1L, -1L),
      dec_convex  = c(-1L, +1L)
    )
    
    if (!is.null(map[[key]])) {
      d1 <- map[[key]][1]
      d2 <- map[[key]][2]
    }
  }
  
  list(d1 = d1, d2 = d2)
}

compute_shape_accuracy <- function(post, true_shapes, sim) {
  sim_m <- adapt_sim_for_metrics(sim)
  
  C <- sim_m$C
  
  total <- 0L
  c1 <- 0L
  c2 <- 0L
  cb <- 0L
  
  per_group <- vector("list", C)
  
  if (is.null(post$shape_class)) {
    return(list(
      per_group = lapply(seq_len(C), function(g) {
        list(K_match = FALSE, d1 = NA_real_, d2 = NA_real_, both = NA_real_)
      }),
      overall_d1 = NA_real_,
      overall_d2 = NA_real_,
      overall = NA_real_
    ))
  }
  
  for (g in seq_len(C)) {
    Kt <- sim_m$K_true[g]
    Kh <- post$K_hat[g]
    
    if (is.na(Kh) || Kh != Kt || is.null(post$shape_class[[g]])) {
      per_group[[g]] <- list(
        K_match = FALSE,
        d1 = NA_real_,
        d2 = NA_real_,
        both = NA_real_
      )
      next
    }
    
    n1 <- 0L
    n2 <- 0L
    nb <- 0L
    
    for (k in seq_len(Kt)) {
      sc <- .get_post_shape_sign(post$shape_class[[g]][[k]])
      ts <- true_shapes[[g]][[k]]
      
      m1 <- !is.na(sc$d1) && !is.na(ts$d1) && sc$d1 == ts$d1
      m2 <- !is.na(sc$d2) && !is.na(ts$d2) && sc$d2 == ts$d2
      
      if (m1) {
        n1 <- n1 + 1L
        c1 <- c1 + 1L
      }
      
      if (m2) {
        n2 <- n2 + 1L
        c2 <- c2 + 1L
      }
      
      if (m1 && m2) {
        nb <- nb + 1L
        cb <- cb + 1L
      }
      
      total <- total + 1L
    }
    
    per_group[[g]] <- list(
      K_match = TRUE,
      d1 = n1 / Kt,
      d2 = n2 / Kt,
      both = nb / Kt
    )
  }
  
  list(
    per_group = per_group,
    overall_d1 = if (total > 0L) c1 / total else NA_real_,
    overall_d2 = if (total > 0L) c2 / total else NA_real_,
    overall = if (total > 0L) cb / total else NA_real_
  )
}

# # Sign-pair sequence of a cluster's full-atom shape_class (length = full-atom K_g),
# # each element c(delta1, delta2) (NA when unknown). Feeds the sign-collapse helpers.
# .post_sign_pairs <- function(post, g) {
#   sc <- if (!is.null(post$shape_class)) post$shape_class[[g]] else NULL
#   if (is.null(sc) || length(sc) == 0L) return(list())
#   lapply(sc, function(s) {
#     sg <- .get_post_shape_sign(s)
#     c(sg$d1, sg$d2)
#   })
# }

compute_all_metrics <- function(
    post,
    sim,
    K_true = NULL,
    J = NULL,
    T_len = NULL,
    f1_margin = 10,
    verbose = TRUE
) {
  sim_m <- adapt_sim_for_metrics(sim)
  
  if (!is.null(J)) sim_m$J <- J
  if (!is.null(T_len)) sim_m$T_len <- T_len
  if (!is.null(K_true)) sim_m$K_true <- as.integer(K_true)
  
  J <- sim_m$J
  T_len <- sim_m$T_len
  C <- sim_m$C
  
  
  
  m <- list()
  
  # --- HierSCP estimate: FULL-ATOM segmentation (NO sign-collapse) ---
  # Every estimated atom boundary is scored as a changepoint, and K_hat is the
  # full-atom count. ARI / Cover / NVI / F1 below are all computed on this same
  # raw segmentation. (NVI and the other new metrics are kept; only the
  # changepoint sign-collapse is removed.)
  m$K_true <- sim_m$K_true
  m$K_hat <- as.integer(post$K_hat)
  m$K_exact_match <- all(m$K_hat == m$K_true)
  m$K_match_per_group <- m$K_hat == m$K_true
  
  ari_g <- haus_g <- cov_g <- nvi_g <- numeric(C)
  f1_g <- vector("list", C)
  
  for (g in seq_len(C)) {
    est_cps <- extract_point_estimates(post$cp_summary[[g]], T_len)
    true_cps <- .normalize_cps(sim_m$tau_group[[g]], T_len)
    
    seg_est <- cps_to_segvec(est_cps, T_len)
    seg_true <- sim_m$S_upper_true[g, ]
    
    ari_g[g] <- metric_ari(seg_true, seg_est)
    haus_g[g] <- metric_hausdorff(true_cps, est_cps)
    cov_g[g] <- metric_covering(seg_true, seg_est)
    nvi_g[g] <- metric_nvi(seg_true, seg_est)
    f1_g[[g]] <- metric_f1(true_cps, est_cps, margin = f1_margin)
  }
  
  m$ARI_group <- ari_g
  m$ARI_group_mean <- mean(ari_g)
  
  m$Haus_group <- haus_g
  m$Haus_group_mean <- mean(haus_g[is.finite(haus_g)])
  if (all(!is.finite(haus_g))) m$Haus_group_mean <- Inf
  
  m$NVI_group <- nvi_g
  m$NVI_group_mean <- mean(nvi_g, na.rm = TRUE)
  
  m$Cov_group <- cov_g
  m$Cov_group_mean <- mean(cov_g)
  
  m$F1_group <- vapply(f1_g, function(x) x$f1, numeric(1))
  m$F1_group_mean <- mean(m$F1_group)
  m$F1_group_detail <- f1_g
  
  ari_s <- haus_s <- cov_s <- nvi_s <- numeric(J)
  f1_s <- vector("list", J)
  
  for (j in seq_len(J)) {
    est_cps <- extract_point_estimates(post$lower_cp_summary[[j]], T_len)
    true_cps <- .normalize_cps(sim_m$tau_series[[j]], T_len)
    
    seg_est <- cps_to_segvec(est_cps, T_len)
    seg_true <- sim_m$S_lower_true[j, ]
    
    ari_s[j] <- metric_ari(seg_true, seg_est)
    haus_s[j] <- metric_hausdorff(true_cps, est_cps)
    cov_s[j] <- metric_covering(seg_true, seg_est)
    nvi_s[j] <- metric_nvi(seg_true, seg_est)
    f1_s[[j]] <- metric_f1(true_cps, est_cps, margin = f1_margin)
  }
  
  m$ARI_series <- ari_s
  m$ARI_series_mean <- mean(ari_s)
  
  m$Haus_series <- haus_s
  m$Haus_series_mean <- mean(haus_s[is.finite(haus_s)])
  if (all(!is.finite(haus_s))) m$Haus_series_mean <- Inf
  
  m$NVI_series <- nvi_s
  m$NVI_series_mean <- mean(nvi_s, na.rm = TRUE)
  
  m$Cov_series <- cov_s
  m$Cov_series_mean <- mean(cov_s)
  
  m$F1_series <- vapply(f1_s, function(x) x$f1, numeric(1))
  m$F1_series_mean <- mean(m$F1_series)
  m$F1_series_detail <- f1_s
  
  if (is.null(post$mu_lo) || is.null(post$mu_hi) || is.null(post$mu_mean)) {
    warning("post$mu_lo, post$mu_hi, or post$mu_mean is missing. ECP and RISE are set to NA.")
    
    m$ECP <- rep(NA_real_, J)
    m$ECP_mean <- NA_real_
    
    m$RISE <- rep(NA_real_, J)
    m$RISE_mean <- NA_real_
  } else {
    ecp <- numeric(J)
    
    for (j in seq_len(J)) {
      covered <- sim_m$mu_true[j, ] >= post$mu_lo[j, ] &
        sim_m$mu_true[j, ] <= post$mu_hi[j, ]
      
      ecp[j] <- mean(covered, na.rm = TRUE)
    }
    
    m$ECP <- ecp
    m$ECP_mean <- mean(ecp, na.rm = TRUE)
    
    rise <- numeric(J)
    
    for (j in seq_len(J)) {
      rise[j] <- sqrt(mean((post$mu_mean[j, ] - sim_m$mu_true[j, ])^2, na.rm = TRUE))
    }
    
    m$RISE <- rise
    m$RISE_mean <- mean(rise, na.rm = TRUE)
  }
  
  true_shapes <- infer_true_shapes(sim_m)
  # Shape accuracy on the FULL-atom shape_class (no sign-collapse), conditional
  # on K_hat == K_true (handled inside compute_shape_accuracy).
  m$SCA <- compute_shape_accuracy(post, true_shapes, sim_m)
  
  if (verbose) {
    print_all_metrics(m, sim_m, f1_margin = f1_margin)
  }
  
  m
}

print_all_metrics <- function(m, sim, f1_margin = 10) {
  sim_m <- adapt_sim_for_metrics(sim)
  
  C <- sim_m$C
  J <- sim_m$J
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  EVALUATION METRICS\n")
  cat(strrep("=", 70), "\n", sep = "")
  
  if (!is.null(sim_m$scenario$name) && !is.na(sim_m$scenario$name)) {
    cat(sprintf("\n  Scenario: %s\n", sim_m$scenario$name))
  }
  
  cat(sprintf(
    "\n  K: true=[%s]  hat=[%s]  %s\n",
    paste(m$K_true, collapse = ","),
    paste(m$K_hat, collapse = ","),
    ifelse(m$K_exact_match, "EXACT MATCH", "MISMATCH")
  ))
  
  for (g in seq_len(C)) {
    cat(sprintf(
      "    G%d: K_true=%d  K_hat=%d  %s\n",
      g,
      m$K_true[g],
      m$K_hat[g],
      ifelse(m$K_match_per_group[g], "OK", "MISS")
    ))
  }
  
  cat(sprintf("\n  %-12s %8s %8s %8s %8s %12s\n", "GROUP", "ARI", "NVI", "Cover", "F1", "Match"))
  cat("  ", strrep("-", 62), "\n", sep = "")
  
  for (g in seq_len(C)) {
    f1d <- m$F1_group_detail[[g]]
    
    tp_str <- sprintf(
      "TP=%d FP=%d FN=%d",
      f1d$TP,
      f1d$FP,
      f1d$FN
    )
    
    cat(sprintf(
      "  G%-11d %7.3f %7.3f %7.3f %7.3f %12s\n",
      g,
      m$ARI_group[g],
      m$NVI_group[g],
      m$Cov_group[g],
      m$F1_group[g],
      tp_str
    ))
  }
  
  cat(sprintf(
    "  %-12s %7.3f %7.3f %7.3f %7.3f\n",
    "MEAN",
    m$ARI_group_mean,
    m$NVI_group_mean,
    m$Cov_group_mean,
    m$F1_group_mean
  ))
  
  cat(sprintf("\n  %-12s %8s %8s %8s %8s\n", "SERIES", "ARI", "NVI", "Cover", "F1"))
  cat("  ", strrep("-", 46), "\n", sep = "")
  
  for (j in seq_len(J)) {
    g <- if (!is.null(sim_m$cluster_assign)) {
      sim_m$cluster_assign[j]
    } else {
      which(sapply(sim_m$group_members, function(idx) j %in% idx))[1]
    }
    
    cat(sprintf(
      "  S%-2d (G%d)    %7.3f %7.3f %7.3f %7.3f\n",
      j,
      g,
      m$ARI_series[j],
      m$NVI_series[j],
      m$Cov_series[j],
      m$F1_series[j]
    ))
  }
  
  cat(sprintf(
    "  %-12s %7.3f %7.3f %7.3f %7.3f\n",
    "MEAN",
    m$ARI_series_mean,
    m$NVI_series_mean,
    m$Cov_series_mean,
    m$F1_series_mean
  ))
  
  cat(sprintf("\n  %-12s %8s %8s\n", "SERIES", "ECP", "RISE"))
  cat("  ", strrep("-", 30), "\n", sep = "")
  
  for (j in seq_len(J)) {
    cat(sprintf(
      "  S%-11d %7.3f %7.3f\n",
      j,
      m$ECP[j],
      m$RISE[j]
    ))
  }
  
  cat(sprintf(
    "  %-12s %7.3f %7.3f\n",
    "MEAN",
    m$ECP_mean,
    m$RISE_mean
  ))
  
  cat("\n  Shape Classification Accuracy:\n")
  
  if (!is.null(m$SCA$overall)) {
    cat(sprintf(
      "    Overall: d1=%.3f  d2=%.3f  both=%.3f\n",
      ifelse(is.na(m$SCA$overall_d1), NA_real_, m$SCA$overall_d1),
      ifelse(is.na(m$SCA$overall_d2), NA_real_, m$SCA$overall_d2),
      ifelse(is.na(m$SCA$overall), NA_real_, m$SCA$overall)
    ))
    
    for (g in seq_len(C)) {
      pg <- m$SCA$per_group[[g]]
      
      if (isTRUE(pg$K_match)) {
        cat(sprintf(
          "    G%d: d1=%.3f  d2=%.3f  both=%.3f\n",
          g,
          pg$d1,
          pg$d2,
          pg$both
        ))
      } else {
        cat(sprintf(
          "    G%d: K mismatch (K_true=%d, K_hat=%d) — skipped\n",
          g,
          m$K_true[g],
          m$K_hat[g]
        ))
      }
    }
  }
  
  cat(sprintf("\n  F1 margin = %d time points\n", f1_margin))
  cat(strrep("=", 70), "\n", sep = "")
}


ppm_prepare_ydata <- function(sim) {
  if (!is.null(sim$Y)) {
    Y <- sim$Y
  } else if (!is.null(sim$truth$Y_std)) {
    Y <- sim$truth$Y_std
  } else {
    stop("Cannot find standardized data: expected sim$Y or sim$truth$Y_std.")
  }
  
  Y <- as.matrix(Y)
  
  T_len <- if (!is.null(sim$truth$T_len)) sim$truth$T_len else ncol(Y)
  
  # Latest simulate_shape_driven returns J x T.
  # If a transposed matrix is accidentally supplied, fix it.
  if (nrow(Y) == T_len && ncol(Y) != T_len) {
    Y <- t(Y)
  }
  
  if (ncol(Y) != T_len) {
    stop("Y does not appear to have T_len columns.")
  }
  
  Y
}

ppm_group_members_from_sim <- function(sim) {
  tr <- sim$truth
  
  if (!is.null(tr$cluster_assign)) {
    C <- tr$C
    return(lapply(seq_len(C), function(g) which(tr$cluster_assign == g)))
  }
  
  J <- tr$J
  C <- tr$C
  base_size <- J %/% C
  
  lapply(seq_len(C), function(g) {
    start <- (g - 1L) * base_size + 1L
    end <- if (g < C) g * base_size else J
    start:end
  })
}

ppm_make_eb_thetas <- function(Y) {
  L <- nrow(Y)
  n <- ncol(Y)
  
  mltypes <- rep(1L, L)
  thetas <- matrix(NA_real_, nrow = L, ncol = 4)
  
  for (i in seq_len(L)) {
    yi <- as.numeric(Y[i, ])
    ni <- length(yi)
    
    m_hat <- mean(yi, na.rm = TRUE)
    s2_hat <- var(yi, na.rm = TRUE)
    
    if (!is.finite(s2_hat) || s2_hat < 1e-8) {
      s2_hat <- 1
    }
    
    acf_i <- tryCatch(
      acf(yi, lag.max = ni - 1L, plot = FALSE)$acf[-1],
      error = function(e) rep(NA_real_, ni - 1L)
    )
    
    pos_lags <- which(is.finite(acf_i) & acf_i > 0)
    
    if (length(pos_lags) == 0L) {
      c_i_ell <- 0.01
    } else {
      ell <- pos_lags[1]
      c_i_ell <- acf_i[ell]
      c_i_ell <- min(max(c_i_ell, 0.01), 0.99)
    }
    
    kappa_i0 <- (1 - c_i_ell) / c_i_ell
    kappa_i0 <- min(max(kappa_i0, 1e-4), 1e4)
    
    d_hat <- mean((yi - m_hat)^4, na.rm = TRUE) / (s2_hat^2)
    
    if (!is.finite(d_hat)) {
      d_hat <- 3
    }
    
    # Pearson kurtosis. Keep it proper and numerically stable.
    d_hat <- min(max(d_hat, 2.1), 50)
    
    mu_i0    <- m_hat
    alpha_i0 <- 0.5 * d_hat
    beta_i0  <- 0.5 * d_hat * (1 - c_i_ell) * s2_hat
    beta_i0  <- max(beta_i0, 1e-8)
    
    thetas[i, ] <- c(mu_i0, kappa_i0, alpha_i0, beta_i0)
  }
  
  list(
    mltypes = mltypes,
    thetas = thetas
  )
}

ppm_make_ccp_prior <- function(L, n, nu0 = 3, rho = 0.5) {
  m0 <- rep(1 / n, L)
  mu0 <- log(m0 / (1 - m0))
  
  sigma0_var <- (1 / n * (1 - 1 / n)) / n
  
  S0 <- matrix(sigma0_var * rho, nrow = L, ncol = L)
  diag(S0) <- sigma0_var
  
  D_inv <- diag(1 / (m0 * (1 - m0)), L)
  
  Sigma0_adj <- ((nu0 - 2) / nu0) * (D_inv %*% S0 %*% D_inv)
  
  list(
    nu0 = nu0,
    mu0 = mu0,
    sigma0 = Sigma0_adj
  )
}

ppm_extract_prob_list <- function(fit, L, n) {
  if (is.null(fit$C)) {
    stop("ppmSuite fit object does not contain $C.")
  }
  
  Cmat <- as.matrix(fit$C)
  expected <- L * (n - 1L)
  
  if (ncol(Cmat) == expected) {
    prob <- colMeans(Cmat, na.rm = TRUE)
  } else if (nrow(Cmat) == expected) {
    prob <- rowMeans(Cmat, na.rm = TRUE)
  } else {
    stop(
      sprintf(
        "Unexpected dimension of fit$C: dim = %s, expected one dimension to equal L*(n-1) = %d.",
        paste(dim(Cmat), collapse = " x "),
        expected
      )
    )
  }
  
  split(as.numeric(prob), rep(seq_len(L), each = n - 1L))
}

ppm_est_cps_from_prob <- function(prob, T_len, threshold = 0.5,
                                  method = c("threshold", "topK"),
                                  K_truth = NULL) {
  method <- match.arg(method)
  
  prob <- as.numeric(prob)
  
  if (method == "threshold") {
    cps <- which(prob > threshold) + 1L
  } else {
    if (is.null(K_truth)) {
      stop("K_truth is required when method = 'topK'.")
    }
    
    n_cp <- max(as.integer(K_truth) - 1L, 0L)
    
    if (n_cp == 0L) {
      cps <- integer(0)
    } else {
      idx <- order(prob, decreasing = TRUE)[seq_len(min(n_cp, length(prob)))]
      cps <- sort(idx + 1L)
    }
  }
  
  ppm_normalize_cps(cps, T_len)
}

ppm_fit_icp <- function(Y, nburn = 10000, nskip = 1, nsave = 10000,
                        a0 = 1, b0 = 20, verbose = TRUE) {
  L <- nrow(Y)
  n <- ncol(Y)
  
  eb <- ppm_make_eb_thetas(Y)
  
  fit <- icp_ppm(
    ydata = Y,
    a0 = a0,
    b0 = b0,
    mltypes = eb$mltypes,
    thetas = eb$thetas,
    nburn = nburn,
    nskip = nskip,
    nsave = nsave,
    verbose = verbose
  )
  
  prob_list <- ppm_extract_prob_list(fit, L = L, n = n)
  
  list(
    fit = fit,
    prob_list = prob_list,
    eb = eb
  )
}

ppm_fit_ccp_by_group <- function(Y, group_members,
                                 nburn = 10000, nskip = 1, nsave = 10000,
                                 nu0 = 3, rho = 0.5, dev_sd = 0.005,
                                 verbose = TRUE) {
  L_all <- nrow(Y)
  n <- ncol(Y)
  
  prob_list <- vector("list", L_all)
  fits <- vector("list", length(group_members))
  priors <- vector("list", length(group_members))
  ebs <- vector("list", length(group_members))
  
  for (g in seq_along(group_members)) {
    idx <- group_members[[g]]
    
    Y_g <- Y[idx, , drop = FALSE]
    L_g <- nrow(Y_g)
    
    message(sprintf("  CCP-PPM group %d: L = %d, n = %d", g, L_g, n))
    
    eb_g <- ppm_make_eb_thetas(Y_g)
    prior_g <- ppm_make_ccp_prior(L = L_g, n = n, nu0 = nu0, rho = rho)
    
    devs_g <- matrix(dev_sd, nrow = L_g, ncol = n - 1L)
    
    fit_g <- ccp_ppm(
      ydata = Y_g,
      model = 1,
      nu0 = prior_g$nu0,
      mu0 = prior_g$mu0,
      sigma0 = prior_g$sigma0,
      mltypes = eb_g$mltypes,
      thetas = eb_g$thetas,
      nburn = nburn,
      nskip = nskip,
      nsave = nsave,
      devs = devs_g,
      verbose = verbose
    )
    
    prob_g <- ppm_extract_prob_list(fit_g, L = L_g, n = n)
    
    for (ii in seq_along(idx)) {
      prob_list[[idx[ii]]] <- prob_g[[ii]]
    }
    
    fits[[g]] <- fit_g
    priors[[g]] <- prior_g
    ebs[[g]] <- eb_g
  }
  
  list(
    fit_by_group = fits,
    prob_list = prob_list,
    prior_by_group = priors,
    eb_by_group = ebs
  )
}

ppm_evaluate_prob_list <- function(prob_list,
                                   truth_lower,
                                   truth_upper,
                                   group_members,
                                   T_len,
                                   threshold = 0.5,
                                   margin = 12,
                                   cp_method = c("threshold", "topK"),
                                   K_true_group = NULL) {
  cp_method <- match.arg(cp_method)
  
  J <- length(prob_list)
  C <- length(group_members)
  
  if (length(truth_lower) != J) {
    stop("truth_lower must be a list of length J.")
  }
  
  if (length(truth_upper) != C) {
    stop("truth_upper must be a list of length C.")
  }
  
  est_lower <- vector("list", J)
  
  for (j in seq_len(J)) {
    K_truth_j <- length(truth_lower[[j]]) + 1L
    
    est_lower[[j]] <- ppm_est_cps_from_prob(
      prob = prob_list[[j]],
      T_len = T_len,
      threshold = threshold,
      method = cp_method,
      K_truth = K_truth_j
    )
  }
  
  series_df <- data.frame(
    series = seq_len(J),
    group = NA_integer_,
    K_true = NA_integer_,
    K_hat = NA_integer_,
    ARI = NA_real_,
    Haus = NA_real_,
    Cover = NA_real_,
    NVI = NA_real_,
    F1 = NA_real_,
    Precision = NA_real_,
    Recall = NA_real_,
    TP = NA_integer_,
    FP = NA_integer_,
    FN = NA_integer_
  )
  
  for (j in seq_len(J)) {
    g <- which(vapply(group_members, function(idx) j %in% idx, logical(1)))
    if (length(g) == 0L) g <- NA_integer_ else g <- g[1]
    
    true_cps <- ppm_normalize_cps(truth_lower[[j]], T_len)
    est_cps  <- ppm_normalize_cps(est_lower[[j]], T_len)
    
    seg_true <- ppm_cps_to_segvec(true_cps, T_len)
    seg_est  <- ppm_cps_to_segvec(est_cps, T_len)
    
    f1 <- ppm_metric_f1(true_cps, est_cps, margin = margin)
    
    series_df$group[j] <- g
    series_df$K_true[j] <- length(true_cps) + 1L
    series_df$K_hat[j] <- length(est_cps) + 1L
    series_df$ARI[j] <- ppm_metric_ari(seg_true, seg_est)
    series_df$Haus[j] <- ppm_metric_hausdorff(true_cps, est_cps)
    series_df$Cover[j] <- ppm_metric_covering(seg_true, seg_est)
    series_df$NVI[j] <- ppm_metric_nvi(seg_true, seg_est)
    series_df$F1[j] <- f1$f1
    series_df$Precision[j] <- f1$precision
    series_df$Recall[j] <- f1$recall
    series_df$TP[j] <- f1$TP
    series_df$FP[j] <- f1$FP
    series_df$FN[j] <- f1$FN
  }
  
  est_upper <- vector("list", C)
  
  for (g in seq_len(C)) {
    idx <- group_members[[g]]
    
    prob_mat <- do.call(
      rbind,
      lapply(idx, function(j) prob_list[[j]])
    )
    
    prob_g <- colMeans(prob_mat, na.rm = TRUE)
    
    K_truth_g <- if (is.null(K_true_group)) {
      length(truth_upper[[g]]) + 1L
    } else {
      K_true_group[g]
    }
    
    est_upper[[g]] <- ppm_est_cps_from_prob(
      prob = prob_g,
      T_len = T_len,
      threshold = threshold,
      method = cp_method,
      K_truth = K_truth_g
    )
  }
  
  group_df <- data.frame(
    group = seq_len(C),
    K_true = NA_integer_,
    K_hat = NA_integer_,
    ARI = NA_real_,
    Haus = NA_real_,
    Cover = NA_real_,
    NVI = NA_real_,
    F1 = NA_real_,
    Precision = NA_real_,
    Recall = NA_real_,
    TP = NA_integer_,
    FP = NA_integer_,
    FN = NA_integer_
  )
  
  for (g in seq_len(C)) {
    true_cps <- ppm_normalize_cps(truth_upper[[g]], T_len)
    est_cps  <- ppm_normalize_cps(est_upper[[g]], T_len)
    
    seg_true <- ppm_cps_to_segvec(true_cps, T_len)
    seg_est  <- ppm_cps_to_segvec(est_cps, T_len)
    
    f1 <- ppm_metric_f1(true_cps, est_cps, margin = margin)
    
    group_df$K_true[g] <- length(true_cps) + 1L
    group_df$K_hat[g] <- length(est_cps) + 1L
    group_df$ARI[g] <- ppm_metric_ari(seg_true, seg_est)
    group_df$Haus[g] <- ppm_metric_hausdorff(true_cps, est_cps)
    group_df$Cover[g] <- ppm_metric_covering(seg_true, seg_est)
    group_df$NVI[g] <- ppm_metric_nvi(seg_true, seg_est)
    group_df$F1[g] <- f1$f1
    group_df$Precision[g] <- f1$precision
    group_df$Recall[g] <- f1$recall
    group_df$TP[g] <- f1$TP
    group_df$FP[g] <- f1$FP
    group_df$FN[g] <- f1$FN
  }
  
  summary_df <- data.frame(
    level = c("group", "series"),
    ARI_mean = c(mean(group_df$ARI), mean(series_df$ARI)),
    Haus_mean = c(ppm_mean_finite(group_df$Haus), ppm_mean_finite(series_df$Haus)),
    Cover_mean = c(mean(group_df$Cover), mean(series_df$Cover)),
    NVI_mean = c(mean(group_df$NVI), mean(series_df$NVI)),
    F1_mean = c(mean(group_df$F1), mean(series_df$F1)),
    Precision_mean = c(mean(group_df$Precision), mean(series_df$Precision)),
    Recall_mean = c(mean(group_df$Recall), mean(series_df$Recall)),
    K_exact_rate = c(
      mean(group_df$K_true == group_df$K_hat),
      mean(series_df$K_true == series_df$K_hat)
    )
  )
  
  list(
    series = series_df,
    group = group_df,
    summary = summary_df,
    est_lower = est_lower,
    est_upper = est_upper
  )
}

ppm_print_eval <- function(method_name, eval_obj) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(sprintf("  %s evaluation\n", method_name))
  cat(strrep("=", 72), "\n", sep = "")
  
  cat("\n[GROUP LEVEL]\n")
  print(eval_obj$group, row.names = FALSE)
  
  cat("\n[SERIES LEVEL]\n")
  print(eval_obj$series, row.names = FALSE)
  
  cat("\n[SUMMARY]\n")
  print(eval_obj$summary, row.names = FALSE)
  
  invisible(NULL)
}


bcp_fit <- function(Y,
                    burnin = 10000,
                    mcmc   = 20000,
                    w0     = NULL,
                    p0     = 0.2,
                    d      = 10,
                    verbose = TRUE) {
  Y <- as.matrix(Y)
  L <- nrow(Y)          # number of series (J)
  n <- ncol(Y)          # series length (T_len)
  
  prob_list <- vector("list", L)
  fits      <- vector("list", L)
  
  for (j in seq_len(L)) {
    if (verbose) message(sprintf("  BCP series %d / %d", j, L))
    
    yj <- as.numeric(Y[j, ])
    
    fit_j <- bcp::bcp(
      y      = yj,
      w0     = w0,
      p0     = p0,
      d      = d,
      burnin = burnin,
      mcmc   = mcmc
    )
    
    pp <- as.numeric(fit_j$posterior.prob)   # length n, last entry NA
    
    # Keep the (n - 1) inter-point boundary probabilities.
    prob_j <- pp[seq_len(n - 1L)]
    prob_j[!is.finite(prob_j)] <- 0
    
    prob_list[[j]] <- prob_j
    fits[[j]]      <- fit_j
  }
  
  list(fit = fits, prob_list = prob_list)
}


.SUMMARY_COLS <- c("method","scenario","seed","level","group",
                   "ARI","Haus","Cover","NVI","F1","RISE","ECP",
                   "dK","K_hat","K_true","K_exact_rate",
                   "SA_d1","SA_d2","SA_both","SA_Kmatch")
.mk_summary_row <- function(method, scenario, seed, level, group = NA_integer_,
                            ARI = NA_real_, Haus = NA_real_, Cover = NA_real_,
                            NVI = NA_real_, F1 = NA_real_,
                            RISE = NA_real_, ECP = NA_real_,
                            dK = NA_real_, K_hat = NA_real_, K_true = NA_real_,
                            K_exact_rate = NA_real_,
                            SA_d1 = NA_real_, SA_d2 = NA_real_, SA_both = NA_real_,
                            SA_Kmatch = NA) {
  data.frame(method = method, scenario = scenario, seed = seed, level = level,
             group = group, ARI = ARI, Haus = Haus, Cover = Cover, NVI = NVI, F1 = F1,
             RISE = RISE, ECP = ECP, dK = dK, K_hat = K_hat, K_true = K_true,
             K_exact_rate = K_exact_rate, SA_d1 = SA_d1, SA_d2 = SA_d2,
             SA_both = SA_both, SA_Kmatch = SA_Kmatch,
             stringsAsFactors = FALSE)
}

extract_hierscp_summary <- function(m, sim_m, scenario, seed) {
  C  <- sim_m$C
  cl <- sim_m$cluster_assign
  Kj_true <- vapply(sim_m$tau_series, function(z) length(z) + 1L, integer(1))   
  Kj_hat  <- as.integer(m$K_hat)[cl]                                            
  dK_j    <- Kj_hat - Kj_true
  
  # ---- upper (cluster) block ----
  upper <- .mk_summary_row(
    "HierSCP", scenario, seed, level = "upper",
    ARI = m$ARI_group_mean, Haus = m$Haus_group_mean, Cover = m$Cov_group_mean,
    NVI = m$NVI_group_mean, F1 = m$F1_group_mean, RISE = mean(m$RISE, na.rm = TRUE), ECP = mean(m$ECP, na.rm = TRUE),
    dK = mean(m$K_hat - m$K_true), K_hat = mean(m$K_hat), K_true = mean(m$K_true),
    K_exact_rate = mean(m$K_hat == m$K_true)
  )
  
  # ---- series_all block (joint grain; carries per-rep regime-pooled SA) ----
  series_all <- .mk_summary_row(
    "HierSCP", scenario, seed, level = "series_all",
    ARI = m$ARI_series_mean, Haus = m$Haus_series_mean, Cover = m$Cov_series_mean,
    NVI = m$NVI_series_mean, F1 = m$F1_series_mean, RISE = m$RISE_mean, ECP = m$ECP_mean,
    dK = mean(dK_j), K_hat = mean(Kj_hat), K_true = mean(Kj_true),
    K_exact_rate = mean(Kj_hat == Kj_true),
    SA_d1 = if (!is.null(m$SCA)) m$SCA$overall_d1 else NA_real_,
    SA_d2 = if (!is.null(m$SCA)) m$SCA$overall_d2 else NA_real_,
    SA_both = if (!is.null(m$SCA)) m$SCA$overall else NA_real_
  )
  
  # ---- per-group lower-series block (carries per-group SA, K-match conditional) ----
  per_group <- do.call(rbind, lapply(seq_len(C), function(g) {
    idx   <- which(cl == g)
    sca_g <- if (!is.null(m$SCA)) m$SCA$per_group[[g]] else NULL
    km    <- isTRUE(sca_g$K_match)
    .mk_summary_row(
      "HierSCP", scenario, seed, level = "series_group", group = g,
      ARI = mean(m$ARI_series[idx]), Haus = .mean_finite(m$Haus_series[idx]),
      Cover = mean(m$Cov_series[idx]), NVI = mean(m$NVI_series[idx]), F1 = mean(m$F1_series[idx]),
      RISE = mean(m$RISE[idx], na.rm = TRUE), ECP = mean(m$ECP[idx], na.rm = TRUE),
      dK = mean(dK_j[idx]), K_hat = as.numeric(m$K_hat[g]),
      K_true = as.numeric(sim_m$K_true[g]),
      K_exact_rate = mean(Kj_hat[idx] == Kj_true[idx]),
      SA_d1 = if (km) sca_g$d1 else NA_real_,
      SA_d2 = if (km) sca_g$d2 else NA_real_,
      SA_both = if (km) sca_g$both else NA_real_,
      SA_Kmatch = km
    )
  }))
  
  rbind(upper, series_all, per_group)[, .SUMMARY_COLS]
}

extract_ppm_summary <- function(eval_obj, method_name, scenario, seed) {
  group_df  <- eval_obj$group
  series_df <- eval_obj$series
  if (is.null(group_df) || is.null(series_df))
    stop("extract_ppm_summary: eval_obj must contain $group and $series.")
  C <- nrow(group_df)
  
  upper <- .mk_summary_row(
    method_name, scenario, seed, level = "upper",
    ARI = mean(group_df$ARI, na.rm = TRUE), Haus = .mean_finite(group_df$Haus),
    Cover = mean(group_df$Cover, na.rm = TRUE),
    NVI = if (!is.null(group_df$NVI)) mean(group_df$NVI, na.rm = TRUE) else NA_real_,
    F1 = mean(group_df$F1, na.rm = TRUE),
    dK = mean(group_df$K_hat - group_df$K_true), K_hat = mean(group_df$K_hat),
    K_true = mean(group_df$K_true), K_exact_rate = mean(group_df$K_hat == group_df$K_true)
  )
  
  series_all <- .mk_summary_row(
    method_name, scenario, seed, level = "series_all",
    ARI = mean(series_df$ARI, na.rm = TRUE), Haus = .mean_finite(series_df$Haus),
    Cover = mean(series_df$Cover, na.rm = TRUE),
    NVI = if (!is.null(series_df$NVI)) mean(series_df$NVI, na.rm = TRUE) else NA_real_,
    F1 = mean(series_df$F1, na.rm = TRUE),
    dK = mean(series_df$K_hat - series_df$K_true), K_hat = mean(series_df$K_hat),
    K_true = mean(series_df$K_true), K_exact_rate = mean(series_df$K_hat == series_df$K_true)
  )
  
  per_group <- do.call(rbind, lapply(seq_len(C), function(g) {
    idx <- which(series_df$group == g)
    if (length(idx) == 0L)
      return(.mk_summary_row(method_name, scenario, seed, level = "series_group", group = g))
    .mk_summary_row(
      method_name, scenario, seed, level = "series_group", group = g,
      ARI = mean(series_df$ARI[idx], na.rm = TRUE), Haus = .mean_finite(series_df$Haus[idx]),
      Cover = mean(series_df$Cover[idx], na.rm = TRUE),
      NVI = if (!is.null(series_df$NVI)) mean(series_df$NVI[idx], na.rm = TRUE) else NA_real_,
      F1 = mean(series_df$F1[idx], na.rm = TRUE),
      dK = mean(series_df$K_hat[idx] - series_df$K_true[idx]),
      K_hat = mean(series_df$K_hat[idx]), K_true = mean(series_df$K_true[idx]),
      K_exact_rate = mean(series_df$K_hat[idx] == series_df$K_true[idx])
    )
  }))
  
  rbind(upper, series_all, per_group)[, .SUMMARY_COLS]
}
ppm_results_to_combined <- function(results_list,
                                    method_map = list(metrics_icp = "ICP-PPM",
                                                      metrics_ccp = "CCP-PPM")) {
  do.call(rbind, lapply(results_list, function(res) {
    do.call(rbind, lapply(names(method_map), function(mk) {
      ev <- res[[mk]]; if (is.null(ev)) return(NULL)
      extract_ppm_summary(ev, method_map[[mk]],
                          scenario = res$scenario, seed = res$seed)
    }))
  }))
}
summarize_dK <- function(combined_df, level = "series_all") {
  d <- combined_df[combined_df$level == level, , drop = FALSE]
  do.call(rbind, lapply(sort(unique(d$method)), function(mm) {
    v  <- d$dK[d$method == mm]
    er <- d$K_exact_rate[d$method == mm]
    data.frame(method = mm,
               dK_mean = mean(v, na.rm = TRUE),
               dK_sd   = sd(v, na.rm = TRUE),
               dK_median = median(v, na.rm = TRUE),
               exact_recovery_rate = mean(er, na.rm = TRUE),
               n_rep = sum(is.finite(v)),
               stringsAsFactors = FALSE)
  }))
}
summarize_hierscp_shape_rise <- function(df_hierscp) {
  sg <- df_hierscp[df_hierscp$level == "series_group", , drop = FALSE]
  groups <- sort(unique(sg$group))
  
  matched <- sg[which(sg$SA_Kmatch %in% TRUE), , drop = FALSE]   # K-matched (g, rep)
  
  # ---- per-group SA (mean over matched reps; K_g fixed in-group => pooling = mean) ----
  sa_by_group <- do.call(rbind, lapply(groups, function(g) {
    r <- matched[matched$group == g, , drop = FALSE]
    n_total <- sum(sg$group == g)
    data.frame(group = paste0("G", g),
               K = if (nrow(r) > 0) r$K_true[1] else NA_real_,
               K_match_rate = nrow(r) / max(n_total, 1L),
               n_Kmatch = nrow(r),
               SA_d1   = if (nrow(r) > 0) mean(r$SA_d1)   else NA_real_,
               SA_d2   = if (nrow(r) > 0) mean(r$SA_d2)   else NA_real_,
               SA_both = if (nrow(r) > 0) mean(r$SA_both) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  
  # ---- regime-pooled overall (weight by K_g => sum matched regimes / sum K_g) ----
  if (nrow(matched) > 0) {
    w <- matched$K_true
    overall <- data.frame(group = "pooled", K = NA_real_,
                          K_match_rate = nrow(matched) / max(nrow(sg), 1L),
                          n_Kmatch = nrow(matched),
                          SA_d1   = sum(matched$SA_d1   * w) / sum(w),
                          SA_d2   = sum(matched$SA_d2   * w) / sum(w),
                          SA_both = sum(matched$SA_both * w) / sum(w),
                          stringsAsFactors = FALSE)
  } else {
    overall <- data.frame(group = "pooled", K = NA_real_, K_match_rate = 0,
                          n_Kmatch = 0L, SA_d1 = NA_real_, SA_d2 = NA_real_,
                          SA_both = NA_real_, stringsAsFactors = FALSE)
  }
  shape_accuracy <- rbind(sa_by_group, overall)
  
  # ---- RISE: per-group + overall series-mean (across reps) ----
  rise_by_group <- do.call(rbind, lapply(groups, function(g) {
    data.frame(group = paste0("G", g),
               RISE = mean(sg$RISE[sg$group == g], na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rise_series_mean <- mean(df_hierscp$RISE[df_hierscp$level == "series_all"], na.rm = TRUE)
  
  list(shape_accuracy   = shape_accuracy,
       rise_by_group     = rise_by_group,
       rise_series_mean  = rise_series_mean)
}
plot_method_boxplots <- function(combined_df,
                                 level        = "series_all",
                                 group_filter = NULL,
                                 metrics      = c("dK", "ARI", "Cover", "NVI", "F1"),
                                 method_order = c("HierSCP", "CCP-PPM", "ICP-PPM", "BCP"),
                                 cols         = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                                 drop_inf     = TRUE,
                                 jitter_amt   = 0.12,
                                 point_cex    = 1.0,
                                 main_prefix  = NULL,
                                 save_pdf     = NULL,
                                 width        = 10,
                                 height       = 8) {
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  d <- combined_df[combined_df$level == level, , drop = FALSE]
  if (!is.null(group_filter)) d <- d[d$group %in% group_filter, , drop = FALSE]
  d <- d[d$method %in% method_order, , drop = FALSE]
  d$method <- factor(d$method, levels = method_order)
  
  metric_labels <- c(
    dK    = "Signed segment-count error  (K_hat - K)",
    ARI   = "Adjusted Rand Index",
    Haus  = "Hausdorff distance",
    Cover = "Covering metric",
    NVI   = "Normalized variation of information  (lower = better)",
    F1    = "F1 score (margin = 10 timepoints)"
  )
  
  if (!is.null(save_pdf)) pdf(save_pdf, width = width, height = height)
  
  n_panel <- length(metrics)
  nr <- ceiling(sqrt(n_panel)); nc <- ceiling(n_panel / nr)
  op <- par(mfrow = c(nr, nc), mar = c(4.5, 4.4, 3.2, 1.2), mgp = c(2.5, 0.8, 0))
  on.exit({ par(op); if (!is.null(save_pdf)) dev.off() }, add = TRUE)
  
  for (met in metrics) {
    vals <- split(d[[met]], d$method)[method_order]
    names(vals) <- method_order
    if (drop_inf) vals <- lapply(vals, function(v) v[is.finite(v)])
    
    finite_all <- unlist(vals)
    yr <- if (length(finite_all) > 0L) range(finite_all, na.rm = TRUE) else c(0, 1)
    if (!all(is.finite(yr))) yr <- c(0, 1)
    if (met == "dK") yr <- range(c(yr, 0))            # always show the 0 line
    pad <- 0.06 * diff(yr); if (pad == 0) pad <- 0.05
    
    ttl <- metric_labels[[met]] %||% met
    if (!is.null(main_prefix)) ttl <- paste0(main_prefix, ": ", ttl)
    
    boxplot(vals, names = method_order, outline = FALSE,
            col = adjustcolor(cols, alpha.f = 0.30), border = cols,
            lwd = 1.8, boxwex = 0.55, ylim = c(yr[1] - pad, yr[2] + pad),
            main = ttl, ylab = met, xlab = "", las = 1, cex.axis = 0.95)
    
    if (met == "dK")
      abline(h = 0, lty = 2, lwd = 1.6, col = "grey35")   # over/under-segmentation reference
    
    for (i in seq_along(method_order)) {
      v <- vals[[method_order[i]]]
      if (length(v) > 0L) {
        xx <- jitter(rep(i, length(v)), amount = jitter_amt)
        points(xx, v, pch = 21, bg = adjustcolor(cols[i], alpha.f = 0.9),
               col = "grey20", cex = point_cex)
      }
    }
  }
  invisible(NULL)
}


run_one_hierscp_rep <- function(seed, scenario,
                                # simulation settings
                                sim_J = 9L, sim_T_len = 120L, sim_C = 3L,
                                sim_cp_upper = list(c(50L), c(42L), c(40L, 85L)),
                                sim_M_basis = 8L, sim_m_min = 10L,
                                # model-spec settings
                                model_M = 8L, model_m_min = 10L,
                                model_K_init = 2, model_K_min = 2, model_K_max = 6,
                                # MCMC settings
                                n_iter = 20000L, n_burnin = 10000L,
                                n_thin = 1L, n_warmup = 100L, ia_every = 20L,
                                N_launch = 2L, verbose_mcmc = 1L,
                                f1_margin = 10, raw_path = NULL) {
  
  sim <- simulate_shape_driven(
    scenario  = scenario,
    J         = sim_J,
    T_len     = sim_T_len,
    C         = sim_C,
    cp_upper  = sim_cp_upper,
    M_basis   = sim_M_basis,
    m_min     = sim_m_min,
    seed      = seed,
    return_scale = "standardized",
    attach_model_objects = FALSE
  )
  
  sim_m <- adapt_sim_for_metrics(sim)
  gm    <- sim_m$group_members
  
  model <- create_model_spec(
    sim$Y, C = sim_C,
    K_init = model_K_init, K_max = model_K_max, K_min = model_K_min,
    M = model_M, m_min = model_m_min, type = "continuous",
    fixed_clusters = TRUE, fixed_cluster_sizes = TRUE,
    group_members = gm
  )
  
  result <- run_mcmc(
    sim$Y, model,
    n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin,
    n_warmup = n_warmup, ia_every = ia_every,
    N_launch = N_launch, verbose = verbose_mcmc
  )
  
  post <- compute_conditional_posterior_summary(result, sim$Y, verbose = FALSE)
  if (!is.null(raw_path))
    saveRDS(list(sim = sim, sim_m = sim_m, post = post, model = model,
                 scenario = scenario, seed = seed), raw_path)
  
  m <- compute_all_metrics(post, sim, f1_margin = f1_margin, verbose = FALSE)
  extract_hierscp_summary(m, sim_m, scenario = scenario, seed = seed)
}



run_hierscp_replications <- function(scenario, seeds,
                                     workers    = 3L,
                                     out_folder = "S3_replicate_folder",
                                     ...) {
  
  stopifnot(length(seeds) >= 1L)
  workers <- max(1L, min(as.integer(workers), length(seeds)))
  dots    <- list(...)
  go_par  <- workers > 1L
  
  ## ── 0) 마스터 점검: 구버전/누락이면 모호한 에러 대신 즉시 중단 ──────────
  need <- c("run_one_hierscp_rep", "simulate_shape_driven",
            "adapt_sim_for_metrics", "create_model_spec", "run_mcmc",
            "compute_conditional_posterior_summary", "compute_all_metrics",
            "extract_hierscp_summary")
  miss <- need[!vapply(need, exists, logical(1), mode = "function")]
  if (length(miss)) stop("마스터에 함수 없음: ", paste(miss, collapse = ", "))
  if (!"scenario" %in% names(formals(simulate_shape_driven)))
    stop("simulate_shape_driven 이 구버전(scenario 인자 없음).")
  
  ## ── 1) 경로를 절대경로로 고정 (워커 wd 가 달라도 안전) ──────────────────
  ##     setwd() 는 이미 됐다는 전제. 이름이 틀리면 여기서 바로 에러난다.
  cpp_abs    <- normalizePath(.HIERSCP_CPP,    mustWork = TRUE)
  bridge_abs <- normalizePath(.HIERSCP_BRIDGE, mustWork = TRUE)
  out_dir    <- file.path(normalizePath(getwd(), winslash = "/"), out_folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## ── 2) 외부포인터 직렬화 에러 차단 ──────────────────────────────────────
  old <- options(future.globals.onReference = "ignore")
  on.exit(options(old), add = TRUE)
  
  ## ── 3) 마스터의 모든 함수(점-헬퍼 포함) + 필요 값들을 워커로 ────────────
  allnm  <- ls(globalenv(), all.names = TRUE)
  is_fun <- vapply(allnm, function(nm)
    is.function(get0(nm, globalenv(), inherits = FALSE)), logical(1))
  ## NON-function top-level constants that exported functions reference by lexical
  ## scope (e.g. .SUMMARY_COLS). WITHOUT exporting these, every worker throws
  ## "object '.SUMMARY_COLS' not found" AFTER the (expensive) MCMC — losing the run.
  is_const <- vapply(allnm, function(nm) {
    obj <- get0(nm, globalenv(), inherits = FALSE)
    !is.null(obj) && is.atomic(obj) && !is.function(obj)
  }, logical(1))
  G <- mget(allnm[is_fun | is_const], envir = globalenv(), inherits = FALSE)
  G[["scenario"]]    <- scenario
  G[[".dots"]]       <- dots
  G[["cpp_abs"]]     <- cpp_abs
  G[["bridge_abs"]]  <- bridge_abs
  G[["out_dir"]]     <- out_dir
  G[[".HIER_READY"]] <- NULL    # 워커가 C++ 셋업을 건너뛰지 않도록 제거
  
  if (go_par) future::plan(future::multisession, workers = workers)
  else        future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  message(sprintf("[HierSCP] scenario=%s | seeds=%d | %s | out=%s",
                  scenario, length(seeds),
                  if (go_par) sprintf("multisession · %d workers", workers) else "sequential",
                  out_folder))
  
  run_seed <- function(sd) {
    ## (워커 1회) C++ 컴파일 + 브리지 연결 — 컴파일 잡음은 로그에 안 남긴다
    if (!exists(".HIER_READY", envir = globalenv())) {
      Rcpp::sourceCpp(cpp_abs)
      source(bridge_abs, local = FALSE)
      assign(".HIER_READY", TRUE, envir = globalenv())
    }
    
    log_path <- file.path(out_dir, sprintf("%s_seed%d.log", scenario, sd))
    out_path <- file.path(out_dir, sprintf("%s_seed%d.rds", scenario, sd))
    out_path_raw <- file.path(out_dir, sprintf("%s_seed%d_raw.rds", scenario, sd))
    
    con   <- file(log_path, open = "wt")
    stamp <- function(txt)
    { cat(sprintf("[%s] seed %d | %s\n",
                  format(Sys.time(), "%H:%M:%S"), sd, txt), file = con); flush(con) }
    
    t0 <- Sys.time()
    stamp(sprintf("START  scenario=%s", scenario))
    
    ## MCMC 등 stdout 출력을 같은 로그로 캡처(verbose_mcmc>=1 이면 반복 진행 보임).
    ## 출력 스트림은 '스택'이라 내가 추가한 1개만 정확히 되돌린다(future 와 안 충돌).
    sink(con, type = "output")
    on.exit({ sink(); flush(con); close(con) }, add = TRUE)
    
    out <- tryCatch(
      do.call(run_one_hierscp_rep,
              c(list(seed = sd, scenario = scenario, raw_path = out_path_raw), .dots)),
      error = function(e) { stamp(sprintf("FAILED: %s", conditionMessage(e))); NULL }
    )
    
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    if (!is.null(out)) {
      saveRDS(out, out_path)
      stamp(sprintf("DONE  %.1fs  ->  %s", el, basename(out_path)))
    } else if (file.exists(out_path_raw)) {
      stamp(sprintf("SUMMARY FAILED but RAW FIT SAVED -> %s (%.1fs); rebuild via resummarize_hierscp_folder(out_dir)",
                    basename(out_path_raw), el))
    } else {
      stamp(sprintf("ABORTED before MCMC finished  %.1fs", el))
    }
    out
  }
  
  res <- future.apply::future_lapply(
    seeds, run_seed,
    future.seed    = TRUE,
    future.globals = G
  )
  
  do.call(rbind, Filter(Negate(is.null), res))
}

resummarize_hierscp_folder <- function(folder, pattern = "_raw\\.rds$", f1_margin = 10) {
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  rows <- lapply(files, function(fp) {
  r <- readRDS(fp)
    
  m <- compute_all_metrics(r$post, r$sim, f1_margin = f1_margin, verbose = FALSE)
  extract_hierscp_summary(m, r$sim_m, scenario = r$scenario, seed = r$seed)
  })
  
  do.call(rbind, Filter(Negate(is.null), rows))
}


run_ppmSuite_hierscp_benchmark <- function(
    scenario = "S1_in_model",
    n_rep = 3L,
    sim_seeds = 2026:(2026 + n_rep - 1L),
    
    J = 9L,
    T_len = 120L,
    C = 3L,
    M_basis = 8L,
    m_min = 10L,
    cp_upper = list(c(50L), c(42L), c(40L, 85L)),
    
    nburn = 20000,
    nskip = 1,
    nsave = 10000,
    
    prob_threshold = 0.5,
    f1_margin = 10,
    cp_method = c("threshold", "topK"),
    
    icp_a0 = 1,
    icp_b0 = 20,
    
    ccp_nu0 = 3,
    ccp_rho = 0.5,
    ccp_dev_sd = 0.005,
    
    verbose_ppm = TRUE,
    print_each_rep = TRUE
) {
  cp_method <- match.arg(cp_method)
  
  if (length(sim_seeds) < n_rep) {
    stop("sim_seeds must have length at least n_rep.")
  }
  
  results <- vector("list", n_rep)
  
  for (r in seq_len(n_rep)) {
    current_seed <- sim_seeds[r]
    
    cat("\n", strrep("#", 80), "\n", sep = "")
    cat(sprintf("Replication %d/%d | scenario = %s | seed = %d\n",
                r, n_rep, scenario, current_seed))
    cat(strrep("#", 80), "\n", sep = "")
    
    sim <- simulate_shape_driven(
      scenario = scenario,
      J = J,
      T_len = T_len,
      C = C,
      cp_upper = cp_upper,
      M_basis = M_basis,
      m_min = m_min,
      seed = current_seed,
      return_scale = "standardized",
      attach_model_objects = FALSE
    )
    
    Y <- ppm_prepare_ydata(sim)
    truth <- sim$truth
    
    group_members <- ppm_group_members_from_sim(sim)
    
    truth_lower <- truth$cp_lower
    truth_upper <- truth$cp_upper
    K_true_group <- truth$K_by_cluster
    
    cat(sprintf("Data dimension for ppmSuite: L = %d, n = %d\n", nrow(Y), ncol(Y)))
    cat(sprintf("True upper K: [%s]\n", paste(K_true_group, collapse = ", ")))
    
    t0 <- Sys.time()
    
    cat("\n--- Fitting ICP-PPM ---\n")
    fit_icp <- ppm_fit_icp(
      Y = Y,
      nburn = nburn,
      nskip = nskip,
      nsave = nsave,
      a0 = icp_a0,
      b0 = icp_b0,
      verbose = verbose_ppm
    )
    
    cat("\n--- Fitting CCP-PPM by true group ---\n")
    fit_ccp <- ppm_fit_ccp_by_group(
      Y = Y,
      group_members = group_members,
      nburn = nburn,
      nskip = nskip,
      nsave = nsave,
      nu0 = ccp_nu0,
      rho = ccp_rho,
      dev_sd = ccp_dev_sd,
      verbose = verbose_ppm
    )
    
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    
    cat(sprintf("\nppmSuite fits completed in %.1f min\n", elapsed))
    
    metrics_icp <- ppm_evaluate_prob_list(
      prob_list = fit_icp$prob_list,
      truth_lower = truth_lower,
      truth_upper = truth_upper,
      group_members = group_members,
      T_len = ncol(Y),
      threshold = prob_threshold,
      margin = f1_margin,
      cp_method = cp_method,
      K_true_group = K_true_group
    )
    
    metrics_ccp <- ppm_evaluate_prob_list(
      prob_list = fit_ccp$prob_list,
      truth_lower = truth_lower,
      truth_upper = truth_upper,
      group_members = group_members,
      T_len = ncol(Y),
      threshold = prob_threshold,
      margin = f1_margin,
      cp_method = cp_method,
      K_true_group = K_true_group
    )
    
    if (isTRUE(print_each_rep)) {
      ppm_print_eval("ICP-PPM", metrics_icp)
      ppm_print_eval("CCP-PPM group-wise", metrics_ccp)
    }
    
    results[[r]] <- list(
      rep = r,
      seed = current_seed,
      scenario = scenario,
      sim = sim,
      fit_icp = fit_icp,
      fit_ccp = fit_ccp,
      metrics_icp = metrics_icp,
      metrics_ccp = metrics_ccp,
      elapsed_min = elapsed
    )
  }
  
  summary_all <- do.call(
    rbind,
    lapply(results, function(res) {
      icp_sum <- res$metrics_icp$summary
      ccp_sum <- res$metrics_ccp$summary
      
      icp_sum$method <- "ICP-PPM"
      ccp_sum$method <- "CCP-PPM"
      
      icp_sum$rep <- res$rep
      ccp_sum$rep <- res$rep
      
      icp_sum$seed <- res$seed
      ccp_sum$seed <- res$seed
      
      rbind(icp_sum, ccp_sum)
    })
  )
  
  summary_all <- summary_all[
    ,
    c(
      "rep", "seed", "method", "level",
      "ARI_mean", "Haus_mean", "Cover_mean", "F1_mean",
      "Precision_mean", "Recall_mean", "K_exact_rate"
    )
  ]
  
  aggregate_summary <- aggregate(
    cbind(
      ARI_mean,
      Haus_mean,
      Cover_mean,
      F1_mean,
      Precision_mean,
      Recall_mean,
      K_exact_rate
    ) ~ method + level,
    data = summary_all,
    FUN = mean
  )
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("Final aggregate summary across replications\n")
  cat(strrep("=", 80), "\n", sep = "")
  print(aggregate_summary, row.names = FALSE)
  
  list(
    results = results,
    summary_by_rep = summary_all,
    aggregate_summary = aggregate_summary
  )
}

run_one_ppm_rep <- function(seed, scenario,
                            methods = c("BCP"),
                            # simulation settings (must match the HierSCP run)
                            sim_J = 9L, sim_T_len = 120L, sim_C = 3L,
                            sim_cp_upper = list(c(50L), c(42L), c(40L, 85L)),
                            sim_M_basis = 8L, sim_m_min = 10L,
                            # shared evaluation settings
                            threshold = 0.5, margin = 10,
                            cp_method = c("threshold", "topK"),
                            # ppmSuite MCMC settings
                            nburn = 10000, nskip = 1L, nsave = 10000,
                            icp_a0 = 1, icp_b0 = 20,
                            ccp_nu0 = 3, ccp_rho = 0.5, ccp_dev_sd = 0.005,
                            # BCP settings
                            bcp_burnin = 10000, bcp_mcmc = 20000,
                            bcp_w0 = NULL, bcp_p0 = 0.2, bcp_d = 10,
                            verbose_fit = FALSE) {
  cp_method <- match.arg(cp_method)
  
  sim <- simulate_shape_driven(
    scenario  = scenario,
    J         = sim_J,
    T_len     = sim_T_len,
    C         = sim_C,
    cp_upper  = sim_cp_upper,
    M_basis   = sim_M_basis,
    m_min     = sim_m_min,
    seed      = seed,
    return_scale = "standardized",
    attach_model_objects = FALSE
  )
  
  Y             <- ppm_prepare_ydata(sim)
  gm            <- ppm_group_members_from_sim(sim)
  truth         <- sim$truth
  truth_lower   <- truth$cp_lower
  truth_upper   <- truth$cp_upper
  K_true_group  <- truth$K_by_cluster
  T_len         <- ncol(Y)
  
  eval_one <- function(prob_list) {
    ppm_evaluate_prob_list(
      prob_list     = prob_list,
      truth_lower   = truth_lower,
      truth_upper   = truth_upper,
      group_members = gm,
      T_len         = T_len,
      threshold     = threshold,
      margin        = margin,
      cp_method     = cp_method,
      K_true_group  = K_true_group
    )
  }
  
  out <- list()
  
  if ("ICP" %in% methods) {
    fit <- ppm_fit_icp(Y, nburn = nburn, nskip = nskip, nsave = nsave,
                       a0 = icp_a0, b0 = icp_b0, verbose = verbose_fit)
    out$ICP <- extract_ppm_summary(eval_one(fit$prob_list),
                                   "ICP-PPM", scenario, seed)
  }
  
  if ("CCP" %in% methods) {
    fit <- ppm_fit_ccp_by_group(Y, group_members = gm,
                                nburn = nburn, nskip = nskip, nsave = nsave,
                                nu0 = ccp_nu0, rho = ccp_rho,
                                dev_sd = ccp_dev_sd, verbose = verbose_fit)
    out$CCP <- extract_ppm_summary(eval_one(fit$prob_list),
                                   "CCP-PPM", scenario, seed)
  }
  
  if ("BCP" %in% methods) {
    fit <- bcp_fit(Y, burnin = bcp_burnin, mcmc = bcp_mcmc,
                   w0 = bcp_w0, p0 = bcp_p0, d = bcp_d, verbose = verbose_fit)
    out$BCP <- extract_ppm_summary(eval_one(fit$prob_list),
                                   "BCP", scenario, seed)
  }
  
  do.call(rbind, out)
}


seeds <- 2026:2075 # 50 replication seeds
df_hierscp <- run_hierscp_replications("S3", seeds = seeds, workers = 3,
                                       n_iter = 20000, n_burnin = 10000,
                                       n_thin = 1, verbose_mcmc = 1,
                                       f1_margin = 10,out_folder = "S3_replicate_50rep")

df_hierscp$scenario <- "S3"
save(df_hierscp,file="S3_Hierscp_50rep_res.RData")


## (b) ICP-PPM + CCP-PPM (Important: CCP-PPM is fitted to series in the same cluster, not the whole data(9 series))

seeds <- 2026:2075
results_ppm <- run_ppmSuite_hierscp_benchmark(
  scenario = "S3", n_rep = 50, sim_seeds = seeds,
  nburn = 10000, nskip = 1, nsave = 10000,
  prob_threshold = 0.5, f1_margin = 10, cp_method = "threshold")
df_ppm <- ppm_results_to_combined(results_ppm$results)   # ICP-PPM + CCP-PPM, incl. dK

df_ppm

## (c) BCP replications (extract_ppm_summary carries dK)
df_bcp <- do.call(rbind, lapply(seeds, function(s)
  run_one_ppm_rep(seed = s, scenario = "S3", methods = "BCP", margin = 10)))

## (d) designated cols for binding
cols <- c("method", "scenario", "seed", "level",
          "ARI", "Cover", "NVI", "F1", "dK", "K_exact_rate")

#rbind(df_ppm[,cols], df_bcp[,cols])

save(rbind(df_ppm[,cols], df_bcp[,cols]),file="S3_PPM_BCP.RData")

combined_df <- rbind(df_hierscp[, cols], df_ppm[, cols], df_bcp[, cols])

plot_method_boxplots(combined_df, level = "series_all",
                     metrics = c("ARI", "Cover", "F1","dK","NVI"),  # settled set
                     save_pdf = "S3_50reps_boxplot.pdf")















