###############################################################################
# 06_local_adjusting.R
# Step 2: Local Adjusting — Direct Independence Sampling (DIS)
#
# CRITICAL INVARIANT: Upper and Lower CP counts are ALWAYS equal (= K-1).
# Local Adjusting only SLIDES the boundary position within the free
# observation set. It never creates or destroys boundaries.
#
# Therefore: NO "no transition" option. Every pair (k, k+1) MUST have
# a transition somewhere in G_L (or G_U). If G_L is empty, the boundary
# is already determined by the constraints and cannot move.
#
# DIS (2026-04-12): Replaces sequential hazard scan with direct categorical
# sampling from the exact full conditional (Lemma 2).
#
# For each candidate boundary position g in the free observation set:
#   log p(g) = [Truncated Geometric prior] + [Segment likelihood]
#
# Prior (Lemma 2):
#   (π*_k)^{g - g^-} × (1-π*_k)^{I(g ≠ g^+)}
#
# The indicator I(g ≠ g^+) means: at the rightmost position g^+,
# the transition is forced by the m_min constraint, so no (1-π*) penalty.
#
# Benefits over sequential scan:
#   - No path dependency: CP can jump to any position in one step
#   - Exact Gibbs: samples from true full conditional
#   - N_launch=1 sufficient (scan needed N_launch≥2 for adequate mixing)
#
# Two sub-steps:
#   (A) Upper-level: slide cluster-level CPs using ALL series in cluster
#   (B) Lower-level: slide series-level CPs within each series
###############################################################################

local_adjusting <- function(state, params, Y, precomp, model) {
  # (A) Upper-level
  for (cc in 1:model$C) {
    state <- upper_local_adjusting(cc, state, params, Y, precomp, model)
  }
  # (B) Lower-level
  for (j in 1:model$J) {
    state <- lower_local_adjusting_series(j, state, params, Y, precomp, model)
  }
  state
}

#' Lower-level Local Adjusting ONLY — upper CPs are fixed.
#' Called between IA rounds: only individual series CPs slide within G_L.
lower_only_local_adjusting <- function(state, params, Y, precomp, model) {
  for (j in 1:model$J) {
    state <- lower_local_adjusting_series(j, state, params, Y, precomp, model)
  }
  state
}

#' Upper-level Local Adjusting ONLY — refine upper CPs after IA.
#' Each cluster's upper boundaries slide pointwise within G_U.
#' Likelihood at each time point t aggregates ALL series in the cluster.
#' S_lower within G_U is synchronized to maintain structural consistency
#' (this is NOT the same as resetting all S_lower = S_upper).
upper_only_local_adjusting <- function(state, params, Y, precomp, model) {
  for (cc in 1:model$C) {
    state <- upper_local_adjusting(cc, state, params, Y, precomp, model)
  }
  state
}

#' [LEGACY] Sequential scan sampler — kept for reference.
#' Replaced by direct categorical sampling (DIS) via sample_categorical_log.
#' The scan was distributionally equivalent but path-dependent.
scan_sample_from_log_probs <- function(log_probs) {
  # DIS: direct categorical sampling replaces sequential scan
  sample_categorical_log(log_probs)
}

###############################################################################
# (A) UPPER-LEVEL LOCAL ADJUSTING
###############################################################################

upper_local_adjusting <- function(cc, state, params, Y, precomp, model) {
  K_c <- state$K[cc]
  if (K_c <= 1) return(state)
  for (k in 1:(K_c - 1)) {
    state <- upper_local_adjusting_pair(cc, k, state, params, Y, precomp, model)
  }
  state
}

#' Slide upper CP between states k and k+1 for cluster cc.
#' Boundary MUST land somewhere in G_U — no "no transition" option.
upper_local_adjusting_pair <- function(cc, k, state, params, Y, precomp, model) {
  m_min <- model$m_min; T_len <- model$T_len

  G_U <- compute_free_obs_upper(k, cc, state, m_min)
  if (length(G_U) == 0) return(state)
  n_free <- length(G_U)

  ps_vec <- params$pi_star[[cc]]
  pi_star_k <- if (k <= length(ps_vec)) ps_vec[k] else 0.5
  pi_star_k <- min(max(pi_star_k, 1e-10), 1 - 1e-10)
  log_pi   <- log(pi_star_k)
  log_1mpi <- log(1 - pi_star_k)

  # [ANCHOR] IA가 확정한 상위 변화점 τ^U_{c,k+1}
  #   Upper LA는 이 anchor 주변에서 상위 변화점을 재조정.
  tau_U_vec <- state$tau_upper[[cc]]
  anchor    <- if ((k + 1) <= length(tau_U_vec)) tau_U_vec[k + 1] else G_U[n_free]

  atoms_c <- params$atoms[[cc]]
  atom_k   <- if (k   <= length(atoms_c) && !is.null(atoms_c[[k]]))   atoms_c[[k]]   else list(gamma1=1, gamma2=1, shape_beta=1.0, shape_gamma=0.0)
  kp1_safe <- min(k + 1, length(atoms_c))
  atom_kp1 <- if (kp1_safe >= 1 && !is.null(atoms_c[[kp1_safe]])) atoms_c[[kp1_safe]] else list(gamma1=-1, gamma2=-1, shape_beta=1.0, shape_gamma=0.0)

  j_in_c <- which(state$cluster == cc)
  if (length(j_in_c) == 0) return(state)

  # ── Segment likelihood across ALL series in the cluster ──
  #   Normal log-likelihood with the series-level scale v_j INSIDE the
  #   variance component:  σ²_{j,t} = v_j · φ_{j,t}^{ξ_{j,t}}.
  #   This is the pre-collapse Normal kernel (the v_j-collapsed Student-t
  #   form previously used here was incorrect for this update).
  #   Fast path delegates to upper_la_pair_cpp(); the exact R fallback below
  #   uses the identical kernel.
  use_cpp <- isTRUE(get0(".hiercpd_rcpp_enabled", ifnotfound = FALSE)) &&
             exists("upper_la_pair_cpp")

  chosen <- NA_integer_
  if (use_cpp) {
    if (exists(".hiercpd_ensure_datx", mode = "function")) .hiercpd_ensure_datx(precomp)
    theta_mat <- params$theta[j_in_c, , drop = FALSE]
    # γ* is segment-invariant, so the k and k+1 intercept vectors coincide.
    gamma_vec <- vapply(j_in_c, function(jj) as.numeric(get_state_intercept(jj, k, params)), numeric(1))
    v_vec <- vapply(j_in_c, function(jj) {
      vj <- params$v[[jj]]
      if (is.null(vj) || length(vj) == 0 || !is.finite(vj[1]) || vj[1] <= 0) 1.0 else as.numeric(vj[1])
    }, numeric(1))
    Y_sub   <- Y[j_in_c, , drop = FALSE]
    phi_sub <- params$phi[j_in_c, , drop = FALSE]
    xi_sub  <- params$xi[j_in_c, , drop = FALSE]
    xi_sub[is.na(xi_sub)] <- 0L
    storage.mode(xi_sub) <- "integer"
    chosen <- tryCatch(
      upper_la_pair_cpp(
        theta_mat,
        gamma_vec, atom_k,
        gamma_vec, atom_kp1,
        as.integer(G_U),
        Y_sub, params$alpha, params$beta[cc],
        v_vec, v_vec,
        phi_sub, xi_sub,
        pi_star_k, as.integer(anchor),
        as.numeric(precomp$x_global)
      ),
      error = function(e) NA_integer_
    )
    if (length(chosen) != 1L || is.na(chosen) || chosen < 1L || chosen > n_free) chosen <- NA_integer_
  }

  if (is.na(chosen)) {
    # ── Exact R fallback: identical Normal kernel + anchor prior ──
    agg_ll_k   <- numeric(n_free)
    agg_ll_kp1 <- numeric(n_free)
    for (j in j_in_c) {
      theta_j <- params$theta[j, ]
      f_k_j   <- eval_shape_at_times(theta_j, atom_k,   G_U, precomp)
      f_kp1_j <- eval_shape_at_times(theta_j, atom_kp1, G_U, precomp)
      gamma_j <- get_state_intercept(j, k, params)   # γ* segment-invariant
      mu_k_j   <- params$alpha + params$beta[cc] + gamma_j + f_k_j
      mu_kp1_j <- params$alpha + params$beta[cc] + gamma_j + f_kp1_j
      v_j <- params$v[[j]]
      v_j <- if (is.null(v_j) || length(v_j) == 0 || !is.finite(v_j[1]) || v_j[1] <= 0) 1 else v_j[1]
      for (idx in 1:n_free) {
        t <- G_U[idx]
        phi_t <- params$phi[j, t]; if (!is.finite(phi_t) || phi_t < 1) phi_t <- 1
        xi_t  <- params$xi[j, t];  if (is.na(xi_t)) xi_t <- 0
        sig2 <- v_j * phi_t^xi_t; if (!is.finite(sig2) || sig2 <= 0) sig2 <- 1
        agg_ll_k[idx]   <- agg_ll_k[idx]   + dnorm_log_var(Y[j, t], mu_k_j[idx],   sig2)
        agg_ll_kp1[idx] <- agg_ll_kp1[idx] + dnorm_log_var(Y[j, t], mu_kp1_j[idx], sig2)
      }
    }

    cum_ll_k   <- cumsum(agg_ll_k)
    cum_ll_kp1 <- cumsum(agg_ll_kp1)
    total_ll_kp1 <- cum_ll_kp1[n_free]

    # [DIS] anchor 기준 prior over G_U 전체 (slack 없음)
    g_plus <- G_U[n_free]
    log_probs <- numeric(n_free)
    for (i in 1:n_free) {
      tau <- G_U[i]
      ll_before <- if (i > 1) cum_ll_k[i - 1] else 0
      ll_after  <- total_ll_kp1 - (if (i > 1) cum_ll_kp1[i - 1] else 0)
      log_prior <- (anchor - tau) * log_pi
      if (tau != (g_plus + 1L)) log_prior <- log_prior + log_1mpi
      log_probs[i] <- log_prior + ll_before + ll_after
    }
    chosen <- sample_categorical_log(log_probs)
  }

  # Update upper state sequence (G_U 전체)
  if (chosen > 1) state$S_upper[cc, G_U[1:(chosen - 1)]] <- as.integer(k)
  state$S_upper[cc, G_U[chosen:n_free]] <- as.integer(k + 1)
  state$tau_upper[[cc]] <- extract_changepoints(state$S_upper[cc, ])

  # [계층 구조] Upper LA 는 상위 변화점만 조정한다. 하위 변화점 S_lower 는 바로 뒤
  #   단계인 Lower LA 가 방금 갱신된 τ^U_{c,k+1} 을 anchor 로 삼아 series 별로 직접
  #   정한다(lower_local_adjusting_pair 의 anchor 사용).  이전 코드는 여기서 S_lower
  #   를 upper 경계로 강제 동기화했는데, 그것이 (i) Lower LA 의 역할을 침범하고
  #   (ii) series 별 m_min 을 보장하지 못해 C2 위반 → enforce_all_invariants 의
  #   강제병합(MCMC 틀 밖 개입)을 유발했다.  따라서 lower 동기화를 제거한다.

  state
}

###############################################################################
# (B) LOWER-LEVEL LOCAL ADJUSTING
# Slide τ^L_{j,k+1} within G^L_{(k,k+1),j}
# Boundary MUST exist — no "no transition" option.
###############################################################################

lower_local_adjusting_series <- function(j, state, params, Y, precomp, model) {
  c_j <- state$cluster[j]
  K_c <- state$K[c_j]
  if (K_c <= 1) return(state)
  for (k in 1:(K_c - 1)) {
    state <- lower_local_adjusting_pair(j, k, state, params, Y, precomp, model)
  }
  state
}

lower_local_adjusting_pair <- function(j, k, state, params, Y, precomp, model) {
  c_j <- state$cluster[j]
  m_min <- model$m_min

  G_L <- compute_free_obs_lower(j, k, state, m_min)
  if (length(G_L) == 0) return(state)
  n_free <- length(G_L)

  # [HDP-HMM] Lower-level π* (upper와 독립)
  ps_vec <- if (!is.null(params$pi_star_lower[[c_j]])) params$pi_star_lower[[c_j]] else params$pi_star[[c_j]]
  pi_star_k <- if (k <= length(ps_vec)) ps_vec[k] else 0.5
  pi_star_k <- min(max(pi_star_k, 1e-10), 1 - 1e-10)
  log_pi   <- log(pi_star_k)
  log_1mpi <- log(1 - pi_star_k)

  # [ANCHOR] 상위 변화점 τ^U_{c,k+1} (IA→Upper LA 결과로 확정)
  tau_U_vec <- state$tau_upper[[c_j]]
  anchor    <- if ((k + 1) <= length(tau_U_vec)) tau_U_vec[k + 1] else G_L[n_free]

  atoms_c <- params$atoms[[c_j]]
  atom_k   <- if (k <= length(atoms_c) && !is.null(atoms_c[[k]])) atoms_c[[k]] else
              list(gamma1=1, gamma2=1, shape_beta=1.0, shape_gamma=0.0)
  kp1_safe <- min(k + 1, length(atoms_c))
  atom_kp1 <- if (kp1_safe >= 1 && !is.null(atoms_c[[kp1_safe]])) atoms_c[[kp1_safe]] else
              list(gamma1=-1, gamma2=-1, shape_beta=1.0, shape_gamma=0.0)

  theta_j   <- params$theta[j, ]
  gamma_k   <- get_state_intercept(j, k, params)
  gamma_kp1 <- get_state_intercept(j, k + 1, params)

  # ── Segment likelihood over G_L ──
  #   Normal log-likelihood with the series-level scale v_j INSIDE the
  #   variance component:  σ²_{j,t} = v_j · φ_{j,t}^{ξ_{j,t}}.
  #   This is the pre-collapse Normal kernel (the v_j-collapsed Student-t
  #   form previously used here was incorrect for this update).
  #   Fast path delegates to lower_la_pair_cpp(); the exact R fallback below
  #   uses the identical kernel.
  v_j <- params$v[[j]]
  v_j <- if (is.null(v_j) || length(v_j) == 0 || !is.finite(v_j[1]) || v_j[1] <= 0) 1 else v_j[1]

  use_cpp <- isTRUE(get0(".hiercpd_rcpp_enabled", ifnotfound = FALSE)) &&
             exists("lower_la_pair_cpp")

  chosen <- NA_integer_
  if (use_cpp) {
    if (exists(".hiercpd_ensure_datx", mode = "function")) .hiercpd_ensure_datx(precomp)
    phi_row <- params$phi[j, ]
    xi_row  <- params$xi[j, ]; xi_row[is.na(xi_row)] <- 0L
    chosen <- tryCatch(
      lower_la_pair_cpp(
        as.numeric(theta_j),
        as.numeric(gamma_k),   atom_k,
        as.numeric(gamma_kp1), atom_kp1,
        as.integer(G_L), as.numeric(Y[j, ]),
        params$alpha, params$beta[c_j],
        as.numeric(v_j), as.numeric(v_j),    # v_k = v_{k+1} = series scalar
        as.numeric(phi_row), as.integer(xi_row),
        as.numeric(pi_star_k), as.integer(anchor),
        as.numeric(precomp$x_global)
      ),
      error = function(e) NA_integer_
    )
    if (length(chosen) != 1L || is.na(chosen) || chosen < 1L || chosen > n_free) chosen <- NA_integer_
  }

  if (is.na(chosen)) {
    # ── Exact R fallback: identical Normal kernel + anchor prior ──
    f_k_vals   <- eval_shape_at_times(theta_j, atom_k,   G_L, precomp)
    f_kp1_vals <- eval_shape_at_times(theta_j, atom_kp1, G_L, precomp)
    mu_k_vals   <- params$alpha + params$beta[c_j] + gamma_k   + f_k_vals
    mu_kp1_vals <- params$alpha + params$beta[c_j] + gamma_kp1 + f_kp1_vals

    ll_k   <- numeric(n_free)
    ll_kp1 <- numeric(n_free)
    for (idx in 1:n_free) {
      t <- G_L[idx]
      phi_t <- params$phi[j, t]; if (!is.finite(phi_t) || phi_t < 1) phi_t <- 1
      xi_t  <- params$xi[j, t];  if (is.na(xi_t)) xi_t <- 0
      sig2  <- v_j * phi_t^xi_t; if (!is.finite(sig2) || sig2 <= 0) sig2 <- 1
      ll_k[idx]   <- dnorm_log_var(Y[j, t], mu_k_vals[idx],   sig2)
      ll_kp1[idx] <- dnorm_log_var(Y[j, t], mu_kp1_vals[idx], sig2)
    }
    cum_ll_k   <- cumsum(ll_k)
    cum_ll_kp1 <- cumsum(ll_kp1)
    total_ll_kp1 <- cum_ll_kp1[n_free]

    # [DIS] anchor 기준 prior over G_L 전체:
    #   p(τ=G_L[i]) ∝ (π*)^{anchor - G_L[i]} (1-π*)^{I(τ ≠ g^+ +1)} × 우도
    g_plus <- G_L[n_free]
    log_probs <- numeric(n_free)
    for (i in 1:n_free) {
      tau <- G_L[i]
      ll_before <- if (i > 1) cum_ll_k[i - 1] else 0
      ll_after  <- total_ll_kp1 - (if (i > 1) cum_ll_kp1[i - 1] else 0)
      log_prior <- (anchor - tau) * log_pi
      if (tau != (g_plus + 1L)) log_prior <- log_prior + log_1mpi
      log_probs[i] <- log_prior + ll_before + ll_after
    }
    chosen <- sample_categorical_log(log_probs)
  }

  # State assignment: G_L 전체에 적용
  if (chosen > 1) state$S_lower[j, G_L[1:(chosen - 1)]] <- as.integer(k)
  state$S_lower[j, G_L[chosen:n_free]] <- as.integer(k + 1)

  state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])
  state
}
cat("06_local_adjusting.R loaded (HDP-HMM + G*_k margin, DIS Lemma 2, Normal v_j-in-variance kernel via lower/upper_la_pair_cpp).\n")
