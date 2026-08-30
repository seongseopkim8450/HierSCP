###############################################################################
# 07_interval_adjusting.R  (v10 — sequential merge scan with p⁻/p/p⁺)
#
# Design implemented here:
#   For each cluster, a single IA call performs exactly two phases
#     split-pass -> merge-pass
#   and EACH phase scans its stage-start targets exactly once.
#
# Split semantics requested by user:
#   - We scan admissible split candidates g from left to right inside current
#     state k = [tau_k, tau_{k+1}-1].
#   - A single NEW LEFT atom is pre-proposed for this whole split attempt.
#   - Candidate actions at nonterminal g:
#       p^- : absorb [tau_k, g-1] into left neighbor (k-1), keep right part as old k
#       p   : confirm split, [tau_k, g-1] uses NEW LEFT state k,
#                              [g, ...] uses old state shifted to k+1
#       p^+ : defer one step right, [tau_k, g] uses NEW LEFT state k,
#                                   [g+1, ...] uses old state shifted to k+1
#   - At the terminal candidate, p^+ is unavailable, so we compare only
#       p^- vs p   (or only p if k is the first state).
#
# Merge semantics remain replacement-based in this file:
#   - left merge deletes current state k and replaces it with k-1 parameters
#   - right merge deletes current state k and replaces it with k+1 parameters
#
# NOTE: The `mode` argument is kept only for backward compatibility and is
# intentionally ignored. IA ALWAYS runs split -> merge in this file.
###############################################################################

interval_adjusting <- function(state, params, Y, precomp, model,
                               mode = c("merge", "split")) {
  if (length(mode) == 0) mode <- "merge"
  for (cc in 1:model$C) {
    res <- interval_adjusting_cluster(cc, state, params, Y, precomp, model, mode = mode)
    state <- res$state
    params <- res$params
  }
  list(state = state, params = params)
}

interval_adjusting_cluster <- function(cc, state, params, Y, precomp, model,
                                       mode = "merge") {
  # v10: split → sequential merge (p⁻/p/p⁺) → same-shape sweep

  res <- split_pass_once(cc, state, params, Y, precomp, model)
  state <- res$state; params <- res$params

  res <- merge_pass_once(cc, state, params, Y, precomp, model)
  state <- res$state; params <- res$params

  # 동일 형상 인접 구간 전용 merge (가중 평균 atom, |J|=1)
  #res <- same_shape_merge_sweep(cc, state, params, Y, precomp, model)
  #state <- res$state; params <- res$params
  
  list(state = state, params = params)
}

###############################################################################
# Helpers for EXACT one-pass scanning per phase
###############################################################################

###############################################################################
# MERGE PASS — Sequential Scan with p⁻/p/p⁺ decisions
#
# split_step의 mirror image로 설계된 sequential merge:
#   k=2(두 번째 상태)부터 K(마지막 상태)까지 순차 스캔.
#   각 위치에서 세 가지 결정을 우도 기반으로 sampling:
#
#     p⁻: 현재 상태 k를 왼쪽 이웃 k-1에 병합 (k 소멸)
#          → old k+1이 position k로 재라벨 → k 유지, 재평가
#     p  : 현재 상태 유지
#          → k ← k+1 (전진)
#     p⁺: 현재 상태 k를 오른쪽 이웃 k+1에 병합 (k 소멸)
#          → 병합된 상태가 position k에 위치 → k 유지, 재평가
#
#   단말 (k = K): p⁺ 불가, p⁻ vs p만 비교
#
# 연쇄 merge: p⁻ 또는 p⁺ 선택 시 K가 감소하고, 새 위치의 상태가
#   즉시 재평가되어 연속 merge가 자연스럽게 발생할 수 있음.
###############################################################################

.merge_eval <- function(cc, k, I_k, j_in_c, state, params, Y, precomp, model) {
  # Returns c(ll_km1, ll_k, ll_kp1) : aggregated interval likelihood under
  #   labels k-1, k, k+1.  Uses merge_eval_cpp when Rcpp active; else R fallback.
  K_now <- state$K[cc]
  has_km1 <- (k - 1L) >= 1L
  has_kp1 <- (k + 1L) <= K_now
  n_k <- length(I_k)
  atoms_c <- params$atoms[[cc]]

  use_cpp <- isTRUE(get0(".hiercpd_rcpp_enabled", ifnotfound = FALSE)) &&
             exists("merge_eval_cpp")

  if (use_cpp) {
    # [PERF] push the D-cube into the module-level cache once (no-op if already
    #        current); the merge kernel now reads G_DATX instead of receiving
    #        the cube as an argument on every call.
    if (exists(".hiercpd_ensure_datx", mode = "function")) .hiercpd_ensure_datx(precomp)
    theta_mat <- params$theta[j_in_c, , drop = FALSE]
    gamma_vec <- vapply(j_in_c, function(jj) {
      gj <- params$gamma[[jj]]
      if (is.null(gj) || length(gj) == 0 || !is.finite(gj[1])) 0.0 else as.numeric(gj[1])
    }, numeric(1))
    v_vec <- vapply(j_in_c, function(jj) {
      vj <- params$v[[jj]]
      if (is.null(vj) || !is.finite(vj[1]) || vj[1] <= 0) 1.0 else as.numeric(vj[1])
    }, numeric(1))
    Y_sub   <- Y[j_in_c, , drop = FALSE]
    phi_sub <- params$phi[j_in_c, , drop = FALSE]
    xi_sub  <- params$xi[j_in_c, , drop = FALSE]
    storage.mode(xi_sub) <- "integer"
    safe_atom <- list(gamma1=1,gamma2=1,shape_beta=1.0,shape_gamma=0.0)
    atom_km1 <- if (has_km1 && !is.null(atoms_c[[k-1L]])) atoms_c[[k-1L]] else safe_atom
    atom_k   <- if (k <= length(atoms_c) && !is.null(atoms_c[[k]])) atoms_c[[k]] else safe_atom
    atom_kp1 <- if (has_kp1 && !is.null(atoms_c[[k+1L]])) atoms_c[[k+1L]] else safe_atom

    res <- tryCatch(
      merge_eval_cpp(
        theta_mat, as.integer(I_k),
        gamma_vec, atom_km1, v_vec,
        gamma_vec, atom_k,   v_vec,
        gamma_vec, atom_kp1, v_vec,
        Y_sub, params$alpha, params$beta[cc],
        phi_sub, xi_sub,
        as.integer(if (has_km1) 1L else 0L),
        as.integer(if (has_kp1) 1L else 0L),
        as.numeric(precomp$x_global)
      ),
      error = function(e) NULL
    )
    if (!is.null(res) && length(res) == 3L && all(is.finite(res) | res == -Inf)) {
      return(as.numeric(res))
    }
    # else fall through
  }

  # ── R fallback ──
  ll_km1 <- if (has_km1) {
    s <- 0; for (j in j_in_c) s <- s + log_lik_interval(j, I_k, rep(k - 1L, n_k), Y, state, params, precomp, model); s
  } else -Inf
  ll_k <- {
    s <- 0; for (j in j_in_c) s <- s + log_lik_interval(j, I_k, rep(k, n_k), Y, state, params, precomp, model); s
  }
  ll_kp1 <- if (has_kp1) {
    s <- 0; for (j in j_in_c) s <- s + log_lik_interval(j, I_k, rep(k + 1L, n_k), Y, state, params, precomp, model); s
  } else -Inf
  c(ll_km1, ll_k, ll_kp1)
}

merge_pass_once <- function(cc, state, params, Y, precomp, model) {
  if (state$K[cc] < 2) return(list(state = state, params = params))

  T_len <- model$T_len; m_min <- model$m_min
  j_in_c <- which(state$cluster == cc)
  if (length(j_in_c) == 0) return(list(state = state, params = params))

  lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0

  k <- 2L
  max_steps <- state$K[cc] * 2L

  for (step in seq_len(max_steps)) {
    K_now <- state$K[cc]
    K_min <- if (!is.null(model$K_min)) model$K_min else 1L
    if (K_now <= K_min || k < 2L || k > K_now) break

    tau_U <- state$tau_upper[[cc]]
    tau_k   <- tau_U[k]
    tau_kp1 <- if (k + 1 <= length(tau_U)) tau_U[k + 1] else T_len + 1
    I_k     <- tau_k:(tau_kp1 - 1)
    n_k     <- length(I_k)

    ps <- params$pi_star[[cc]]
    ps_k   <- min(max(ps[k], 1e-10), 1 - 1e-10)
    ps_km1 <- if (k - 1 >= 1 && k - 1 <= length(ps)) {
      min(max(ps[k - 1], 1e-10), 1 - 1e-10)
    } else 0.5

    log_poisson_merge_bonus <- log(max(K_now, 1)) - log(lambda0)

    # (B-i 대칭) merge 는 split 의 역연산: 상태 k 를 흡수하면 그 atom a_{z_k} 가
    #   풀로 되돌아간다. detailed balance 를 위해 split 의 제안 q(z^new) 에 대응하는
    #   역제안 q(z_k) = π_{z_k} / Σ_{m∈A'} π_m 를 merge(p-, p+) 확률에 더한다.
    #   [Sequential Gibbs] split 의 prior 가중 제거와 짝을 맞춰, merge 의 역제안
    #   항(log_q_rev)도 제거한다.  atom 배정은 PU 의 B2 가 책임지므로 merge 가
    #   사라지는 atom 의 prior 가중을 따로 더하면 이중계산이다.
    terminal <- (k >= K_now)

    # ── 우도 3-가설을 cpp 한 번에 (p-, p, p+) ──
    ll3 <- .merge_eval(cc, k, I_k, j_in_c, state, params, Y, precomp, model)

    if (!terminal) {
      log_prior <- c(
        n_k * log(ps_km1) + log(1 - ps_km1) + log_poisson_merge_bonus,   # p-
        max(n_k - m_min, 0) * log(ps_k) + log(1 - ps_k),                             # p (유지)
        log(1 - ps_k) + log_poisson_merge_bonus                          # p+
      )
      log_lik <- c(ll3[1], ll3[2], ll3[3])
      decision <- sample_categorical_log(log_prior + log_lik)

      if (decision == 1L) {
        res <- execute_merge(cc, k, "left", state, params, model)
        state <- res$state; params <- res$params
      } else if (decision == 2L) {
        k <- k + 1L
      } else {
        res <- execute_merge(cc, k, "right", state, params, model)
        state <- res$state; params <- res$params
      }

    } else {
      # 단말: p-, p 만 (p+ = -Inf)
      log_prior <- c(
        n_k * log(ps_km1) + log(1 - ps_km1) + log_poisson_merge_bonus,   # p-
        max(n_k - m_min, 0) * log(ps_k) + log(1 - ps_k)                              # p (유지)
      )
      log_lik <- c(ll3[1], ll3[2])
      decision <- sample_categorical_log(log_prior + log_lik)

      if (decision == 1L) {
        res <- execute_merge(cc, k, "left", state, params, model)
        state <- res$state; params <- res$params
      }
      break
    }
  }

  list(state = state, params = params)
}


split_pass_once <- function(cc, state, params, Y, precomp, model) {
  tau_start <- state$tau_upper[[cc]]
  if (length(tau_start) < 1) {
    return(list(state = state, params = params))
  }

  # Stage-start state list only, exactly once each.
  scan_targets <- tau_start

  for (tau_ref in scan_targets) {
    if (state$K[cc] >= model$K_max) break

    tau_cur <- state$tau_upper[[cc]]
    k_cur <- match(tau_ref, tau_cur)

    # If not found, skip rather than revisiting anything new.
    if (is.na(k_cur)) next

    res <- split_step(k_cur, cc, state, params, Y, precomp, model)
    state <- res$state
    params <- res$params
  }

  list(state = state, params = params)
}

insert_after_state <- function(vec, k, new_value = NULL) {
  if (is.null(vec)) return(NULL)
  if (length(vec) == 0) return(new_value)
  if (is.null(new_value)) new_value <- vec[min(k, length(vec))]
  c(vec[seq_len(k)], new_value, if (k < length(vec)) vec[(k + 1):length(vec)] else NULL)
}

insert_before_state <- function(vec, k, new_value = NULL) {
  if (is.null(vec)) return(NULL)
  if (length(vec) == 0) return(new_value)
  k <- max(1L, min(as.integer(k), length(vec) + 1L))
  if (is.null(new_value)) {
    ref_idx <- min(max(k, 1L), length(vec))
    new_value <- vec[ref_idx]
  }
  left <- if (k > 1L) vec[1:(k - 1L)] else NULL
  right <- if (k <= length(vec)) vec[k:length(vec)] else NULL
  c(left, new_value, right)
}

build_split_left_proposal_params <- function(cc, k, proposed_left_atom, state, params) {
  split_params <- params
  K_c <- state$K[cc]

  # atoms: insert the NEW LEFT atom at label k, shift old k:K to k+1:(K+1)
  old_atoms <- params$atoms[[cc]]
  new_atoms <- vector("list", K_c + 1L)
  if (k > 1L) {
    for (kk in 1:(k - 1L)) new_atoms[[kk]] <- old_atoms[[min(kk, length(old_atoms))]]
  }
  new_atoms[[k]] <- proposed_left_atom
  for (kk in k:K_c) new_atoms[[kk + 1L]] <- old_atoms[[min(kk, length(old_atoms))]]
  split_params$atoms[[cc]] <- new_atoms

  # [γ*(c,j)] γ*는 시계열 단위 스칼라이므로 split 시 조작 불필요
  j_in_c <- which(state$cluster == cc)

  # provisional stick weights / conditional persistence
  if (!is.null(params$pi_weights[[cc]])) {
    op <- params$pi_weights[[cc]]
    np <- numeric(K_c + 1L)
    if (k > 1L) {
      for (kk in 1:(k - 1L)) np[kk] <- op[min(kk, length(op))]
    }
    np[k] <- op[min(k, length(op))] / 2
    np[k + 1L] <- np[k]
    if (k < K_c) {
      for (kk in (k + 1L):K_c) np[kk + 1L] <- op[min(kk, length(op))]
    }
    np <- np / sum(np)
    split_params$pi_weights[[cc]] <- np
    split_params$pi_star[[cc]] <- compute_pi_star(np)
  }

  split_params
}

###############################################################################
# 2-STAGE DATA-INFORMED ATOM PROPOSAL FOR SPLIT
#
# Stage 1: 2차 다항식 적합 → δ₁, δ₂ 결정
#   - 잔차 r_{j,t} = Y - α - β_c, 평균 차감으로 γ* 흡수
#   - Centered quadratic: r̃ ~ (x-x̄) + (x-x̄)²
#   - 정밀도 가중 a₁, a₂ 집계 → sign(a₁)→δ₁, sign(a₂)→δ₂
#   - 부호 매핑: f'' = -δ₂·b·H''(x), H'' < 0 전형적
#     → sign(f'') = sign(δ₂) → sign(a₂) = sign(δ₂)
#
# Stage 2: (b, g) approximate full conditional 에서 확률적 샘플링
#   - δ₁, δ₂ 고정 시 f^shape = g·A_j(x) + b·B_j(x) (선형)
#   - WLS 정규방정식 → 사후 mode + Σ → truncated Normal 샘플링
#
# 기존 대비 이점:
#   - δ₁,δ₂: 4-way 우도 비교(×4 f_shape 재계산) → O(n) 회귀 1회
#   - b,g: 사전분포 blind 추출 → full conditional mode 근방 제안
#     → "tailored proposal" (Brooks et al. 2003)
###############################################################################

propose_atom_spline_conditional <- function(cc, k, tau_k, g_left, g_right,
                                             j_in_c, state, params, Y, precomp, model) {
  hyper <- params$hyper
  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 3.0
  r_beta  <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 0.5
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 2.5
  r_gamma <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 1.0

  Mp1 <- ncol(params$theta)

  # ── Fallback: prior에서 추출 ──
  fallback <- function(d1 = NULL, d2 = NULL) {
    sb <- max(rgamma(1, shape = a_beta, rate = r_beta), 0.01)
    sg <- max(rgamma(1, shape = a_gamma, rate = r_gamma), 0)
    if (is.null(d1)) d1 <- sample(c(-1L, 1L), 1)
    if (is.null(d2)) d2 <- sample(c(-1L, 1L), 1)
    canonicalize_atom(list(gamma1 = d1, gamma2 = d2,
                           shape_beta = sb, shape_gamma = sg))
  }

  # ── 참조 구간 ──
  g_mid <- as.integer((g_left + g_right) / 2)
  I_left <- tau_k:(g_mid - 1L)
  if (length(I_left) < 3L) return(fallback())
  x_left <- precomp$x_global[I_left]
  n_left <- length(I_left)

  # ==================================================================
  # Step 1: 시리즈별 H_j(x_t), 잔차, 가중치 사전계산
  # ==================================================================
  series_data <- vector("list", length(j_in_c))
  any_valid <- FALSE

  for (idx in seq_along(j_in_c)) {
    j <- j_in_c[idx]
    theta_j <- params$theta[j, ]
    theta_shape <- theta_j[2:Mp1]
    SZ_safe <- max(sum(theta_shape^2), 1e-12)

    theta_outer <- tcrossprod(theta_shape)
    D_sub <- precomp$D_at_x[2:Mp1, 2:Mp1, , drop = FALSE]
    H_vals <- numeric(n_left)
    for (i in seq_along(I_left)) {
      H_vals[i] <- sum(theta_outer * D_sub[, , I_left[i]]) / SZ_safe
    }

    # [COLLAPSED] γ* 미차감. r = Y - α - β_c. 가중 centering 으로 γ* 적분 소거.
    #   (update_atoms 와 일치; f 의 segment 절편↔γ* 흡수 차단)
    r_j <- Y[j, I_left] - params$alpha - params$beta[cc]

    v_j <- get_state_variance(j, k, params)
    w_t <- rep(1 / v_j, n_left)
    for (ti in seq_along(I_left)) {
      t_idx <- I_left[ti]
      xi_t <- if (!is.null(params$xi)) params$xi[j, t_idx] else 0
      if (is.na(xi_t)) xi_t <- 0
      if (xi_t == 1L) {
        phi_t <- if (!is.null(params$phi)) params$phi[j, t_idx] else 1
        if (!is.finite(phi_t) || phi_t < 1) phi_t <- 1
        w_t[ti] <- 1 / (v_j * phi_t)
      }
    }

    # [COLLAPSED] γ* 소거를 위한 가중 centering (분모 W_j = Σw + 1/σ²_c)
    sigma2_c_val <- if (!is.null(params$sigma2_gamma_c)) {
      params$sigma2_gamma_c[cc]
    } else if (!is.null(params$sigma2_gamma)) {
      params$sigma2_gamma
    } else 0.01
    if (!is.finite(sigma2_c_val) || sigma2_c_val <= 0) sigma2_c_val <- 0.01
    W_j <- sum(w_t) + 1.0 / sigma2_c_val
    r_c <- r_j - sum(w_t * r_j) / W_j

    series_data[[idx]] <- list(x = x_left, H = H_vals, r_c = r_c, w = w_t, W_j = W_j)
    any_valid <- TRUE
  }

  if (!any_valid) return(fallback())

  # ==================================================================
  # Step 2: 완전 prior 추출 + 이산 Gibbs 선택
  #
  # δ₁, δ₂ ~ Uniform({+1,-1})  (매 제안마다 25% 확률로 모든 형상 탐색)
  # (b, g) ~ Base H: b ~ Gamma(a_β, r_β), g ~ Gamma(a_γ, r_γ)
  # N_cand개 후보 생성 → sub-interval 우도 비례 이산 Gibbs 선택
  #
  # Ishwaran & James (2001) 응용: prior에서 뽑고 likelihood로 선별
  # 틀린 (δ₁,δ₂,b,g) → 우도 낮음 → MH 기각
  # ==================================================================
  delta1 <- sample(c(-1L, 1L), 1)
  delta2 <- sample(c(-1L, 1L), 1)

  N_cand <- 50L
  cand_b <- pmax(rgamma(N_cand, shape = a_beta,  rate = r_beta),  0.01)
  cand_g <- pmax(rgamma(N_cand, shape = a_gamma, rate = r_gamma), 0.00)

  # ── 후보별 sub-interval 우도 스코어 ──
  #   [선택 B] 무거운 N_cand × series × n_left 루프만 C++(profile_atom_scores_cpp).
  #   난수(δ, Gamma draw, categorical)는 R 유지 → R RNG 재현성·검증 보존.
  use_cpp_ps <- isTRUE(get0(".hiercpd_rcpp_enabled", ifnotfound = FALSE)) &&
                exists("profile_atom_scores_cpp")
  if (use_cpp_ps) {
    # series_data 평탄화: x,H,r_c,w 이어붙이고 시리즈 경계·W_j 전달
    x_all <- numeric(0); H_all <- numeric(0); rc_all <- numeric(0); w_all <- numeric(0)
    seg_start <- integer(length(series_data)); seg_len <- integer(length(series_data))
    W_j_vec <- numeric(length(series_data))
    off <- 0L
    for (si in seq_along(series_data)) {
      s <- series_data[[si]]
      Ls <- length(s$x)
      x_all <- c(x_all, s$x); H_all <- c(H_all, s$H)
      rc_all <- c(rc_all, s$r_c); w_all <- c(w_all, s$w)
      seg_start[si] <- off; seg_len[si] <- Ls; W_j_vec[si] <- s$W_j
      off <- off + Ls
    }
    log_scores <- tryCatch(
      profile_atom_scores_cpp(x_all, H_all, rc_all, w_all,
                              as.integer(seg_start), as.integer(seg_len),
                              W_j_vec, as.numeric(delta1), as.numeric(delta2),
                              cand_b, cand_g),
      error = function(e) NULL)
    if (is.null(log_scores)) use_cpp_ps <- FALSE   # 실패 시 R 폴백
  }
  if (!use_cpp_ps) {
    # ── R 폴백 (C++ 와 동일 수학) ──
    log_scores <- numeric(N_cand)
    for (ci in 1:N_cand) {
      g_i <- cand_g[ci];  b_i <- cand_b[ci]
      ll <- 0
      for (s in series_data) {
        z_g_raw <- delta1 * s$x
        z_b_raw <- (delta1 + delta2) / 2 * s$x - delta2 * s$H
        A_c <- z_g_raw - sum(s$w * z_g_raw) / s$W_j
        B_c <- z_b_raw - sum(s$w * z_b_raw) / s$W_j
        f_i <- A_c * g_i + B_c * b_i
        ll <- ll - 0.5 * sum(s$w * (s$r_c - f_i)^2)
      }
      log_scores[ci] <- ll
    }
  }

  # 우도 비례 이산 Gibbs 선택 (난수 R)
  chosen <- sample_categorical_log(log_scores)
  sb <- cand_b[chosen]
  sg <- cand_g[chosen]

  canonicalize_atom(list(gamma1 = delta1, gamma2 = delta2,
                         shape_beta = sb, shape_gamma = sg))
}

###############################################################################
# SAME-SHAPE MERGE SWEEP
# 동일 (δ₁,δ₂) 인접 구간을 대상으로 한 전용 RJMCMC merge.
# 일반 merge는 형상이 다른 구간도 시도하여 기각이 많지만,
# 이 sweep은 동일 형상만 겨냥하여 높은 수용률을 달성합니다.
# 상세 균형: 표준 MH 제안이므로 정확히 성립.
###############################################################################

###############################################################################
# SAME-SHAPE MERGE SWEEP (가중 평균 atom + 차분 매개변수화)
#
# 수용비:
#   α = L(merged) / L(current) × p(K-1)/p(K) × n_same/|G_split|
#
# 야코비안: 차분 u = atom_L - atom_R 매개변수화에서 |J| = 1
# Prior/proposal 상쇄: 역방향 split 제안을 prior와 일치시켜 H, p(γ*), p(v) 항 소거
###############################################################################

same_shape_merge_sweep <- function(cc, state, params, Y, precomp, model) {
  K_c <- state$K[cc]
  if (K_c < 2) return(list(state = state, params = params))

  T_len <- model$T_len; m_min <- model$m_min
  j_in_c <- which(state$cluster == cc)
  if (length(j_in_c) == 0) return(list(state = state, params = params))

  for (k in seq(K_c, 2, by = -1)) {
    if (k > state$K[cc]) next
    K_now <- state$K[cc]
    K_min <- if (!is.null(model$K_min)) model$K_min else 1L
    atoms_now <- params$atoms[[cc]]
    if (k > K_now || K_now <= K_min) next

    atom_left  <- atoms_now[[k - 1]]
    atom_right <- atoms_now[[k]]
    if (is.null(atom_left) || is.null(atom_right)) next

    if (sign(atom_left$gamma1) != sign(atom_right$gamma1) ||
        sign(atom_left$gamma2) != sign(atom_right$gamma2)) next

    # ── 구간 정보 ──
    tau_U <- state$tau_upper[[cc]]
    tau_k <- tau_U[k]
    tau_kp1 <- if (k + 1 <= length(tau_U)) tau_U[k + 1] else T_len + 1
    tau_km1 <- tau_U[k - 1]
    I_left <- tau_km1:(tau_k - 1); I_right <- tau_k:(tau_kp1 - 1)
    I_merged <- tau_km1:(tau_kp1 - 1)
    n_l <- length(I_left); n_r <- length(I_right)
    w_l <- n_l / (n_l + n_r); w_r <- 1 - w_l

    # ── 현재 우도 ──
    ll_current <- 0
    for (j in j_in_c) {
      ll_current <- ll_current +
        log_lik_interval(j, I_left,  rep(k - 1L, length(I_left)),  Y, state, params, precomp, model) +
        log_lik_interval(j, I_right, rep(k,      length(I_right)), Y, state, params, precomp, model)
    }

    # ── 가중 평균 파라미터 생성 ──
    avg_atom <- list(
      gamma1      = atom_left$gamma1,
      gamma2      = atom_left$gamma2,
      shape_beta  = w_l * atom_left$shape_beta  + w_r * atom_right$shape_beta,
      shape_gamma = w_l * atom_left$shape_gamma + w_r * atom_right$shape_gamma
    )

    surv_k <- k - 1
    old_atom <- params$atoms[[cc]][[surv_k]]
    params$atoms[[cc]][[surv_k]] <- avg_atom

    # [γ*(c,j)] γ*는 시계열 단위 스칼라이므로 merge 시 조작 불필요

    # ── Merged 우도 ──
    ll_merged <- 0
    for (j in j_in_c) {
      ll_merged <- ll_merged +
        log_lik_interval(j, I_merged, rep(surv_k, length(I_merged)), Y, state, params, precomp, model)
    }

    # ── 원상 복구 ──
    params$atoms[[cc]][[surv_k]] <- old_atom

    # ── MH 수용비 (야코비안=1, prior/proposal 상쇄) ──
    # [PATCH] Poisson(lambda0_K) prior: merge bonus = log(K/lambda0)
    lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0
    log_prior_ratio <- log(max(K_now, 1)) - log(lambda0)

    n_same <- 0
    for (kk in 2:K_now) {
      a_l <- atoms_now[[kk-1]]; a_r <- atoms_now[[kk]]
      if (!is.null(a_l) && !is.null(a_r) &&
          sign(a_l$gamma1)==sign(a_r$gamma1) && sign(a_l$gamma2)==sign(a_r$gamma2))
        n_same <- n_same + 1
    }
    if (n_same < 1) n_same <- 1
    G_split_size <- max(length(I_merged) - 2 * m_min + 1, 1)
    log_proposal_ratio <- log(n_same) - log(G_split_size)

    log_alpha <- (ll_merged - ll_current) + log_prior_ratio + log_proposal_ratio

    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      params$atoms[[cc]][[surv_k]] <- avg_atom
      res <- execute_merge(cc, k, "left", state, params, model)
      state <- res$state
      params <- res$params
      print("same-shape-merged!!!!!!")
    }
  }

  list(state = state, params = params)
}

###############################################################################
# MERGE (replacement-based)
###############################################################################

merge_step <- function(k, cc, state, params, Y, precomp, model) {
  K_c <- state$K[cc]; T_len <- model$T_len; m_min <- model$m_min
  tau_U <- state$tau_upper[[cc]]
  j_in_c <- which(state$cluster == cc)
  if (length(j_in_c) == 0) return(list(state = state, params = params, merged = FALSE))

  tau_k   <- tau_U[k]
  tau_kp1 <- if (k + 1 <= length(tau_U)) tau_U[k + 1] else T_len + 1
  I_k  <- tau_k:(tau_kp1 - 1)
  n_k  <- length(I_k)

  ps <- params$pi_star[[cc]]
  ps_k   <- min(max(ps[k], 1e-10), 1 - 1e-10)
  ps_km1 <- if (k >= 2 && k - 1 <= length(ps)) min(max(ps[k - 1], 1e-10), 1 - 1e-10) else 0.5

  # [PATCH] Poisson(lambda0_K) prior on K_c
  # merge로 K→K-1: log P(K-1)/P(K) = log(K_c / lambda0)
  lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0
  log_poisson_merge_bonus <- log(max(K_c, 1)) - log(lambda0)

  log_prior <- c(
    n_k * log(ps_km1) + log(1 - ps_km1) + log_poisson_merge_bonus,  # [PATCH] left merge
    max(n_k - m_min, 0) * log(ps_k) + log(1 - ps_k),                # maintain (no K change)
    log(1 - ps_k) + log_poisson_merge_bonus                           # [PATCH] right merge
  )
  if (k <= 1)   log_prior[1] <- -Inf
  if (k >= K_c) log_prior[3] <- -Inf
  # [K_min] merge로 K < K_min이면 merge 선택지 차단
  K_min_val <- if (!is.null(model$K_min)) model$K_min else 1L
  if (K_c <= K_min_val) { log_prior[1] <- -Inf; log_prior[3] <- -Inf }

  log_lik <- c(0, 0, 0)
  for (j in j_in_c) {
    if (is.finite(log_prior[1])) {
      log_lik[1] <- log_lik[1] +
        log_lik_interval(j, I_k, rep(k - 1L, n_k), Y, state, params, precomp, model)
    }

    log_lik[2] <- log_lik[2] +
      log_lik_interval(j, I_k, rep(k, n_k), Y, state, params, precomp, model)

    if (is.finite(log_prior[3])) {
      log_lik[3] <- log_lik[3] +
        log_lik_interval(j, I_k, rep(k + 1L, n_k), Y, state, params, precomp, model)
    }
  }

  decision <- sample_categorical_log(log_prior + log_lik)
  merged <- FALSE

  if (decision == 1) {
    res <- execute_merge(cc, k, "left", state, params, model)
    state <- res$state; params <- res$params; merged <- TRUE
  } else if (decision == 3) {
    res <- execute_merge(cc, k, "right", state, params, model)
    state <- res$state; params <- res$params; merged <- TRUE
  }

  list(state = state, params = params, merged = merged)
}
###############################################################################
# SPLIT via DIS (Direct Independence Sampling on changepoint position tau)
#
#   Replaces the sequential p^-/p/p^+ scan with a single global position draw.
#   segment k = [tau_k, tau_kp1 - 1]; admissible right-block starts
#       g in {g_left, ..., g_right},  g_left = tau_k + m_min, g_right = tau_kp1 - m_min
#   Split at g: left block [tau_k, g-1] -> NEW LEFT atom, right block [g, ...] -> old k.
#   Plus a NOSPLIT option (whole segment stays under old atom k).
#
#   Likelihood: cumulative sums over the segment => O(n_k * J_c) (cf. O(n_k^2) scan).
#   Prior (image / Lemma, base = pre-split pi*_k):
#       split at g : (pi*_k)^{(g - tau_k) - m_min} (1 - pi*_k)
#       nosplit    : (pi*_k)^{(n_k)        - m_min}
#
#   Merge compatibility: split-DIS only INCREASES K (split vs nosplit).  The
#   K-decreasing direction (left/right absorption) is owned by merge_pass_once.
#   On a split decision we call execute_split(), which inserts the new atom,
#   syncs z_state (IJ), and relabels.  Identical execute_split() as merge path.
#
#   Rcpp: split_scan_upper_cpp does the position scan in C++ when available;
#   otherwise an exact R fallback (same cumulative-sum math) is used.
###############################################################################

split_step <- function(k, cc, state, params, Y, precomp, model) {
  K_c   <- state$K[cc]
  T_len <- model$T_len
  m_min <- model$m_min
  no_op <- list(state = state, params = params, did_split = FALSE)

  tau_U <- state$tau_upper[[cc]]
  if (k > length(tau_U)) return(no_op)

  tau_k   <- tau_U[k]
  tau_kp1 <- if (k + 1 <= length(tau_U)) tau_U[k + 1] else T_len + 1
  if (is.na(tau_k) || is.na(tau_kp1)) return(no_op)

  seg_len <- tau_kp1 - tau_k
  m_split <- m_min
  if (seg_len < 2L * m_split) return(no_op)

  g_left  <- tau_k + m_split
  g_right <- tau_kp1 - m_split
  if (g_left > g_right) return(no_op)
  G_split <- g_left:g_right
  n_cand  <- length(G_split)

  j_in_c <- which(state$cluster == cc)
  if (length(j_in_c) == 0) return(no_op)

  # ── pre-split persistence (position-prior base; IMAGE-consistent) ──
  ps_cur <- params$pi_star[[cc]]
  ps_k   <- if (k <= length(ps_cur)) min(max(ps_cur[k], 1e-10), 1 - 1e-10) else 0.5
  # k-1 persistence for the ABSORB (p-) option
  has_left <- (k > 1L)
  ps_km1 <- if (has_left && (k - 1L) <= length(ps_cur)) {
    min(max(ps_cur[k - 1L], 1e-10), 1 - 1e-10)
  } else 0.5

  # ── NEW LEFT atom 을 profile(데이터 적합) 로 제안 ─────────────────────────
  #   propose_atom_spline_conditional: 좌블록 참조구간의 잔차에 (δ1,δ2,b,g) 를
  #   적합(50 후보 우도가중 선택) → 그 토막에 맞는 atom.  이걸로 p⁰ 좌블록 우도
  #   LL_left^new(=cum_new) 가 강해져 split 이 적절히 일어난다.
  #   [double-fitting 회피] 이 atom 은 execute_split 에서 atoms 뷰에만 잠정 삽입되고
  #   pool(atom_pool) 에는 안 들어간다.  다음 PU 가 B1(pool 사후갱신)+B2(z 재배정)
  #   후 atoms_view ← pool[z] 로 재구성하며 이 제안 atom 을 덮어쓴다(휘발).
  #   따라서 데이터로 제안해도 모델에 영구 커밋되지 않아 사후 왜곡이 없다.
  #   z_state 동기화용 인덱스는 pool 에서 별도로 하나 뽑는다(인접 제약).
  proposed_left_atom <- propose_atom_spline_conditional(
    cc, k, tau_k, g_left, g_right, j_in_c, state, params, Y, precomp, model)
  N_pool <- if (!is.null(model$N_atom)) model$N_atom else 20L
  z_cur  <- if (!is.null(params$z_state) && !is.null(params$z_state[[cc]])) params$z_state[[cc]] else NULL
  left_nb  <- if (!is.null(z_cur) && k > 1L && (k-1L) <= length(z_cur)) z_cur[k-1L] else NA_integer_
  right_nb <- if (!is.null(z_cur) && k <= length(z_cur)) z_cur[k] else NA_integer_
  forbid_split <- c(left_nb, right_nb); forbid_split <- forbid_split[!is.na(forbid_split)]
  atom_pi_c <- if (!is.null(params$atom_pi)) params$atom_pi[[cc]] else NULL
  z_new_k <- .propose_atom_index(N_pool, forbid_split, atom_pi_c)
  # [Sequential Gibbs] atom 배정 prior 가중은 PU(B2) 책임 → split p⁰ 에 더하지 않음.
  log_q_split <- 0
  old_atom_k <- params$atoms[[cc]][[min(k, length(params$atoms[[cc]]))]]
  # k-1 atom for ABSORB (p-): the previous segment's shape
  atom_km1 <- if (has_left) params$atoms[[cc]][[k - 1L]] else NULL

  lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0
  log_poisson_split_penalty <- log(lambda0) - log(K_c + 1)

  # ── manuscript Proposition 3.1 순차 스캔 {p⁻/p⁰/p⁺} ──
  scan <- .split_seq_scan(
    cc, k, tau_k, tau_kp1, G_split, j_in_c,
    proposed_left_atom, old_atom_k, atom_km1,
    ps_k, ps_km1, has_left, m_min, log_poisson_split_penalty,
    state, params, Y, precomp, model,
    log_q_atom = log_q_split)

  if (scan$action == "nosplit") {
    return(no_op)                                  # 끝까지 이연 → 분할 안 함
  }

  if (scan$action == "split") {
    # ── p⁰ SPLIT: 새 좌 atom, K -> K+1 ──
    g_star <- scan$g
    res <- execute_split(cc, k, g_star, state, params, model,
                         proposed_left_atom = proposed_left_atom,
                         z_new_index = z_new_k)
    return(list(state = res$state, params = res$params, did_split = TRUE))
  }

  # ── p⁻ ABSORB: 좌블록 [tau_k, g-1] 을 k-1 에 흡수 (K 불변) ──
  g_star <- scan$g
  res <- .absorb_left_into_prev(cc, k, g_star, state, params, model)
  list(state = res$state, params = res$params, did_split = FALSE)
}

###############################################################################
# .absorb_left_into_prev — ABSORB (p-) state move
#   Left block [tau_k, g-1] of state k is reassigned to state k-1 (K unchanged).
#   Upper: S_upper[cc, tau_k:(g-1)] <- k-1.  Lower: synchronize S_lower for series
#   whose current label is k within that window.  tau recomputed; NO relabel
#   (label set unchanged), NO atom/pi insertion (K constant).
###############################################################################

.absorb_left_into_prev <- function(cc, k, g, state, params, model) {
  if (k <= 1L) return(list(state = state, params = params))   # no previous state
  T_len <- model$T_len
  tau_k <- state$tau_upper[[cc]][k]
  if (is.na(tau_k) || g <= tau_k) return(list(state = state, params = params))

  # Upper: reassign [tau_k, g-1] to k-1 (vectorized)
  win <- tau_k:(g - 1L)
  su <- state$S_upper[cc, ]
  sel <- win[su[win] == k]
  if (length(sel) > 0L) state$S_upper[cc, sel] <- as.integer(k - 1L)
  state$tau_upper[[cc]] <- extract_changepoints(state$S_upper[cc, ])

  # Lower: upper 경계를 통째로 이식 (split 과 동일 원칙).
  #   absorb 후 각 series 의 S_lower 를 새 S_upper 로 그대로 복사 → lower 세그먼트
  #   = upper 세그먼트가 되어 C2 가 보장되고, series 별 하위 변화점의 미세조정은
  #   같은 iteration 의 Lower LA 가 자유관측집합 G^(k,k+1) 에서 순차적으로 수행한다.
  #   (이전엔 win 안의 시점만 k-1 로 바꿔 series 별 upper 와 어긋나 C2 위반 →
  #    enforce 강제병합을 유발했으므로 제거.)
  j_in_c <- which(state$cluster == cc)
  for (j in j_in_c) {
    state$S_lower[j, ] <- state$S_upper[cc, ]
    state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])
  }
  list(state = state, params = params)
}

###############################################################################
# .split_seq_scan — manuscript Proposition 3.1 의 순차 split 스캔 (p⁻/p⁰/p⁺)
#   후보 g 를 좌→우로 스캔. 각 g 에서 세 사건의 (사전 × 우도)를 비교해 하나 추출:
#     p⁻ : 좌블록 [τ_k, g-1] 을 k-1 에 흡수 (K 불변)
#            (π*_{k-1})^{g-τ_k} (1-π*_{k-1}) · L_absorb(g, k-1)
#     p⁰ : g 에서 split 확정 (K→K+1)
#            (π*_k)^{g-τ_k-m_min} (1-π*_k) · (λ0/(K+1)) · q(z^new) · L_split(g, k)
#     p⁺ : 이연. manuscript 식(3): g 이후 G^U 에서 나올 수 있는 모든 경우의
#          (사전 × 우도) 합 =
#            Σ_{h=g+1}^{g⁺} [ (π*_k)^{n_seg-m_min}(1-π*_k) L_split(h,k)
#                            + (π*_{k-1})^{h-τ_k}(1-π*_{k-1}) L_absorb(h,k-1) ]
#            + (π*_k)^{n_seg-m_min} L_nosplit
#   p⁰ 면 split·종료, p⁻ 면 absorb·종료, p⁺ 면 g+1 로 이동. terminal 에선 p⁺ 없음.
#   끝까지 p⁺ 면 nosplit.
#   반환: list(action = "split"|"absorb"|"nosplit", g = g_star 또는 NA)
#   참고: manuscript 식(3) 의 p⁺ 합 구조는 C++ split_scan_upper_cpp 가 모르므로
#         이 함수는 R 전용(split 은 IA round 에만 발생 → 비용 감당 가능).
###############################################################################

###############################################################################
# .marginal_loglik_left_cd — (ㄷ) 좌블록 주변우도 (atom 적분)
#   split 결정의 p⁰ 좌블록 우도를, "잠정 atom 하나의 조건부 우도"가 아니라
#   atom (δ1,δ2,b,g) 를 사전에서 적분해 없앤 주변우도(model evidence)로 계산한다.
#
#   모델: r̃_t = g·z_g(t) + b·z_b(t) + ε_t,  ε~N(0,1/w_t)
#     z_g = δ1·x,  z_b = (δ1+δ2)/2·x − δ2·H
#   Gaussian working prior: g~N(m_g,1/p_g), b~N(m_b,1/p_b)  (Gamma 를 mode-match)
#   (b,g) 선형-가우시안 → 주변우도 닫힌형:
#     log∫L·H dβ = ½logdetP0 − ½logdetP + ½(β̂ᵀPβ̂ − m0ᵀP0 m0) − ½Σw r̃²
#     P = ZᵀWZ + P0,  β̂ = P⁻¹(ZᵀW r̃ + P0 m0)
#   + 양수 절단보정 log P(b,g>0 | 사후)  (mvtnorm::pmvnorm, .biv_pos_mass)
#   δ 4조합은 logsumexp 로 합산(주변화).
#
#   입력은 좌블록 충분통계량(누적합):  조합 di 별
#     Sgg=Σw z_g², Sgb=Σw z_g z_b, Sbb=Σw z_b², Sgr=Σw z_g r̃, Sbr=Σw z_b r̃,
#     Srr=Σw r̃²  (Srr 은 조합 무관)
#   → 각 g 에서 prefix-sum 누적값으로 O(1) 평가.
###############################################################################
.marg_evidence_one_combo <- function(Sgg, Sgb, Sbb, Sgr, Sbr, Srr,
                                      m_g, m_b, p_g, p_b) {
  P11 <- Sgg + p_g
  P12 <- Sgb
  P22 <- Sbb + p_b
  detP <- P11 * P22 - P12 * P12
  if (!is.finite(detP) || detP < 1e-20) return(list(lev = -Inf, mu = c(m_g, m_b),
                                                    Sig = diag(c(1/p_g, 1/p_b))))
  invd <- 1 / detP
  S11 <-  P22 * invd; S12 <- -P12 * invd; S22 <- P11 * invd
  rhs_g <- Sgr + p_g * m_g
  rhs_b <- Sbr + p_b * m_b
  g_m <- S11 * rhs_g + S12 * rhs_b
  b_m <- S12 * rhs_g + S22 * rhs_b
  quad_post <- g_m * g_m * P11 + 2 * g_m * b_m * P12 + b_m * b_m * P22
  quad_0    <- m_g * m_g * p_g + m_b * m_b * p_b
  lev_untrunc <- 0.5 * log(p_g * p_b) - 0.5 * log(detP) +
                 0.5 * (quad_post - quad_0) - 0.5 * Srr
  list(lev = lev_untrunc, mu = c(g_m, b_m),
       Sig = matrix(c(S11, S12, S12, S22), 2, 2))
}

# 좌블록 [seg_start, g-1] 의 4조합 충분통계량 누적(prefix)을 받아 주변 loglik 반환.
#   ss: 길이 n_seg 의 list, 각 원소 = 시점별 (per-combo z_g,z_b,r̃,w) 가 아니라,
#       아래 .split_seq_scan 에서 누적합 행렬로 만들어 전달.
.marginal_loglik_left_cd <- function(cumSgg, cumSgb, cumSbb, cumSgr, cumSbr, cumSrr,
                                     sl, m_g, m_b, p_g, p_b, use_trunc) {
  # sl = 좌블록 점 개수 (prefix index)
  levs <- numeric(4)
  for (di in 1:4) {
    ev <- .marg_evidence_one_combo(
      cumSgg[di, sl], cumSgb[di, sl], cumSbb[di, sl],
      cumSgr[di, sl], cumSbr[di, sl], cumSrr[sl],
      m_g, m_b, p_g, p_b)
    lev <- ev$lev
    if (use_trunc && is.finite(lev)) {
      pmass <- tryCatch(.biv_pos_mass(ev$mu, ev$Sig), error = function(e) NA_real_)
      if (is.finite(pmass) && pmass > 0) lev <- lev + log(pmass)
    }
    levs[di] <- if (is.finite(lev)) lev else -Inf
  }
  m <- max(levs)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(levs - m)))   # logsumexp over 4 δ-combos
}

###############################################################################
# .profile_loglik_left_cd — (일관-B) 좌블록 profile 조건부 우도
#   p⁰ 좌블록을 atom 적분(주변우도)이 아니라, 사후 mode atom â 로 고정한 *조건부*
#   가우시안 우도로 계산한다.  → p⁻(기존 atom 조건부)·nosplit 과 같은 차원·같은
#   정규화상수 스케일에서 공정 비교.  적분의 logdet 부풀림이 없어 split 폭발 방지.
#
#   각 (δ1,δ2) 조합에서 사후 mode β̂=(ĝ,b̂) (= .marg_evidence_one_combo$mu, 양수
#   clamp) 로 조건부 로그우도:
#     log L(β̂) = −½[ Srr − 2(ĝ·Sgr+b̂·Sbr) + (ĝ²Sgg+2ĝb̂Sgb+b̂²Sbb) ] + C
#   C = −½Σ_t log(2π/w_t) = −½(n·log(2π) − Σlog w)  (좌블록 정규화상수)
#   4조합 중 (절단보정) evidence 최대 조합을 argmax 선택, 그 조건부 우도 반환.
#   (atom 은 평가용; execute_split 가 잠정 삽입 후 다음 PU 가 사후 갱신 → 안전)
###############################################################################
.profile_loglik_left_cd <- function(cumSgg, cumSgb, cumSbb, cumSgr, cumSbr, cumSrr,
                                     cumLogW, sl, m_g, m_b, p_g, p_b, use_trunc) {
  best_di <- 1L; best_ev <- -Inf; best_g <- m_g; best_b <- m_b
  for (di in 1:4) {
    ev <- .marg_evidence_one_combo(
      cumSgg[di, sl], cumSgb[di, sl], cumSbb[di, sl],
      cumSgr[di, sl], cumSbr[di, sl], cumSrr[sl],
      m_g, m_b, p_g, p_b)
    lev <- ev$lev
    if (use_trunc && is.finite(lev)) {
      pmass <- tryCatch(.biv_pos_mass(ev$mu, ev$Sig), error = function(e) NA_real_)
      if (is.finite(pmass) && pmass > 0) lev <- lev + log(pmass)
    }
    if (is.finite(lev) && lev > best_ev) {
      best_ev <- lev; best_di <- di
      best_g <- ev$mu[1]; best_b <- ev$mu[2]
    }
  }
  # 양수 clamp (점추정이므로 절단질량 불필요)
  gh <- max(best_g, 1e-6); bh <- max(best_b, 1e-6)
  # mode atom 고정 조건부 로그우도 (충분통계량으로)
  Sgg <- cumSgg[best_di, sl]; Sgb <- cumSgb[best_di, sl]; Sbb <- cumSbb[best_di, sl]
  Sgr <- cumSgr[best_di, sl]; Sbr <- cumSbr[best_di, sl]; Srr <- cumSrr[sl]
  quad <- Srr - 2*(gh*Sgr + bh*Sbr) + (gh*gh*Sgg + 2*gh*bh*Sgb + bh*bh*Sbb)
  n_pts <- sl
  C <- -0.5 * (n_pts * log(2*pi) - cumLogW[sl])   # Σ log w 누적
  val <- -0.5 * quad + C
  if (is.finite(val)) val else -Inf
}

.split_seq_scan <- function(cc, k, tau_k, tau_kp1, G_split, j_in_c,
                            atom_left, atom_right, atom_km1,
                            ps_k, ps_km1, has_left, m_min,
                            log_poisson_split_penalty,
                            state, params, Y, precomp, model,
                            log_q_atom = 0) {
  n_cand <- length(G_split)
  if (n_cand == 0) return(list(action = "nosplit", g = NA_integer_))

  # 반환 정수 인코딩 → list(action, g) 변환 헬퍼 (C++/R 공통)
  .decode <- function(dec) {
    if (is.na(dec) || dec < 1L) return(NULL)
    if (dec == (2L * n_cand + 1L)) return(list(action = "nosplit", g = NA_integer_))
    if (dec <= n_cand)            return(list(action = "split",  g = G_split[dec]))
    if (dec <= 2L * n_cand)       return(list(action = "absorb", g = G_split[dec - n_cand]))
    NULL
  }

  # ── C++ 우선 경로: split_seq_scan_upper_cpp ──
  #   [(ㄷ) 주변우도] split p⁰/p⁺ 의 좌블록 우도가 atom 적분 주변우도로 바뀌었는데
  #   C++ 커널(split_seq_scan_upper_cpp)은 이를 모르고 잠정 atom 조건부 우도를
  #   쓴다 → 결과 불일치.  따라서 (ㄷ) 가 켜진 동안 split 결정은 R 경로를 강제한다.
  #   (split 은 IA round 에만 발생 → 비용 감당 가능.  C++ 에 주변우도 구현 시 해제.)
  use_cpp <- FALSE
  if (use_cpp) {
    if (exists(".hiercpd_ensure_datx", mode = "function")) .hiercpd_ensure_datx(precomp)
    theta_mat <- params$theta[j_in_c, , drop = FALSE]
    gamma_vec <- vapply(j_in_c, function(jj) {
      gj <- params$gamma[[jj]]
      if (is.null(gj) || length(gj) == 0 || !is.finite(gj[1])) 0.0 else as.numeric(gj[1])
    }, numeric(1))
    v_vec <- vapply(j_in_c, function(jj) {
      vj <- params$v[[jj]]
      if (is.null(vj) || !is.finite(vj[1]) || vj[1] <= 0) 1.0 else as.numeric(vj[1])
    }, numeric(1))
    Y_sub   <- Y[j_in_c, , drop = FALSE]
    phi_sub <- params$phi[j_in_c, , drop = FALSE]
    xi_sub  <- params$xi[j_in_c, , drop = FALSE]
    xi_sub[is.na(xi_sub)] <- 0L
    storage.mode(xi_sub) <- "integer"
    atom_km1_arg <- if (has_left && !is.null(atom_km1)) atom_km1 else old_atom_safe(atom_right)
    dec <- tryCatch(
      split_seq_scan_upper_cpp(
        theta_mat,
        gamma_vec, atom_left,
        gamma_vec, atom_right,
        atom_km1_arg, gamma_vec,
        as.integer(G_split),
        as.integer(tau_k), as.integer(tau_kp1),
        Y_sub, params$alpha, params$beta[cc],
        v_vec, v_vec, v_vec,
        phi_sub, xi_sub,
        as.numeric(ps_k), as.numeric(ps_km1),
        as.integer(m_min),
        as.numeric(log_poisson_split_penalty),
        as.integer(if (has_left) 1L else 0L),
        as.numeric(log_q_atom),
        as.numeric(precomp$x_global)
      ),
      error = function(e) NA_integer_
    )
    res_cpp <- .decode(dec)
    if (!is.null(res_cpp)) return(res_cpp)
    # else fall through to R
  }

  # ── 시점별 우도 재료 (좌:새 atom / 좌:k-1 atom / 블록:old k) ──
  split_params <- build_split_left_proposal_params(cc, k, atom_left, state, params)
  km1_params <- params
  if (has_left && !is.null(atom_km1)) {
    km1_params$atoms[[cc]] <- params$atoms[[cc]]
    km1_params$atoms[[cc]][[k]] <- atom_km1
  }
  seg_start <- tau_k
  n_seg     <- (tau_kp1 - 1L) - seg_start + 1L

  ll_new   <- numeric(n_seg)   # 좌블록 시점별: 새 atom
  ll_km1   <- numeric(n_seg)   # 좌블록 시점별: k-1 atom
  ll_right <- numeric(n_seg)   # 블록 시점별:  old k
  for (loc in seq_len(n_seg)) {
    t1 <- seg_start + loc - 1L
    sn <- 0; sm <- 0; sr <- 0
    for (j in j_in_c) {
      sn <- sn + log_lik_interval(j, t1, k, Y, state, split_params, precomp, model)
      sr <- sr + log_lik_interval(j, t1, k, Y, state, params,       precomp, model)
      if (has_left) sm <- sm + log_lik_interval(j, t1, k, Y, state, km1_params, precomp, model)
    }
    ll_new[loc]   <- sn
    ll_right[loc] <- sr
    ll_km1[loc]   <- sm
  }
  cum_new   <- cumsum(ll_new)
  cum_km1   <- cumsum(ll_km1)
  cum_right <- cumsum(ll_right)
  total_right <- cum_right[n_seg]

  # ── (ㄷ) 좌블록 주변우도용 충분통계량 prefix-sum ──────────────────────────
  #   4개 (δ1,δ2) 조합별로 시점 t 에서 z_g=δ1·x, z_b=(δ1+δ2)/2·x − δ2·H 를 만들고
  #   Σw z_g², Σw z_g z_b, Σw z_b², Σw z_g r̃, Σw z_b r̃ 를 누적.  Σw r̃² 는 조합 무관.
  #   r̃: γ* 가중 centering 으로 소거(아래 W_j 분모).  좌블록은 새 atom 가설이므로
  #   atom 자유 → 적분.  (PU 가 사후 갱신하므로 split 은 atom 미커밋 → double-fitting X)
  combos <- list(c(1L,1L), c(1L,-1L), c(-1L,1L), c(-1L,-1L))
  Sgg <- matrix(0, 4, n_seg); Sgb <- matrix(0, 4, n_seg); Sbb <- matrix(0, 4, n_seg)
  Sgr <- matrix(0, 4, n_seg); Sbr <- matrix(0, 4, n_seg); Srr <- numeric(n_seg)
  SlogW <- numeric(n_seg)   # Σ log w 누적 (profile 조건부 정규화상수 C 용)
  Mp1_cd <- ncol(params$theta)
  sigma2_c_cd <- if (!is.null(params$sigma2_gamma_c)) params$sigma2_gamma_c[cc]
                 else if (!is.null(params$sigma2_gamma)) params$sigma2_gamma else 0.01
  if (!is.finite(sigma2_c_cd) || sigma2_c_cd <= 0) sigma2_c_cd <- 0.01
  # 시점별 (계열 합산된) per-combo 통계 누적
  # 먼저 각 계열의 좌블록 r,w,x,H,centering 분모를 미리 모은다 (좌블록 전체 = seg)
  cd_series <- vector("list", length(j_in_c))
  for (ji in seq_along(j_in_c)) {
    j <- j_in_c[ji]
    tt <- seg_start:(tau_kp1 - 1L)
    rj <- as.numeric(Y[j, tt] - params$alpha - params$beta[cc])   # γ* 미차감
    s2 <- get_effective_variance_all(j, state, params)[tt]
    s2[!is.finite(s2) | s2 <= 0] <- 1
    wj <- 1.0 / s2
    xj <- precomp$x_global[tt]
    theta_j <- params$theta[j, ]
    if (Mp1_cd > 1) {
      th <- theta_j[2:Mp1_cd]; SZ <- max(sum(th^2), 1e-12); tho <- tcrossprod(th)
      Hj <- vapply(seq_along(tt), function(ii)
        sum(tho * precomp$D_at_x[2:Mp1_cd, 2:Mp1_cd, tt[ii]]) / SZ, numeric(1))
    } else Hj <- rep(0, length(tt))
    Wj <- sum(wj) + 1.0 / sigma2_c_cd        # γ* 가중 centering 분모
    cd_series[[ji]] <- list(r = rj, w = wj, x = xj, H = Hj, Wj = Wj)
  }
  # 시점 loc 별 per-combo 증분 (계열 합산) → 누적
  for (loc in seq_len(n_seg)) {
    sgg <- numeric(4); sgb <- numeric(4); sbb <- numeric(4)
    sgr <- numeric(4); sbr <- numeric(4); srr <- 0; slogw <- 0
    for (s in cd_series) {
      w <- s$w[loc]; x <- s$x[loc]; H <- s$H[loc]
      slogw <- slogw + log(max(w, 1e-300))
      # γ* 가중 centering: r̃, 설계열도 centering (update_atoms COLLAPSED 기준과 일치)
      rc <- s$r[loc] - sum(s$w * s$r) / s$Wj
      srr <- srr + w * rc * rc
      for (di in 1:4) {
        d1 <- combos[[di]][1]; d2 <- combos[[di]][2]
        zg_raw <- d1 * x
        zb_raw <- (d1 + d2) / 2 * x - d2 * H
        zg <- zg_raw - sum(s$w * (d1 * s$x)) / s$Wj
        zb <- zb_raw - sum(s$w * ((d1 + d2) / 2 * s$x - d2 * s$H)) / s$Wj
        sgg[di] <- sgg[di] + w * zg * zg
        sgb[di] <- sgb[di] + w * zg * zb
        sbb[di] <- sbb[di] + w * zb * zb
        sgr[di] <- sgr[di] + w * zg * rc
        sbr[di] <- sbr[di] + w * zb * rc
      }
    }
    if (loc == 1L) {
      Sgg[, 1] <- sgg; Sgb[, 1] <- sgb; Sbb[, 1] <- sbb
      Sgr[, 1] <- sgr; Sbr[, 1] <- sbr; Srr[1] <- srr; SlogW[1] <- slogw
    } else {
      Sgg[, loc] <- Sgg[, loc-1] + sgg; Sgb[, loc] <- Sgb[, loc-1] + sgb
      Sbb[, loc] <- Sbb[, loc-1] + sbb; Sgr[, loc] <- Sgr[, loc-1] + sgr
      Sbr[, loc] <- Sbr[, loc-1] + sbr; Srr[loc] <- Srr[loc-1] + srr
      SlogW[loc] <- SlogW[loc-1] + slogw
    }
  }
  # Gaussian working prior mode/precision (Gamma(a,r): mode=(a-1)/r, prec=r²/(a-1))
  hyp <- params$hyper
  a_b_cd <- if (!is.null(hyp$a_shape_beta))  hyp$a_shape_beta  else 3.0
  r_b_cd <- if (!is.null(hyp$b_shape_beta))  hyp$b_shape_beta  else 0.5
  a_g_cd <- if (!is.null(hyp$a_shape_gamma)) hyp$a_shape_gamma else 2.5
  r_g_cd <- if (!is.null(hyp$b_shape_gamma)) hyp$b_shape_gamma else 1.0
  .mp <- function(a, r) {
    if (is.finite(a) && a > 1 && is.finite(r) && r > 0) {
      md <- (a-1)/r; pr <- (a-1)/(md*md)
    } else { md <- 1.0; pr <- 1.0 }
    if (!is.finite(md) || md <= 0) md <- 1.0
    if (!is.finite(pr) || pr <= 0) pr <- 1.0
    c(md, pr)
  }
  wp_b_cd <- .mp(a_b_cd, r_b_cd); wp_g_cd <- .mp(a_g_cd, r_g_cd)
  m_b_cd <- wp_b_cd[1]; p_b_cd <- wp_b_cd[2]
  m_g_cd <- wp_g_cd[1]; p_g_cd <- wp_g_cd[2]
  use_trunc_cd <- requireNamespace("mvtnorm", quietly = TRUE)

  pk   <- min(max(ps_k,   1e-10), 1 - 1e-10)
  pkm1 <- min(max(ps_km1, 1e-10), 1 - 1e-10)
  log_pi_k <- log(pk);  log_1mpi_k <- log(1 - pk)
  log_pi_m <- log(pkm1); log_1mpi_m <- log(1 - pkm1)
  eff_ns <- max(n_seg - m_min, 0L)   # split 항 지수 (n_seg - m_min; h 무관)

  # ── 후보별 사건 log(사전×우도) 미리 계산 ──
  #   split_loc(g) = g - seg_start = 좌블록 점 개수;  유효 범위 1..n_seg-1
  loc_of <- function(g) g - seg_start
  # p⁰(g): split at g
  #   [Sequential Gibbs] log_q_atom 은 호출부에서 0 으로 전달된다 (atom 배정은
  #   PU 의 B2 가 책임 → split 에서 prior 가중 중복 금지).  atom 자체는 우도
  #   LL_left(=cum_new) 계산에만 쓰인다.
  logp0 <- function(g) {
    sl <- loc_of(g)
    if (sl < 1L || sl > n_seg - 1L) return(-Inf)
    # 좌블록 우도 = 제안된 profile atom 고정 조건부 우도 (cum_new).
    #   atom 은 propose_atom_spline_conditional 로 데이터 적합 → 우도 강함.
    #   p⁻(cum_km1, 기존 atom 조건부)와 같은 조건부 스케일·같은 γ* 처리.
    LL_left <- cum_new[sl]; LL_right <- total_right - cum_right[sl]
    val <- max((g - tau_k) - m_min, 0L) * log_pi_k + log_1mpi_k +
           log_poisson_split_penalty + log_q_atom + LL_left + LL_right
    if (is.finite(val)) val else -Inf
  }
  # p⁻(g): absorb [τ_k,g-1] into k-1
  logpm <- function(g) {
    if (!has_left) return(-Inf)
    sl <- loc_of(g)
    if (sl < 1L || sl > n_seg - 1L) return(-Inf)
    LL_left <- cum_km1[sl]; LL_right <- total_right - cum_right[sl]
    val <- (g - tau_k) * log_pi_m + log_1mpi_m + LL_left + LL_right   # 지수: g-τ_k (m_min 할인 없음)
    if (is.finite(val)) val else -Inf
  }
  # nosplit: 전체 세그먼트 old k
  logp_nosplit <- eff_ns * log_pi_k + total_right

  # p⁺(g): manuscript 식(3) — h=g+1..g_right 의 (split(h)+absorb(h)) 합 + nosplit
  #   split(h) 항 지수 = eff_ns (= n_seg-m_min, h 무관),  absorb(h) 지수 = h-τ_k
  logp_plus <- function(g) {
    terms <- c()
    hs <- G_split[G_split > g]
    for (h in hs) {
      sl <- loc_of(h)
      if (sl < 1L || sl > n_seg - 1L) next
      LL_right_h <- total_right - cum_right[sl]
      # split(h): 제안된 profile atom 고정 조건부 우도 (logp0 과 동일 기준)
      t_split  <- eff_ns * log_pi_k + log_1mpi_k + cum_new[sl] + LL_right_h
      terms <- c(terms, t_split)
      if (has_left) {
        t_absorb <- (h - tau_k) * log_pi_m + log_1mpi_m + cum_km1[sl] + LL_right_h
        terms <- c(terms, t_absorb)
      }
    }
    terms <- c(terms, logp_nosplit)   # 끝까지 안 쪼갬
    terms <- terms[is.finite(terms)]
    if (length(terms) == 0) return(-Inf)
    M <- max(terms); M + log(sum(exp(terms - M)))   # logsumexp
  }

  # ── 순차 스캔 ──
  for (ci in seq_len(n_cand)) {
    g <- G_split[ci]
    terminal <- (ci == n_cand)
    lp0 <- logp0(g)
    lpm <- logpm(g)
    if (terminal) {
      # 마지막 후보: p⁺(이연) 의 종착점은 nosplit.  manuscript 식(3) 의 L_nosplit
      #   항이 여기서 '안 쪼갬' 출구가 된다.  {p⁻, p⁰, nosplit} 비교
      #   (has_left=F 면 p⁻ 제외 → {p⁰, nosplit}).
      cand_lp  <- c(lpm, lp0, logp_nosplit)
      cand_act <- c("absorb", "split", "nosplit")
    } else {
      lpp <- logp_plus(g)
      cand_lp  <- c(lpm, lp0, lpp)
      cand_act <- c("absorb", "split", "defer")
    }
    keep <- is.finite(cand_lp)
    if (!any(keep)) {            # 이 g 에서 아무 사건도 불가 → 다음 후보로
      if (terminal) return(list(action = "nosplit", g = NA_integer_)) else next
    }
    cand_lp  <- cand_lp[keep]; cand_act <- cand_act[keep]
    sel <- sample_categorical_log(cand_lp)
    act <- cand_act[sel]
    if (act == "split")  return(list(action = "split",  g = g))
    if (act == "absorb") return(list(action = "absorb", g = g))
    # act == "defer" → 다음 후보 g 로 (재귀 스캔)
  }
  list(action = "nosplit", g = NA_integer_)   # 끝까지 이연 → nosplit
}


# helper: safe atom when k-1 missing (never used in score because has_left=0 path)
old_atom_safe <- function(a) if (is.null(a)) list(gamma1=1,gamma2=1,shape_beta=1.0,shape_gamma=0.0) else a



###############################################################################
# execute_merge, execute_split, relabel_after_change
###############################################################################

execute_merge <- function(cc, k, direction, state, params, model) {
  # [K_min] merge 후 K < K_min이면 merge 거부
  K_min <- if (!is.null(model$K_min)) model$K_min else 1L
  if (state$K[cc] <= K_min) return(list(state = state, params = params))

  target <- if (direction == "left") k - 1 else k + 1

  t_in_k <- which(state$S_upper[cc, ] == k)
  state$S_upper[cc, t_in_k] <- target

  j_in_c <- which(state$cluster == cc)
  for (j in j_in_c) {
    t_j <- which(state$S_lower[j, ] == k)
    state$S_lower[j, t_j] <- target
  }

  relabel_after_change(cc, state, params, model)
}

# ── 인접 제약을 지키며 pool 에서 atom 인덱스 하나 제안 추출 ──
#   forbid: 피해야 할 인덱스(좌/우 이웃). pi 가중이 있으면 비례, 없으면 균등.
.propose_atom_index <- function(N, forbid, atom_pi = NULL) {
  allowed <- setdiff(seq_len(N), forbid)
  if (length(allowed) == 0) allowed <- seq_len(N)  # 안전장치
  if (!is.null(atom_pi) && length(atom_pi) == N) {
    w <- atom_pi[allowed]; w <- pmax(w, 1e-300); w <- w / sum(w)
    return(sample(allowed, 1, prob = w))
  }
  sample(allowed, 1)
}

# ── (B-i) 풀 인덱스 제안의 로그확률 ──────────────────────────────────────────
#   q(n) = π_n · 1(n ∈ A) / Σ_{m∈A} π_m ,  A = {1..N} \ forbid
#   detailed balance 를 위해 split 은 +log q(z^new), merge 는 +log q(z_k) 를 더한다.
#   (atom_pi 없거나 길이 불일치 시 균등 제안: 1/|allowed|)
.atom_proposal_logprob <- function(idx, N, forbid, atom_pi = NULL) {
  allowed <- setdiff(seq_len(N), forbid)
  if (length(allowed) == 0) allowed <- seq_len(N)
  if (!(idx %in% allowed)) return(-Inf)
  if (!is.null(atom_pi) && length(atom_pi) == N) {
    w <- pmax(atom_pi[allowed], 1e-300)
    return(log(pmax(atom_pi[idx], 1e-300)) - log(sum(w)))
  }
  -log(length(allowed))  # 균등
}

# ===========================================================================
# execute_split (z_state 동기화 추가)
# ===========================================================================
execute_split <- function(cc, k, g, state, params, model, proposed_left_atom = NULL,
                          z_new_index = NULL) {
  T_len <- model$T_len; K_c <- state$K[cc]; m_min <- model$m_min

  # ----- upper states (vectorized) -----
  su <- state$S_upper[cc, ]
  su[su >= k] <- su[su >= k] + 1L                       # shift labels >= k up by 1
  if (g > 1L) {
    idxL <- 1:(g - 1L)
    sel <- idxL[su[idxL] == (k + 1L)]                   # left piece -> NEW LEFT state k
    if (length(sel) > 0L) su[sel] <- k
  }
  state$S_upper[cc, ] <- su
  state$K[cc] <- K_c + 1L
  state$tau_upper[[cc]] <- extract_changepoints(state$S_upper[cc, ])

  # ----- lower states: upper 경계를 통째로 이식 -----
  #   split 직후 각 series 의 하위 상태열을 새 upper 상태열로 그대로 복사한다.
  #   → lower 세그먼트 = upper 세그먼트이므로 G_U 의 m_min 마진이 lower 에도 그대로
  #     적용되어 C2 가 보장된다(어떤 series 도 짧은 세그먼트가 생기지 않음).
  #   이후 같은 iteration 의 Lower LA 가 이 split 변화점들을 anchor 로 삼아,
  #   자유관측집합 G^(k,k+1) 에서 series 별 하위 변화점을 순차적으로 조정한다.
  #   (이전: split 대상 세그먼트만 g 로 나누고 나머지 lower 경계는 series 고유값을
  #    남겨둠 → lower 가 upper 와 다른 series 에서 m_min 위반 → 강제병합 유발. 제거.)
  j_in_c <- which(state$cluster == cc)
  for (j in j_in_c) {
    state$S_lower[j, ] <- state$S_upper[cc, ]
    state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])
  }

  # ----- params: atoms (펼친 뷰) insert NEW LEFT at k -----
  new_left_atom <- if (is.null(proposed_left_atom)) sample_atom_from_base(params) else proposed_left_atom
  old_atoms <- params$atoms[[cc]]
  new_K <- K_c + 1L
  new_atoms <- vector("list", new_K)
  if (k > 1L) for (kk in 1:(k - 1L)) new_atoms[[kk]] <- old_atoms[[min(kk, length(old_atoms))]]
  new_atoms[[k]] <- new_left_atom
  for (kk in k:K_c) new_atoms[[kk + 1L]] <- old_atoms[[min(kk, length(old_atoms))]]
  params$atoms[[cc]] <- new_atoms

  # ----- [동기화] z_state: k 위치에 새 atom 인덱스 삽입, 나머지 시프트 -----
  #   (B-i 해석 B) split_step 이 풀에서 뽑아 우도 평가에 사용한 바로 그 인덱스
  #   z_new_index 를 z[k] 로 박는다. → atoms[[cc]][[k]] = pool[z[k]] 가 보장되어
  #   다음 PU 의 atoms_view <- pool[z] 재구성과 완전히 일치 (split↔PU 풀 단일화).
  if (!is.null(params$z_state) && !is.null(params$z_state[[cc]])) {
    N <- if (!is.null(model$N_atom)) model$N_atom else 20L
    oz <- params$z_state[[cc]]
    if (length(oz) != K_c) {
      # 길이 불일치 안전장치: 재초기화(1,2,1,2,..)
      oz <- integer(K_c)
      if (K_c >= 1) oz[1] <- 1L
      if (K_c >= 2) for (kk in 2:K_c) oz[kk] <- if (oz[kk - 1] == 1L) 2L else 1L
    }
    if (!is.null(z_new_index) && z_new_index >= 1L && z_new_index <= N) {
      z_new_k <- as.integer(z_new_index)              # split_step 이 뽑은 인덱스 재사용
    } else {
      # 안전장치(인덱스 미전달): 이웃과 겹치지 않게 풀에서 제안
      left_nb  <- if (k > 1L) oz[k - 1L] else NA_integer_
      right_nb <- if (k <= K_c) oz[k] else NA_integer_
      forbid <- c(left_nb, right_nb); forbid <- forbid[!is.na(forbid)]
      z_new_k <- .propose_atom_index(N, forbid, params$atom_pi[[cc]])
    }
    nz <- integer(new_K)
    if (k > 1L) nz[1:(k - 1L)] <- oz[1:(k - 1L)]
    nz[k] <- z_new_k
    if (k <= K_c) nz[(k + 1L):new_K] <- oz[k:K_c]
    params$z_state[[cc]] <- nz
  }

  # ----- pi_weights / pi_star (변화점 persistence; 원본과 동일) -----
  op <- params$pi_weights[[cc]]
  np <- numeric(new_K)
  if (k > 1L) for (kk in 1:(k - 1L)) np[kk] <- op[min(kk, length(op))]
  np[k] <- op[min(k, length(op))] / 2
  np[k + 1L] <- np[k]
  if (k < K_c) for (kk in (k + 1L):K_c) np[kk + 1L] <- op[min(kk, length(op))]
  np <- np / sum(np)
  params$pi_weights[[cc]] <- np
  params$pi_star[[cc]] <- compute_pi_star(np)

  list(state = state, params = params)
}

# ===========================================================================
# relabel_after_change (z_state 동기화 추가)
#   execute_merge 가 S_upper 의 라벨 k 를 target 으로 바꾼 뒤 호출됨.
#   여기서 old_labels(살아남은 라벨) 기준으로 모든 것을 1..K_new 로 재라벨.
#   z_state 도 old_labels 순서대로 보존.
# ===========================================================================
relabel_after_change <- function(cc, state, params, model) {
  T_len <- model$T_len

  old_labels <- sort(unique(state$S_upper[cc, ]))
  K_new      <- length(old_labels)
  max_old <- max(old_labels)
  lut <- integer(max_old)                 # lut[old_label] = new_label
  lut[old_labels] <- seq_len(K_new)

  # upper (vectorized)
  su <- state$S_upper[cc, ]
  state$S_upper[cc, ] <- lut[su]
  state$K[cc]           <- K_new
  state$tau_upper[[cc]] <- extract_changepoints(state$S_upper[cc, ])

  j_in_c <- which(state$cluster == cc)
  for (j in j_in_c) {
    sl <- state$S_lower[j, ]
    mapped <- integer(length(sl))
    in_range <- sl >= 1L & sl <= max_old
    present <- in_range; present[in_range] <- lut[sl[in_range]] > 0L
    mapped[present] <- lut[sl[present]]
    if (any(!present)) mapped[!present] <- pmin(pmax(sl[!present], 1L), K_new)
    sl <- cummax(mapped)                              # monotone nondecreasing
    state$S_lower[j, ] <- sl
    state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])
  }

  # atoms (펼친 뷰): old_labels 순서대로 보존
  old_atoms <- params$atoms[[cc]]
  new_atoms <- vector("list", K_new)
  for (i in seq_along(old_labels)) {
    ok <- old_labels[i]
    new_atoms[[i]] <- if (ok <= length(old_atoms) && !is.null(old_atoms[[ok]])) {
      old_atoms[[ok]]
    } else sample_atom_from_base(params)
  }
  params$atoms[[cc]] <- new_atoms

  # [동기화] z_state: old_labels 순서대로 보존 (살아남은 상태의 atom 인덱스 유지)
  if (!is.null(params$z_state) && !is.null(params$z_state[[cc]])) {
    oz <- params$z_state[[cc]]
    nz <- integer(K_new)
    for (i in seq_along(old_labels)) {
      ok <- old_labels[i]
      nz[i] <- if (ok <= length(oz)) oz[ok] else 1L
    }
    # merge 로 인접이 바뀌어 제약(z_k≠z_{k±1})이 깨질 수 있으나,
    # 갈래 2: 직후 update_atoms(B2)가 인접 제약 하에 재배출하여 교정.
    params$z_state[[cc]] <- nz
  }

  # pi_weights / pi_star (원본과 동일)
  op <- params$pi_weights[[cc]]
  np <- numeric(K_new)
  for (i in seq_along(old_labels)) {
    np[i] <- if (old_labels[i] <= length(op)) op[old_labels[i]] else 1 / K_new
  }
  np <- np / sum(np)
  params$pi_weights[[cc]] <- np
  params$pi_star[[cc]]    <- compute_pi_star(np)

  list(state = state, params = params)
}

cat("07_interval_adjusting.R loaded (v11: 2-stage spline+conditional atom proposal + sequential merge p-/p/p+ + Poisson prior on K).\n")
