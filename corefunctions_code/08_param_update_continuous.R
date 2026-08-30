###############################################################################
# 08_param_update_continuous.R
# Step 3: Parameter Updating — Continuous Response Model
# Order: alpha -> beta -> lambda2 -> C_j,p -> tau2,eta2,theta0 -> 
#        r_j,kr,xi_r -> theta_j -> gamma_tilde -> v -> xi -> theta_out -> phi -> nu
###############################################################################

#' One complete Parameter Updating sweep (continuous model)
#'
#' @param state State list
#' @param params Parameter list
#' @param Y Data matrix (J x T)
#' @param precomp Precomputed basis
#' @param model Model specification
#' @return Updated params list
param_update_continuous <- function(state, params, Y, precomp, model) {
  
  # Shared cache for alpha/beta updates:
  #   base_j[t]   = gamma_jt + f_jt  (does not depend on alpha/beta)
  #   sigma2_j[t] = v_j,S[j,t] * phi_jt^xi_jt
  # This preserves the exact Bayesian updates while avoiding repeated inner-loop
  # recomputation across alpha/beta.
  cache_ab <- build_continuous_ab_cache(state, params, Y, precomp, model)

  # (i) Global intercept alpha
  params <- update_alpha_continuous(state, params, Y, precomp, model, cache_ab)
  
  # (i) Cluster effects beta_j
  params <- update_beta_continuous(state, params, Y, precomp, model, cache_ab)
  
  # (i) Cluster shrinkage variance lambda2_c
  params <- update_lambda2(state, params, model)
  
  # (ii) Cluster allocation C_j
  if (isTRUE(model$fixed_clusters)) {
    # [PATCH] ── 클러스터 완전 고정 모드 ──
    # C_j가 분석가에 의해 사전 지정됨. 할당 업데이트 없음.
    # α, β, θ0, τ², η² 등 클러스터 수준 파라미터만 갱신됨.
    # 클러스터 간 정보 공유는 α(전역 절편)과 hyperparameter를 통해 유지.
  } else if (isTRUE(model$fixed_cluster_sizes)) {
    # ── 클러스터 크기 고정 모드 (Group Swap) ──
    # [PATCH v2] PU 안에서 swap하지 않음.
    # Group swap은 IA round에서만 수행 (11_mcmc_main.R에서 호출).
    # 이유: 매 PU마다 swap하면 K 진동 + 파라미터 적응 파괴.
  } else {
    # ── 클러스터 크기 자유 모드 ──
    # 단일 사이트 Gibbs + 가중치 갱신 (swap은 IA round에서)
    result <- update_cluster_allocation(state, params, Y, precomp, model)
    state <- result$state
    params <- result$params
    params <- update_cluster_weights(state, params, model)
  }
  
  # (iii) Smoothing variances tau2_{m,c}
  params <- update_tau2(state, params, model)
  
  # (iii) Upper GP variance eta2_c
  params <- update_eta2(state, params, model)
  
  # (iv) Upper GP coefficients theta_{0,m,c}
  params <- update_theta0(state, params, model)
  
  
  
  # (iii) Hyperparameters kr, xi_r (slice sampling)
  params <- update_kr_xir(state, params, model)
  
  # (iv) Lower GP coefficients theta_j (ESS)
  params <- update_theta_j(state, params, Y, precomp, model)
  
  # (iii) Smoothness parameter r_j (slice sampling)
  params <- update_rj(state, params, model)

  # (vi) Regime-specific DP atoms (gamma1, gamma2, shape_beta, shape_gamma)
  params <- update_atoms(state, params, Y, precomp, model)

  # (vi-b) Shape strength rate hyperparameters r_β, r_γ (conjugate Gibbs)
  #   atom 의 곡률/선형 강도 스케일을 데이터에 맞춰 추론 (그룹별 이질성 흡수).
  #   shape a_β=a_γ 는 고정, rate 만 업데이트. update_atoms 직후 = 최신 atom 값 사용.
  params <- update_shape_rate_hyper(state, params, model)
  
  # (va) Series-specific state intercepts gamma^*_{j,k}
  params <- update_state_intercepts_continuous(state, params, Y, precomp, model)
  
  # (va-b) Post-sweep centering: α+β_c+γ*_j ridge 제거
  # Papaspiliopoulos et al. (2007, Stat. Sci.) partially centered parameterization.
  # μ_{j,t} 불변, β₁=0 제약 유지, σ²_c 양의 피드백 루프 차단.
  params <- center_gamma_sweep(state, params, model)
  
  # (vb) RW shrinkage variance sigma2_gamma (conjugate Gibbs)
  params <- update_sigma2_gamma(state, params, model)
  
  # (vi) Robust error block: v, xi, theta_out, phi, nu
  base_cache <- build_continuous_base_cache(state, params, Y, precomp, model)
  params <- update_v_continuous(state, params, Y, precomp, model, base_cache)
  params <- update_xi_continuous(state, params, Y, precomp, model, base_cache)
  params <- update_theta_out(state, params, model)
  params <- update_phi_continuous(state, params, Y, precomp, model, base_cache)
  params <- update_nu(state, params, model)
  
  # DP stick-breaking weights
  params <- update_stick_breaking(state, params, model)
  
  # 클러스터 할당 변경이 반영된 state와 params를 모두 반환
  list(state = state, params = params)
}


###############################################################################
# Shared continuous-model caches (Bayesian logic unchanged; speed only)
###############################################################################

compute_series_base_sigma_entry <- function(j, state, params, Y, precomp, model) {
  c_j <- state$cluster[j]
  K_c <- state$K[c_j]

  if (exists("compute_series_base_sigma_fast", mode = "function", inherits = TRUE)) {
    out <- compute_series_base_sigma_fast(j, state, params, precomp, model)
    base_j <- as.numeric(out$base)
    sigma2_j <- as.numeric(out$sigma2)
  } else {
    f_j <- compute_f_all_timepoints(j, state, params, precomp, model)
    gamma_j <- compute_gamma_all_timepoints(j, state, params, model)
    base_j <- as.numeric(gamma_j + f_j)
    # [v_j 스칼라화] get_effective_variance_all이 내부적으로 스칼라 v_j 사용
    sigma2_j <- as.numeric(get_effective_variance_all(j, state, params))
  }

  if (length(base_j) != model$T_len) base_j <- rep(0, model$T_len)
  if (length(sigma2_j) != model$T_len) sigma2_j <- rep(1, model$T_len)
  base_j[!is.finite(base_j)] <- 0
  sigma2_j[!is.finite(sigma2_j) | sigma2_j <= 0] <- 1
  list(base = base_j, sigma2 = sigma2_j)
}

build_continuous_ab_cache <- function(state, params, Y, precomp, model) {
  J <- model$J
  cache <- vector("list", J)
  for (j in 1:J) cache[[j]] <- compute_series_base_sigma_entry(j, state, params, Y, precomp, model)
  cache
}

build_continuous_base_cache <- function(state, params, Y, precomp, model) {
  lapply(build_continuous_ab_cache(state, params, Y, precomp, model), `[[`, "base")
}

###############################################################################
# (i) Fixed effects
###############################################################################

update_alpha_continuous <- function(state, params, Y, precomp, model, cache_ab = NULL) {
  J <- model$J
  hyper <- params$hyper
  if (is.null(cache_ab)) cache_ab <- build_continuous_ab_cache(state, params, Y, precomp, model)

  sum_inv_sig2 <- 0
  sum_resid <- 0

  for (j in 1:J) {
    c_j <- state$cluster[j]
    base_j <- cache_ab[[j]]$base
    sigma2_j <- cache_ab[[j]]$sigma2
    resid_j <- Y[j, ] - params$beta[c_j] - base_j

    sum_inv_sig2 <- sum_inv_sig2 + sum(1 / sigma2_j)
    sum_resid <- sum_resid + sum(resid_j / sigma2_j)
  }

  V_alpha <- 1 / (1 / hyper$sigma2_alpha + sum_inv_sig2)
  m_star <- V_alpha * (hyper$m_alpha / hyper$sigma2_alpha + sum_resid)

  if (!is.finite(V_alpha) || V_alpha <= 0) {
    V_alpha <- hyper$sigma2_alpha
  }
  if (!is.finite(m_star)) {
    m_star <- hyper$m_alpha
  }

  params$alpha <- rnorm(1, m_star, sqrt(V_alpha))
  if (!is.finite(params$alpha)) params$alpha <- m_star
  params
}

update_beta_continuous <- function(state, params, Y, precomp, model, cache_ab = NULL) {
  C <- model$C
  if (is.null(cache_ab)) cache_ab <- build_continuous_ab_cache(state, params, Y, precomp, model)

  for (cc in 1:C) {
    # [식별 제약] β₁ = 0: cluster 1이 기준 수준, α가 흡수
    if (cc == 1) {
      params$beta[1] <- 0
      next
    }
    
    j_in_c <- which(state$cluster == cc)
    if (length(j_in_c) == 0) next

    lambda2_c <- params$lambda2[cc]
    sum_inv_sig2 <- 0
    sum_resid <- 0

    for (j in j_in_c) {
      base_j <- cache_ab[[j]]$base
      sigma2_j <- cache_ab[[j]]$sigma2
      resid_j <- Y[j, ] - params$alpha - base_j

      sum_inv_sig2 <- sum_inv_sig2 + sum(1 / sigma2_j)
      sum_resid <- sum_resid + sum(resid_j / sigma2_j)
    }

    V_beta <- 1 / (1 / lambda2_c + sum_inv_sig2)
    m_star <- V_beta * sum_resid

    if (!is.finite(V_beta) || V_beta <= 0) V_beta <- lambda2_c
    if (!is.finite(m_star)) m_star <- 0

    params$beta[cc] <- rnorm(1, m_star, sqrt(V_beta))
    if (!is.finite(params$beta[cc])) params$beta[cc] <- m_star
  }
  params
}

update_lambda2 <- function(state, params, model) {
  C <- model$C
  hyper <- params$hyper

  for (cc in 1:C) {
    shape <- hyper$a0_beta + 0.5
    scale <- hyper$b0_beta + params$beta[cc]^2 / 2
    params$lambda2[cc] <- rinvgamma(1, shape, scale)
  }
  params
}

###############################################################################
# (ii) Cluster allocation and weights
###############################################################################

update_cluster_allocation <- function(state, params, Y, precomp, model) {
  J <- model$J; C <- model$C; M <- model$M

  for (j in 1:J) {
    log_w <- numeric(C)
    
    # [백업] 현재 상태의 파라미터들을 모두 안전하게 저장
    old_cluster <- state$cluster[j]
    old_S_lower <- state$S_lower[j, ]
    old_gamma <- params$gamma[[j]][1]
    old_theta <- params$theta[j, ]

    for (cc in 1:C) {
      # 1. DP 사전확률 가중치
      log_w[cc] <- log(params$p[cc])

      # 2. 뼈대(Theta) 교체 및 공정한 Prior 평가
      # 새 클러스터 cc의 표준 뼈대(theta0)를 임시로 차용
      params$theta[j, ] <- params$theta0[cc, ]

      # theta_j를 theta0_cc로 두었으므로, P(theta_j | theta0_cc)의 
      # 이차항(거리)은 0이 되고 정규화 상수만 남음.
      # [θ_{j,0}=0 고정] m=0 제외
      for (mm in 1:M) {
        var_m <- params$tau2[cc, mm + 1] * exp(-mm * params$r[j])
        if (var_m <= 0 || !is.finite(var_m)) var_m <- 1e-10
        log_w[cc] <- log_w[cc] - 0.5 * log(2 * pi * var_m)
      }

      # 3. 임시 할당
      state$cluster[j] <- cc
      state$S_lower[j, ] <- state$S_upper[cc, ]

      # 4. [γ*(c,j)] 스칼라 임시 절편: 전체 잔차 평균
      f_j_temp <- compute_f_all_timepoints(j, state, params, precomp, model)
      resid_temp <- Y[j, ] - params$alpha - params$beta[cc] - f_j_temp
      new_gamma <- mean(resid_temp)
      if (!is.finite(new_gamma)) new_gamma <- 0
      params$gamma[[j]] <- new_gamma
      # [v_j 스칼라화] v는 시리즈 단위 스칼라이므로 K 변화에 따른 resize 불필요

      # 5. 공정한 우도(Likelihood) 채점
      tryCatch({
        if (model$type == "continuous") {
          ll_val <- log_lik_series_continuous(j, Y, state, params, precomp, model)
          if (is.finite(ll_val)) log_w[cc] <- log_w[cc] + ll_val else log_w[cc] <- -Inf
        } else {
          ll_val <- log_lik_series_count(j, Y, state, params, precomp, model)
          if (is.finite(ll_val)) log_w[cc] <- log_w[cc] + ll_val else log_w[cc] <- -Inf
        }
      }, error = function(e) {
        log_w[cc] <<- -Inf
      })
    }

    # 6. 평가 종료 후 원상 복구 (샘플링 전 필수)
    params$theta[j, ] <- old_theta
    params$gamma[[j]] <- old_gamma
    state$cluster[j] <- old_cluster
    state$S_lower[j, ] <- old_S_lower

    # 7. 새로운 클러스터 확률적 샘플링
    new_c <- sample_categorical_log(log_w)

    # 8. [핵심] 실제로 이사를 가는 경우 파라미터 완전 갱신
    if (new_c != state$cluster[j]) {
      state$cluster[j] <- new_c
      state$S_lower[j, ] <- state$S_upper[new_c, ]
      state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])

      # 새집의 뼈대와 가구를 완전히 새로 부여 (과거의 잔재 삭제)
      params$theta[j, ] <- params$theta0[new_c, ]

      K_new <- state$K[new_c]

      f_j_temp <- compute_f_all_timepoints(j, state, params, precomp, model)
      resid_temp <- Y[j, ] - params$alpha - params$beta[new_c] - f_j_temp
      new_gamma <- mean(resid_temp)
      if (!is.finite(new_gamma)) new_gamma <- 0
      params$gamma[[j]] <- new_gamma
      # [v_j 스칼라화] v는 시리즈 단위 스칼라이므로 K 변화에 따른 resize 불필요
    }
  }

  list(state = state, params = params)
}

###############################################################################
# (ii-b) Pair Swap Block Update for Cluster Allocation
#
# 단일 사이트 Gibbs의 자유 에너지 장벽을 극복하기 위한 MH 블록 업데이트.
# 서로 다른 클러스터의 두 시계열을 동시에 교환하여,
# 중간 상태(한 시계열만 빈 클러스터에 홀로 이동)를 건너뜁니다.
###############################################################################

#' 시계열 j를 클러스터 cc에 임시 배치했을 때의 log-weight 계산
#' (update_cluster_allocation 내부 로직과 동일)
#'
#' @return log-weight (log-prior + log-likelihood)
eval_series_in_cluster_logw <- function(j, cc, state, params, Y, precomp, model) {
  M <- model$M

  # 1. theta prior 정규화 상수 (theta=theta0이므로 이차항=0)
  # [θ_{j,0}=0 고정] m=0 제외
  lw <- 0
  for (mm in 1:M) {
    var_m <- params$tau2[cc, mm + 1] * exp(-mm * params$r[j])
    if (var_m <= 0 || !is.finite(var_m)) var_m <- 1e-10
    lw <- lw - 0.5 * log(2 * pi * var_m)
  }

  # 2. 임시 파라미터 설정
  old_theta <- params$theta[j, ]
  old_gamma <- params$gamma[[j]][1]
  old_cl    <- state$cluster[j]
  old_sl    <- state$S_lower[j, ]

  params$theta[j, ] <- params$theta0[cc, ]
  state$cluster[j]  <- cc
  state$S_lower[j, ] <- state$S_upper[cc, ]

  # 3. [γ*(c,j)] 스칼라 임시 gamma 계산
  f_j <- compute_f_all_timepoints(j, state, params, precomp, model)
  resid <- Y[j, ] - params$alpha - params$beta[cc] - f_j
  new_gamma <- mean(resid)
  if (!is.finite(new_gamma)) new_gamma <- 0
  params$gamma[[j]] <- new_gamma
  # [v_j 스칼라화] v는 시리즈 단위 스칼라이므로 K 변화에 따른 resize 불필요

  # 4. 우도 계산
  ll <- tryCatch({
    if (model$type == "continuous") {
      log_lik_series_continuous(j, Y, state, params, precomp, model)
    } else {
      log_lik_series_count(j, Y, state, params, precomp, model)
    }
  }, error = function(e) -Inf)
  if (!is.finite(ll)) ll <- -Inf

  # 5. 원상 복구
  params$theta[j, ]  <- old_theta
  params$gamma[[j]]  <- old_gamma
  state$cluster[j]   <- old_cl
  state$S_lower[j, ] <- old_sl

  lw + ll
}

#' 시계열 j를 클러스터 new_c로 실제 전환 (파라미터 갱신)
apply_cluster_switch <- function(j, new_c, state, params, Y, precomp, model) {
  state$cluster[j]   <- new_c
  state$S_lower[j, ] <- state$S_upper[new_c, ]
  state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])

  K_new <- state$K[new_c]

  params$theta[j, ] <- params$theta0[new_c, ]

  f_j <- compute_f_all_timepoints(j, state, params, precomp, model)
  resid <- Y[j, ] - params$alpha - params$beta[new_c] - f_j
  new_gamma <- mean(resid)
  if (!is.finite(new_gamma)) new_gamma <- 0
  params$gamma[[j]] <- new_gamma
  # [v_j 스칼라화] v는 시리즈 단위 스칼라이므로 K 변화에 따른 resize 불필요

  list(state = state, params = params)
}

#' Pair Swap Block Update
#'
#' 서로 다른 클러스터의 두 시계열을 동시에 교환하는 MH 제안.
#' 제안이 대칭(i↔j)이므로 수용비는 순수 우도비로 결정됩니다.
#'
#' @param n_attempts 시도 횟수 (기본: J)
update_cluster_pair_swap <- function(state, params, Y, precomp, model,
                                     n_attempts = NULL) {
  J <- model$J; C <- model$C
  if (C < 2) return(list(state = state, params = params))
  if (is.null(n_attempts)) n_attempts <- J

  for (att in 1:n_attempts) {
    # 1. 서로 다른 클러스터의 두 시계열 무작위 선택
    i <- sample(J, 1)
    ci <- state$cluster[i]

    other_idx <- which(state$cluster != ci)
    if (length(other_idx) == 0) next  # 모든 시계열이 같은 클러스터
    j <- if (length(other_idx) == 1) other_idx else sample(other_idx, 1)
    cj <- state$cluster[j]

    # 2. 현재 배치의 log-weight
    lw_i_ci <- eval_series_in_cluster_logw(i, ci, state, params, Y, precomp, model)
    lw_j_cj <- eval_series_in_cluster_logw(j, cj, state, params, Y, precomp, model)

    # 3. 교환 배치의 log-weight
    lw_i_cj <- eval_series_in_cluster_logw(i, cj, state, params, Y, precomp, model)
    lw_j_ci <- eval_series_in_cluster_logw(j, ci, state, params, Y, precomp, model)

    # 4. MH 수용비 (대칭 제안이므로 proposal ratio = 1)
    log_alpha <- (lw_i_cj + lw_j_ci) - (lw_i_ci + lw_j_cj)

    # 5. 수용/기각
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      res <- apply_cluster_switch(i, cj, state, params, Y, precomp, model)
      state <- res$state; params <- res$params
      res <- apply_cluster_switch(j, ci, state, params, Y, precomp, model)
      state <- res$state; params <- res$params
    }
  }

  list(state = state, params = params)
}

###############################################################################
# Group-Level Cluster Swap (IA round 전용)
#
# 시리즈 단위 pair swap 대신, 분석가가 지정한 Group 전체를 swap.
# 
# 장점:
#   1. Group coherence 구조적 보장 → coherence 필터 불필요
#   2. K 진동 방지: Group 내 시리즈가 분열되지 않음
#   3. K 고착 탈출: 과분할된 cluster에서 다른 K 구조로 전체 리셋
#
# 설계:
#   - IA round에서만 호출 (매 PU가 아님)
#   - C(C-1)/2 개의 Group 쌍 중 n_attempts 쌍을 무작위 시도
#   - 대칭 MH: log α = Σ_{j∈g1} lw(j→c2) + Σ_{j∈g2} lw(j→c1)
#                     - Σ_{j∈g1} lw(j→c1) - Σ_{j∈g2} lw(j→c2)
###############################################################################

#' Group-level cluster swap
#'
#' @param state State list
#' @param params Parameter list
#' @param Y Data matrix (J x T)
#' @param precomp Precomputed basis
#' @param model Model specification (must have model$group_members)
#' @param n_attempts Number of group-pair swap attempts (default: C)
#' @return list(state, params) with updated cluster assignments
update_cluster_group_swap <- function(state, params, Y, precomp, model,
                                      n_attempts = NULL) {
  C <- model$C
  if (C < 2) return(list(state = state, params = params))
  
  # Group 정보: model$group_members[[g]] = 그룹 g의 시리즈 인덱스 벡터
  group_members <- model$group_members
  if (is.null(group_members)) {
    # Fallback: 균등 분할
    n_per <- model$J %/% C
    group_members <- lapply(1:C, function(g) {
      start <- (g - 1) * n_per + 1
      end <- if (g < C) g * n_per else model$J
      start:end
    })
  }
  n_groups <- length(group_members)
  if (n_groups < 2) return(list(state = state, params = params))
  
  if (is.null(n_attempts)) n_attempts <- n_groups  # 기본: 그룹 수만큼 시도
  
  for (att in 1:n_attempts) {
    # 1. 서로 다른 클러스터에 속한 두 그룹 선택
    g1 <- sample(n_groups, 1)
    c1 <- state$cluster[group_members[[g1]][1]]  # g1의 현재 클러스터
    
    # g1과 다른 클러스터에 속한 그룹 후보
    other_groups <- which(sapply(1:n_groups, function(g) {
      state$cluster[group_members[[g]][1]] != c1
    }))
    if (length(other_groups) == 0) next
    g2 <- if (length(other_groups) == 1) other_groups else sample(other_groups, 1)
    c2 <- state$cluster[group_members[[g2]][1]]
    
    # 2. 현재 배치의 log-weight: g1→c1, g2→c2
    lw_current <- 0
    for (j in group_members[[g1]]) {
      lw_current <- lw_current + eval_series_in_cluster_logw(j, c1, state, params, Y, precomp, model)
    }
    for (j in group_members[[g2]]) {
      lw_current <- lw_current + eval_series_in_cluster_logw(j, c2, state, params, Y, precomp, model)
    }
    
    # 3. 교환 배치의 log-weight: g1→c2, g2→c1
    lw_swap <- 0
    for (j in group_members[[g1]]) {
      lw_swap <- lw_swap + eval_series_in_cluster_logw(j, c2, state, params, Y, precomp, model)
    }
    for (j in group_members[[g2]]) {
      lw_swap <- lw_swap + eval_series_in_cluster_logw(j, c1, state, params, Y, precomp, model)
    }
    
    # 4. MH 수용비 (대칭 제안)
    log_alpha <- lw_swap - lw_current
    
    # 5. 수용/기각
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      # g1의 모든 시리즈를 c2로
      for (j in group_members[[g1]]) {
        res <- apply_cluster_switch(j, c2, state, params, Y, precomp, model)
        state <- res$state; params <- res$params
      }
      # g2의 모든 시리즈를 c1으로
      for (j in group_members[[g2]]) {
        res <- apply_cluster_switch(j, c1, state, params, Y, precomp, model)
        state <- res$state; params <- res$params
      }
    }
  }
  
  list(state = state, params = params)
}

update_cluster_weights <- function(state, params, model) {
  C <- model$C
  hyper <- params$hyper
  
  alpha0 <- hyper$alpha0_dirichlet
  n_c <- tabulate(state$cluster, nbins = C)
  
  params$p <- rdirichlet(alpha0 + n_c)
  params
}

###############################################################################
# (iii) Smoothing variances and GP hyperparameters
###############################################################################

update_tau2 <- function(state, params, model) {
  C <- model$C; M <- model$M
  hyper <- params$hyper
  
  for (cc in 1:C) {
    j_in_c <- which(state$cluster == cc)
    nc <- length(j_in_c)
    
    # [θ_{j,0}=0 고정] m=0은 θ=θ₀=0이므로 τ² 업데이트 불필요, m=1부터 시작
    for (mm in 1:M) {
      ss <- 0
      for (j in j_in_c) {
        ss <- ss + exp(mm * params$r[j]) * 
              (params$theta[j, mm + 1] - params$theta0[cc, mm + 1])^2
      }
      
      shape <- hyper$a_tau0 + nc / 2
      scale <- hyper$b_tau0 + ss / 2
      
      params$tau2[cc, mm + 1] <- rinvgamma(1, shape, scale)
      if (!is.finite(params$tau2[cc, mm + 1]) || params$tau2[cc, mm + 1] <= 1e-10) {
        params$tau2[cc, mm + 1] <- 1e-10
      }
    }
  }
  params
}

update_eta2 <- function(state, params, model) {
  C <- model$C; M <- model$M
  hyper <- params$hyper
  
  for (cc in 1:C) {
    ss <- 0
    # [θ_{j,0}=0 고정] m=0은 θ₀=0이므로 기여 없음, m=1부터
    for (mm in 1:M) {
      prior_var_scale <- (1 + mm / params$xi_r)^params$kr
      ss <- ss + params$theta0[cc, mm + 1]^2 * prior_var_scale
    }
    
    shape <- hyper$a_eta0 + M / 2  # M+1 → M (m=0 제외)
    scale <- hyper$b_eta0 + ss / 2
    
    params$eta2[cc] <- rinvgamma(1, shape, scale)
  }
  params
}

update_theta0 <- function(state, params, model) {
  C <- model$C; M <- model$M
  
  for (cc in 1:C) {
    j_in_c <- which(state$cluster == cc)
    
    # [θ_{j,0}=0 고정] m=0 절편은 0으로 고정
    params$theta0[cc, 1] <- 0
    
    for (mm in 1:M) {
      # Prior variance for theta_{0,m,c}
      prior_var <- params$eta2[cc] * (1 + mm / params$xi_r)^(-params$kr)
      
      # Data contribution
      sum_data <- 0
      sum_prec <- 0
      for (j in j_in_c) {
        var_jm <- params$tau2[cc, mm + 1] * exp(-mm * params$r[j])
        
        # [방어벽 1] var_jm이 0으로 언더플로우 되어 1/0 = Inf가 되는 것을 원천 차단!
        var_jm <- max(var_jm, 1e-10) 
        
        sum_prec <- sum_prec + 1 / var_jm
        sum_data <- sum_data + params$theta[j, mm + 1] / var_jm
      }
      
      V_post <- 1 / (1 / prior_var + sum_prec)
      m_post <- V_post * sum_data
      
      # [방어벽 2] NaN 전염 차단
      if (!is.finite(V_post) || V_post <= 0) V_post <- prior_var
      if (!is.finite(m_post)) m_post <- 0
      
      params$theta0[cc, mm + 1] <- rnorm(1, m_post, sqrt(V_post))
    }
  }
  params
}

update_rj <- function(state, params, model) {
  J <- model$J; M <- model$M
  
  for (j in 1:J) {
    c_j <- state$cluster[j]
    
    log_target <- function(rj) {
      if (rj <= 0) return(-Inf)
      
      # Prior: r_j ~ Gamma(shape = kr, rate = xi_r)   [BUGFIX: was scale]
      #   Lenk 정합성: E[exp(-m r_j)] = (1 + m/xi_r)^(-kr) 가 theta0 의
      #   algebraic decay (1 + m/xi_r)^(-kr) 와 일치하려면 xi_r 은 RATE 여야 함.
      #   (rate 모수화에서 -rj * xi_r;  이전의 -rj / xi_r 은 scale 해석으로 불일치)
      log_prior <- (params$kr - 1) * log(rj) - rj * params$xi_r
      
      # Likelihood of theta_{j,m} given r_j
      # [θ_{j,0}=0 고정] m=0은 상수항이므로 r_j에 의존하지 않음, m=1부터
      log_lik <- 0
      for (mm in 1:M) {
        var_m <- params$tau2[c_j, mm + 1] * exp(-mm * rj)
        log_lik <- log_lik + dnorm_log_var(
          params$theta[j, mm + 1], params$theta0[c_j, mm + 1], var_m)
      }
      
      log_prior + log_lik
    }
    
    params$r[j] <- slice_sample_positive(params$r[j], log_target, w = 0.5, m = 10)
  }
  params
}

update_kr_xir <- function(state, params, model) {
  J <- model$J; C <- model$C; M <- model$M
  
  # Update kr
  log_target_kr <- function(kr) {
    if (kr <= 0) return(-Inf)
    
    log_val <- dgamma(kr, shape = 2, rate = 1, log = TRUE)  # prior
    
    for (j in 1:J) {
      log_val <- log_val + dgamma(params$r[j], shape = kr,
                                   rate = params$xi_r, log = TRUE)   # [BUGFIX] scale→rate
    }
    
    for (cc in 1:C) {
      # [θ_{j,0}=0 고정] m=0 제외
      for (mm in 1:M) {
        prior_var <- params$eta2[cc] * (1 + mm / params$xi_r)^(-kr)
        log_val <- log_val + dnorm_log_var(params$theta0[cc, mm + 1], 0, prior_var)
      }
    }
    
    log_val
  }

  # Update xi_r
  log_target_xir <- function(xir) {
    if (xir <= 0) return(-Inf)

    log_val <- dgamma(xir, shape = 2, rate = 1, log = TRUE)  # prior

    for (j in 1:J) {
      log_val <- log_val + dgamma(params$r[j], shape = params$kr,
                                   rate = xir, log = TRUE)   # [BUGFIX] scale→rate
    }

    for (cc in 1:C) {
      # [θ_{j,0}=0 고정] m=0 제외
      for (mm in 1:M) {
        prior_var <- params$eta2[cc] * (1 + mm / xir)^(-params$kr)
        log_val <- log_val + dnorm_log_var(params$theta0[cc, mm + 1], 0, prior_var)
      }
    }

    log_val
  }

  # ── ξ_r 만 학습 (k_r = 1.5 고정) ──
  #   k_r, ξ_r 동시학습은 강한 능선(상관~0.89, 축 = kr/ξ_r = E[r_j])을 이뤄
  #   혼합이 느림. 우선 ξ_r(스케일, 데이터 정보 많음)만 풀어 능선을 회피한다.
  #   k_r(다항 차수)은 약식별이라 1.5 고정. ξ_r 학습이 안정되면 추후 k_r 도
  #   풀되 그때는 블록 MH(둘 동시) 권장.
  #   target 은 매 호출 params$kr / params$xi_r 를 새로 읽으므로 정합 유지.
  params$kr <- 1.5
  params$xi_r <- slice_sample_positive(params$xi_r, log_target_xir, w = 0.5, m = 10)
  params
}

###############################################################################
# (iv) GP coefficients
###############################################################################

update_theta_j <- function(state, params, Y, precomp, model) {
  J <- model$J; M <- model$M
  
  for (j in 1:J) {
    c_j <- state$cluster[j]
    
    # [θ_{j,0}=0 고정] ESS는 m=1:M 성분에서만 수행 (M차원)
    # m=0 절편은 α + β + γ*가 담당하므로 GP에서 제거
    params$theta[j, 1] <- 0  # 명시적 고정
    
    prior_mean_sub <- params$theta0[c_j, 2:(M + 1)]
    prior_cov_diag_sub <- numeric(M)
    for (mm in 1:M) {
      prior_cov_diag_sub[mm] <- params$tau2[c_j, mm + 1] * exp(-mm * params$r[j])
      if (!is.finite(prior_cov_diag_sub[mm]) || prior_cov_diag_sub[mm] <= 0) {
        prior_cov_diag_sub[mm] <- 1e-10
      }
    }
    
    # Log-likelihood function: ESS 후보 = m=1:M 부분벡터
    log_lik_fn <- function(theta_sub_new) {
      old_theta <- params$theta[j, ]
      params$theta[j, ] <- c(0, theta_sub_new)  # m=0은 항상 0
      
      ll <- 0
      if (model$type == "continuous") {
        ll <- log_lik_series_continuous(j, Y, state, params, precomp, model)
      } else {
        ll <- log_lik_series_count(j, Y, state, params, precomp, model)
      }
      
      params$theta[j, ] <- old_theta  # restore
      ll
    }
    
    theta_sub_result <- elliptical_slice_sample(
      theta_current = params$theta[j, 2:(M + 1)],
      prior_mean = prior_mean_sub,
      prior_cov_diag = prior_cov_diag_sub,
      log_likelihood = log_lik_fn
    )
    params$theta[j, ] <- c(0, theta_sub_result)
  }
  params
}


###############################################################################
# (v) Regime-specific atoms — Truncated Blocked Gibbs (Ishwaran & James, 2001)
#
# Decoupled base distribution H:
#   (δ1, δ2)    ~ Uniform({-1,+1}^2)
#   shape_beta  ~ Gamma(a_shape_beta, b_shape_beta)   [curvature-favoring]
#   shape_gamma ~ Gamma(a_shape_gamma, b_shape_gamma)  [sparsity-inducing]
#
# The state intercept gamma^*_{j,k} is NOT in the DP atom.
# It is updated separately with a Gaussian random-walk shrinkage prior.
#
# For truncated DP with K finite atoms:
#   Empty state:    atom_{c,k} ~ H  (exact Gibbs)
#   Occupied state: π(atom_{c,k}|rest) ∝ H(atom) × ∏ L(y|atom)
#
# Discrete Gibbs (Ishwaran & James 2001):
#   1. Draw N_cand candidate atoms from H (+ include current atom)
#   2. All candidates come from H → H cancels in posterior ratio

#   3. Weight_i ∝ L(data_in_state_k | candidate_i) × H(candidate_i)
#   4. Sample one candidate from this discrete distribution
#
# This is exact Gibbs on a finite support set. Always moves.
# Sign flips in (γ1, γ2) occur naturally since H = N_2(μ0, Σ0)
# generates candidates with all sign combinations.
###############################################################################

###############################################################################
# (vi) DP Atom Update — Ishwaran & James (2001) Blocked Gibbs
#
# Occupied state k에 대해:
#   1. 4개 (δ₁,δ₂) 조합에 대해 (b,g)의 조건부 사후를 WLS로 계산
#   2. (b,g)를 Laplace 근사로 적분 소거 → P(δ₁,δ₂ | data_k) 산출
#   3. 이산 분포에서 (δ₁,δ₂) 샘플링
#   4. 선택된 (δ₁,δ₂) 하에서 (b,g) ~ truncated Normal 사후 샘플링
#
# f^shape(x) = [δ₁·g + b·(δ₁+δ₂)/2]·x - δ₂·b·H_j(x)
#   = g·(δ₁·x) + b·[(δ₁+δ₂)/2·x - δ₂·H_j(x)]
# → (b,g) 선형 → WLS conjugate
###############################################################################

###############################################################################
# update_atoms — Collapsed Gibbs with γ*(c,j) integrated out
#
# 변경 핵심:
#   기존: r_{j,t} = Y - α - β_c - γ*(j)  → γ*에 조건부 고정
#   변경: r_{j,t} = Y - α - β_c           → γ* 미차감
#         시리즈별 가중 centering으로 γ* 해석적 적분 소거
#
# 수학적 근거:
#   γ*_j ~ N(0, σ²_c) 하에서 적분:
#     p(Y_j | δ,b,g) = ∫ ∏_t N(Y_t | α+β_c+γ*+f_t, σ²_t) · N(γ*|0,σ²_c) dγ*
#   = profiled WLS with intercept, 등가적으로:
#     W_j = Σ_t w_t + 1/σ²_c,  w̄(·) = Σ_t w_t(·) / W_j
#     r̃_t = r_t - w̄(r),  z̃_g = z_g - w̄(z_g),  z̃_b = z_b - w̄(z_b)
#     → centered WLS on (r̃, z̃_g, z̃_b)
#
# 효과:
#   δ₁ 선택이 γ* 수준에 묶이지 않고, 구간 내 "형태"만으로 결정됨.
#   split 제안(propose_atom_spline_conditional)의 centering과 일치.
###############################################################################

###############################################################################
# update_atoms — Ishwaran & James (2001) Truncated Blocked Gibbs
#
# 구조 (해석 1: distinct shape atom 풀을 상태들이 공유):
#   - 그룹별 DP: atom_pool[[cc]] 는 길이 N (=truncation level) 고정 풀
#   - 상태(segment) k 는 atom 인덱스 z_state[[cc]][k] ∈ {1..N} 에 배정
#   - 인접 제약: z_k ≠ z_{k-1}, z_k ≠ z_{k+1}  (변화점 = 인접 형상 전환)
#     → 비인접 tie 허용 (G3의 성장↔회복 형상 재사용)
#   - 펼치기: 끝에서 atoms[[cc]][[k]] <- atom_pool[[cc]][[ z_state[[cc]][k] ]]
#     → 외부 코드(eval_shape, split/merge, 시각화)는 atoms[[cc]][[k]] 그대로 사용
#
# 세 블록 (그룹마다):
#   (B1) atom pool 갱신: atom n 에 배정된 모든 상태의 (시점,계열) 데이터를 모아
#        Gaussian evidence 로 (δ1,δ2,b,g) 추출.
#        γ* 는 conditional 로 차감(collapse 안 함). center_gamma_sweep
#        (Papaspiliopoulos 재매개화) 가 ridge 를 이미 처리하므로 중복 방지.
#        부호 evidence = 비절단 Gaussian evidence + log P(b,g>0 | 사후)  [양수 절단 반영]
#        (b,g) = 양수 절단 이변량 정규.  미배정(m_n=0) atom 은 base H.
#   (B2) 배정 z_k 재추출 (순차 in-place):
#        P(z_k=n) ∝ 1[n≠z_{k-1}] 1[n≠z_{k+1}] · π_n · L_k(n)
#        L_k(n) = 고정된 atom_pool[[n]] 의 plug-in 우도 (γ* 조건부, compute_atom_loglik)
#        π_n = V_n ∏_{l<n}(1-V_l)
#   (B3) atom DP stick-breaking:
#        V_n ~ Beta(1+m_n, ς^atom_c + Σ_{l>n} m_l),  V_N=1,  m_n=#{k:z_k=n}
#        ς^atom_c slice 갱신 (변화점 varsigma 와 분리)
#
# 의존: canonicalize_atom, sample_atom_from_base, sample_categorical_log,
#       get_effective_variance_all, compute_atom_loglik, slice_sample_positive
#       mvtnorm (pmvnorm)
###############################################################################

# ── 이변량 정규 양사분면 질량 P(X>0, Y>0), N(mu, Sig) ──
# mvtnorm::pmvnorm 사용. lower=c(0,0), upper=c(Inf,Inf) 가 양사분면 적분.
#   mu = c(g_m, b_m), Sig = 2x2 사후 공분산 (첫 성분 g=shape_gamma, 둘째 b=shape_beta).
.biv_pos_mass <- function(mu, Sig) {
  # 수치 안정: 공분산 대칭화 + 미소 jitter (pmvnorm 의 비양정치 회피)
  Sig <- (Sig + t(Sig)) / 2
  d <- sqrt(pmax(diag(Sig), 1e-12))
  # 양정치 보정: 최소 고유값이 음수면 ridge 추가
  ev <- tryCatch(min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values),
                 error = function(e) -1)
  if (!is.finite(ev) || ev <= 1e-10) {
    Sig <- Sig + diag(rep(1e-8 * max(d)^2 + 1e-12, 2))
  }
  p <- tryCatch(
    as.numeric(mvtnorm::pmvnorm(
      lower = c(0, 0), upper = c(Inf, Inf),
      mean  = as.numeric(mu), sigma = Sig)),
    error = function(e) NA_real_)
  if (!is.finite(p)) {
    # fallback: 독립 근사 (상관 무시) — pmvnorm 실패 시에만
    s1 <- sqrt(max(Sig[1, 1], 1e-12)); s2 <- sqrt(max(Sig[2, 2], 1e-12))
    p <- pnorm(mu[1] / s1) * pnorm(mu[2] / s2)
  }
  max(min(p, 1), 1e-300)
}

# ── 한 데이터 묶음(여러 상태의 시점·계열)에서 atom 갱신 ──
#   series_data: list of list(r, w, x, H, W_j)  — γ* centering 분모 W_j 포함
#   반환: list(atom = canonicalized atom)
.update_one_atom_from_data <- function(series_data, params,
                                        a_b, r_b, a_g, r_g,
                                        b_mode0, g_mode0, b_prec0, g_prec0) {
  if (length(series_data) == 0) {
    return(canonicalize_atom(sample_atom_from_base(params)))
  }
  delta_combos <- list(c(1L, 1L), c(1L, -1L), c(-1L, 1L), c(-1L, -1L))
  log_ev   <- numeric(4)
  bg_modes <- vector("list", 4)
  bg_covs  <- vector("list", 4)

  for (di in 1:4) {
    d1 <- delta_combos[[di]][1]; d2 <- delta_combos[[di]][2]
    r_pool <- c(); w_pool <- c(); zg_pool <- c(); zb_pool <- c()
    for (s in series_data) {
      # [CONDITIONAL] γ* 이미 차감됨 → centering 불필요. raw 설계열 사용.
      z_g_raw <- d1 * s$x
      z_b_raw <- (d1 + d2) / 2 * s$x - d2 * s$H
      r_pool  <- c(r_pool,  s$r)
      w_pool  <- c(w_pool,  s$w)
      zg_pool <- c(zg_pool, z_g_raw)
      zb_pool <- c(zb_pool, z_b_raw)
    }
    wzg <- w_pool * zg_pool; wzb <- w_pool * zb_pool
    P11 <- sum(wzg * zg_pool) + g_prec0
    P12 <- sum(wzg * zb_pool)
    P22 <- sum(wzb * zb_pool) + b_prec0
    detP <- P11 * P22 - P12 * P12
    if (!is.finite(detP) || detP < 1e-20) {
      log_ev[di] <- -1e300
      bg_modes[[di]] <- c(g_mode0, b_mode0)
      bg_covs[[di]]  <- diag(c(1 / g_prec0, 1 / b_prec0))
      next
    }
    invd <- 1 / detP
    S11 <-  P22 * invd; S12 <- -P12 * invd; S22 <- P11 * invd
    wr <- w_pool * r_pool
    rhs_g <- sum(zg_pool * wr) + g_prec0 * g_mode0
    rhs_b <- sum(zb_pool * wr) + b_prec0 * b_mode0
    g_m <- S11 * rhs_g + S12 * rhs_b
    b_m <- S12 * rhs_g + S22 * rhs_b
    Sig <- matrix(c(S11, S12, S12, S22), 2, 2)
    bg_modes[[di]] <- c(g_m, b_m)
    bg_covs[[di]]  <- Sig

    # 정확한 비절단 Gaussian evidence (Gaussian prior N(m0, P0^{-1}) 가정)
    quad_post <- g_m * g_m * P11 + 2 * g_m * b_m * P12 + b_m * b_m * P22
    quad_0    <- g_mode0 * g_mode0 * g_prec0 + b_mode0 * b_mode0 * b_prec0
    log_ev_untrunc <- 0.5 * log(g_prec0 * b_prec0) - 0.5 * log(detP) +
                      0.5 * (quad_post - quad_0) - 0.5 * sum(w_pool * r_pool^2)
    # 양수 절단 보정: + log P(b,g>0 | 사후)
    p_pos <- .biv_pos_mass(c(g_m, b_m), Sig)
    log_ev[di] <- log_ev_untrunc + log(p_pos)
    if (!is.finite(log_ev[di])) log_ev[di] <- -1e300
  }

  # (δ1,δ2) 범주형 추출 (절단보정 evidence)
  chosen <- sample_categorical_log(log_ev)
  d1 <- delta_combos[[chosen]][1]; d2 <- delta_combos[[chosen]][2]
  mu <- bg_modes[[chosen]]; S <- bg_covs[[chosen]]

  # (b,g) 양수 절단 이변량 정규 (2D Cholesky + rejection)
  L11 <- sqrt(max(S[1, 1], 1e-10))
  L21 <- S[2, 1] / L11
  L22 <- sqrt(max(S[2, 2] - L21^2, 1e-10))
  g_new <- max(mu[1], 0.01); b_new <- max(mu[2], 0.01)
  for (try_i in 1:200) {
    z1 <- rnorm(1); z2 <- rnorm(1)
    g_cand <- mu[1] + L11 * z1
    b_cand <- mu[2] + L21 * z1 + L22 * z2
    if (g_cand > 0 && b_cand > 0) { g_new <- g_cand; b_new <- b_cand; break }
  }
  canonicalize_atom(list(gamma1 = d1, gamma2 = d2,
                         shape_beta = b_new, shape_gamma = g_new))
}

###############################################################################
# .update_one_atom_candidate_IJ   [DEAD — 호출되지 않음]
#   (가) 정통 I-J 사후갱신 채택으로 update_atoms 의 점유 slot 은
#   .update_one_atom_from_data (데이터 사후 draw) 를 쓴다.  이 후보기반 변형은
#   점유 slot 을 데이터 사후로 정제하지 못해 split 제안 우도가 낮아지는 문제가
#   있었으므로 더 이상 사용하지 않는다.  참고용으로 정의만 보존.
#
#   Ishwaran & James (2001) truncated blocked-Gibbs ATOM update, candidate form.
#   Replaces the closed-form (Gaussian-working-prior) draw of
#   .update_one_atom_from_data with the canonical I-J atom full-conditional
#
#     atom*_n | rest  ∝  H(atom*_n) · ∏_{k: z_k=n} L(y_k | atom*_n)
#
#   sampled exactly on a finite base-drawn support set:
#     candidate set = { current atom } ∪ { N_cand-1 draws from base H } .
#   Every candidate comes from the SAME base H, so H cancels in the posterior
#   ratio and the categorical weight is the pure segment likelihood
#       w_i ∝ L(merged_data | candidate_i) ,   i = 1, ..., N_cand .
#   Including the current atom makes the move reversible (a high-likelihood
#   incumbent can be retained) and preserves detailed balance.
#   The Gaussian segment log-likelihood uses the SAME conditional residual and
#   raw (un-centered) design columns as .collect_state_series_data supplies
#   (gamma^*_{c,j} already subtracted into s$r), so the likelihood kernel is
#   identical to the old WLS path -- only the way the atom is DRAWN changes.
###############################################################################
.update_one_atom_candidate_IJ <- function(series_data, params, current_atom,
                                          N_cand = 50L) {
  if (length(series_data) == 0) {
    return(canonicalize_atom(sample_atom_from_base(params)))
  }
  N_cand <- max(as.integer(N_cand), 2L)

  # Candidate set: { current } U { N_cand-1 base draws }
  cand <- vector("list", N_cand)
  start <- 1L
  if (!is.null(current_atom)) {
    cand[[1]] <- canonicalize_atom(current_atom)
    start <- 2L
  }
  for (i in start:N_cand) {
    cand[[i]] <- canonicalize_atom(sample_atom_from_base(params))
  }

  # Pooled Gaussian segment log-likelihood per candidate (raw design, no centering)
  log_scores <- numeric(N_cand)
  for (ci in seq_len(N_cand)) {
    at <- cand[[ci]]
    d1 <- at$gamma1; d2 <- at$gamma2
    g_i <- at$shape_gamma; b_i <- at$shape_beta
    ll <- 0
    for (s in series_data) {
      z_g <- d1 * s$x
      z_b <- (d1 + d2) / 2 * s$x - d2 * s$H
      f_i <- g_i * z_g + b_i * z_b
      ll <- ll - 0.5 * sum(s$w * (s$r - f_i)^2)
    }
    log_scores[ci] <- if (is.finite(ll)) ll else -Inf
  }

  chosen <- sample_categorical_log(log_scores)
  canonicalize_atom(cand[[chosen]])
}

# ── 한 상태 k 의 (시점·계열) 데이터 수집 (γ* centering 재료) ──
#   특정 atom 으로 본 게 아니라, 상태 k 의 raw 재료(r, w, x, H, W_j)
.collect_state_series_data <- function(cc, k, j_in_c, state, params, Y, precomp) {
  series_data <- list()
  for (j in j_in_c) {
    t_idx <- which(state$S_lower[j, ] == k)
    if (length(t_idx) == 0) next
    # [CONDITIONAL] γ* 를 조건부로 차감 (collapse 안 함).
    #   center_gamma_sweep (Papaspiliopoulos 재매개화) 가 ridge 를 이미 처리하므로
    #   atom 갱신은 현재 γ*_j 값을 고정해 빼는 conditional 우도를 사용.
    gamma_j <- params$gamma[[j]]
    if (is.null(gamma_j) || length(gamma_j) == 0 || !is.finite(gamma_j[1])) gamma_j <- 0
    r_jt <- as.numeric(Y[j, t_idx] - params$alpha - params$beta[cc] - gamma_j[1])
    sigma2_jt <- get_effective_variance_all(j, state, params)[t_idx]
    sigma2_jt[!is.finite(sigma2_jt) | sigma2_jt <= 0] <- 1
    w_jt <- 1.0 / sigma2_jt
    x_jt <- precomp$x_global[t_idx]
    theta_j <- params$theta[j, ]; Mp1 <- length(theta_j)
    if (Mp1 > 1) {
      theta_shape <- theta_j[2:Mp1]
      SZ_safe <- max(sum(theta_shape^2), 1e-12)
      theta_outer <- tcrossprod(theta_shape)
      H_jt <- numeric(length(t_idx))
      for (ti in seq_along(t_idx)) {
        H_jt[ti] <- sum(theta_outer * precomp$D_at_x[2:Mp1, 2:Mp1, t_idx[ti]]) / SZ_safe
      }
    } else H_jt <- rep(0, length(t_idx))
    series_data[[length(series_data) + 1]] <- list(
      r = r_jt, w = w_jt, x = x_jt, H = H_jt)
  }
  series_data
}

###############################################################################
# update_shape_rate_hyper — Gamma rate hyperparameters (conjugate Gibbs)
#
#   atom 의 곡률/선형 강도:
#     b_{c,k} = shape_beta  ~ Gamma(a_β, r_β)
#     g_{c,k} = shape_gamma ~ Gamma(a_γ, r_γ)
#   shape a_β, a_γ 는 고정. rate r_β, r_γ 만 데이터에 맞춰 업데이트.
#
#   conjugate hyperprior:  r_β ~ Gamma(a0_β, b0_β),  r_γ ~ Gamma(a0_γ, b0_γ)
#   full conditional:
#     r_β | {b} ~ Gamma(a0_β + N_b·a_β,  b0_β + Σ b_{c,k})
#     r_γ | {g} ~ Gamma(a0_γ + N_g·a_γ,  b0_γ + Σ g_{c,k})
#   (N_b = 배정된 distinct atom 수, Σ = 그 atom 들의 shape 합)
#
#   IJ 구조: distinct shape atom 은 atom_pool 의 배정된 슬롯(m_n>0).
#     - 펼친 뷰(atoms[[cc]])는 tie 시 중복 → 사용하지 않음
#     - 미배정 슬롯(m_n=0)은 base H 에서 매 sweep 재추출되는 정보 없는 atom → 제외
#   atom_pool 이 없으면(IJ 미적용) 펼친 뷰 atoms[[cc]] 로 fallback.
###############################################################################
update_shape_rate_hyper <- function(state, params, model) {
  hyper <- params$hyper
  C <- model$C

  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 2.0
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 2.0

  # rate 의 hyperprior 모수 (약정보 Gamma(2,1) 기본)
  a0_beta  <- if (!is.null(hyper$a0_rbeta))  hyper$a0_rbeta  else 2.0
  b0_beta  <- if (!is.null(hyper$b0_rbeta))  hyper$b0_rbeta  else 1.0
  a0_gamma <- if (!is.null(hyper$a0_rgamma)) hyper$a0_rgamma else 2.0
  b0_gamma <- if (!is.null(hyper$b0_rgamma)) hyper$b0_rgamma else 1.0

  # ── distinct atom 의 shape_beta / shape_gamma 수집 ──
  b_vals <- numeric(0)
  g_vals <- numeric(0)

  use_pool <- !is.null(params$atom_pool) && !is.null(params$z_state)

  for (cc in 1:C) {
    if (use_pool && !is.null(params$atom_pool[[cc]]) && !is.null(params$z_state[[cc]])) {
      # IJ: 배정된 distinct 슬롯만 (m_n > 0)
      z <- params$z_state[[cc]]
      assigned <- unique(z[z >= 1 & z <= length(params$atom_pool[[cc]])])
      pool <- params$atom_pool[[cc]]
      for (n in assigned) {
        at <- pool[[n]]
        if (is.null(at)) next
        if (!is.null(at$shape_beta)  && is.finite(at$shape_beta)  && at$shape_beta  > 0)
          b_vals <- c(b_vals, at$shape_beta)
        if (!is.null(at$shape_gamma) && is.finite(at$shape_gamma) && at$shape_gamma > 0)
          g_vals <- c(g_vals, at$shape_gamma)
      }
    } else {
      # fallback: 펼친 뷰 (IJ 미적용 시)
      atoms_c <- params$atoms[[cc]]
      if (is.null(atoms_c)) next
      for (k in seq_along(atoms_c)) {
        at <- atoms_c[[k]]
        if (is.null(at)) next
        if (!is.null(at$shape_beta)  && is.finite(at$shape_beta)  && at$shape_beta  > 0)
          b_vals <- c(b_vals, at$shape_beta)
        if (!is.null(at$shape_gamma) && is.finite(at$shape_gamma) && at$shape_gamma > 0)
          g_vals <- c(g_vals, at$shape_gamma)
      }
    }
  }

  # ── r_β conjugate Gibbs ──
  N_b <- length(b_vals)
  if (N_b > 0) {
    post_shape <- a0_beta + N_b * a_beta
    post_rate  <- b0_beta + sum(b_vals)
    if (is.finite(post_shape) && is.finite(post_rate) && post_rate > 0) {
      r_beta_new <- rgamma(1, shape = post_shape, rate = post_rate)
      if (is.finite(r_beta_new) && r_beta_new > 1e-8) {
        params$hyper$b_shape_beta <- r_beta_new
      }
    }
  }

  # ── r_γ conjugate Gibbs ──
  N_g <- length(g_vals)
  if (N_g > 0) {
    post_shape <- a0_gamma + N_g * a_gamma
    post_rate  <- b0_gamma + sum(g_vals)
    if (is.finite(post_shape) && is.finite(post_rate) && post_rate > 0) {
      r_gamma_new <- rgamma(1, shape = post_shape, rate = post_rate)
      if (is.finite(r_gamma_new) && r_gamma_new > 1e-8) {
        params$hyper$b_shape_gamma <- r_gamma_new
      }
    }
  }

  params
}

update_atoms <- function(state, params, Y, precomp, model) {
  C <- model$C
  hyper <- params$hyper
  a_b <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 5
  r_b <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 0.5
  a_g <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 3
  r_g <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 1.0
  N <- if (!is.null(model$N_atom)) model$N_atom else 20L

  # [(가) 정통 I-J 사후갱신] 점유 slot 은 데이터 사후에서 draw 한다
  #   (.update_one_atom_from_data).  그 함수는 Gamma prior b~Gamma(a_b,r_b),
  #   g~Gamma(a_g,r_g) 를 mode-matched Gaussian working prior 로 근사해 닫힌형
  #   사후를 만든다.  Gamma(a, rate r) 의 mode = (a-1)/r (a>1),  mode 근처 곡률
  #   매칭 precision = (a-1)/mode^2 = r^2/(a-1).  (a<=1 이면 안전한 기본값.)
  .wp_mode_prec <- function(a, r) {
    if (is.finite(a) && a > 1 && is.finite(r) && r > 0) {
      mode <- (a - 1) / r
      prec <- (a - 1) / (mode * mode)        # = r^2/(a-1)
    } else {
      mode <- 1.0; prec <- 1.0
    }
    if (!is.finite(mode) || mode <= 0) mode <- 1.0
    if (!is.finite(prec) || prec <= 0) prec <- 1.0
    c(mode = mode, prec = prec)
  }
  .wp_b <- .wp_mode_prec(a_b, r_b)
  .wp_g <- .wp_mode_prec(a_g, r_g)
  b_mode0 <- .wp_b[["mode"]]; b_prec0 <- .wp_b[["prec"]]
  g_mode0 <- .wp_g[["mode"]]; g_prec0 <- .wp_g[["prec"]]

  # 컨테이너 보장
  if (is.null(params$atom_pool))     params$atom_pool     <- vector("list", C)
  if (is.null(params$z_state))       params$z_state       <- vector("list", C)
  if (is.null(params$atom_V))        params$atom_V        <- vector("list", C)
  if (is.null(params$atom_pi))       params$atom_pi       <- vector("list", C)
  if (is.null(params$atom_varsigma)) params$atom_varsigma <- rep(1.0, C)

  for (cc in 1:C) {
    K_c <- state$K[cc]
    j_in_c <- which(state$cluster == cc)

    # [CONDITIONAL] γ* collapse 제거 → σ²_c 불필요 (atom 갱신에서)

    # ── 풀/배정/막대 초기화 (없거나 길이 안 맞으면) ──
    pool <- params$atom_pool[[cc]]
    if (is.null(pool) || length(pool) != N) {
      pool <- vector("list", N)
      for (n in 1:N) pool[[n]] <- canonicalize_atom(sample_atom_from_base(params))
    }
    z <- params$z_state[[cc]]
    if (is.null(z) || length(z) != K_c) {
      # 인접 제약 만족하는 초기 배정: 1,2,1,2,... 를 풀 앞쪽에 매핑
      z <- integer(K_c)
      if (K_c >= 1) z[1] <- 1L
      if (K_c >= 2) for (k in 2:K_c) z[k] <- if (z[k - 1] == 1L) 2L else 1L
    }
    Vv <- params$atom_V[[cc]]
    if (is.null(Vv) || length(Vv) != N) {
      Vv <- c(rbeta(N - 1, 1, params$atom_varsigma[cc]), 1)
      Vv <- pmin(pmax(Vv, 1e-10), 1 - 1e-10); Vv[N] <- 1
    }
    vs_c <- params$atom_varsigma[cc]
    if (!is.finite(vs_c) || vs_c <= 0) vs_c <- 1.0

    # 상태별 데이터 재료 1회 수집 (B1, B2 공용)
    state_data <- vector("list", K_c)
    for (k in 1:K_c) {
      state_data[[k]] <- .collect_state_series_data(
        cc, k, j_in_c, state, params, Y, precomp)
    }

    # ══════════════ (B1) atom pool 갱신 (정통 I-J 사후갱신) ══════════════
    # Ishwaran & James (2001) truncated blocked Gibbs:
    #   atom_n | rest ∝ H(atom_n) ∏_{k:z_k=n} L_k(atom_n)
    #   · 점유 slot (배정된 상태 있음): 데이터 사후에서 draw
    #       → .update_one_atom_from_data: (δ1,δ2) 4조합의 절단보정 evidence 로
    #         categorical, 선택 조합에서 (b,g) truncated bivariate Gaussian draw.
    #         (Gamma prior 를 mode-matched Gaussian working prior 로 근사한 닫힌형)
    #   · 빈 slot (미배정): 우도항 없음 → 사후=prior → base H 에서 draw
    # [변경] 이전의 후보기반(.update_one_atom_candidate_IJ: {현재}∪{base 50개}에서
    #   고르기)은 점유 slot 을 데이터 사후로 정제하지 못해, pool atom 이 세그먼트에
    #   덜 적합 → split 후보 atom 우도 L_split 가 낮음 → 분할이 잘 안 됨.  정통
    #   사후갱신으로 점유 slot 을 데이터에 적합시켜 split 제안 우도를 높인다.
    for (n in 1:N) {
      states_n <- which(z == n)
      if (length(states_n) == 0) {
        # 미배정: base H 에서 직접 추출
        pool[[n]] <- canonicalize_atom(sample_atom_from_base(params))
      } else {
        merged_data <- list()
        for (k in states_n) {
          sd_k <- state_data[[k]]
          for (s in sd_k) merged_data[[length(merged_data) + 1]] <- s
        }
        pool[[n]] <- .update_one_atom_from_data(
          merged_data, params,
          a_b, r_b, a_g, r_g,
          b_mode0, g_mode0, b_prec0, g_prec0)
      }
    }

    # ══════════════ (B2) 배정 z_k 재추출 (순차 in-place) ══════════════
    # π_n = V_n ∏_{l<n}(1-V_l)
    log_pi <- numeric(N)
    log_one_minus <- log(pmax(1 - Vv, 1e-300))
    cum <- 0
    for (n in 1:N) {
      log_pi[n] <- log(max(Vv[n], 1e-300)) + cum
      cum <- cum + log_one_minus[n]
    }
    if (K_c >= 1) {
      # 펼친 atoms 임시 구성 (plug-in 우도 계산용; 매 z 변경 후 갱신 불필요—
      #  L_k(n)은 atom_pool[[n]] 고정 plug-in 이므로 z 와 무관)
      for (k in 1:K_c) {
        z_prev <- if (k > 1) z[k - 1] else NA_integer_
        z_next <- if (k < K_c) z[k + 1] else NA_integer_   # 옛 값(아직 미갱신)
        # plug-in 우도 L_k(n): 상태 k 데이터를 고정 atom_pool[[n]] 로 본 우도
        # compute_atom_loglik 은 (cc,k,atom,...) 시그니처. 상위구간 t_in_k_upper 필요.
        t_in_k_upper <- which(state$S_upper[cc, ] == k)
        log_post <- numeric(N)
        for (n in 1:N) {
          if (!is.na(z_prev) && n == z_prev) { log_post[n] <- -Inf; next }
          if (!is.na(z_next) && n == z_next) { log_post[n] <- -Inf; next }
          Ln <- compute_atom_loglik(cc, k, pool[[n]], j_in_c, t_in_k_upper,
                                    state, params, Y, precomp, model)
          log_post[n] <- log_pi[n] + Ln
        }
        z[k] <- sample_categorical_log(log_post)
      }
    }

    # ══════════════ (B3) atom DP stick-breaking ══════════════
    m_n <- tabulate(z, nbins = N)
    Vnew <- numeric(N)
    tail_sum <- sum(m_n)   # Σ_{l>=1} m_l, 아래서 차감
    for (n in 1:N) {
      tail_after <- tail_sum - cumsum(m_n)[n]  # Σ_{l>n} m_l
      if (n < N) {
        Vnew[n] <- rbeta(1, 1 + m_n[n], vs_c + tail_after)
        Vnew[n] <- min(max(Vnew[n], 1e-10), 1 - 1e-10)
      } else Vnew[n] <- 1
    }
    Vv <- Vnew

    # ── ς^atom_c slice 갱신 (변화점 varsigma 와 분리) ──
    a_vs <- if (!is.null(hyper$a_atom_varsigma)) hyper$a_atom_varsigma else 2.0
    b_vs <- if (!is.null(hyper$b_atom_varsigma)) hyper$b_atom_varsigma else 1.0
    log_target_vs <- function(vv) {
      if (vv <= 0) return(-Inf)
      lp <- dgamma(vv, shape = a_vs, rate = b_vs, log = TRUE)
      ll <- 0
      for (n in 1:(N - 1)) ll <- ll + dbeta(Vv[n], 1, vv, log = TRUE)
      lp + ll
    }
    vs_c <- slice_sample_positive(vs_c, log_target_vs, w = 1, m = 5)

    # ── 펼치기: 외부 코드용 atoms[[cc]][[k]] ──
    atoms_view <- vector("list", K_c)
    for (k in 1:K_c) atoms_view[[k]] <- pool[[ z[k] ]]

    # 저장
    params$atom_pool[[cc]]     <- pool
    params$z_state[[cc]]       <- z
    params$atom_V[[cc]]        <- Vv
    params$atom_pi[[cc]]       <- exp(log_pi)
    params$atom_varsigma[cc]   <- vs_c
    params$atoms[[cc]]         <- atoms_view
  }
  params
}


###############################################################################
# (va) γ*(c,j): 시계열별 스칼라 intercept (구간 불변, segment-invariant)
#
# [γ*(c,j) 구조 변경] γ*(j,k) → γ*(c,j)
#   기존: 시계열 j, 구간 k마다 별도 intercept → 감소를 γ*가 흡수 → δ₁=-1 억제
#   변경: 시계열 j 전체에 하나의 intercept → 구간 간 변화는 반드시 f^shape가 설명
#   v*(j,k)→v_j 교체와 동일 논리: 구간별 자유도 제거 → 식별성 개선
#
# Prior:    γ*(c,j) ~ N(0, σ²_c)
# Likelihood: Y_{j,t} = α + β_c + γ*(c,j) + f^shape_{j,k(t)}(x_t) + ε_{j,t}
#
# Full conditional (conjugate Normal):
#   r_{j,t} = Y_{j,t} - α - β_c - f^shape_{j,k(t)}(x_t)
#   w_{j,t} = 1 / (v_j · φ_{j,t}^{ξ_{j,t}})
#   σ²_post = (Σ_t w_{j,t} + 1/σ²_c)^{-1}
#   μ_post  = σ²_post · Σ_t w_{j,t} · r_{j,t}
#   γ*(c,j) | rest ~ N(μ_post, σ²_post)
###############################################################################

update_state_intercepts_continuous <- function(state, params, Y, precomp, model) {
  for (j in 1:model$J) {
    c_j <- state$cluster[j]
    
    # σ²_c for this cluster
    sigma2_c <- params$sigma2_gamma_c[c_j]
    if (!is.finite(sigma2_c) || sigma2_c <= 0) sigma2_c <- 0.01
    
    # Residuals excluding γ* (γ*는 스칼라이므로 f_j 계산 시 현재 γ*가 포함되지 않음)
    f_j <- compute_f_all_timepoints(j, state, params, precomp, model)
    r_j <- Y[j, ] - params$alpha - params$beta[c_j] - f_j
    
    # Effective variance weights
    sigma2_j <- get_effective_variance_all(j, state, params)
    sigma2_j[!is.finite(sigma2_j) | sigma2_j <= 0] <- 1
    w_j <- 1 / sigma2_j
    
    # Conjugate Normal posterior
    prec_post <- sum(w_j) + 1 / sigma2_c
    sigma2_post <- 1 / prec_post
    mu_post <- sigma2_post * sum(w_j * r_j)
    
    if (!is.finite(mu_post)) mu_post <- 0
    if (!is.finite(sigma2_post) || sigma2_post <= 0) sigma2_post <- sigma2_c
    
    params$gamma[[j]] <- rnorm(1, mu_post, sqrt(sigma2_post))
    if (!is.finite(params$gamma[[j]])) params$gamma[[j]] <- mu_post
  }
  params
}

###############################################################################
# (vb) σ²_c: 클러스터별 γ* 수축 분산
#
# Prior:    σ²_c ~ InvGamma(c_a0, c_b0)
# Full conditional:
#   σ²_c | rest ~ InvGamma(c_a0 + n_c/2, c_b0 + Σ_{j∈c} γ*(c,j)² / 2)
#
# 여기서 n_c = 클러스터 c에 속한 시계열 수
###############################################################################

update_sigma2_gamma <- function(state, params, model) {
  hyper <- params$hyper
  c_a0 <- if (!is.null(hyper$c_a0)) hyper$c_a0 else 10
  c_b0 <- if (!is.null(hyper$c_b0)) hyper$c_b0 else 0.02
  
  for (cc in 1:model$C) {
    j_in_c <- which(state$cluster == cc)
    n_c <- length(j_in_c)
    
    SS_c <- 0
    for (j in j_in_c) {
      g_j <- params$gamma[[j]][1]
      if (is.finite(g_j)) SS_c <- SS_c + g_j^2
    }
    
    shape_post <- c_a0 + n_c / 2
    scale_post <- c_b0 + SS_c / 2
    
    params$sigma2_gamma_c[cc] <- rinvgamma(1, shape_post, scale_post)
    if (!is.finite(params$sigma2_gamma_c[cc]) || params$sigma2_gamma_c[cc] <= 0) {
      params$sigma2_gamma_c[cc] <- c_b0 / max(c_a0 - 1, 0.5)
    }
  }
  
  # 하위 호환: sigma2_gamma에 전체 평균 저장 (진단용)
  params$sigma2_gamma <- mean(params$sigma2_gamma_c)
  
  params
}

# Bivariate normal log-density helper
dmvnorm_log <- function(x, mu, Sigma) {
  d <- length(x)
  diff <- x - mu
  -0.5 * d * log(2 * pi) - 0.5 * log(det(Sigma)) - 
    0.5 * as.numeric(t(diff) %*% solve(Sigma) %*% diff)
}

###############################################################################
# (vi) Robust error block
###############################################################################

update_v_continuous <- function(state, params, Y, precomp, model, base_cache = NULL) {
  J <- model$J; T_len <- model$T_len
  hyper <- params$hyper
  
  for (j in 1:J) {
    if (is.null(base_cache)) {
      base_j <- compute_series_base_sigma_entry(j, state, params, Y, precomp, model)$base
    } else {
      base_j <- base_cache[[j]]
    }
    mu_j <- params$alpha + params$beta[state$cluster[j]] + base_j
    e_j <- Y[j, ] - mu_j
    
    # [v_j 스칼라화] 시리즈 전체의 잔차를 모아 하나의 InvGamma 업데이트
    phi_xi <- numeric(T_len)
    for (t in 1:T_len) {
      xi_val <- params$xi[j, t]; if (is.na(xi_val)) xi_val <- 0
      phi_val <- params$phi[j, t]; if (!is.finite(phi_val) || phi_val < 1) phi_val <- 1
      phi_xi[t] <- phi_val^xi_val
    }
    phi_xi[phi_xi <= 0 | !is.finite(phi_xi)] <- 1
    ss <- sum(e_j^2 / phi_xi)
    if (!is.finite(ss)) ss <- 0
    
    shape <- hyper$v01 + T_len / 2
    scale <- hyper$v02 + ss / 2
    params$v[[j]] <- rinvgamma(1, shape, scale)
    if (!is.finite(params$v[[j]]) || params$v[[j]] <= 0) params$v[[j]] <- 1
  }
  params
}

update_xi_continuous <- function(state, params, Y, precomp, model, base_cache = NULL) {
  J <- model$J; T_len <- model$T_len
  
  for (j in 1:J) {
    if (is.null(base_cache)) {
      base_j <- compute_series_base_sigma_entry(j, state, params, Y, precomp, model)$base
    } else {
      base_j <- base_cache[[j]]
    }
    mu_j <- params$alpha + params$beta[state$cluster[j]] + base_j
    e_j <- Y[j, ] - mu_j
    
    for (t in 1:T_len) {
      # [v_j 스칼라화] regime 인덱싱 제거
      v_j <- params$v[[j]]
      if (!is.finite(v_j) || v_j <= 0) v_j <- 1
      phi_jt <- params$phi[j, t]
      if (!is.finite(phi_jt) || phi_jt < 1) phi_jt <- 1
      
      log_p1 <- log(max(params$theta_out, 1e-10)) + dnorm_log_var(e_j[t], 0, v_j * phi_jt)
      log_p0 <- log(max(1 - params$theta_out, 1e-10)) + dnorm_log_var(e_j[t], 0, v_j)
      
      if (!is.finite(log_p1)) log_p1 <- -1e300
      if (!is.finite(log_p0)) log_p0 <- -1e300
      
      prob1 <- 1 / (1 + exp(log_p0 - log_p1))
      if (!is.finite(prob1)) prob1 <- 0.5
      params$xi[j, t] <- as.integer(runif(1) < prob1)
    }
  }
  params
}

update_theta_out <- function(state, params, model) {
  hyper <- params$hyper
  J <- model$J; T_len <- model$T_len
  
  sum_xi <- sum(params$xi, na.rm = TRUE)
  total <- J * T_len
  
  a_post <- max(hyper$k0 * hyper$m0 + sum_xi, 0.01)
  b_post <- max(hyper$k0 * (1 - hyper$m0) + total - sum_xi, 0.01)
  
  params$theta_out <- rbeta(1, a_post, b_post)
  if (!is.finite(params$theta_out)) params$theta_out <- hyper$m0
  params
}

update_phi_continuous <- function(state, params, Y, precomp, model, base_cache = NULL) {
  J <- model$J; T_len <- model$T_len
  
  for (j in 1:J) {
    if (is.null(base_cache)) {
      base_j <- compute_series_base_sigma_entry(j, state, params, Y, precomp, model)$base
    } else {
      base_j <- base_cache[[j]]
    }
    mu_j <- params$alpha + params$beta[state$cluster[j]] + base_j
    e_j <- Y[j, ] - mu_j
    
    for (t in 1:T_len) {
      # [v_j 스칼라화] regime 인덱싱 제거
      v_j <- params$v[[j]]
      if (!is.finite(v_j) || v_j <= 0) v_j <- 1
      nu <- params$nu
      xi_jt <- params$xi[j, t]
      if (is.na(xi_jt)) xi_jt <- 0L
      
      if (xi_jt == 0) {
        params$phi[j, t] <- 1.0
      } else {
        e2v <- e_j[t]^2 / v_j
        if (!is.finite(e2v)) e2v <- 0
        shape <- (nu + 1) / 2
        scale <- (nu + e2v) / 2
        # Truncated InvGamma(shape, scale) restricted to [1, Inf)
        # 기각 샘플링으로 phi >= 1 보장
        max_tries <- 200
        accepted <- FALSE
        for (try_i in 1:max_tries) {
          proposal <- rinvgamma(1, shape, scale)
          if (is.finite(proposal) && proposal >= 1.0) {
            params$phi[j, t] <- proposal
            accepted <- TRUE
            break
          }
        }
        if (!accepted) {
          # 극단적으로 드문 경우: scale/shape >> 1이면 mode < 1일 수 있음
          # 이 경우 mode를 1로 설정 (안전한 fallback)
          params$phi[j, t] <- 1.0
        }
      }
      
      # 안전장치: NaN/Inf 방어
      if (!is.finite(params$phi[j, t])) {
        params$phi[j, t] <- 1.0
      }
    }
  }
  params
}

update_nu <- function(state, params, model) {
  J <- model$J; T_len <- model$T_len

  # Keep current value inside the declared support.
  if (is.null(params$nu) || !is.finite(params$nu)) params$nu <- 5
  params$nu <- min(max(as.numeric(params$nu)[1], 1.001), 40)

  log_target <- function(nu) {
    if (length(nu) != 1 || is.na(nu) || !is.finite(nu)) return(-Inf)
    if (nu < 1 || nu > 40) return(-Inf)

    log_val <- 0  # uniform prior on [1, 40]
    for (j in 1:J) {
      for (t in 1:T_len) {
        # phi는 xi=1(이상치)일 때만 InvGamma에서 추출됨.
        # xi=0(정상)이면 phi=1로 결정적 고정이므로 nu의 우도에 포함하면 안 됨.
        if (params$xi[j, t] == 0L) next
        phi_jt <- params$phi[j, t]
        if (is.na(phi_jt) || !is.finite(phi_jt) || phi_jt <= 0) return(-Inf)
        dens <- dinvgamma_log(phi_jt, nu / 2, nu / 2)
        if (is.na(dens) || !is.finite(dens)) return(-Inf)
        log_val <- log_val + dens
      }
    }
    log_val
  }

  # Support must match log_target(). The previous lower=8 caused invalid
  # brackets when the chain started below 8 (e.g. nu=5 at initialization).
  params$nu <- slice_sample(params$nu, log_target, w = 2, m = 5,
                             lower = 1, upper = 40)
  #params$nu <- 20

  params
}

###############################################################################
# DP stick-breaking weights
###############################################################################

update_stick_breaking <- function(state, params, model) {
  C <- model$C
  
  for (cc in 1:C) {
    K_c <- state$K[cc]
    T_len <- model$T_len
    
    # Count transitions for each state from the upper-level state sequence
    n_stay <- numeric(K_c)
    n_trans <- numeric(K_c)
    
    S_upper_c <- state$S_upper[cc, ]
    for (t in 1:(T_len - 1)) {
      k <- S_upper_c[t]
      if (k >= 1 && k <= K_c) {
        if (S_upper_c[t + 1] == k) {
          n_stay[k] <- n_stay[k] + 1
        } else if (S_upper_c[t + 1] == k + 1) {
          n_trans[k] <- n_trans[k] + 1
        }
      }
    }
    
    # Sample v_k from posterior
    # π*_k = v_k (conditional persistence = stay probability)
    # Likelihood: (π*_k)^{n_stay} * (1-π*_k)^{n_trans}
    # Prior: v_k ~ Beta(1, varsigma)
    # Posterior: v_k ~ Beta(1 + n_stay_k, varsigma + n_trans_k)
    varsigma_c <- params$varsigma[cc]
    v_k <- numeric(K_c)
    for (k in 1:K_c) {
      if (k < K_c) {
        v_k[k] <- rbeta(1, 1 + n_stay[k], varsigma_c + n_trans[k])
        v_k[k] <- min(max(v_k[k], 1e-10), 1 - 1e-10)
      } else {
        v_k[k] <- 1  # last stick gets everything
      }
    }
    
    # Construct pi_weights from stick-breaking: pi_k = v_k * prod_{l<k}(1-v_l)
    pi_weights <- numeric(K_c)
    pi_weights[1] <- v_k[1]
    if (K_c >= 2) {
      cum_one_minus_v <- cumprod(1 - v_k)
      for (k in 2:K_c) {
        pi_weights[k] <- v_k[k] * cum_one_minus_v[k - 1]
      }
    }
    
    # Normalize for safety
    pi_weights <- pmax(pi_weights, 1e-10)
    pi_weights <- pi_weights / sum(pi_weights)
    
    params$pi_weights[[cc]] <- pi_weights
    params$pi_star[[cc]] <- compute_pi_star(pi_weights)
    
    # Update concentration parameter varsigma_c
    # Target: p(varsigma) * prod_{k<K} Beta(v_k; 1, varsigma)
    log_target_varsigma <- function(vs) {
      if (vs <= 0) return(-Inf)
      hyper <- params$hyper
      log_prior <- dgamma(vs, shape = hyper$a_varsigma, rate = hyper$b_varsigma, log = TRUE)
      # Likelihood of v_k under Beta(1, varsigma)
      log_lik <- 0
      for (k in 1:(K_c - 1)) {
        log_lik <- log_lik + dbeta(v_k[k], 1, vs, log = TRUE)
      }
      log_prior + log_lik
    }
    
    params$varsigma[cc] <- slice_sample_positive(
      params$varsigma[cc], log_target_varsigma, w = 1, m = 5)
  }
  
  # ═══════════════════════════════════════════════════════════════
  # [HDP-HMM] Lower-level stick-breaking: π*_lower
  #
  # v^L_{c,k} ~ Beta(a_pi_lower + Σ_j n^L_{j,stay,k},
  #                   b_pi_lower + Σ_j n^L_{j,trans,k})
  #
  # a_pi_lower >> 1 → π*_lower 를 1 근방으로 유도
  # → 하위 CP geometric prior가 약정보 (경계 끌림 억제)
  # ═══════════════════════════════════════════════════════════════
  hyper <- params$hyper
  a_pi_L <- if (!is.null(hyper$a_pi_lower)) hyper$a_pi_lower else 200
  b_pi_L <- if (!is.null(hyper$b_pi_lower)) hyper$b_pi_lower else 1
  
  if (is.null(params$pi_star_lower)) params$pi_star_lower <- vector("list", C)
  
  for (cc in 1:C) {
    K_c <- state$K[cc]
    T_len <- model$T_len
    j_in_c <- which(state$cluster == cc)
    
    # Aggregate lower-level stay/trans counts across series in cluster
    n_stay_L <- numeric(K_c)
    n_trans_L <- numeric(K_c)
    for (j in j_in_c) {
      S_L <- state$S_lower[j, ]
      for (t in 1:(T_len - 1)) {
        k <- S_L[t]
        if (k >= 1 && k <= K_c) {
          if (S_L[t + 1] == k) {
            n_stay_L[k] <- n_stay_L[k] + 1
          } else if (S_L[t + 1] == k + 1) {
            n_trans_L[k] <- n_trans_L[k] + 1
          }
        }
      }
    }
    
    # Sample v^L_k from posterior
    v_k_L <- numeric(K_c)
    for (k in 1:K_c) {
      if (k < K_c) {
        v_k_L[k] <- rbeta(1, a_pi_L + n_stay_L[k], b_pi_L + n_trans_L[k])
        v_k_L[k] <- min(max(v_k_L[k], 1e-10), 1 - 1e-10)
      } else {
        v_k_L[k] <- 1
      }
    }
    
    # Stick-breaking → pi_star_lower
    pi_w_L <- numeric(K_c)
    pi_w_L[1] <- v_k_L[1]
    if (K_c >= 2) {
      cum_1mv <- cumprod(1 - v_k_L)
      for (k in 2:K_c) pi_w_L[k] <- v_k_L[k] * cum_1mv[k - 1]
    }
    pi_w_L <- pmax(pi_w_L, 1e-10)
    pi_w_L <- pi_w_L / sum(pi_w_L)
    params$pi_star_lower[[cc]] <- compute_pi_star(pi_w_L)
  }
  
  params
}


###############################################################################
# γ* Post-sweep Centering
#
# α + β_c + γ*_j 간의 ridge (비식별성)를 제거하여 MCMC 혼합 개선.
# μ_{j,t} = α + β_c + γ*_j + f_j(t) 는 불변.
#
# 알고리즘:
#   1. 각 클러스터 c에서 γ̄_c = mean(γ*_j : j ∈ c) 계산
#   2. γ*_j ← γ*_j - γ̄_c  (클러스터 내 센터링)
#   3. β_c ← β_c + γ̄_c    (평균 흡수)
#   4. β₁ = 0 제약 복원:
#      Δ = β₁ (현재값), α ← α + Δ, β_c ← β_c - Δ (∀c)
#
# Ref: Papaspiliopoulos et al. (2007, Stat. Sci.)
###############################################################################
center_gamma_sweep <- function(state, params, model) {
  C <- model$C
  
  # Step 1-3: 클러스터별 γ* 센터링 + β 흡수
  for (cc in 1:C) {
    j_in_c <- which(state$cluster == cc)
    if (length(j_in_c) == 0) next
    
    gamma_vals <- sapply(j_in_c, function(j) {
      g <- params$gamma[[j]]
      if (is.null(g) || length(g) == 0 || !is.finite(g[1])) 0 else g[1]
    })
    
    gamma_bar <- mean(gamma_vals)
    
    # γ*_j ← γ*_j - γ̄_c
    for (j in j_in_c) {
      params$gamma[[j]] <- params$gamma[[j]] - gamma_bar
    }
    
    # β_c ← β_c + γ̄_c
    params$beta[cc] <- params$beta[cc] + gamma_bar
  }
  
  # Step 4: β₁ = 0 제약 복원
  delta <- params$beta[1]
  params$alpha <- params$alpha + delta
  params$beta <- params$beta - delta
  
  params
}
