###############################################################################
# 10_initialization.R (REVISED: Anti-Overfitting Version)
###############################################################################

initialize_model <- function(Y, model) {
  J <- model$J; T_len <- model$T_len; C <- model$C
  K_init <- model$K_init; M <- model$M; m_min <- model$m_min
  
  # === Data standardization ===
  Y_means <- rowMeans(Y)
  Y_sds <- apply(Y, 1, sd)
  Y_sds[Y_sds < 1e-10] <- 1
  Y_std <- (Y - Y_means) / Y_sds
  std_info <- list(means = Y_means, sds = Y_sds, Y_original = Y)
  
  # === Precompute K-L basis ===
  precomp <- precompute_basis(M, n_grid = 500)
  precomp <- precompute_global_x(T_len, M, precomp)
  precomp$std_info <- std_info
  
  # === Hyperparameters ===
  hyper <- list(
    m_alpha = 0, sigma2_alpha = 25,
    a0_beta = 2, b0_beta = 0.5,
    alpha0_dirichlet = J / C,  # 균등 클러스터 사전분포 (상세 균형 보존)
    a_tau0 = 2, b_tau0 = 0.5,
    a_eta0 = 2, b_eta0 = 0.5,
    
    # Decoupled Atom Hyperpriors
    
    a_shape_beta = 2.0,  b_shape_beta = 1.5,   # rate 는 초기값, Gibbs 업데이트됨
    a_shape_gamma = 1.5 , b_shape_gamma = 2.0,
    a0_rbeta = 3.0,  b0_rbeta = 2.0,         # [신규] r_β ~ Gamma(2,1)
    a0_rgamma = 4.0, b0_rgamma = 2.0,        # [신규] r_γ ~ Gamma(2,1)

    a_gamma = 20, b_gamma = 0.04,
    
    # [γ*(c,j)] 클러스터별 γ* 수축 분산의 InvGamma 사전분포
    # σ²_c ~ IG(c_a0, c_b0), mode = c_b0/(c_a0+1)
    # [v2 수정] c_a0=10, c_b0=0.02 (mode≈0.002)는 γ*를 사실상 0으로 고정하여
    #   α가 δ₁ offset을 흡수하는 피드백 루프를 유발.
    #   c_a0=5, c_b0=0.5 (mode≈0.083)으로 완화하여 γ*가 부분 흡수 가능하게 함.
    #   Collapsed centering이 atom update 내에서 δ₁ 식별을 보장하므로,
    #   σ²_c 완화가 Mode B를 유발하지 않음.
    c_a0 = 5, c_b0 = 0.5,
    
    v01 = 2, v02 = 1.0,
    k0 = 50, m0 = 0.03,
    
    # DP Concentration (Upper): 0.1로 낮추어 과분할 억제
    a_varsigma = 2, b_varsigma = 20, 
    
    # [HDP-HMM] Lower-level stick-breaking hyperparameters
    # v^L_k ~ Beta(a_pi_lower + n_stay, b_pi_lower + n_trans)
    # a_pi_lower >> 1 → π*_lower ≈ 1 → 약정보 geometric prior (경계 끌림 억제)
    # (0.99)^27 ≈ 0.76 vs (0.97)^27 ≈ 0.44
    a_pi_lower = 50, b_pi_lower = 1,
    # ς_lower 농도 파라미터의 사전분포
    a_varsigma_lower = 1, b_varsigma_lower = 1000,
    
    as_count = 1, bs_count = 9, r_od = 5,
    
    # [PATCH] Poisson prior on K_c: K_c ~ TruncPoisson(lambda0_K, {1,...,K_max})
    lambda0_K = if (!is.null(model$lambda0_K)) model$lambda0_K else 1.2
  )
  
  # === State & Params ===
  fixed_cs <- isTRUE(model$fixed_cluster_sizes)
  state <- init_state(J, T_len, C, K_init, m_min, fixed_cluster_sizes = fixed_cs)
  state$cluster <- integer(J)
for (cc in seq_len(C)) {
  state$cluster[model$group_members[[cc]]] <- cc
}
  params <- list(hyper = hyper)
  params$alpha <- 0
  params$beta <- rep(0, C)
  params$beta[1] <- 0  # [식별 제약] β₁ = 0: cluster 1이 기준 수준
  params$lambda2 <- rep(hyper$b0_beta/max(hyper$a0_beta-1,0.5),C)  # beta를 위한 수축 파라미터 초기화
  params$p <- rep(1/C, C)
  params$sigma2_gamma <- hyper$c_b0 / max(hyper$c_a0 - 1, 0.5)  # backward compat

  # [θ_{j,0}=0 고정] m=0 절편은 α+β+γ*가 담당하므로 GP에서 제거
  # m=1,2만 데이터 투영으로 초기화, m=0은 0 고정, m≥3은 0
  params$theta <- matrix(0, J, M + 1)
  for (j in 1:J) {
    params$theta[j, 1] <- 0  # m=0: 항상 0
    for (mm in 1:M) {
      if (mm <= 2) {
        basis_vals <- sqrt(2) * cos(mm * pi * precomp$x_global)
        params$theta[j, mm + 1] <- (sum(Y_std[j, ] * basis_vals) / T_len) * 0.5
      } else {
        params$theta[j, mm + 1] <- 0.0
      }
    }
  }
  
  # [θ_{j,0}=0 고정] theta0도 m=0은 0 고정
  params$theta0 <- matrix(0, C, M + 1)
  for (cc in 1:C) {
    j_in_c <- which(state$cluster == cc)
    if (length(j_in_c) > 0) params$theta0[cc, ] <- colMeans(params$theta[j_in_c, , drop = FALSE])
    params$theta0[cc, 1] <- 0  # m=0 고정
  }
  
  params$tau2 <- matrix(hyper$b_tau0/max(hyper$a_tau0-1,0.5), C, M + 1)
  params$eta2 <- rep(hyper$b_eta0/max(hyper$a_eta0-1,0.5), C)
  
  params$r <- rep(1.0, J)
  params$kr <- 1.5; params$xi_r <- 2
  
  # 원자 초기화 — spline+conditional 데이터 기반 초기화
  # 각 구간의 데이터에 2차 다항식을 적합하여 δ₁,δ₂를 결정하고
  # WLS로 b,g를 추정. H에서 blind 추출 대신 데이터 적합 atom으로 시작.
  params$atoms <- vector("list", C)
  for (cc in 1:C) {
    params$atoms[[cc]] <- vector("list", K_init)
    j_in_c <- which(state$cluster == cc)
    for (k in 1:K_init) {
      params$atoms[[cc]][[k]] <- init_atom_from_data(
        Y_std, state, k, cc, j_in_c, params, precomp, hyper)
    }
  }
  
  # [γ*(c,j)] 시계열별 스칼라 intercept (구간 불변)
  # 기존 γ*(j,k) → γ*(c,j): 구간별 자유도 제거 → 식별성 개선
  params$gamma <- vector("list", J)
  for (j in 1:J) {
    params$gamma[[j]] <- 0  # scalar, NOT K-length vector
  }
  
  # σ²_c: 클러스터별 γ* 수축 분산
  params$sigma2_gamma_c <- rep(hyper$c_b0 / max(hyper$c_a0 - 1, 0.5), C)
  
  # [수정] DP 집중도 파라미터 0.1로 하향 (과분할 방어벽)
  params$varsigma <- rep(hyper$a_varsigma/hyper$b_varsigma,C)
  params$pi_weights <- vector("list", C)
  params$pi_star <- vector("list", C)
  for (cc in 1:C) {
    params$pi_weights[[cc]] <- sample_stick_breaking(K_init, params$varsigma[cc])
    params$pi_star[[cc]] <- compute_pi_star(params$pi_weights[[cc]])
  }
  
  # [HDP-HMM] Lower-level stick-breaking: 별도의 π*_lower (경계 끌림 억제)
  params$varsigma_lower <- rep(hyper$b_varsigma_lower / max(hyper$a_varsigma_lower, 1), C)
  params$pi_star_lower <- vector("list", C)
  for (cc in 1:C) {
    # 초기값: a_pi_lower / (a_pi_lower + b_pi_lower) ≈ 0.995
    pi_lower_init <- hyper$a_pi_lower / (hyper$a_pi_lower + hyper$b_pi_lower)
    params$pi_star_lower[[cc]] <- rep(pi_lower_init, K_init)
  }
  
  # Robust error & Count model (기존 로직 유지)
  # [v_j 스칼라화] 시리즈당 하나의 기저 분산 (regime-invariant)
  params$v <- vector("list", J); for (j in 1:J) params$v[[j]] <- 1.0
  params$xi <- matrix(0L, J, T_len); params$theta_out <- 0.005
  params$phi <- matrix(1, J, T_len); params$nu <- 8
  
  if (model$type == "count") {
    params$z_count <- matrix(0L, J, T_len); params$eta <- matrix(1, J, T_len)
    params$s_count <- 0.1; params$r_od <- hyper$r_od
  }
  
  list(state = state, params = params, precomp = precomp, Y_std = Y_std)
}

###############################################################################
# DATA-INFORMED ATOM INITIALIZATION (Spline+Conditional)
#
# 07_interval_adjusting.R의 propose_atom_spline_conditional과 동일 원리:
#   Stage 1: 2차 다항식 적합 → δ₁, δ₂ 결정
#   Stage 2: δ₁,δ₂ 조건부 WLS → b, g 추정
#
# 초기화 시 θ_j가 이미 설정되어 있으므로 H_j(x) 계산 가능.
# 이를 통해 K_init개 구간 각각이 데이터에 맞는 atom으로 시작.
# → merge가 "진짜 불필요한" 구간만 제거, "필요한" 구간은 보존.
###############################################################################

init_atom_from_data <- function(Y, state, k, cc, j_in_c, params, precomp, hyper) {
  a_sb <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 5.0
  b_sb <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 0.5
  a_sg <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 0.5
  b_sg <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 1.0

  # Fallback
  fallback <- function(d1 = NULL, d2 = NULL) {
    sb <- max(rgamma(1, shape = a_sb, rate = b_sb), 0.01)
    sg <- max(rgamma(1, shape = a_sg, rate = b_sg), 0)
    if (is.null(d1)) d1 <- sample(c(-1L, 1L), 1)
    if (is.null(d2)) d2 <- sample(c(-1L, 1L), 1)
    canonicalize_atom(list(gamma1 = d1, gamma2 = d2,
                           shape_beta = sb, shape_gamma = sg))
  }

  t_in_k <- which(state$S_upper[cc, ] == k)
  if (length(t_in_k) < 3 || length(j_in_c) == 0) return(fallback())

  x_seg <- precomp$x_global[t_in_k]
  x_c   <- x_seg - mean(x_seg)
  alpha  <- params$alpha
  beta_c <- params$beta[cc]
  Mp1    <- ncol(params$theta)

  # ==================================================================
  # Stage 1: 2차 다항식 적합 → δ₁, δ₂
  # ==================================================================
  a1_agg <- 0; a2_agg <- 0; w_total <- 0
  X_quad <- cbind(1, x_c, x_c^2)

  for (j in j_in_c) {
    r_j <- Y[j, t_in_k] - alpha - beta_c
    r_j <- r_j - mean(r_j)

    v_j <- var(Y[j, t_in_k])
    if (!is.finite(v_j) || v_j <= 0) v_j <- 1
    w_j <- length(t_in_k) / v_j

    fit <- tryCatch(.lm.fit(X_quad, r_j), error = function(e) NULL)
    if (is.null(fit) || length(fit$coefficients) < 3L) next

    a1_agg <- a1_agg + w_j * fit$coefficients[2]
    a2_agg <- a2_agg + w_j * fit$coefficients[3]
    w_total <- w_total + w_j
  }

  if (w_total < 1e-12) return(fallback())

  delta1 <- if (a1_agg >= 0) +1L else -1L
  delta2 <- if (a2_agg >= 0) +1L else -1L

  # ==================================================================
  # Stage 2: WLS → b, g (δ₁, δ₂ 조건부)
  # ==================================================================
  n_seg <- length(t_in_k)
  XtWX <- matrix(0, 2, 2)
  XtWr <- numeric(2)

  for (j in j_in_c) {
    theta_j <- params$theta[j, ]
    theta_shape <- theta_j[2:Mp1]
    SZ_safe <- max(sum(theta_shape^2), 1e-12)

    # H_j(x_t)
    theta_outer <- tcrossprod(theta_shape)
    D_sub <- precomp$D_at_x[2:Mp1, 2:Mp1, , drop = FALSE]
    H_vals <- numeric(n_seg)
    for (i in seq_along(t_in_k)) {
      H_vals[i] <- sum(theta_outer * D_sub[, , t_in_k[i]]) / SZ_safe
    }

    A_j <- delta1 * x_seg
    B_j <- (delta1 + delta2) / 2 * x_seg - delta2 * H_vals
    r_j <- Y[j, t_in_k] - alpha - beta_c

    # Per-series centering → absorb γ*
    A_c <- A_j - mean(A_j); B_c <- B_j - mean(B_j); r_c <- r_j - mean(r_j)

    # Weight: use empirical variance at init (no v_j params yet)
    v_j <- var(Y[j, t_in_k])
    if (!is.finite(v_j) || v_j <= 0) v_j <- 1
    w_t <- rep(1 / v_j, n_seg)

    XtWX[1,1] <- XtWX[1,1] + sum(w_t * A_c^2)
    XtWX[1,2] <- XtWX[1,2] + sum(w_t * A_c * B_c)
    XtWX[2,1] <- XtWX[2,1] + sum(w_t * A_c * B_c)
    XtWX[2,2] <- XtWX[2,2] + sum(w_t * B_c^2)
    XtWr[1]   <- XtWr[1]   + sum(w_t * A_c * r_c)
    XtWr[2]   <- XtWr[2]   + sum(w_t * B_c * r_c)
  }

  # Gamma prior regularization
  g_mode <- max((a_sg - 1) / b_sg, 0.01)
  b_mode <- max((a_sb - 1) / b_sb, 0.01)
  prec_g <- b_sg^2 / max(a_sg - 1, 0.5)
  prec_b <- b_sb^2 / max(a_sb - 1, 0.5)

  XtWX[1,1] <- XtWX[1,1] + prec_g
  XtWX[2,2] <- XtWX[2,2] + prec_b
  XtWr[1]   <- XtWr[1]   + prec_g * g_mode
  XtWr[2]   <- XtWr[2]   + prec_b * b_mode

  # Solve 2×2
  det_val <- XtWX[1,1] * XtWX[2,2] - XtWX[1,2]^2
  if (abs(det_val) < 1e-20 || !is.finite(det_val)) return(fallback(delta1, delta2))

  Sigma <- matrix(c(XtWX[2,2], -XtWX[1,2], -XtWX[2,1], XtWX[1,1]), 2, 2) / det_val
  mode_gb <- as.numeric(Sigma %*% XtWr)
  if (!all(is.finite(mode_gb))) return(fallback(delta1, delta2))

  # 초기화에서는 mode 사용 (확률적 샘플링 대신 안정적 시작점)
  sg <- max(mode_gb[1], 0)
  sb <- max(mode_gb[2], 0.01)

  canonicalize_atom(list(gamma1 = delta1, gamma2 = delta2,
                         shape_beta = sb, shape_gamma = sg))
}

# [PATCH] fixed_clusters: 클러스터 할당 완전 고정 (pair swap 포함 제거)
# [PATCH] lambda0_K: K_c ~ TruncPoisson(lambda0_K, {1,...,K_max}) 사전분포의 모수
# [PATCH] group_members: 분석가 지정 그룹 (고정, pair swap에 불변)
#   Group = 분석가가 사전 지정한 시계열 집합 (MCMC 전체에서 불변)
#   Cluster = MCMC iteration마다 pair swap으로 변할 수 있는 라벨
#   Posterior averaging은 Group 중심으로 수행됨
create_model_spec <- function(Y,
                              C = 2,
                              K_init = 3,
                              K_max = 5,
                              K_min = 2,
                              M = 8,
                              m_min = 12,
                              type = "continuous",
                              fixed_cluster_sizes = TRUE,
                              fixed_clusters = TRUE,
                              lambda0_K = 1.0,
                              N_atom = NULL,
                              group_members=list(c(1,2,3),c(4,5,6),c(7,8,9)),           # <- NO DEFAULT (required)
                              virtual_augment = TRUE) {
  
  J <- nrow(Y)
  
  
  
  # Coerce each element to integer and validate
  group_members <- lapply(seq_along(group_members), function(g) {
    gm <- group_members[[g]]
    if (length(gm) == 0L) {
      stop(sprintf("group_members[[%d]] is empty. Each cluster must contain at least one series.", g),
           call. = FALSE)
    }
    if (!is.numeric(gm) || any(gm != as.integer(gm))) {
      stop(sprintf("group_members[[%d]] must be integer-valued.", g), call. = FALSE)
    }
    as.integer(gm)
  })
  
  # All series indices must be in {1, ..., J}, no duplicates, all J series covered
  all_idx <- unlist(group_members)
  
  if (any(all_idx < 1L) || any(all_idx > J)) {
    stop(sprintf(
      "group_members contains indices outside 1:%d. Found range [%d, %d].",
      J, min(all_idx), max(all_idx)
    ), call. = FALSE)
  }
  
  if (anyDuplicated(all_idx)) {
    dup <- unique(all_idx[duplicated(all_idx)])
    stop(sprintf(
      "group_members has duplicated series indices: %s. Each series must belong to exactly one cluster.",
      paste(dup, collapse = ", ")
    ), call. = FALSE)
  }
  
  if (length(all_idx) != J) {
    missing_idx <- setdiff(seq_len(J), all_idx)
    stop(sprintf(
      "group_members covers %d series but J = %d. Missing indices: %s",
      length(all_idx), J, paste(missing_idx, collapse = ", ")
    ), call. = FALSE)
  }
  


init_atom_dp_fields <- function(params, state, model) {
  C <- model$C
  N <- if (!is.null(model$N_atom)) model$N_atom else 20L

  params$atom_pool     <- vector("list", C)  # 길이 N 풀
  params$z_state       <- vector("list", C)  # 길이 K_c 배정
  params$atom_V        <- vector("list", C)  # 길이 N stick 막대
  params$atom_pi       <- vector("list", C)  # 길이 N stick 가중
  params$atom_varsigma <- rep(1.0, C)        # atom DP 농도 (slice 갱신됨)

  for (cc in 1:C) {
    K_c <- state$K[cc]
    # 풀: base H 에서 N개 추출
    pool <- vector("list", N)
    for (n in 1:N) pool[[n]] <- canonicalize_atom(sample_atom_from_base(params))
    params$atom_pool[[cc]] <- pool

    # 배정: 인접 제약 만족하는 1,2,1,2,... (앞쪽 풀에 매핑)
    z <- integer(K_c)
    if (K_c >= 1) z[1] <- 1L
    if (K_c >= 2) for (k in 2:K_c) z[k] <- if (z[k-1]==1L) 2L else 1L
    params$z_state[[cc]] <- z

    # stick 막대: Beta(1, ς) prior 추출
    Vv <- c(rbeta(N-1, 1, params$atom_varsigma[cc]), 1)
    Vv <- pmin(pmax(Vv, 1e-10), 1-1e-10); Vv[N] <- 1
    params$atom_V[[cc]] <- Vv

    # 가중
    log_pi <- numeric(N); cum <- 0
    lom <- log(pmax(1-Vv,1e-300))
    for (n in 1:N) { log_pi[n] <- log(max(Vv[n],1e-300))+cum; cum <- cum+lom[n] }
    params$atom_pi[[cc]] <- exp(log_pi)

    # 초기 펼친 뷰
    av <- vector("list", K_c)
    for (k in 1:K_c) av[[k]] <- pool[[ z[k] ]]
    params$atoms[[cc]] <- av
  }
  params
}

  # ---------------------------------------------------------------------------
  # Virtual augmentation offset (Mastrantonio et al. 2022)
  # ---------------------------------------------------------------------------
  virtual_offset <- if (isTRUE(virtual_augment)) as.integer(m_min) else 0L

  # ---------------------------------------------------------------------------
  # Atom-pool truncation level N (single source of truth).
  #   Ishwaran-James truncation must satisfy N >= floor(T / m_min); we also keep
  #   a floor of 20.  ALL code paths (init pool, update_atoms, split, merge,
  #   .propose_atom_index) must read model$N_atom so the pool index space is
  #   identical everywhere (otherwise split can propose slots PU never created).
  N_atom_val <- if (!is.null(N_atom)) as.integer(N_atom) else
                max(20L, as.integer(floor(ncol(Y) / max(m_min, 1L))) + 5L)

  # ---------------------------------------------------------------------------
  # Return spec
  # ---------------------------------------------------------------------------
  list(
    type                = type,
    J                   = J,
    T_len               = ncol(Y),
    C                   = C,
    K_init              = K_init,
    K_max               = K_max,
    K_min               = K_min,
    M                   = M,
    m_min               = m_min,
    N_atom              = N_atom_val,
    fixed_cluster_sizes = fixed_cluster_sizes,
    fixed_clusters      = fixed_clusters,
    lambda0_K           = lambda0_K,
    group_members       = group_members,
    virtual_offset      = virtual_offset
  )
}