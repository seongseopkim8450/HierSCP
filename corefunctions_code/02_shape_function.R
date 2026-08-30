###############################################################################
# 02_shape_function.R
# Shape-restricted smoothing function f_j(x_t) — GLOBAL time normalization
#
# x_t = (t-1)/(T-1) ∈ [0, 1] is a FIXED global coordinate.
# When changepoints move, only the atom (δ1, δ2, shape_beta, shape_gamma) governing each
# time point changes — the x_t values themselves never change.
#
# IMPORTANT REFACTOR:
#   The state intercept gamma^*_{j,k} is no longer part of the DP atom.
#   It now lives in params$gamma[[j]][k] with a Gaussian random-walk
#   shrinkage prior across states.
###############################################################################

#' Precompute global time grid and all x-dependent quantities
precompute_global_x <- function(T_len, M, precomp) {
  x_global <- (0:(T_len - 1)) / max(T_len - 1, 1)
  
  phi_at_x <- matrix(0, nrow = M + 1, ncol = T_len)
  phi_at_x[1, ] <- 1
  for (mm in 1:M) phi_at_x[mm + 1, ] <- sqrt(2) * cos(mm * pi * x_global)
  
  D_at_x <- array(0, dim = c(M + 1, M + 1, T_len))
  for (mi in 0:M) for (mj in mi:M) {
    D_vals <- compute_D_analytical(mi, mj, x_global)
    D_at_x[mi + 1, mj + 1, ] <- D_vals
    if (mi != mj) D_at_x[mj + 1, mi + 1, ] <- D_vals
  }
  
  precomp$x_global <- x_global
  precomp$phi_at_x <- phi_at_x
  precomp$D_at_x <- D_at_x
  precomp$T_len <- T_len
  precomp
}

#' Evaluate SHAPE-ONLY function at specified time indices using global x_t
#'
#' Decoupled shape function (slack variable reparameterization):
#'   f_j^shape(x_t) = L(δ1,δ2) * (shape_beta * SZ + shape_gamma) * x_t
#'                     - δ2 * shape_beta * DI(x_t)
#' where
#'   L(δ1,δ2) = δ1 { |δ2+1|/2 + |δ2-1||δ1+1|/4 + |δ1δ2+1|/2 }
#'   SZ       = Σ θ²   over m=1,...,M
#'   DI(x)    = θᵀ D(x) θ
#'
#' Guaranteed shape classes for all (shape_beta > 0, shape_gamma ≥ 0):
#'   (+,+) increasing convex, (+,-) increasing concave,
#'   (-,+) decreasing convex, (-,-) decreasing concave.
#'
#' The state intercept gamma^*_{j,k} is handled separately in compute_mu_all().
shape_delta1 <- function(atom) normalize_sign_carrier(atom$gamma1, default = 1)
shape_delta2 <- function(atom) normalize_sign_carrier(atom$gamma2, default = 1)
shape_linear_coef <- function(delta1, delta2) {
  a1 <- abs(delta2 + 1) / 2
  a2 <- abs(delta2 - 1) * abs(delta1 + 1) / 4
  a3 <- abs(delta1 * delta2 + 1) / 2
  delta1 * (a1 + a2 + a3)
}

eval_shape_at_times <- function(theta_j, atom, t_indices, precomp) {
  if (is.null(atom) || length(t_indices) == 0) return(rep(0, length(t_indices)))
  atom <- canonicalize_atom(atom)
  delta1 <- shape_delta1(atom)
  delta2 <- shape_delta2(atom)
  shape_beta  <- atom$shape_beta
  shape_gamma <- atom$shape_gamma

  Mp1 <- length(theta_j)
  if (Mp1 <= 1) {
    return(rep(0, length(t_indices)))
  }
  theta_shape <- theta_j[2:Mp1]

  SZ <- sum(theta_shape^2)
  x_vals <- precomp$x_global[t_indices]

  theta_outer <- tcrossprod(theta_shape)
  DI <- numeric(length(t_indices))
  D_sub <- precomp$D_at_x[2:Mp1, 2:Mp1, , drop = FALSE]
  for (i in seq_along(t_indices)) {
    DI[i] <- sum(theta_outer * D_sub[, , t_indices[i]])
  }

  # Decoupled formula (C++ hiercpd_core.cpp와 동일):
  #   f(x) = linear_coef * x - delta2 * shape_beta * H(x)
  # 여기서 linear_coef = delta1 * shape_gamma + shape_beta * (delta1 + delta2) / 2
  #
  # 이 수식은 네 가지 형상 클래스 모두에서 단조성을 보장한다:
  #   f'(위험지점) = ±sg > 0  (sg > 0이면 항상 성립)

  # 1. 0으로 나누는 것을 방지
  SZ_safe <- ifelse(SZ < 1e-12, 1e-12, SZ)

  # 2. 정규화된 형태 프로파일
  H_vals <- DI / SZ_safe

  # 3. C++ shape_value_single과 동일한 수식
  linear_coef <- delta1 * shape_gamma + shape_beta * (delta1 + delta2) / 2
  term_linear <- linear_coef * x_vals
  term_curve  <- delta2 * shape_beta * H_vals

  return(term_linear - term_curve)
}


#' Normalized shape profile H_θ(x_t) = (θᵀ D(x_t) θ) / max(SZ, 1e-12) at given times
#'
#' SZ = Σ_{m≥1} θ²_m.  This is exactly the per-time quantity used to build the
#' (shape_beta, shape_gamma) design columns in .collect_state_series_data().
#' Pure-R reference; 15_rcpp_bridge.R overrides it with compute_H_at_times_cpp
#' (which reads the cached D-cube) when Rcpp acceleration is active.
compute_H_at_times <- function(theta_j, t_indices, precomp) {
  if (length(t_indices) == 0) return(numeric(0))
  Mp1 <- length(theta_j)
  if (Mp1 <= 1) return(rep(0, length(t_indices)))
  theta_shape <- theta_j[2:Mp1]
  SZ_safe <- max(sum(theta_shape^2), 1e-12)
  theta_outer <- tcrossprod(theta_shape)
  D_sub <- precomp$D_at_x[2:Mp1, 2:Mp1, , drop = FALSE]
  H <- numeric(length(t_indices))
  for (i in seq_along(t_indices)) {
    H[i] <- sum(theta_outer * D_sub[, , t_indices[i]]) / SZ_safe
  }
  H
}

#' Get state intercept gamma^*_{c,j} (scalar, segment-invariant)
get_state_intercept <- function(j, k, params) {
  g <- params$gamma[[j]]
  if (is.null(g) || length(g) == 0 || !is.finite(g[1])) return(0)
  g[1]  # [γ*(c,j)] 스칼라이므로 k 무관
}

#' Compute gamma^*_{j,t} for all time points of series j
compute_gamma_all_timepoints <- function(j, state, params, model) {
  # [γ*(c,j)] 시계열별 스칼라 intercept → 전 시점 동일값
  g <- params$gamma[[j]]
  if (is.null(g) || length(g) == 0 || !is.finite(g[1])) return(numeric(model$T_len))
  rep(g[1], model$T_len)
}

#' Compute gamma^*_{j,t} for a hypothetical state sequence
compute_gamma_hypothetical_interval <- function(j, state_seq, params) {
  # [γ*(c,j)] 스칼라이므로 state_seq 무관하게 동일값 반복
  g <- params$gamma[[j]]
  val <- if (!is.null(g) && length(g) > 0 && is.finite(g[1])) g[1] else 0
  rep(val, length(state_seq))
}

#' Compute f_j for ALL time points of series j
compute_f_all_timepoints <- function(j, state, params, precomp, model) {
  T_len <- model$T_len; c_j <- state$cluster[j]; K_c <- state$K[c_j]
  f_vals <- numeric(T_len); theta_j <- params$theta[j, ]
  
  for (k in 1:K_c) {
    t_in_k <- which(state$S_lower[j, ] == k)
    if (length(t_in_k) == 0) next
    if (k > length(params$atoms[[c_j]]) || is.null(params$atoms[[c_j]][[k]])) {
      f_vals[t_in_k] <- 0; next
    }
    f_vals[t_in_k] <- eval_shape_at_times(theta_j, params$atoms[[c_j]][[k]], t_in_k, precomp)
  }
  f_vals
}

#' Compute mu_{j,t} = alpha + beta_{c_j} + gamma^*_{j,t} + f_j^shape(x_t)
compute_mu_all <- function(j, state, params, precomp, model) {
  c_j <- state$cluster[j]
  params$alpha + params$beta[c_j] +
    compute_gamma_all_timepoints(j, state, params, model) +
    compute_f_all_timepoints(j, state, params, precomp, model)
}

#' Compute f for a hypothetical state sequence
compute_f_hypothetical_interval <- function(j, t_indices, state_seq, state,
                                             params, precomp, model) {
  c_j <- state$cluster[j]; theta_j <- params$theta[j, ]
  f_vals <- numeric(length(t_indices))
  for (k in unique(state_seq)) {
    idx_k <- which(state_seq == k); t_k <- t_indices[idx_k]
    if (length(t_k) == 0) next
    if (k < 1 || k > length(params$atoms[[c_j]]) || is.null(params$atoms[[c_j]][[k]])) {
      f_vals[idx_k] <- 0; next
    }
    f_vals[idx_k] <- eval_shape_at_times(theta_j, params$atoms[[c_j]][[k]], t_k, precomp)
  }
  f_vals
}

cat("02_shape_function.R loaded (global x_t version).\n")
