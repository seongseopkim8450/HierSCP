###############################################################################
# 05_likelihood.R
# Log-likelihood computation — global x_t version
###############################################################################

#' Log-likelihood for a single observation (continuous)
log_lik_continuous <- function(y, mu, sigma2) {
  dnorm_log_var(y, mu, sigma2)
}

#' Log-likelihood for all observations of series j (continuous)
log_lik_series_continuous <- function(j, Y, state, params, precomp, model) {
  T_len <- model$T_len
  c_j <- state$cluster[j]
  mu_j <- tryCatch(
    compute_mu_all(j, state, params, precomp, model),
    error = function(e) rep(params$alpha + params$beta[c_j], T_len)
  )
  sigma2_j <- get_effective_variance_all(j, state, params)
  sigma2_j[!is.finite(sigma2_j) | sigma2_j <= 0] <- 1
  ll <- sum(dnorm_log_var(Y[j, ], mu_j, sigma2_j))
  if (!is.finite(ll)) ll <- -1e300
  ll
}

#' Log-likelihood for a time interval under hypothetical states (continuous)
log_lik_interval_continuous <- function(j, t_indices, state_seq, Y, state,
                                        params, precomp, model,
                                        var_state_override = NULL) {
  if (length(t_indices) == 0) return(0)

  c_j <- state$cluster[j]
  f_vals <- tryCatch(
    compute_f_hypothetical_interval(j, t_indices, state_seq, state,
                                     params, precomp, model),
    error = function(e) rep(0, length(t_indices))
  )
  gamma_vals <- compute_gamma_hypothetical_interval(j, state_seq, params)
  mu_vals <- params$alpha + params$beta[c_j] + gamma_vals + f_vals

  sigma2_vals <- numeric(length(t_indices))
  # [v_j 스칼라화] regime-invariant 분산
  v_j <- params$v[[j]]
  if (!is.finite(v_j) || v_j <= 0) v_j <- 1
  for (i in seq_along(t_indices)) {
    t_i <- t_indices[i]
    phi_val <- if (t_i >= 1 && t_i <= ncol(params$phi)) params$phi[j, t_i] else 1
    xi_val  <- if (t_i >= 1 && t_i <= ncol(params$xi))  params$xi[j, t_i]  else 0
    if (is.na(xi_val)) xi_val <- 0
    sigma2_vals[i] <- v_j * phi_val^xi_val
    if (!is.finite(sigma2_vals[i]) || sigma2_vals[i] <= 0) sigma2_vals[i] <- 1
  }

  ll <- sum(dnorm_log_var(Y[j, t_indices], mu_vals, sigma2_vals))
  if (!is.finite(ll)) ll <- -1e300
  ll
}

#' Log-likelihood for all observations of series j (count)
log_lik_series_count <- function(j, Y, state, params, precomp, model) {
  T_len <- model$T_len
  c_j <- state$cluster[j]
  mu_j <- tryCatch(
    compute_mu_all(j, state, params, precomp, model),
    error = function(e) rep(params$alpha + params$beta[c_j], T_len)
  )
  lambda_j <- exp(mu_j) * params$eta[j, ] * params$gamma_od[j, ]
  lambda_j[lambda_j <= 0 | !is.finite(lambda_j)] <- 1e-10
  ll <- sum(dpois_log(Y[j, ], lambda_j))
  if (!is.finite(ll)) ll <- -1e300
  ll
}

#' Log-likelihood for a time interval (count)
###############################################################################
# Count likelihood — Poisson and NB (collapsed γ_od)
###############################################################################

log_lik_interval_count_poisson <- function(j, t_indices, state_seq, Y, state,
                                            params, precomp, model) {
  if (length(t_indices) == 0) return(0)
  c_j <- state$cluster[j]
  f_vals <- tryCatch(
    compute_f_hypothetical_interval(j, t_indices, state_seq, state,
                                     params, precomp, model),
    error = function(e) rep(0, length(t_indices)))
  gamma_vals <- compute_gamma_hypothetical_interval(j, state_seq, params)
  mu_vals <- params$alpha + params$beta[c_j] + gamma_vals + f_vals
  lambda_vals <- exp(mu_vals)
  lambda_vals[lambda_vals <= 0 | !is.finite(lambda_vals)] <- 1e-10
  ll <- sum(dpois(Y[j, t_indices], lambda_vals, log = TRUE))
  if (!is.finite(ll)) ll <- -1e300
  ll
}

log_lik_interval_count_nb <- function(j, t_indices, state_seq, Y, state,
                                       params, precomp, model) {
  if (length(t_indices) == 0) return(0)
  c_j <- state$cluster[j]
  f_vals <- tryCatch(
    compute_f_hypothetical_interval(j, t_indices, state_seq, state,
                                     params, precomp, model),
    error = function(e) rep(0, length(t_indices)))
  gamma_vals <- compute_gamma_hypothetical_interval(j, state_seq, params)
  mu_vals <- params$alpha + params$beta[c_j] + gamma_vals + f_vals
  lambda_vals <- exp(mu_vals)
  lambda_vals[lambda_vals <= 0 | !is.finite(lambda_vals)] <- 1e-10
  r_od <- if (!is.null(params$r_seg)) {
    k_seg <- state_seq[1]  # segment index for this interval
    params$r_seg[c_j, min(k_seg, ncol(params$r_seg))]
  } else if (!is.null(params$r_od)) params$r_od else 5
  ll <- sum(dnbinom(Y[j, t_indices], size = r_od, mu = lambda_vals, log = TRUE))
  if (!is.finite(ll)) ll <- -1e300
  ll
}

# Default: Poisson (split/merge — symmetric, fair K decision)
log_lik_interval_count <- log_lik_interval_count_poisson

#' Dispatch: default — Poisson for count (split/merge)
log_lik_interval <- function(j, t_indices, state_seq, Y, state, params,
                              precomp, model, var_state_override = NULL) {
  if (model$type == "continuous") {
    log_lik_interval_continuous(j, t_indices, state_seq, Y, state, params,
                                precomp, model, var_state_override)
  } else {
    log_lik_interval_count_poisson(j, t_indices, state_seq, Y, state, params,
                                    precomp, model)
  }
}

#' Dispatch: LA — NB (stable positioning)
log_lik_interval_la <- function(j, t_indices, state_seq, Y, state, params,
                                 precomp, model, var_state_override = NULL) {
  if (model$type == "continuous") {
    log_lik_interval_continuous(j, t_indices, state_seq, Y, state, params,
                                precomp, model, var_state_override)
  } else {
    log_lik_interval_count_nb(j, t_indices, state_seq, Y, state, params,
                               precomp, model)
  }
}

cat("05_likelihood.R loaded (global x_t version).\n")

###############################################################################
# Segmentation quality score — Full Log-Likelihood
#
# Conditional MAP (Fearnhead 2006) Stage 2 requires maximizing
#   P(τ, θ | K_hat, y) ∝ L(y | τ, θ, K_hat) × π(τ, θ | K_hat)
#
# The log-likelihood properly accounts for heterogeneous variance
# structure (state-specific v_{j,k} and outlier inflation φ^ξ),
# unlike MSE which ignores the variance weighting entirely.
#
# Score = Σ_j ℓ_j  (higher is better)
#   where ℓ_j = Σ_t [ -½ log(2π) - ½ log(σ²_{j,t})
#                      - ½(y_{j,t} - μ_{j,t})² / σ²_{j,t} ]
#
# When Rcpp is active, log_lik_series_continuous dispatches to
# log_lik_series_cpp automatically — no performance regression.
###############################################################################

###############################################################################
# Segmentation quality score — Unnormalized Log-Posterior
#
# Conditional MAP (Fearnhead 2006) Stage 2 requires maximizing
#   P(τ, θ | K_hat, y) ∝ L(y | τ, θ, K_hat) × π(τ, θ | K_hat)
#
# 과거의 잘못된 MISE(단순 오차 최소화)를 버리고 진정한 Log-Posterior로 교체.
# Log-Likelihood에 전역 파라미터들의 Log-Prior 패널티를 합산하여 계산합니다.
###############################################################################

compute_log_posterior <- function(state, params, Y, precomp, model, ...) {
  J <- model$J; C <- model$C; M <- model$M
  hyper <- params$hyper

  # 1. Log-Likelihood (분산과 이상치가 완벽히 고려된 우도)
  total_ll <- 0
  for (j in 1:J) {
    total_ll <- total_ll +
      log_lik_series_continuous(j, Y, state, params, precomp, model)
  }
  if (!is.finite(total_ll)) total_ll <- -1e300

  # 2. Log-Prior (Conditional MAP에 필요한 주요 파라미터들의 사전분포)
  total_lp <- 0

  # (a) alpha ~ N(m_alpha, sigma2_alpha)
  total_lp <- total_lp + dnorm_log_var(params$alpha, hyper$m_alpha, hyper$sigma2_alpha)

  # (b) beta_c ~ N(0, lambda2_c) — 클러스터별 실제 수축 분산 사용
  for (cc in 1:C) {
    lam2 <- params$lambda2[cc]
    if (!is.finite(lam2) || lam2 <= 0) lam2 <- 1
    total_lp <- total_lp + dnorm_log_var(params$beta[cc], 0, lam2)
  }

  # (c) theta_{j,m} ~ N(theta0_{c,m}, tau2_{c,m} * exp(-m*r_j))
  # [θ_{j,0}=0 고정] m=0 제외 (상수 기여만 하므로)
  for (j in 1:J) {
    c_j <- state$cluster[j]
    for (mm in 1:M) {
      var_m <- params$tau2[c_j, mm + 1] * exp(-mm * params$r[j])
      if (!is.finite(var_m) || var_m <= 0) var_m <- 1
      total_lp <- total_lp + dnorm_log_var(params$theta[j, mm + 1],
                                            params$theta0[c_j, mm + 1], var_m)
    }
  }

  # (d) gamma*_{j,k} ~ N(0, sigma2_gamma)
  sigma2_g <- params$sigma2_gamma
  if (is.null(sigma2_g) || !is.finite(sigma2_g) || sigma2_g <= 0) sigma2_g <- 0.05
  for (j in 1:J) {
    gamma_j <- params$gamma[[j]]
    if (!is.null(gamma_j) && length(gamma_j) > 0) {
      total_lp <- total_lp + sum(dnorm_log_var(gamma_j, 0, sigma2_g))
    }
  }

  # (e) atom prior: H(atom) for each occupied atom
  for (cc in 1:C) {
    K_c <- state$K[cc]
    for (k in 1:min(K_c, length(params$atoms[[cc]]))) {
      atom_k <- params$atoms[[cc]][[k]]
      if (!is.null(atom_k)) {
        lp_atom <- atom_log_prior(atom_k, params)
        if (is.finite(lp_atom)) total_lp <- total_lp + lp_atom
      }
    }
  }

  # (f) v_j ~ InvGamma(v01, v02) — [v_j 스칼라화] 시리즈당 하나
  for (j in 1:J) {
    v_j <- params$v[[j]]
    if (is.finite(v_j) && v_j > 0) {
      total_lp <- total_lp + dinvgamma_log(v_j, hyper$v01, hyper$v02)
    }
  }
  

  # 3. Unnormalized Log-Posterior
  total_log_posterior <- total_ll + total_lp
  return(total_log_posterior)
}

