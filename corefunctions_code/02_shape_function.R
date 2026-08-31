###############################################################################
# 02_shape_function.R
# Shape-restricted smoothing function f_j(x_t) — GLOBAL time normalization
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


  SZ_safe <- ifelse(SZ < 1e-12, 1e-12, SZ)

  H_vals <- DI / SZ_safe

  linear_coef <- delta1 * shape_gamma + shape_beta * (delta1 + delta2) / 2
  term_linear <- linear_coef * x_vals
  term_curve  <- delta2 * shape_beta * H_vals

  return(term_linear - term_curve)
}

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

get_state_intercept <- function(j, k, params) {
  g <- params$gamma[[j]]
  if (is.null(g) || length(g) == 0 || !is.finite(g[1])) return(0)
  g[1]  
}

compute_gamma_all_timepoints <- function(j, state, params, model) {
  g <- params$gamma[[j]]
  if (is.null(g) || length(g) == 0 || !is.finite(g[1])) return(numeric(model$T_len))
  rep(g[1], model$T_len)
}


compute_gamma_hypothetical_interval <- function(j, state_seq, params) {
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
