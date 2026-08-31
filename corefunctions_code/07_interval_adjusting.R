###############################################################################
# 07_interval_adjusting.R  (v10 — sequential merge scan with p⁻/p/p⁺)

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

  res <- split_pass_once(cc, state, params, Y, precomp, model)
  state <- res$state; params <- res$params

  res <- merge_pass_once(cc, state, params, Y, precomp, model)
  state <- res$state; params <- res$params

  
  list(state = state, params = params)
}


.merge_eval <- function(cc, k, I_k, j_in_c, state, params, Y, precomp, model) {
  
  K_now <- state$K[cc]
  has_km1 <- (k - 1L) >= 1L
  has_kp1 <- (k + 1L) <= K_now
  n_k <- length(I_k)
  atoms_c <- params$atoms[[cc]]

  use_cpp <- isTRUE(get0(".hiercpd_rcpp_enabled", ifnotfound = FALSE)) &&
             exists("merge_eval_cpp")

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

    
    terminal <- (k >= K_now)

    
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


propose_atom_spline_conditional <- function(cc, k, tau_k, g_left, g_right,
                                             j_in_c, state, params, Y, precomp, model) {
  hyper <- params$hyper
  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 3.0
  r_beta  <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 0.5
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 2.5
  r_gamma <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 1.0

  Mp1 <- ncol(params$theta)

  
  fallback <- function(d1 = NULL, d2 = NULL) {
    sb <- max(rgamma(1, shape = a_beta, rate = r_beta), 0.01)
    sg <- max(rgamma(1, shape = a_gamma, rate = r_gamma), 0)
    if (is.null(d1)) d1 <- sample(c(-1L, 1L), 1)
    if (is.null(d2)) d2 <- sample(c(-1L, 1L), 1)
    canonicalize_atom(list(gamma1 = d1, gamma2 = d2,
                           shape_beta = sb, shape_gamma = sg))
  }

  g_mid <- as.integer((g_left + g_right) / 2)
  I_left <- tau_k:(g_mid - 1L)
  if (length(I_left) < 3L) return(fallback())
  x_left <- precomp$x_global[I_left]
  n_left <- length(I_left)

  
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

  delta1 <- sample(c(-1L, 1L), 1)
  delta2 <- sample(c(-1L, 1L), 1)

  N_cand <- 50L
  cand_b <- pmax(rgamma(N_cand, shape = a_beta,  rate = r_beta),  0.01)
  cand_g <- pmax(rgamma(N_cand, shape = a_gamma, rate = r_gamma), 0.00)


  use_cpp_ps <- isTRUE(get0(".hiercpd_rcpp_enabled", ifnotfound = FALSE)) &&
                exists("profile_atom_scores_cpp")
  if (use_cpp_ps) {
    
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
    if (is.null(log_scores)) use_cpp_ps <- FALSE   
  }
  if (!use_cpp_ps) {
    
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

### Legacy function(same_shape_merge_sweep)
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

    
    tau_U <- state$tau_upper[[cc]]
    tau_k <- tau_U[k]
    tau_kp1 <- if (k + 1 <= length(tau_U)) tau_U[k + 1] else T_len + 1
    tau_km1 <- tau_U[k - 1]
    I_left <- tau_km1:(tau_k - 1); I_right <- tau_k:(tau_kp1 - 1)
    I_merged <- tau_km1:(tau_kp1 - 1)
    n_l <- length(I_left); n_r <- length(I_right)
    w_l <- n_l / (n_l + n_r); w_r <- 1 - w_l

    ll_current <- 0
    for (j in j_in_c) {
      ll_current <- ll_current +
        log_lik_interval(j, I_left,  rep(k - 1L, length(I_left)),  Y, state, params, precomp, model) +
        log_lik_interval(j, I_right, rep(k,      length(I_right)), Y, state, params, precomp, model)
    }

    avg_atom <- list(
      gamma1      = atom_left$gamma1,
      gamma2      = atom_left$gamma2,
      shape_beta  = w_l * atom_left$shape_beta  + w_r * atom_right$shape_beta,
      shape_gamma = w_l * atom_left$shape_gamma + w_r * atom_right$shape_gamma
    )

    surv_k <- k - 1
    old_atom <- params$atoms[[cc]][[surv_k]]
    params$atoms[[cc]][[surv_k]] <- avg_atom

    ll_merged <- 0
    for (j in j_in_c) {
      ll_merged <- ll_merged +
        log_lik_interval(j, I_merged, rep(surv_k, length(I_merged)), Y, state, params, precomp, model)
    }

    params$atoms[[cc]][[surv_k]] <- old_atom

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

  
  lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0
  log_poisson_merge_bonus <- log(max(K_c, 1)) - log(lambda0)

  log_prior <- c(
    n_k * log(ps_km1) + log(1 - ps_km1) + log_poisson_merge_bonus,  
    max(n_k - m_min, 0) * log(ps_k) + log(1 - ps_k),                
    log(1 - ps_k) + log_poisson_merge_bonus                           
  )
  if (k <= 1)   log_prior[1] <- -Inf
  if (k >= K_c) log_prior[3] <- -Inf
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

  proposed_left_atom <- propose_atom_spline_conditional(
    cc, k, tau_k, g_left, g_right, j_in_c, state, params, Y, precomp, model)
  N_pool <- if (!is.null(model$N_atom)) model$N_atom else 20L
  z_cur  <- if (!is.null(params$z_state) && !is.null(params$z_state[[cc]])) params$z_state[[cc]] else NULL
  left_nb  <- if (!is.null(z_cur) && k > 1L && (k-1L) <= length(z_cur)) z_cur[k-1L] else NA_integer_
  right_nb <- if (!is.null(z_cur) && k <= length(z_cur)) z_cur[k] else NA_integer_
  forbid_split <- c(left_nb, right_nb); forbid_split <- forbid_split[!is.na(forbid_split)]
  atom_pi_c <- if (!is.null(params$atom_pi)) params$atom_pi[[cc]] else NULL
  z_new_k <- .propose_atom_index(N_pool, forbid_split, atom_pi_c)
  # Sequential Gibbs
  log_q_split <- 0
  old_atom_k <- params$atoms[[cc]][[min(k, length(params$atoms[[cc]]))]]
  # k-1 atom for ABSORB: the previous segment's shape
  atom_km1 <- if (has_left) params$atoms[[cc]][[k - 1L]] else NULL

  lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0
  log_poisson_split_penalty <- log(lambda0) - log(K_c + 1)

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

    g_star <- scan$g
    res <- execute_split(cc, k, g_star, state, params, model,
                         proposed_left_atom = proposed_left_atom,
                         z_new_index = z_new_k)
    return(list(state = res$state, params = res$params, did_split = TRUE))
  }

  g_star <- scan$g
  res <- .absorb_left_into_prev(cc, k, g_star, state, params, model)
  list(state = res$state, params = res$params, did_split = FALSE)
}

.absorb_left_into_prev <- function(cc, k, g, state, params, model) {
  if (k <= 1L) return(list(state = state, params = params))   # no previous state
  T_len <- model$T_len
  tau_k <- state$tau_upper[[cc]][k]
  if (is.na(tau_k) || g <= tau_k) return(list(state = state, params = params))

  win <- tau_k:(g - 1L)
  su <- state$S_upper[cc, ]
  sel <- win[su[win] == k]
  if (length(sel) > 0L) state$S_upper[cc, sel] <- as.integer(k - 1L)
  state$tau_upper[[cc]] <- extract_changepoints(state$S_upper[cc, ])

  
  j_in_c <- which(state$cluster == cc)
  for (j in j_in_c) {
    state$S_lower[j, ] <- state$S_upper[cc, ]
    state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])
  }
  list(state = state, params = params)
}

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
.marginal_loglik_left_cd <- function(cumSgg, cumSgb, cumSbb, cumSgr, cumSbr, cumSrr,
                                     sl, m_g, m_b, p_g, p_b, use_trunc) {
  
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
  m + log(sum(exp(levs - m)))  
}

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
  gh <- max(best_g, 1e-6); bh <- max(best_b, 1e-6)
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

  .decode <- function(dec) {
    if (is.na(dec) || dec < 1L) return(NULL)
    if (dec == (2L * n_cand + 1L)) return(list(action = "nosplit", g = NA_integer_))
    if (dec <= n_cand)            return(list(action = "split",  g = G_split[dec]))
    if (dec <= 2L * n_cand)       return(list(action = "absorb", g = G_split[dec - n_cand]))
    NULL
  }

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

  
  split_params <- build_split_left_proposal_params(cc, k, atom_left, state, params)
  km1_params <- params
  if (has_left && !is.null(atom_km1)) {
    km1_params$atoms[[cc]] <- params$atoms[[cc]]
    km1_params$atoms[[cc]][[k]] <- atom_km1
  }
  seg_start <- tau_k
  n_seg     <- (tau_kp1 - 1L) - seg_start + 1L

  ll_new   <- numeric(n_seg)   
  ll_km1   <- numeric(n_seg)   
  ll_right <- numeric(n_seg)  
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
  combos <- list(c(1L,1L), c(1L,-1L), c(-1L,1L), c(-1L,-1L))
  Sgg <- matrix(0, 4, n_seg); Sgb <- matrix(0, 4, n_seg); Sbb <- matrix(0, 4, n_seg)
  Sgr <- matrix(0, 4, n_seg); Sbr <- matrix(0, 4, n_seg); Srr <- numeric(n_seg)
  SlogW <- numeric(n_seg)   # Σ log w 누적 (profile 조건부 정규화상수 C 용)
  Mp1_cd <- ncol(params$theta)
  sigma2_c_cd <- if (!is.null(params$sigma2_gamma_c)) params$sigma2_gamma_c[cc]
                 else if (!is.null(params$sigma2_gamma)) params$sigma2_gamma else 0.01
  if (!is.finite(sigma2_c_cd) || sigma2_c_cd <= 0) sigma2_c_cd <- 0.01
  
  cd_series <- vector("list", length(j_in_c))
  for (ji in seq_along(j_in_c)) {
    j <- j_in_c[ji]
    tt <- seg_start:(tau_kp1 - 1L)
    rj <- as.numeric(Y[j, tt] - params$alpha - params$beta[cc])   
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
    Wj <- sum(wj) + 1.0 / sigma2_c_cd        
    cd_series[[ji]] <- list(r = rj, w = wj, x = xj, H = Hj, Wj = Wj)
  }
  
  for (loc in seq_len(n_seg)) {
    sgg <- numeric(4); sgb <- numeric(4); sbb <- numeric(4)
    sgr <- numeric(4); sbr <- numeric(4); srr <- 0; slogw <- 0
    for (s in cd_series) {
      w <- s$w[loc]; x <- s$x[loc]; H <- s$H[loc]
      slogw <- slogw + log(max(w, 1e-300))
  
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
  eff_ns <- max(n_seg - m_min, 0L)   

  
  loc_of <- function(g) g - seg_start
  
  logp0 <- function(g) {
    sl <- loc_of(g)
    if (sl < 1L || sl > n_seg - 1L) return(-Inf)
    
    LL_left <- cum_new[sl]; LL_right <- total_right - cum_right[sl]
    val <- max((g - tau_k) - m_min, 0L) * log_pi_k + log_1mpi_k +
           log_poisson_split_penalty + log_q_atom + LL_left + LL_right
    if (is.finite(val)) val else -Inf
  }
  
  logpm <- function(g) {
    if (!has_left) return(-Inf)
    sl <- loc_of(g)
    if (sl < 1L || sl > n_seg - 1L) return(-Inf)
    LL_left <- cum_km1[sl]; LL_right <- total_right - cum_right[sl]
    val <- (g - tau_k) * log_pi_m + log_1mpi_m + LL_left + LL_right   
    if (is.finite(val)) val else -Inf
  }
  
  logp_nosplit <- eff_ns * log_pi_k + total_right

  
  logp_plus <- function(g) {
    terms <- c()
    hs <- G_split[G_split > g]
    for (h in hs) {
      sl <- loc_of(h)
      if (sl < 1L || sl > n_seg - 1L) next
      LL_right_h <- total_right - cum_right[sl]
      
      t_split  <- eff_ns * log_pi_k + log_1mpi_k + cum_new[sl] + LL_right_h
      terms <- c(terms, t_split)
      if (has_left) {
        t_absorb <- (h - tau_k) * log_pi_m + log_1mpi_m + cum_km1[sl] + LL_right_h
        terms <- c(terms, t_absorb)
      }
    }
    terms <- c(terms, logp_nosplit)  
    terms <- terms[is.finite(terms)]
    if (length(terms) == 0) return(-Inf)
    M <- max(terms); M + log(sum(exp(terms - M)))   # logsumexp
  }

  
  for (ci in seq_len(n_cand)) {
    g <- G_split[ci]
    terminal <- (ci == n_cand)
    lp0 <- logp0(g)
    lpm <- logpm(g)
    if (terminal) {
    
      cand_lp  <- c(lpm, lp0, logp_nosplit)
      cand_act <- c("absorb", "split", "nosplit")
    } else {
      lpp <- logp_plus(g)
      cand_lp  <- c(lpm, lp0, lpp)
      cand_act <- c("absorb", "split", "defer")
    }
    keep <- is.finite(cand_lp)
    if (!any(keep)) {            
      if (terminal) return(list(action = "nosplit", g = NA_integer_)) else next
    }
    cand_lp  <- cand_lp[keep]; cand_act <- cand_act[keep]
    sel <- sample_categorical_log(cand_lp)
    act <- cand_act[sel]
    if (act == "split")  return(list(action = "split",  g = g))
    if (act == "absorb") return(list(action = "absorb", g = g))
    
  }
  list(action = "nosplit", g = NA_integer_) 
}


# helper: safe atom when k-1 missing (never used in score because has_left=0 path)
old_atom_safe <- function(a) if (is.null(a)) list(gamma1=1,gamma2=1,shape_beta=1.0,shape_gamma=0.0) else a




# execute_merge, execute_split, relabel_after_change

execute_merge <- function(cc, k, direction, state, params, model) {
  
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

.propose_atom_index <- function(N, forbid, atom_pi = NULL) {
  allowed <- setdiff(seq_len(N), forbid)
  if (length(allowed) == 0) allowed <- seq_len(N)  # 안전장치
  if (!is.null(atom_pi) && length(atom_pi) == N) {
    w <- atom_pi[allowed]; w <- pmax(w, 1e-300); w <- w / sum(w)
    return(sample(allowed, 1, prob = w))
  }
  sample(allowed, 1)
}

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

execute_split <- function(cc, k, g, state, params, model, proposed_left_atom = NULL,
                          z_new_index = NULL) {
  T_len <- model$T_len; K_c <- state$K[cc]; m_min <- model$m_min


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

  
  j_in_c <- which(state$cluster == cc)
  for (j in j_in_c) {
    state$S_lower[j, ] <- state$S_upper[cc, ]
    state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])
  }

  
  new_left_atom <- if (is.null(proposed_left_atom)) sample_atom_from_base(params) else proposed_left_atom
  old_atoms <- params$atoms[[cc]]
  new_K <- K_c + 1L
  new_atoms <- vector("list", new_K)
  if (k > 1L) for (kk in 1:(k - 1L)) new_atoms[[kk]] <- old_atoms[[min(kk, length(old_atoms))]]
  new_atoms[[k]] <- new_left_atom
  for (kk in k:K_c) new_atoms[[kk + 1L]] <- old_atoms[[min(kk, length(old_atoms))]]
  params$atoms[[cc]] <- new_atoms

  
  if (!is.null(params$z_state) && !is.null(params$z_state[[cc]])) {
    N <- if (!is.null(model$N_atom)) model$N_atom else 20L
    oz <- params$z_state[[cc]]
    if (length(oz) != K_c) {
      oz <- integer(K_c)
      if (K_c >= 1) oz[1] <- 1L
      if (K_c >= 2) for (kk in 2:K_c) oz[kk] <- if (oz[kk - 1] == 1L) 2L else 1L
    }
    if (!is.null(z_new_index) && z_new_index >= 1L && z_new_index <= N) {
      z_new_k <- as.integer(z_new_index)              
    } else {
      
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

relabel_after_change <- function(cc, state, params, model) {
  T_len <- model$T_len

  old_labels <- sort(unique(state$S_upper[cc, ]))
  K_new      <- length(old_labels)
  max_old <- max(old_labels)
  lut <- integer(max_old)                 
  lut[old_labels] <- seq_len(K_new)

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

  old_atoms <- params$atoms[[cc]]
  new_atoms <- vector("list", K_new)
  for (i in seq_along(old_labels)) {
    ok <- old_labels[i]
    new_atoms[[i]] <- if (ok <= length(old_atoms) && !is.null(old_atoms[[ok]])) {
      old_atoms[[ok]]
    } else sample_atom_from_base(params)
  }
  params$atoms[[cc]] <- new_atoms

  if (!is.null(params$z_state) && !is.null(params$z_state[[cc]])) {
    oz <- params$z_state[[cc]]
    nz <- integer(K_new)
    for (i in seq_along(old_labels)) {
      ok <- old_labels[i]
      nz[i] <- if (ok <= length(oz)) oz[ok] else 1L
    }
    
    params$z_state[[cc]] <- nz
  }

  # pi_weights / pi_star
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
