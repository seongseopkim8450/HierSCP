# Regime paths, changepoints, free observation sets and path bookkeeping

# Equally spaced initial regimes; lower paths start at the cluster path
init_state <- function(J, T_len, C, K_init, m_min, fixed_cluster_sizes = FALSE) {
  if (fixed_cluster_sizes) {
    cluster <- rep(1:C, each = ceiling(J / C))[1:J]
  } else {
    cluster <- sample(1:C, J, replace = TRUE)
  }

  S_upper <- matrix(0L, nrow = C, ncol = T_len)
  tau_upper <- vector("list", C)
  K <- rep(K_init, C)

  for (cc in 1:C) {
    breaks <- round(seq(1, T_len + 1, length.out = K_init + 1))
    tau_upper[[cc]] <- breaks[1:K_init]
    for (k in 1:K_init) {
      t_start <- breaks[k]
      t_end <- breaks[k + 1] - 1
      S_upper[cc, t_start:t_end] <- k
    }
  }

  S_lower <- matrix(0L, nrow = J, ncol = T_len)
  tau_lower <- vector("list", J)

  for (j in 1:J) {
    cc <- cluster[j]
    S_lower[j, ] <- S_upper[cc, ]
    tau_lower[[j]] <- tau_upper[[cc]]
  }

  list(
    S_upper = S_upper,
    S_lower = S_lower,
    K = K,
    cluster = cluster,
    tau_upper = tau_upper,
    tau_lower = tau_lower
  )
}

# Start times of the regimes (the first is t = 1)
extract_changepoints <- function(S_vec) {
  T_len <- length(S_vec)
  cps <- 1L
  for (t in 2:T_len) {
    if (S_vec[t] != S_vec[t - 1]) {
      cps <- c(cps, t)
    }
  }
  cps
}

# Relabel the paths to 1, ..., K and recompute the changepoints
rebuild_changepoints <- function(state, model) {
  C <- model$C; J <- model$J

  for (cc in 1:C) {
    S_U <- state$S_upper[cc, ]

    for (t in 2:length(S_U)) {
      if (S_U[t] < S_U[t-1]) S_U[t] <- S_U[t-1]
    }

    old_labels <- unique(S_U)
    new_map <- setNames(seq_along(old_labels), old_labels)
    S_U_new <- as.integer(new_map[as.character(S_U)])
    state$S_upper[cc, ] <- S_U_new
    state$K[cc] <- max(S_U_new)
    state$tau_upper[[cc]] <- extract_changepoints(S_U_new)

    if (!is.null(state$atom_reorder_map)) rm(state$atom_reorder_map)
    if (!is.null(attr(state, "relabel_map"))) attr(state, "relabel_map") <- NULL
    # old labels, read by ensure_atoms_consistency()
    attr(state, paste0("relabel_map_", cc)) <- old_labels

    j_in_c <- which(state$cluster == cc)
    for (j in j_in_c) {
      S_L <- state$S_lower[j, ]

      for (t in 2:length(S_L)) {
        if (S_L[t] < S_L[t-1]) S_L[t] <- S_L[t-1]
      }

      S_L_new <- integer(length(S_L))
      for (t in seq_along(S_L)) {
        idx <- match(S_L[t], old_labels)
        S_L_new[t] <- if (!is.na(idx)) as.integer(idx) else S_L[t]
      }

      S_L_new <- pmin(pmax(S_L_new, 1L), state$K[cc])

      for (t in 2:length(S_L_new)) {
        if (S_L_new[t] < S_L_new[t-1]) S_L_new[t] <- S_L_new[t-1]
      }

      state$S_lower[j, ] <- S_L_new
      state$tau_lower[[j]] <- extract_changepoints(S_L_new)
    }
  }

  state
}

# Reorder the atoms after relabeling; match atoms and pi_star to K_c
ensure_atoms_consistency <- function(state, params, model) {
  C <- model$C

  for (cc in 1:C) {
    K_c <- state$K[cc]
    old_labels <- attr(state, paste0("relabel_map_", cc))

    if (!is.null(old_labels) && length(old_labels) > 0) {
      old_atoms <- params$atoms[[cc]]
      new_atoms <- vector("list", K_c)
      for (new_k in 1:K_c) {
        old_k <- old_labels[new_k]
        if (!is.na(old_k) && old_k <= length(old_atoms) && !is.null(old_atoms[[old_k]])) {
          new_atoms[[new_k]] <- old_atoms[[old_k]]
        } else {
          new_atoms[[new_k]] <- sample_atom_from_base(params)
        }
      }
      params$atoms[[cc]] <- new_atoms
      attr(state, paste0("relabel_map_", cc)) <- NULL
    }

    while (length(params$atoms[[cc]]) < K_c) {
      params$atoms[[cc]][[length(params$atoms[[cc]]) + 1]] <- sample_atom_from_base(params)
    }
    if (length(params$atoms[[cc]]) > K_c) {
      params$atoms[[cc]] <- params$atoms[[cc]][1:K_c]
    }

    j_in_c <- which(state$cluster == cc)

    while (length(params$pi_star[[cc]]) < K_c) {
      params$pi_star[[cc]] <- c(params$pi_star[[cc]], 1)
    }
    if (length(params$pi_star[[cc]]) > K_c) {
      params$pi_star[[cc]] <- params$pi_star[[cc]][1:K_c]
    }
  }

  params
}

# Upper free set G^U_{k,k+1,c} = {tau_k + m_min, ..., tau_{k+2} - m_min}
compute_free_obs_upper <- function(k, c_idx, state, m_min, virtual_offset = 0L) {
  K_c <- state$K[c_idx]
  tau <- state$tau_upper[[c_idx]]
  T_len <- ncol(state$S_upper)

  if (k >= K_c) return(integer(0))
  if (k > length(tau)) return(integer(0))

  tau_k_eff <- if (k == 1L && virtual_offset > 0L) tau[k] - virtual_offset else tau[k]
  left <- tau_k_eff + m_min
  left <- max(left, tau[k] + 1L)
  if (is.na(left)) return(integer(0))

  if (k + 1 < K_c && (k + 2) <= length(tau)) {
    right <- tau[k + 2] - m_min
  } else {
    right <- T_len - m_min + 1  # tau_{K+1} = T + 1
  }
  if (is.na(right)) return(integer(0))

  if (left > right) return(integer(0))
  left:right
}

# Lower free set G^L_{k,k+1,j}, with the changepoints of series j
compute_free_obs_lower <- function(j, k, state, m_min, virtual_offset = 0L) {
  c_j <- state$cluster[j]
  K_c <- state$K[c_j]
  tau_L <- state$tau_lower[[j]]
  tau_U <- state$tau_upper[[c_j]]
  T_len <- ncol(state$S_lower)

  if (k >= K_c) return(integer(0))
  if (k > length(tau_L)) return(integer(0))

  tau_k_eff <- if (k == 1L && virtual_offset > 0L) tau_L[k] - virtual_offset else tau_L[k]
  left <- tau_k_eff + m_min
  left <- max(left, tau_L[k] + 1L)
  if (is.na(left)) return(integer(0))

  if (k + 1 < K_c && (k + 2) <= length(tau_L)) {
    right <- tau_L[k + 2] - m_min
  } else {
    right <- T_len - m_min + 1  # tau_{K+1} = T + 1
  }
  if (is.na(right)) return(integer(0))

  if (left > right) return(integer(0))
  left:right
}

# Log PAR of a lower-level candidate t; reference tau_0 = upper boundary of k + 1
compute_log_PAR <- function(t, k, j, state, params, G_L) {
  c_j <- state$cluster[j]
  K_c <- state$K[c_j]
  tau_U <- state$tau_upper[[c_j]]

  ps_vec <- if (!is.null(params$pi_star_lower) &&
                !is.null(params$pi_star_lower[[c_j]])) {
    params$pi_star_lower[[c_j]]
  } else {
    params$pi_star[[c_j]]
  }
  pi_star_k <- ps_vec[k]
  pi_star_k <- min(max(pi_star_k, 1e-10), 1 - 1e-10)

  if (k + 1 <= K_c) {
    tau_U_kp1 <- tau_U[k + 1]
  } else {
    tau_U_kp1 <- ncol(state$S_lower) + 1
  }

  g_plus <- max(G_L)

  log_par <- (t - tau_U_kp1) * log(pi_star_k)

  # no leave factor at g+, where (C2) forces the transition
  if (t != g_plus) {
    log_par <- log_par + log(1 - pi_star_k)
  }

  log_par
}

# sigma^2_{j,t} = v_j phi_{j,t}^xi_{j,t}
get_effective_variance <- function(j, t, state, params) {
  v_j <- params$v[[j]]
  phi_jt <- params$phi[j, t]
  xi_jt <- params$xi[j, t]
  v_j * phi_jt^xi_jt
}

get_effective_variance_all <- function(j, state, params) {
  T_len <- ncol(state$S_lower)
  v_j <- params$v[[j]]
  if (!is.finite(v_j) || v_j <= 0) v_j <- 1
  sigma2 <- numeric(T_len)
  for (t in 1:T_len) {
    phi_val <- if (j <= nrow(params$phi) && t <= ncol(params$phi)) params$phi[j, t] else 1
    xi_val <- if (j <= nrow(params$xi) && t <= ncol(params$xi)) params$xi[j, t] else 0
    sigma2[t] <- v_j * phi_val^xi_val
    if (!is.finite(sigma2[t]) || sigma2[t] <= 0) sigma2[t] <- 1
  }
  sigma2
}

# Baseline variance v_j, shared by all regimes (k is not used)
get_state_variance <- function(j, k, params) {
  v_j <- params$v[[j]]
  if (is.finite(v_j) && v_j > 0) return(v_j)
  return(1)
}

# Restore (C1)-(C2) on every path and match parameter lengths to K_c
enforce_all_invariants <- function(state, params, model) {
  C <- model$C; J <- model$J; m_min <- model$m_min; T_len <- model$T_len

  for (cc in 1:C) {
    S_U <- state$S_upper[cc, ]
    S_U <- fix_state_vector(S_U, m_min, T_len)
    state$S_upper[cc, ] <- S_U
    K_c <- max(S_U)
    state$K[cc] <- K_c
    state$tau_upper[[cc]] <- extract_changepoints(S_U)

    j_in_c <- which(state$cluster == cc)
    for (j in j_in_c) {
      S_L <- state$S_lower[j, ]
      S_L <- fix_state_vector(S_L, m_min, T_len)

      K_L <- max(S_L)
      # a series must have as many regimes as its cluster
      if (K_L != K_c) {
        S_L <- S_U
      }

      state$S_lower[j, ] <- S_L
      state$tau_lower[[j]] <- extract_changepoints(S_L)
    }

    while (length(params$atoms[[cc]]) < K_c)
      params$atoms[[cc]][[length(params$atoms[[cc]]) + 1]] <- sample_atom_from_base(params)
    if (length(params$atoms[[cc]]) > K_c)
      params$atoms[[cc]] <- params$atoms[[cc]][1:K_c]

    for (j in j_in_c) {
      # gamma*_j is a scalar: nothing to resize
    }
    while (length(params$pi_star[[cc]]) < K_c)
      params$pi_star[[cc]] <- c(params$pi_star[[cc]], 1)
    if (length(params$pi_star[[cc]]) > K_c)
      params$pi_star[[cc]] <- params$pi_star[[cc]][1:K_c]

    if (!is.null(params$pi_star_lower)) {
      while (length(params$pi_star_lower[[cc]]) < K_c)
        params$pi_star_lower[[cc]] <- c(params$pi_star_lower[[cc]], 0.995)
      if (length(params$pi_star_lower[[cc]]) > K_c)
        params$pi_star_lower[[cc]] <- params$pi_star_lower[[cc]][1:K_c]
    }
  }

  list(state = state, params = params)
}

# Make a path nondecreasing, merge runs shorter than m_min, relabel to 1:K
fix_state_vector <- function(S, m_min, T_len) {
  for (t in 2:T_len) {
    if (is.na(S[t]) || is.na(S[t-1])) { S[t] <- if (is.na(S[t-1])) 1L else S[t-1]; next }
    if (S[t] < S[t-1]) S[t] <- S[t-1]
  }
  if (is.na(S[1])) S[1] <- 1L

  repeat {
    rle_S <- rle(S)
    short_idx <- which(rle_S$lengths < m_min)
    if (length(short_idx) == 0) break

    # the shortest run joins its longer neighbor
    i <- short_idx[which.min(rle_S$lengths[short_idx])]

    end_pos <- cumsum(rle_S$lengths)
    start_pos <- c(1L, end_pos[-length(end_pos)] + 1L)

    left_len  <- if (i > 1) rle_S$lengths[i-1] else 0
    right_len <- if (i < length(rle_S$lengths)) rle_S$lengths[i+1] else 0

    if (left_len >= right_len && i > 1) {
      S[start_pos[i]:end_pos[i]] <- rle_S$values[i-1]
    } else if (i < length(rle_S$lengths)) {
      S[start_pos[i]:end_pos[i]] <- rle_S$values[i+1]
    } else if (i > 1) {
      S[start_pos[i]:end_pos[i]] <- rle_S$values[i-1]
    } else {
      break
    }
  }

  rle_S <- rle(S)
  K_new <- length(rle_S$lengths)
  end_pos <- cumsum(rle_S$lengths)
  start_pos <- c(1L, end_pos[-length(end_pos)] + 1L)
  for (k in 1:K_new) {
    S[start_pos[k]:end_pos[k]] <- as.integer(k)
  }

  S
}

cat("03_state_management.R loaded.\n")
