# Parameter Updating for the continuous model

# One round of Parameter Updating
param_update_continuous <- function(state, params, Y, precomp, model) {
  # level parameters
  cache_ab <- build_continuous_ab_cache(state, params, Y, precomp, model)
  params <- update_alpha_continuous(state, params, Y, precomp, model, cache_ab)
  params <- update_beta_continuous(state, params, Y, precomp, model, cache_ab)
  params <- update_lambda2(state, params, model)

  if (isTRUE(model$fixed_clusters)) {
    # allocation fixed
  } else if (isTRUE(model$fixed_cluster_sizes)) {
    # sizes fixed: group swaps only, in run_mcmc()
  } else {
    result <- update_cluster_allocation(state, params, Y, precomp, model)
    state <- result$state
    params <- result$params
    params <- update_cluster_weights(state, params, model)
  }

  # Gaussian process coefficients and smoothness
  params <- update_tau2(state, params, model)
  params <- update_eta2(state, params, model)
  params <- update_theta0(state, params, model)
  params <- update_kr_xir(state, params, model)
  params <- update_theta_j(state, params, Y, precomp, model)
  params <- update_rj(state, params, model)

  # shape atoms and magnitude rates
  params <- update_atoms(state, params, Y, precomp, model)
  params <- update_shape_rate_hyper(state, params, model)

  # series effects, centering move and their variances
  params <- update_state_intercepts_continuous(state, params, Y, precomp, model)
  params <- center_gamma_sweep(state, params, model)
  params <- update_sigma2_gamma(state, params, model)

  # robust error block
  base_cache <- build_continuous_base_cache(state, params, Y, precomp, model)
  params <- update_v_continuous(state, params, Y, precomp, model, base_cache)
  params <- update_xi_continuous(state, params, Y, precomp, model, base_cache)
  params <- update_theta_out(state, params, model)
  params <- update_nu(state, params, model)
  params <- update_phi_continuous(state, params, Y, precomp, model, base_cache)

  # stay probabilities
  params <- update_stick_breaking(state, params, model)

  list(state = state, params = params)
}

# gamma*_j + f_{j,k}(x_t) and sigma^2_{j,t} of series j
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

    sigma2_j <- as.numeric(get_effective_variance_all(j, state, params))
  }

  if (length(base_j) != model$T_len) base_j <- rep(0, model$T_len)
  if (length(sigma2_j) != model$T_len) sigma2_j <- rep(1, model$T_len)
  base_j[!is.finite(base_j)] <- 0
  sigma2_j[!is.finite(sigma2_j) | sigma2_j <= 0] <- 1
  list(base = base_j, sigma2 = sigma2_j)
}

# The same for all series
build_continuous_ab_cache <- function(state, params, Y, precomp, model) {
  J <- model$J
  cache <- vector("list", J)
  for (j in 1:J) cache[[j]] <- compute_series_base_sigma_entry(j, state, params, Y, precomp, model)
  cache
}

build_continuous_base_cache <- function(state, params, Y, precomp, model) {
  lapply(build_continuous_ab_cache(state, params, Y, precomp, model), `[[`, "base")
}

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
    if (cc == 1) {
      # reference cluster: beta_1 = 0
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

# Cluster labels; each candidate is scored with theta_j = theta0_c and a plug-in gamma*_j
update_cluster_allocation <- function(state, params, Y, precomp, model) {
  J <- model$J; C <- model$C; M <- model$M

  for (j in 1:J) {
    log_w <- numeric(C)

    old_cluster <- state$cluster[j]
    old_S_lower <- state$S_lower[j, ]
    old_gamma <- params$gamma[[j]][1]
    old_theta <- params$theta[j, ]

    for (cc in 1:C) {
      log_w[cc] <- log(params$p[cc])

      params$theta[j, ] <- params$theta0[cc, ]

      for (mm in 1:M) {
        var_m <- params$tau2[cc, mm + 1] * exp(-mm * params$r[j])
        if (var_m <= 0 || !is.finite(var_m)) var_m <- 1e-10
        log_w[cc] <- log_w[cc] - 0.5 * log(2 * pi * var_m)
      }

      state$cluster[j] <- cc
      state$S_lower[j, ] <- state$S_upper[cc, ]

      f_j_temp <- compute_f_all_timepoints(j, state, params, precomp, model)
      resid_temp <- Y[j, ] - params$alpha - params$beta[cc] - f_j_temp
      new_gamma <- mean(resid_temp)
      if (!is.finite(new_gamma)) new_gamma <- 0
      params$gamma[[j]] <- new_gamma

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

    params$theta[j, ] <- old_theta
    params$gamma[[j]] <- old_gamma
    state$cluster[j] <- old_cluster
    state$S_lower[j, ] <- old_S_lower

    new_c <- sample_categorical_log(log_w)

    if (new_c != state$cluster[j]) {
      state$cluster[j] <- new_c
      state$S_lower[j, ] <- state$S_upper[new_c, ]
      state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])

      params$theta[j, ] <- params$theta0[new_c, ]

      K_new <- state$K[new_c]

      f_j_temp <- compute_f_all_timepoints(j, state, params, precomp, model)
      resid_temp <- Y[j, ] - params$alpha - params$beta[new_c] - f_j_temp
      new_gamma <- mean(resid_temp)
      if (!is.finite(new_gamma)) new_gamma <- 0
      params$gamma[[j]] <- new_gamma
    }
  }

  list(state = state, params = params)
}

# Log-weight of series j in cluster cc, as in update_cluster_allocation()
eval_series_in_cluster_logw <- function(j, cc, state, params, Y, precomp, model) {
  M <- model$M

  lw <- 0
  for (mm in 1:M) {
    var_m <- params$tau2[cc, mm + 1] * exp(-mm * params$r[j])
    if (var_m <= 0 || !is.finite(var_m)) var_m <- 1e-10
    lw <- lw - 0.5 * log(2 * pi * var_m)
  }

  old_theta <- params$theta[j, ]
  old_gamma <- params$gamma[[j]][1]
  old_cl    <- state$cluster[j]
  old_sl    <- state$S_lower[j, ]

  params$theta[j, ] <- params$theta0[cc, ]
  state$cluster[j]  <- cc
  state$S_lower[j, ] <- state$S_upper[cc, ]

  f_j <- compute_f_all_timepoints(j, state, params, precomp, model)
  resid <- Y[j, ] - params$alpha - params$beta[cc] - f_j
  new_gamma <- mean(resid)
  if (!is.finite(new_gamma)) new_gamma <- 0
  params$gamma[[j]] <- new_gamma

  ll <- tryCatch({
    if (model$type == "continuous") {
      log_lik_series_continuous(j, Y, state, params, precomp, model)
    } else {
      log_lik_series_count(j, Y, state, params, precomp, model)
    }
  }, error = function(e) -Inf)
  if (!is.finite(ll)) ll <- -Inf

  params$theta[j, ]  <- old_theta
  params$gamma[[j]]  <- old_gamma
  state$cluster[j]   <- old_cl
  state$S_lower[j, ] <- old_sl

  lw + ll
}

# Move series j to cluster new_c and reset theta_j and gamma*_j
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

  list(state = state, params = params)
}

# Metropolis swap of the cluster labels of two series
update_cluster_pair_swap <- function(state, params, Y, precomp, model,
                                     n_attempts = NULL) {
  J <- model$J; C <- model$C
  if (C < 2) return(list(state = state, params = params))
  if (is.null(n_attempts)) n_attempts <- J

  for (att in 1:n_attempts) {
    i <- sample(J, 1)
    ci <- state$cluster[i]

    other_idx <- which(state$cluster != ci)
    if (length(other_idx) == 0) next
    j <- if (length(other_idx) == 1) other_idx else sample(other_idx, 1)
    cj <- state$cluster[j]

    lw_i_ci <- eval_series_in_cluster_logw(i, ci, state, params, Y, precomp, model)
    lw_j_cj <- eval_series_in_cluster_logw(j, cj, state, params, Y, precomp, model)

    lw_i_cj <- eval_series_in_cluster_logw(i, cj, state, params, Y, precomp, model)
    lw_j_ci <- eval_series_in_cluster_logw(j, ci, state, params, Y, precomp, model)

    log_alpha <- (lw_i_cj + lw_j_ci) - (lw_i_ci + lw_j_cj)

    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      res <- apply_cluster_switch(i, cj, state, params, Y, precomp, model)
      state <- res$state; params <- res$params
      res <- apply_cluster_switch(j, ci, state, params, Y, precomp, model)
      state <- res$state; params <- res$params
    }
  }

  list(state = state, params = params)
}

# Metropolis swap of the clusters of two groups in model$group_members
update_cluster_group_swap <- function(state, params, Y, precomp, model,
                                      n_attempts = NULL) {
  C <- model$C
  if (C < 2) return(list(state = state, params = params))

  group_members <- model$group_members
  if (is.null(group_members)) {
    n_per <- model$J %/% C
    group_members <- lapply(1:C, function(g) {
      start <- (g - 1) * n_per + 1
      end <- if (g < C) g * n_per else model$J
      start:end
    })
  }
  n_groups <- length(group_members)
  if (n_groups < 2) return(list(state = state, params = params))

  if (is.null(n_attempts)) n_attempts <- n_groups

  for (att in 1:n_attempts) {
    g1 <- sample(n_groups, 1)
    c1 <- state$cluster[group_members[[g1]][1]]

    other_groups <- which(sapply(1:n_groups, function(g) {
      state$cluster[group_members[[g]][1]] != c1
    }))
    if (length(other_groups) == 0) next
    g2 <- if (length(other_groups) == 1) other_groups else sample(other_groups, 1)
    c2 <- state$cluster[group_members[[g2]][1]]

    lw_current <- 0
    for (j in group_members[[g1]]) {
      lw_current <- lw_current + eval_series_in_cluster_logw(j, c1, state, params, Y, precomp, model)
    }
    for (j in group_members[[g2]]) {
      lw_current <- lw_current + eval_series_in_cluster_logw(j, c2, state, params, Y, precomp, model)
    }

    lw_swap <- 0
    for (j in group_members[[g1]]) {
      lw_swap <- lw_swap + eval_series_in_cluster_logw(j, c2, state, params, Y, precomp, model)
    }
    for (j in group_members[[g2]]) {
      lw_swap <- lw_swap + eval_series_in_cluster_logw(j, c1, state, params, Y, precomp, model)
    }

    log_alpha <- lw_swap - lw_current

    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      for (j in group_members[[g1]]) {
        res <- apply_cluster_switch(j, c2, state, params, Y, precomp, model)
        state <- res$state; params <- res$params
      }

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

update_tau2 <- function(state, params, model) {
  C <- model$C; M <- model$M
  hyper <- params$hyper

  for (cc in 1:C) {
    j_in_c <- which(state$cluster == cc)
    nc <- length(j_in_c)

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

    for (mm in 1:M) {
      prior_var_scale <- (1 + mm / params$xi_r)^params$kr
      ss <- ss + params$theta0[cc, mm + 1]^2 * prior_var_scale
    }

    shape <- hyper$a_eta0 + M / 2
    scale <- hyper$b_eta0 + ss / 2

    params$eta2[cc] <- rinvgamma(1, shape, scale)
  }
  params
}

update_theta0 <- function(state, params, model) {
  C <- model$C; M <- model$M

  for (cc in 1:C) {
    j_in_c <- which(state$cluster == cc)

    params$theta0[cc, 1] <- 0

    for (mm in 1:M) {
      prior_var <- params$eta2[cc] * (1 + mm / params$xi_r)^(-params$kr)

      sum_data <- 0
      sum_prec <- 0
      for (j in j_in_c) {
        var_jm <- params$tau2[cc, mm + 1] * exp(-mm * params$r[j])

        var_jm <- max(var_jm, 1e-10)

        sum_prec <- sum_prec + 1 / var_jm
        sum_data <- sum_data + params$theta[j, mm + 1] / var_jm
      }

      V_post <- 1 / (1 / prior_var + sum_prec)
      m_post <- V_post * sum_data

      if (!is.finite(V_post) || V_post <= 0) V_post <- prior_var
      if (!is.finite(m_post)) m_post <- 0

      params$theta0[cc, mm + 1] <- rnorm(1, m_post, sqrt(V_post))
    }
  }
  params
}

# Slice sampling of the decay rates r_j
update_rj <- function(state, params, model) {
  J <- model$J; M <- model$M

  for (j in 1:J) {
    c_j <- state$cluster[j]

    log_target <- function(rj) {
      if (rj <= 0) return(-Inf)

      log_prior <- (params$kr - 1) * log(rj) - rj * params$xi_r

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

# k_r fixed at 1.5; xi_r by slice sampling
update_kr_xir <- function(state, params, model) {
  J <- model$J; C <- model$C; M <- model$M

  params$kr <- 1.5

  log_target_xir <- function(xir) {
    if (xir <= 0) return(-Inf)

    log_val <- dgamma(xir, shape = 2, rate = 1, log = TRUE)

    for (j in 1:J) {
      log_val <- log_val + dgamma(params$r[j], shape = params$kr,
                                   rate = xir, log = TRUE)
    }

    for (cc in 1:C) {
      for (mm in 1:M) {
        prior_var <- params$eta2[cc] * (1 + mm / xir)^(-params$kr)
        log_val <- log_val + dnorm_log_var(params$theta0[cc, mm + 1], 0, prior_var)
      }
    }

    log_val
  }

  params$xi_r <- slice_sample_positive(params$xi_r, log_target_xir, w = 0.5, m = 10)
  params
}

# Elliptical slice sampling of (theta_{j,1}, ..., theta_{j,M}); theta_{j,0} = 0
update_theta_j <- function(state, params, Y, precomp, model) {
  J <- model$J; M <- model$M

  for (j in 1:J) {
    c_j <- state$cluster[j]

    params$theta[j, 1] <- 0

    prior_mean_sub <- params$theta0[c_j, 2:(M + 1)]
    prior_cov_diag_sub <- numeric(M)
    for (mm in 1:M) {
      prior_cov_diag_sub[mm] <- params$tau2[c_j, mm + 1] * exp(-mm * params$r[j])
      if (!is.finite(prior_cov_diag_sub[mm]) || prior_cov_diag_sub[mm] <= 0) {
        prior_cov_diag_sub[mm] <- 1e-10
      }
    }

    log_lik_fn <- function(theta_sub_new) {
      old_theta <- params$theta[j, ]
      params$theta[j, ] <- c(0, theta_sub_new)

      ll <- 0
      if (model$type == "continuous") {
        ll <- log_lik_series_continuous(j, Y, state, params, precomp, model)
      } else {
        ll <- log_lik_series_count(j, Y, state, params, precomp, model)
      }

      params$theta[j, ] <- old_theta
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

.biv_pos_logmass <- function(mu, Sig) {
  if (!exists("biv_pos_logmass", mode = "function", inherits = TRUE)) {
    stop("biv_pos_logmass helper is not loaded; source 04_sampling.R first.",
         call. = FALSE)
  }
  biv_pos_logmass(mu, Sig)
}

.biv_pos_mass <- function(mu, Sig) {
  lp <- .biv_pos_logmass(mu, Sig)
  if (!is.finite(lp)) return(0)
  if (lp < log(.Machine$double.xmin)) return(0)
  exp(lp)
}

# Atom update from the data of the regimes it serves: sign pair from the four
# integrated evidences, then (g, b) from the positive-truncated bivariate normal
.update_one_atom_from_data <- function(series_data, params,
                                       a_b, r_b, a_g, r_g,
                                       b_mode0, g_mode0, b_prec0, g_prec0) {
  if (length(series_data) == 0) {
    return(canonicalize_atom(sample_atom_from_base(params)))
  }

  if (!is.finite(g_prec0) || g_prec0 <= 0 ||
      !is.finite(b_prec0) || b_prec0 <= 0) {
    stop("Magnitude-prior precisions must be finite and positive.", call. = FALSE)
  }

  delta_combos <- list(c(1L, 1L), c(1L, -1L), c(-1L, 1L), c(-1L, -1L))
  log_ev   <- numeric(4)
  bg_modes <- vector("list", 4)
  bg_covs  <- vector("list", 4)

  for (di in seq_len(4L)) {
    d1 <- delta_combos[[di]][1L]
    d2 <- delta_combos[[di]][2L]

    r_pool <- numeric(0)
    w_pool <- numeric(0)
    zg_pool <- numeric(0)
    zb_pool <- numeric(0)

    for (s in series_data) {
      z_g_raw <- d1 * s$x
      z_b_raw <- (d1 + d2) / 2 * s$x - d2 * s$H
      r_pool  <- c(r_pool, s$r)
      w_pool  <- c(w_pool, s$w)
      zg_pool <- c(zg_pool, z_g_raw)
      zb_pool <- c(zb_pool, z_b_raw)
    }

    if (!length(r_pool) ||
        length(r_pool) != length(w_pool) ||
        length(r_pool) != length(zg_pool) ||
        length(r_pool) != length(zb_pool) ||
        any(!is.finite(r_pool)) || any(!is.finite(w_pool)) ||
        any(!is.finite(zg_pool)) || any(!is.finite(zb_pool)) ||
        any(w_pool <= 0)) {
      stop("Invalid data supplied to the shape-magnitude update.", call. = FALSE)
    }

    wzg <- w_pool * zg_pool
    wzb <- w_pool * zb_pool
    P11 <- sum(wzg * zg_pool) + g_prec0
    P12 <- sum(wzg * zb_pool)
    P22 <- sum(wzb * zb_pool) + b_prec0
    P <- matrix(c(P11, P12, P12, P22), 2L, 2L)

    if (any(!is.finite(P))) {
      stop("Non-finite Gaussian posterior precision.", call. = FALSE)
    }
    RP <- tryCatch(
      chol(P),
      error = function(e) {
        stop("Gaussian posterior precision is not positive definite.", call. = FALSE)
      }
    )
    Sig <- chol2inv(RP)

    wr <- w_pool * r_pool
    rhs <- c(
      sum(zg_pool * wr) + g_prec0 * g_mode0,
      sum(zb_pool * wr) + b_prec0 * b_mode0
    )
    mu_post <- as.numeric(Sig %*% rhs)
    if (any(!is.finite(mu_post))) {
      stop("Non-finite Gaussian posterior mean.", call. = FALSE)
    }

    bg_modes[[di]] <- mu_post
    bg_covs[[di]] <- Sig

    log_detP <- 2 * sum(log(diag(RP)))
    quad_post <- sum(rhs * mu_post)
    quad_0 <- g_mode0^2 * g_prec0 + b_mode0^2 * b_prec0
    log_ev_untrunc <-
      0.5 * (log(g_prec0) + log(b_prec0) - log_detP) +
      0.5 * (quad_post - quad_0) -
      0.5 * sum(w_pool * r_pool^2)

    # the N+ normalizers are common to the four sign pairs and cancel
    log_pos <- .biv_pos_logmass(mu_post, Sig)
    log_ev[di] <- log_ev_untrunc + log_pos
    if (is.na(log_ev[di]) || log_ev[di] == Inf) {
      stop("Invalid integrated sign evidence.", call. = FALSE)
    }
  }

  if (!any(is.finite(log_ev))) {
    stop("All four integrated sign cells have zero/non-finite weight.",
         call. = FALSE)
  }
  chosen <- sample_categorical_log(log_ev)
  d1 <- delta_combos[[chosen]][1L]
  d2 <- delta_combos[[chosen]][2L]
  mu <- bg_modes[[chosen]]
  S <- bg_covs[[chosen]]

  if (length(mu) != 2L || any(!is.finite(mu)) ||
      !is.matrix(S) || !identical(dim(S), c(2L, 2L)) ||
      any(!is.finite(S)) || !isSymmetric(S)) {
    stop("Invalid selected Gaussian parameters.", call. = FALSE)
  }
  if (!exists("rbivnorm_pos_exact", mode = "function", inherits = TRUE)) {
    stop("rbivnorm_pos_exact helper is not loaded; source 04_sampling.R first.",
         call. = FALSE)
  }

  candidate <- rbivnorm_pos_exact(mu, S)
  g_new <- candidate[1L]
  b_new <- candidate[2L]

  return(canonicalize_atom(list(
    gamma1 = d1,
    gamma2 = d2,
    shape_beta = b_new,
    shape_gamma = g_new
  )))
}

# Residuals y - alpha - beta_c - gamma*_j, weights, x_t and H_j(x_t) in regime k
.collect_state_series_data <- function(cc, k, j_in_c, state, params, Y, precomp) {
  series_data <- list()
  for (j in j_in_c) {
    t_idx <- which(state$S_lower[j, ] == k)
    if (length(t_idx) == 0) next

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

# Slice sampling of the magnitude rates r_b and r_g given the occupied atoms
update_shape_rate_hyper <- function(state, params, model) {
  hyper <- params$hyper
  C <- model$C
  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 2.0
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 1.5
  a0_beta  <- if (!is.null(hyper$a0_rbeta))  hyper$a0_rbeta  else 3.0
  b0_beta  <- if (!is.null(hyper$b0_rbeta))  hyper$b0_rbeta  else 2.0
  a0_gamma <- if (!is.null(hyper$a0_rgamma)) hyper$a0_rgamma else 4.0
  b0_gamma <- if (!is.null(hyper$b0_rgamma)) hyper$b0_rgamma else 2.0
  if (a_beta <= 1 || a_gamma <= 1) stop("Magnitude prior requires a_beta,a_gamma > 1.", call. = FALSE)

  b_vals <- numeric(0); g_vals <- numeric(0)
  use_pool <- !is.null(params$atom_pool) && !is.null(params$z_state)
  for (cc in seq_len(C)) {
    if (use_pool && !is.null(params$atom_pool[[cc]]) && !is.null(params$z_state[[cc]])) {
      z <- params$z_state[[cc]]
      assigned <- unique(z[z >= 1L & z <= length(params$atom_pool[[cc]])])
      for (n in assigned) {
        at <- params$atom_pool[[cc]][[n]]
        if (is.null(at)) next
        if (is.finite(at$shape_beta)  && at$shape_beta  > 0) b_vals <- c(b_vals, at$shape_beta)
        if (is.finite(at$shape_gamma) && at$shape_gamma > 0) g_vals <- c(g_vals, at$shape_gamma)
      }
    } else {
      for (at in params$atoms[[cc]]) {
        if (is.null(at)) next
        if (is.finite(at$shape_beta)  && at$shape_beta  > 0) b_vals <- c(b_vals, at$shape_beta)
        if (is.finite(at$shape_gamma) && at$shape_gamma > 0) g_vals <- c(g_vals, at$shape_gamma)
      }
    }
  }

  # log p(r | x) = log Gamma(r; a0, b0) + N log r - sum_i {r x_i - (a - 1)}^2 / {2(a - 1)}
  update_one_rate <- function(current, x, a, a0, b0) {
    if (!length(x)) return(rgamma(1L, shape = a0, rate = b0))
    log_target <- function(r) {
      if (!is.finite(r) || r <= 0) return(-Inf)
      dgamma(r, shape = a0, rate = b0, log = TRUE) +
        length(x) * log(r) -
        sum((r * x - (a - 1))^2) / (2 * (a - 1))
    }
    slice_sample_positive(current, log_target, w = 0.35, m = 20L)
  }

  params$hyper$b_shape_beta <- update_one_rate(
    params$hyper$b_shape_beta, b_vals, a_beta, a0_beta, b0_beta)
  params$hyper$b_shape_gamma <- update_one_rate(
    params$hyper$b_shape_gamma, g_vals, a_gamma, a0_gamma, b0_gamma)
  params
}

# Blocked Gibbs update of the atom pool (Ishwaran and James, 2001)
update_atoms <- function(state, params, Y, precomp, model) {
  C <- model$C
  hyper <- params$hyper
  a_b <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 2.0
  r_b <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 1.5
  a_g <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 1.5
  r_g <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 2.0
  N <- if (!is.null(model$N_atom)) model$N_atom else 20L

 .prior_mode_prec <- function(a, r) {
    if (!is.finite(a) || a <= 1 || !is.finite(r) || r <= 0) {
      stop("Actual positive-truncated Gaussian magnitude prior requires a>1 and r>0.", call. = FALSE)
    }
    c(mode = (a - 1) / r, prec = r^2 / (a - 1))
  }
  .pr_b <- .prior_mode_prec(a_b, r_b)
  .pr_g <- .prior_mode_prec(a_g, r_g)
  b_mode0 <- .pr_b[["mode"]]; b_prec0 <- .pr_b[["prec"]]
  g_mode0 <- .pr_g[["mode"]]; g_prec0 <- .pr_g[["prec"]]

  if (is.null(params$atom_pool))     params$atom_pool     <- vector("list", C)
  if (is.null(params$z_state))       params$z_state       <- vector("list", C)
  if (is.null(params$atom_V))        params$atom_V        <- vector("list", C)
  if (is.null(params$atom_pi))       params$atom_pi       <- vector("list", C)
  if (is.null(params$atom_varsigma)) params$atom_varsigma <- rep(1.0, C)

  for (cc in 1:C) {
    K_c <- state$K[cc]
    j_in_c <- which(state$cluster == cc)

    pool <- params$atom_pool[[cc]]
    if (is.null(pool) || length(pool) != N) {
      pool <- vector("list", N)
      for (n in 1:N) pool[[n]] <- canonicalize_atom(sample_atom_from_base(params))
    }
    z <- params$z_state[[cc]]
    if (is.null(z) || length(z) != K_c) {
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

    # occupied entries from their data, empty entries from F0
    state_data <- vector("list", K_c)
    for (k in 1:K_c) {
      state_data[[k]] <- .collect_state_series_data(
        cc, k, j_in_c, state, params, Y, precomp)
    }
    for (n in 1:N) {
      states_n <- which(z == n)
      if (length(states_n) == 0) {
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

    log_pi <- numeric(N)
    log_one_minus <- log(pmax(1 - Vv, 1e-300))
    cum <- 0
    for (n in 1:N) {
      log_pi[n] <- log(max(Vv[n], 1e-300)) + cum
      cum <- cum + log_one_minus[n]
    }
    # allocation z_k of each regime; neighboring regimes may not share an entry
    if (K_c >= 1) {
      for (k in 1:K_c) {
        z_prev <- if (k > 1) z[k - 1] else NA_integer_
        z_next <- if (k < K_c) z[k + 1] else NA_integer_
        t_in_k_upper <- seq_len(ncol(state$S_lower))
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

    # sticks V_n | z ~ Beta(1 + m_n, varsigma + sum_{l > n} m_l)
    m_n <- tabulate(z, nbins = N)
    Vnew <- numeric(N)
    tail_sum <- sum(m_n)
    for (n in 1:N) {
      tail_after <- tail_sum - cumsum(m_n)[n]
      if (n < N) {
        Vnew[n] <- rbeta(1, 1 + m_n[n], vs_c + tail_after)
        Vnew[n] <- min(max(Vnew[n], 1e-10), 1 - 1e-10)
      } else Vnew[n] <- 1
    }
    Vv <- Vnew
    if (length(Vv) != N || any(!is.finite(Vv)) ||
        any(Vv <= 0) || any(Vv > 1) || Vv[N] != 1 ||
        (N > 1L && any(Vv[seq_len(N - 1L)] >= 1))) {
      stop("Invalid finite atom stick-breaking state.", call. = FALSE)
    }

    log_pi <- numeric(N)
    log_remaining <- 0
    for (n in seq_len(N)) {
      log_pi[n] <- log(Vv[n]) + log_remaining
      if (n < N) {
        log_remaining <- log_remaining + log1p(-Vv[n])
      }
    }

    # pool concentration by slice sampling
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

    atoms_view <- vector("list", K_c)
    for (k in 1:K_c) atoms_view[[k]] <- pool[[ z[k] ]]

    params$atom_pool[[cc]]     <- pool
    params$z_state[[cc]]       <- z
    params$atom_V[[cc]]        <- Vv
    params$atom_pi[[cc]]       <- exp(log_pi)
    params$atom_varsigma[cc]   <- vs_c
    params$atoms[[cc]]         <- atoms_view
  }
  params
}

# Gibbs update of gamma*_j given the cluster variance sigma2_gamma_c
update_state_intercepts_continuous <- function(state, params, Y, precomp, model) {
  for (j in 1:model$J) {
    c_j <- state$cluster[j]

    sigma2_c <- params$sigma2_gamma_c[c_j]
    if (!is.finite(sigma2_c) || sigma2_c <= 0) sigma2_c <- 0.01

    f_j <- compute_f_all_timepoints(j, state, params, precomp, model)
    r_j <- Y[j, ] - params$alpha - params$beta[c_j] - f_j

    sigma2_j <- get_effective_variance_all(j, state, params)
    sigma2_j[!is.finite(sigma2_j) | sigma2_j <= 0] <- 1
    w_j <- 1 / sigma2_j

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

# sigma2_gamma_c | . ~ InvGamma; sigma2_gamma is their mean
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

  params$sigma2_gamma <- mean(params$sigma2_gamma_c)

  params
}

dmvnorm_log <- function(x, mu, Sigma) {
  d <- length(x)
  diff <- x - mu
  -0.5 * d * log(2 * pi) - 0.5 * log(det(Sigma)) -
    0.5 * as.numeric(t(diff) %*% solve(Sigma) %*% diff)
}

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

# Outlier proportion (vartheta)
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

# Variance multipliers phi_{j,t} from their truncated inverse-gamma conditionals
update_phi_continuous <- function(state, params, Y, precomp, model,
                                  base_cache = NULL) {
  J <- model$J
  T_len <- model$T_len
  nu <- params$nu
  if (length(nu) != 1L || !is.finite(nu) || nu < 1 || nu > 40) {
    stop("update_phi_continuous: nu must lie in [1, 40].")
  }

  for (j in seq_len(J)) {
    if (is.null(base_cache)) {
      base_j <- compute_series_base_sigma_entry(
        j, state, params, Y, precomp, model)$base
    } else {
      base_j <- base_cache[[j]]
    }
    v_j <- params$v[[j]]
    xi_j <- params$xi[j, ]
    if (length(v_j) != 1L || !is.finite(v_j) || v_j <= 0) {
      stop("update_phi_continuous: invalid series variance.")
    }
    if (length(xi_j) != T_len || anyNA(xi_j) ||
        any(!xi_j %in% c(0L, 1L))) {
      stop("update_phi_continuous: xi must contain only 0 and 1.")
    }
    mu_j <- params$alpha + params$beta[state$cluster[j]] + base_j
    e_j <- Y[j, ] - mu_j
    if (length(e_j) != T_len || any(!is.finite(e_j))) {
      stop("update_phi_continuous: invalid residual vector.")
    }

    shape_j <- (nu + xi_j) / 2
    scale_j <- nu / 2 + xi_j * e_j^2 / (2 * v_j)
    if (any(!is.finite(scale_j)) || any(scale_j <= 0)) {
      stop("update_phi_continuous: invalid inverse-gamma scale.")
    }

    params$phi[j, ] <- rtrunc_invgamma_lower1(shape_j, scale_j)
  }
  params
}

# nu by slice sampling with the phi of unflagged points integrated out (partially
# collapsed); update_phi_continuous() then refreshes all phi
update_nu <- function(state, params, model) {
  if (length(params$nu) != 1L || !is.finite(params$nu) ||
      params$nu < 1 || params$nu > 40) {
    stop("update_nu: current nu must lie in [1, 40].")
  }
  if (anyNA(params$xi) || any(!params$xi %in% c(0L, 1L))) {
    stop("update_nu: xi must contain only 0 and 1.")
  }

  phi_active <- params$phi[params$xi == 1L]
  if (any(!is.finite(phi_active)) || any(phi_active < 1)) {
    stop("update_nu: flagged scales must lie in [1, infinity).")
  }
  N1 <- length(phi_active)

  # no flagged points: draw from the Uniform[1, 40] prior
  if (N1 == 0L) {
    params$nu <- runif(1, 1, 40)
    return(params)
  }

  log_target <- function(nu) {
    if (length(nu) != 1L || !is.finite(nu) || nu < 1 || nu > 40) return(-Inf)
    a <- nu / 2
    log_Z <- pgamma(1, shape = a, rate = a,
                    lower.tail = TRUE, log.p = TRUE)
    if (!is.finite(log_Z)) return(-Inf)
    log_density <- dinvgamma_log(phi_active, shape = a, scale = a)
    if (any(!is.finite(log_density))) return(-Inf)
    sum(log_density) - N1 * log_Z
  }

  params$nu <- slice_sample(params$nu, log_target, w = 2, m = 5,
                            lower = 1, upper = 40)
  params
}

# Free stays (n_k - m_min) and free departures of each nonterminal regime; a
# departure is free only if the next regime is longer than m_min
.count_free_sojourn_decisions <- function(S, K_c, m_min) {
  if (!is.numeric(K_c) || length(K_c) != 1L || !is.finite(K_c) ||
      K_c < 1 || K_c != as.integer(K_c)) {
    stop("Free-decision counts require a positive integer regime count.", call. = FALSE)
  }
  if (!is.numeric(m_min) || length(m_min) != 1L || !is.finite(m_min) ||
      m_min < 1 || m_min != as.integer(m_min)) {
    stop("Free-decision counts require a positive integer minimum length.", call. = FALSE)
  }
  if (!is.numeric(S) || !length(S) || anyNA(S) || any(!is.finite(S)) ||
      any(S < 1 | S != as.integer(S))) {
    stop("Free-decision counts require a finite integer state path.", call. = FALSE)
  }
  K_c <- as.integer(K_c)
  m_min <- as.integer(m_min)
  runs <- rle(as.integer(S))
  if (!identical(runs$values, seq_len(K_c))) {
    stop("State path must visit regimes 1:K consecutively without returning.", call. = FALSE)
  }
  if (any(runs$lengths < m_min)) {
    stop("State path violates the minimum-segment-length constraint.", call. = FALSE)
  }

  n_stay <- numeric(K_c)
  n_trans <- numeric(K_c)
  if (K_c > 1L) {
    nonterminal <- seq_len(K_c - 1L)
    n_stay[nonterminal] <- runs$lengths[nonterminal] - m_min
    n_trans[nonterminal] <- as.numeric(runs$lengths[nonterminal + 1L] > m_min)
  }
  list(n_stay = n_stay, n_trans = n_trans)
}

# log of a Gamma(shape, 1) draw; for shape < 1 via G(shape + 1) U^(1/shape)
.draw_log_gamma_unit <- function(shape) {
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0) {
    stop("Log-Gamma draw requires a finite positive shape.", call. = FALSE)
  }
  if (shape < 1) {
    g <- rgamma(1, shape = shape + 1, rate = 1)
    u <- runif(1)
    if (!is.finite(g) || g <= 0 || !is.finite(u) || u <= 0 || u >= 1) {
      stop("Invalid base draw in log-Gamma sampler.", call. = FALSE)
    }
    ans <- log(g) + log(u) / shape
  } else {
    g <- rgamma(1, shape = shape, rate = 1)
    if (!is.finite(g) || g <= 0) {
      stop("Invalid base draw in log-Gamma sampler.", call. = FALSE)
    }
    ans <- log(g)
  }
  if (!is.finite(ans)) {
    stop("Log-Gamma draw exceeds floating-point range.", call. = FALSE)
  }
  ans
}

# Beta(shape1, shape2) draw as A/(A + B); log(1 - v) is kept for the varsigma update
.draw_beta_persistence <- function(shape1, shape2) {
  log_a <- .draw_log_gamma_unit(shape1)
  log_b <- .draw_log_gamma_unit(shape2)
  d <- log_a - log_b
  if (!is.finite(d)) {
    stop("Log-Beta draw exceeds floating-point range.", call. = FALSE)
  }
  if (d >= 0) {
    correction <- log1p(exp(-d))
    log_value <- -correction
    log_complement <- -d - correction
    value <- -expm1(log_complement)
  } else {
    correction <- log1p(exp(d))
    log_value <- d - correction
    log_complement <- -correction
    value <- exp(log_value)
  }
  list(value = min(max(value, 1e-10), 1 - 1e-10),
       log_value = log_value, log_complement = log_complement)
}

# Upper sticks and varsigma_c, and the lower stay probabilities
update_stick_breaking <- function(state, params, model) {
  C <- model$C
  m_min <- model$m_min
  hyper <- params$hyper

  # upper level: v_k ~ Beta(1 + free stays, varsigma_c + free departures)
  for (cc in seq_len(C)) {
    K_c <- state$K[cc]
    counts_U <- .count_free_sojourn_decisions(state$S_upper[cc, ], K_c, m_min)
    n_stay <- counts_U$n_stay
    n_trans <- counts_U$n_trans

    varsigma_c <- params$varsigma[cc]
    if (!is.finite(varsigma_c) || varsigma_c <= 0) {
      stop("Invalid upper stick-breaking concentration.", call. = FALSE)
    }

    v_k <- numeric(K_c)
    log_one_minus_v <- numeric(K_c)
    if (K_c > 1L) {
      for (k in seq_len(K_c - 1L)) {
        draw <- .draw_beta_persistence(1 + n_stay[k],
                                        varsigma_c + n_trans[k])
        v_k[k] <- draw$value
        log_one_minus_v[k] <- draw$log_complement
      }
    }
    v_k[K_c] <- 1

    pi_weights <- numeric(K_c)
    remaining <- 1
    for (k in seq_len(K_c)) {
      pi_weights[k] <- v_k[k] * remaining
      remaining <- remaining * (1 - v_k[k])
    }
    params$pi_weights[[cc]] <- pi_weights
    params$pi_star[[cc]] <- pmin(pmax(v_k, 1e-10), 1 - 1e-10)

    if (K_c > 1L) {
      # varsigma_c ~ Gamma(a + K - 1, b - sum_{k < K} log(1 - v_k))
      post_shape <- hyper$a_varsigma + K_c - 1L
      post_rate <- hyper$b_varsigma - sum(log_one_minus_v[seq_len(K_c - 1L)])
      if (!is.finite(post_rate) || post_rate <= 0) {
        stop("Invalid upper concentration posterior rate.", call. = FALSE)
      }
      params$varsigma[cc] <- rgamma(1, shape = post_shape, rate = post_rate)
    } else {
      params$varsigma[cc] <- rgamma(1, shape = hyper$a_varsigma,
                                   rate = hyper$b_varsigma)
    }
  }

  # lower level: pi*L_k ~ Beta(a + free stays, b + free departures), pooled over the cluster
  a_pi_L <- if (!is.null(hyper$a_pi_lower)) hyper$a_pi_lower else 50
  b_pi_L <- if (!is.null(hyper$b_pi_lower)) hyper$b_pi_lower else 1
  if (is.null(params$pi_star_lower)) params$pi_star_lower <- vector("list", C)

  for (cc in seq_len(C)) {
    K_c <- state$K[cc]
    j_in_c <- which(state$cluster == cc)
    n_stay_L <- numeric(K_c)
    n_trans_L <- numeric(K_c)

    for (j in j_in_c) {
      counts_L <- .count_free_sojourn_decisions(state$S_lower[j, ], K_c, m_min)
      n_stay_L <- n_stay_L + counts_L$n_stay
      n_trans_L <- n_trans_L + counts_L$n_trans
    }

    pi_L <- numeric(K_c)
    if (K_c > 1L) {
      for (k in seq_len(K_c - 1L)) {
        draw <- .draw_beta_persistence(a_pi_L + n_stay_L[k],
                                        b_pi_L + n_trans_L[k])
        pi_L[k] <- draw$value
      }
    }
    pi_L[K_c] <- 1
    params$pi_star_lower[[cc]] <- pi_L
  }

  params
}

# Centering move (Papaspiliopoulos et al., 2007): zero-mean gamma*_j within each
# cluster and beta_1 = 0, with the fitted mean unchanged
center_gamma_sweep <- function(state, params, model) {
  C <- model$C

  for (cc in 1:C) {
    j_in_c <- which(state$cluster == cc)
    if (length(j_in_c) == 0) next

    gamma_vals <- sapply(j_in_c, function(j) {
      g <- params$gamma[[j]]
      if (is.null(g) || length(g) == 0 || !is.finite(g[1])) 0 else g[1]
    })

    gamma_bar <- mean(gamma_vals)

    for (j in j_in_c) {
      params$gamma[[j]] <- params$gamma[[j]] - gamma_bar
    }

    params$beta[cc] <- params$beta[cc] + gamma_bar
  }

  delta <- params$beta[1]
  params$alpha <- params$alpha + delta
  params$beta <- params$beta - delta

  params
}
