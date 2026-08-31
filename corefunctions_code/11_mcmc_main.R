

#' Run the full MCMC sampler
#'
#' @param Y Data matrix (J x T)
#' @param model Model specification
#' @param n_iter Total number of MCMC iterations
#' @param n_burnin Number of burn-in iterations
#' @param n_thin Thinning interval
#' @param N_launch Number of launch state refinement rounds (default 3)
#' @param verbose Print progress every verbose iterations (0 = silent)
#' @return List with posterior samples and diagnostics



init_atom_dp_fields <- function(params, state, model) {
  C <- model$C
  N <- if (!is.null(model$N_atom)) model$N_atom else 20L

  params$atom_pool     <- vector("list", C)  
  params$z_state       <- vector("list", C)  
  params$atom_V        <- vector("list", C)  
  params$atom_pi       <- vector("list", C)  
  params$atom_varsigma <- rep(1.0, C)        

  for (cc in 1:C) {
    K_c <- state$K[cc]
    pool <- vector("list", N)
    for (n in 1:N) pool[[n]] <- canonicalize_atom(sample_atom_from_base(params))
    params$atom_pool[[cc]] <- pool

    z <- integer(K_c)
    if (K_c >= 1) z[1] <- 1L
    if (K_c >= 2) for (k in 2:K_c) z[k] <- if (z[k-1]==1L) 2L else 1L
    params$z_state[[cc]] <- z

    Vv <- c(rbeta(N-1, 1, params$atom_varsigma[cc]), 1)
    Vv <- pmin(pmax(Vv, 1e-10), 1-1e-10); Vv[N] <- 1
    params$atom_V[[cc]] <- Vv

    log_pi <- numeric(N); cum <- 0
    lom <- log(pmax(1-Vv,1e-300))
    for (n in 1:N) { log_pi[n] <- log(max(Vv[n],1e-300))+cum; cum <- cum+lom[n] }
    params$atom_pi[[cc]] <- exp(log_pi)

    av <- vector("list", K_c)
    for (k in 1:K_c) av[[k]] <- pool[[ z[k] ]]
    params$atoms[[cc]] <- av
  }
  params
}


run_mcmc <- function(Y, model, n_iter = 5000, n_burnin = 1000, n_thin = 5,
                     N_launch = 3, n_warmup = 20, ia_every = 5, verbose = 100,
                     diagnostic_callback = NULL) {
  
  cat("=== CPD-HierBSAR MCMC Sampler ===\n")
  cat(sprintf("Model type: %s\n", model$type))
  cat(sprintf("J=%d series, T=%d time points, C=%d clusters\n", 
              model$J, model$T_len, model$C))
  cat(sprintf("K_init=%d, M=%d, m_min=%d\n", model$K_init, model$M, model$m_min))
  cat(sprintf("Iterations: %d (burn-in: %d, thin: %d)\n", n_iter, n_burnin, n_thin))
  cat(sprintf("Warm-up: %d iters (param only), IA every %d iters\n", n_warmup, ia_every))
  cat(sprintf("Launch state refinement rounds: %d\n", N_launch))
  
  
  if (isTRUE(model$fixed_clusters)) {
    cat("Cluster allocation: FIXED (no pair swap, no reallocation)\n")
  } else if (isTRUE(model$fixed_cluster_sizes)) {
    cat("Cluster allocation: fixed sizes (pair swap only)\n")
  } else {
    cat("Cluster allocation: free (Gibbs + pair swap)\n")
  }
  lambda0_K <- if (!is.null(model$lambda0_K)) model$lambda0_K else 2.0
  cat(sprintf("K prior: Poisson(lambda0=%.1f), K_max=%d\n", lambda0_K, model$K_max))
  
  # Initialize
  init <- initialize_model(Y, model)
  state <- init$state
  params <- init$params
  precomp <- init$precomp
  # Use standardized data for all MCMC computations
  Y_work <- if (!is.null(init$Y_std)) init$Y_std else Y
  
  # Storage for posterior samples
  n_save <- floor((n_iter - n_burnin) / n_thin)
  samples <- init_storage(n_save, model)
  
  save_idx <- 0
  start_time <- Sys.time()
  
  # ================================================================
  # Conditional MAP registry (Group-Based Independent):
  #   Stage 1: Track K vector frequency across post-burnin iterations
  #   Stage 2: For each distinct K, track the best log-posterior iteration
  #   After MCMC: Marginalize K_g independently per group, select best state
  # ================================================================
  cmap_registry <- list()
  
  # Choose parameter update function based on model type
  param_update_fn <- if (model$type == "continuous") {
    param_update_continuous
  } else {
    param_update_count
  }
  
  if (model$type == "continuous") {
  params <- init_atom_dp_fields(params, state, model)
 }


  # ================================================================
  # Warm-up — Parameter Updating ONLY (no state changes)
  # This calibrates alpha, beta, theta, atoms to the initial segmentation
  # before any boundaries are allowed to move.
  # ================================================================
  if (n_warmup > 0) {
    cat(sprintf("Warm-up phase: %d iterations of Parameter Updating only...\n", n_warmup))
    for (w in 1:n_warmup) {
      pu_result <- param_update_fn(state, params, Y_work, precomp, model)
      state <- pu_result$state; params <- pu_result$params
    }
    cat(sprintf("Warm-up done. alpha=%.3f\n", params$alpha))
  }
  
  for (iter in 1:n_iter) {
    
    if (iter %% ia_every == 0) {
      
      # ==============================================================
      # IA round: split→merge → LA refinement → Group Swap
      # ==============================================================
      
      res_ia <- tryCatch(interval_adjusting(state, params, Y_work, precomp, model),
        error = function(e) stop(sprintf("[iter %d] interval_adjusting: %s", iter, conditionMessage(e))))
      state <- res_ia$state; params <- res_ia$params
      
      state <- tryCatch(rebuild_changepoints(state, model),
        error = function(e) stop(sprintf("[iter %d] rebuild_changepoints: %s", iter, conditionMessage(e))))
      params <- tryCatch(ensure_atoms_consistency(state, params, model),
        error = function(e) stop(sprintf("[iter %d] ensure_atoms_consistency: %s", iter, conditionMessage(e))))
      res_inv <- tryCatch(enforce_all_invariants(state, params, model),
        error = function(e) stop(sprintf("[iter %d] enforce_all_invariants (post-IA): %s", iter, conditionMessage(e))))
      state <- res_inv$state; params <- res_inv$params
      
      
      #state <- tryCatch(upper_only_local_adjusting(state, params, Y_work, precomp, model),
      #error = function(e) stop(sprintf("[iter %d] upper_LA (IA round, pass %d): %s", iter, n_l, conditionMessage(e))))
      #res_inv <- tryCatch(enforce_all_invariants(state, params, model),
      #error = function(e) stop(sprintf("[iter %d] enforce_inv (post-upperLA, pass %d): %s", iter, n_l, conditionMessage(e))))
      #state <- res_inv$state; params <- res_inv$params
      
      
        state <- tryCatch(lower_only_local_adjusting(state, params, Y_work, precomp, model),
          error = function(e) stop(sprintf("[iter %d] lower_LA (IA round, pass %d): %s", iter, n_l, conditionMessage(e))))
        res_inv <- tryCatch(enforce_all_invariants(state, params, model),
          error = function(e) stop(sprintf("[iter %d] enforce_inv (post-lowerLA, pass %d): %s", iter, n_l, conditionMessage(e))))
        state <- res_inv$state; params <- res_inv$params
      for (n_l in 1:N_launch) {
        pu_result <- tryCatch(param_update_fn(state, params, Y_work, precomp, model),
          error = function(e) stop(sprintf("[iter %d] PU (post-lowerLA, pass %d): %s", iter, n_l, conditionMessage(e))))
        state <- pu_result$state; params <- pu_result$params
      }
      
      
      if (!isTRUE(model$fixed_clusters)) {
        gs_result <- tryCatch(update_cluster_group_swap(state, params, Y_work, precomp, model),
          error = function(e) stop(sprintf("[iter %d] cluster_swap: %s", iter, conditionMessage(e))))
        state <- gs_result$state; params <- gs_result$params
        res_inv <- tryCatch(enforce_all_invariants(state, params, model),
          error = function(e) stop(sprintf("[iter %d] enforce_inv (post-swap): %s", iter, conditionMessage(e))))
        state <- res_inv$state; params <- res_inv$params
      }
      
    } else {
      
        state <- tryCatch(upper_only_local_adjusting(state, params, Y_work, precomp, model),
        error = function(e) stop(sprintf("[iter %d] upper_LA (non-IA, pass %d): %s", iter, n_l, conditionMessage(e))))
        res_inv <- tryCatch(enforce_all_invariants(state, params, model),
        error = function(e) stop(sprintf("[iter %d] enforce_inv (non-IA upper, pass %d): %s", iter, n_l, conditionMessage(e))))
        state <- res_inv$state; params <- res_inv$params
      
        state <- tryCatch(lower_only_local_adjusting(state, params, Y_work, precomp, model),
          error = function(e) stop(sprintf("[iter %d] lower_LA (non-IA, pass %d): %s", iter, n_l, conditionMessage(e))))
        res_inv <- tryCatch(enforce_all_invariants(state, params, model),
          error = function(e) stop(sprintf("[iter %d] enforce_inv (non-IA lower, pass %d): %s", iter, n_l, conditionMessage(e))))
        state <- res_inv$state; params <- res_inv$params
        for (n_l in 1:N_launch) {
        pu_result <- tryCatch(param_update_fn(state, params, Y_work, precomp, model),
          error = function(e) stop(sprintf("[iter %d] PU (non-IA lower, pass %d): %s", iter, n_l, conditionMessage(e))))
        state <- pu_result$state; params <- pu_result$params
      }
    }
    
    # ================================================================
    # Save samples + Conditional MAP tracking
    # ================================================================
    if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0) {
      save_idx <- save_idx + 1
      samples <- save_iteration(samples, save_idx, state, params, model)
      
      log_post <- compute_log_posterior(state, params, Y_work, precomp, model)
      samples$bf_score[save_idx] <- log_post
      
      K_key <- paste(state$K, collapse = ",")
      
      if (is.null(cmap_registry[[K_key]])) {
        cmap_registry[[K_key]] <- list(
          log_post = log_post, state = state, params = params,
          count = 1L, save_idx = save_idx
        )
      } else {
        cmap_registry[[K_key]]$count <- cmap_registry[[K_key]]$count + 1L
        if (log_post > cmap_registry[[K_key]]$log_post) {
          cmap_registry[[K_key]]$log_post <- log_post
          cmap_registry[[K_key]]$state <- state
          cmap_registry[[K_key]]$params <- params
          cmap_registry[[K_key]]$save_idx <- save_idx
        }
      }
    }
    
    # ================================================================
    # Progress report
    # ================================================================
    if (verbose > 0 && iter %% verbose == 0) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      K_vec <- paste(state$K, collapse = ",")
      cat(sprintf("\n--- Iter %d/%d [%.1f min] | K=[%s] | alpha=%.3f | outlier=%.3f ---\n",
                  iter, n_iter, as.numeric(elapsed), K_vec, 
                  params$alpha, 
                  ifelse(model$type == "continuous", params$theta_out, params$s_count)))
      
      for (cc in 1:model$C) {
        tau_str <- paste(state$tau_upper[[cc]], collapse = ",")
        cat(sprintf("  Cluster %d: Upper CPs = [%s]\n", cc, tau_str))
        K_c <- state$K[cc]
        atoms_str <- ""
        for (k in 1:min(K_c, length(params$atoms[[cc]]))) {
          a <- params$atoms[[cc]][[k]]
          if (!is.null(a)) {
            d1 <- ifelse(sign(a$gamma1) >= 0, "+", "-")
            d2 <- ifelse(sign(a$gamma2) >= 0, "+", "-")
            sb <- if (!is.null(a$shape_beta)) sprintf("%.2f", a$shape_beta) else "?"
            sg <- if (!is.null(a$shape_gamma)) sprintf("%.2f", a$shape_gamma) else "?"
            atoms_str <- paste0(atoms_str, sprintf("k%d=(%s,%s)[b=%s,g=%s] ", k, d1, d2, sb, sg))
          }
        }
        cat(sprintf("    Shapes: %s\n", atoms_str))
      }
      
      for (j in 1:model$J) {
        cp_j <- extract_changepoints(state$S_lower[j, ])
        cp_str <- paste(cp_j, collapse = ",")
        cat(sprintf("  S%d (C%d): Lower CPs = [%s]\n", j, state$cluster[j], cp_str))
      }
    }
    
    # ================================================================
    # Diagnostic callback
    # ================================================================
    if (!is.null(diagnostic_callback)) {
      diagnostic_callback(iter, state, params, Y_work, precomp, model)
    }
  }
  
  elapsed_total <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("\n=== MCMC completed in %.1f minutes ===\n", as.numeric(elapsed_total)))
    if (length(cmap_registry) > 0) {
    C <- model$C
    
    
    K_marginal <- vector("list", C)
    for (cc in 1:C) K_marginal[[cc]] <- list()
    
    for (key in names(cmap_registry)) {
      K_vec <- as.integer(strsplit(key, ",")[[1]])
      cnt <- cmap_registry[[key]]$count
      for (cc in 1:C) {
        k_str <- as.character(K_vec[cc])
        if (is.null(K_marginal[[cc]][[k_str]])) {
          K_marginal[[cc]][[k_str]] <- cnt
        } else {
          K_marginal[[cc]][[k_str]] <- K_marginal[[cc]][[k_str]] + cnt
        }
      }
    }
    
    
    K_hat <- integer(C)
    cat("\n--- Group-Based Independent Conditional MAP (Frequency-Based) ---\n")
    total_samples <- sum(sapply(cmap_registry, function(x) x$count))
    
    lambda0 <- if (!is.null(model$lambda0_K)) model$lambda0_K else 2.0
    
    for (cc in 1:C) {
      counts <- unlist(K_marginal[[cc]])
      k_vals <- as.integer(names(counts))
      
      
      K_hat[cc] <- k_vals[which.max(counts)]
      
      grp_members <- if (!is.null(model$group_members)) {
        paste(model$group_members[[cc]], collapse = ",")
      } else {
        paste(((cc-1) * (model$J %/% C) + 1):(cc * (model$J %/% C)), collapse = ",")
      }
      
    
      total_g <- sum(counts)
      post_probs <- counts / total_g  # P(K_g = k | Y)
      sorted_idx <- order(-counts)
      
      k_hat_count <- max(counts)
      k_hat_val   <- K_hat[cc]
      if (length(sorted_idx) >= 2) {
        k_2nd_val   <- k_vals[sorted_idx[2]]
        k_2nd_count <- counts[sorted_idx[2]]
        posterior_odds <- k_hat_count / max(k_2nd_count, 1)
        # Prior odds: Poisson(lambda0) ratio
        prior_odds <- dpois(k_hat_val, lambda0) / max(dpois(k_2nd_val, lambda0), 1e-300)
        bf_val <- posterior_odds / max(prior_odds, 1e-300)
        # Kass-Raftery interpretation
        bf_interp <- if (bf_val > 150) "very strong"
                     else if (bf_val > 20) "strong"
                     else if (bf_val > 3) "substantial"
                     else "anecdotal"
      } else {
        bf_val <- Inf; bf_interp <- "only one K observed"
      }
      
      cat(sprintf("  Group %d (series %s): K_hat=%d, P(K=%d|Y)=%.1f%%, BF=%.1f (%s)\n",
                  cc, grp_members, K_hat[cc], K_hat[cc],
                  100 * k_hat_count / total_g, bf_val, bf_interp))
      
    
      for (i in sorted_idx) {
        cat(sprintf("    K=%d: n=%d, P(K=%d|Y)=%.1f%%\n",
                    k_vals[i], counts[i], k_vals[i],
                    100 * counts[i] / total_g))
      }
    }
    
    
    K_hat_key <- paste(K_hat, collapse = ",")
    if (!is.null(cmap_registry[[K_hat_key]])) {
      best_entry <- cmap_registry[[K_hat_key]]
    } else {
      
      best_log_post <- -Inf
      best_entry <- NULL
      for (key in names(cmap_registry)) {
        K_vec <- as.integer(strsplit(key, ",")[[1]])
        if (all(K_vec == K_hat)) {
          if (cmap_registry[[key]]$log_post > best_log_post) {
            best_log_post <- cmap_registry[[key]]$log_post
            best_entry <- cmap_registry[[key]]
          }
        }
      }
      
      if (is.null(best_entry)) {
        all_lp <- sapply(cmap_registry, function(x) x$log_post)
        best_entry <- cmap_registry[[which.max(all_lp)]]
        cat(sprintf("  Note: K_hat=[%s] not in registry; using best log-posterior state.\n",
                    K_hat_key))
      }
    }
    
    best_state  <- best_entry$state
    best_params <- best_entry$params
    best_idx    <- best_entry$save_idx
    best_log_post <- best_entry$log_post
    
    cat(sprintf("K_hat (Group-independent) = [%s]\n", paste(K_hat, collapse = ",")))
    cat(sprintf("Best log-posterior: %.4f (iteration %d)\n", best_log_post, best_idx))
  } else {
    best_state <- state; best_params <- params; best_idx <- save_idx
    cat("Warning: No post-burnin samples collected.\n")
  }
  
  list(
    samples = samples,
    final_state = state,
    final_params = params,
    best_state = best_state,
    best_params = best_params,
    best_idx = best_idx,
    cmap_registry = cmap_registry,
    model = model,
    precomp = precomp
  )
}

init_storage <- function(n_save, model) {
  J <- model$J; T_len <- model$T_len; C <- model$C
  M <- if (!is.null(model$M)) model$M else 6

  list(
  
    alpha         = numeric(n_save),
    beta          = matrix(NA_real_, n_save, C),
    K             = matrix(NA_integer_, n_save, C),
    S_lower       = array(NA_integer_, dim = c(n_save, J, T_len)),
    S_upper       = array(NA_integer_, dim = c(n_save, C, T_len)),
    cluster       = matrix(NA_integer_, n_save, J),
    theta_out     = numeric(n_save),
    nu            = numeric(n_save),
    sigma2_gamma  = numeric(n_save),
    bf_score      = numeric(n_save),

    
    theta         = array(NA_real_, dim = c(n_save, J, M + 1)),
    v             = matrix(NA_real_, n_save, J),
    gamma_store   = vector("list", n_save),
    atoms_store   = vector("list", n_save),
    phi           = array(NA_real_, dim = c(n_save, J, T_len)),
    xi_store      = array(NA_integer_, dim = c(n_save, J, T_len)),

    n_saved       = 0
  )
}

#' Save one iteration to storage (EXTENDED)
save_iteration <- function(samples, idx, state, params, model) {
  J <- model$J; C <- model$C; T_len <- model$T_len

  
  samples$alpha[idx]       <- params$alpha
  samples$beta[idx, ]      <- params$beta
  samples$K[idx, ]         <- state$K
  samples$S_lower[idx,,]   <- state$S_lower
  samples$S_upper[idx,,]   <- state$S_upper
  samples$cluster[idx, ]   <- state$cluster
  samples$theta_out[idx]   <- if (!is.null(params$theta_out)) params$theta_out else NA
  samples$nu[idx]          <- if (!is.null(params$nu)) params$nu else NA
  samples$sigma2_gamma[idx]<- if (!is.null(params$sigma2_gamma)) params$sigma2_gamma else NA
  samples$n_saved          <- idx

  
  if (!is.null(samples$theta) && !is.null(params$theta)) {
    M_store <- dim(samples$theta)[3]
    M_actual <- ncol(params$theta)
    m_use <- min(M_store, M_actual)
    samples$theta[idx, , 1:m_use] <- params$theta[, 1:m_use]
  }

  
  if (!is.null(samples$v)) {
    for (j in 1:J) {
      samples$v[idx, j] <- if (!is.null(params$v[[j]])) params$v[[j]] else 1
    }
  }

  
  if (!is.null(samples$gamma_store)) {
    g_list <- vector("list", J)
    for (j in 1:J) {
      g_list[[j]] <- if (!is.null(params$gamma[[j]])) as.numeric(params$gamma[[j]]) else 0
    }
    samples$gamma_store[[idx]] <- g_list
  }

  
  if (!is.null(samples$atoms_store)) {
    a_list <- vector("list", C)
    for (cc in 1:C) {
      K_c <- state$K[cc]
      a_list[[cc]] <- vector("list", K_c)
      for (k in 1:K_c) {
        atom <- if (k <= length(params$atoms[[cc]])) params$atoms[[cc]][[k]] else NULL
        if (!is.null(atom)) {
          a_list[[cc]][[k]] <- list(
            gamma1     = as.numeric(atom$gamma1),
            gamma2     = as.numeric(atom$gamma2),
            shape_beta = as.numeric(atom$shape_beta),
            shape_gamma= as.numeric(atom$shape_gamma)
          )
        } else {
          a_list[[cc]][[k]] <- list(gamma1=1, gamma2=1, shape_beta=1, shape_gamma=0)
        }
      }
    }
    samples$atoms_store[[idx]] <- a_list
  }

  
  if (!is.null(samples$phi) && !is.null(params$phi)) {
    samples$phi[idx,,] <- params$phi
  }
  if (!is.null(samples$xi_store) && !is.null(params$xi)) {
    samples$xi_store[idx,,] <- params$xi
  }

  samples
}



cat("11_mcmc_main.R loaded (v3: frequency-based K_hat + BF reporting + extended storage).\n")
