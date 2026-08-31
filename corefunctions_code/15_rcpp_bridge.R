###############################################################################
# 15_rcpp_bridge.R
# Rcpp acceleration bridge 

.hiercpd_rcpp_enabled <- FALSE


.hiercpd_ensure_datx <- function(precomp) invisible(FALSE)

.hiercpd_has_compiled_symbol <- function(sym) {
  exists(sym, mode = "function", inherits = TRUE)
}

.hiercpd_find_cpp_path <- function() {
  candidate_dirs <- character(0)
  for (i in rev(seq_len(sys.nframe()))) {
    fr <- sys.frame(i)
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      candidate_dirs <- c(candidate_dirs, dirname(normalizePath(get("ofile", envir = fr), mustWork = FALSE)))
    }
  }
  candidate_dirs <- unique(c(candidate_dirs, getwd()))
  for (dd in candidate_dirs) {
    cpp_path <- file.path(dd, "hiercpd_core.cpp")
    if (file.exists(cpp_path)) return(cpp_path)
  }
  NULL
}

.hiercpd_try_compile_cpp <- function() {
  if (!requireNamespace("Rcpp", quietly = TRUE) ||
      !requireNamespace("RcppArmadillo", quietly = TRUE)) {
    message("15_rcpp_bridge.R: Rcpp/RcppArmadillo not available; using R reference path.")
    return(FALSE)
  }

  cpp_path <- .hiercpd_find_cpp_path()
  if (is.null(cpp_path) || !file.exists(cpp_path)) {
    message("15_rcpp_bridge.R: hiercpd_core.cpp not found; using R reference path.")
    return(FALSE)
  }

  ok <- tryCatch({
    Rcpp::sourceCpp(cpp_path)
    TRUE
  }, error = function(e) {
    message("15_rcpp_bridge.R: sourceCpp failed; using R reference path.\n  -> ", conditionMessage(e))
    FALSE
  })

  if (ok) {
    compat_ok <- tryCatch({
      test_atom <- list(gamma1 = 1, gamma2 = 1, shape_beta = 1.0, shape_gamma = 0.5)
      test_theta <- c(0, 1.0)            # intercept + 1 basis
      test_x <- c(0, 0.5, 1.0)
      test_D <- array(0, dim = c(2, 2, 3))
      set_precomp_datx(test_D)           # push cube into module cache (no-cube kernels)
      res <- eval_shape_at_times_cpp(test_theta, test_atom, 1:3, test_x)
      
      abs(res[2] - 0.75) < 0.01
    }, error = function(e) FALSE)
    if (!compat_ok) {
      message("15_rcpp_bridge.R: C++ does not support decoupled atoms + cube cache (set_precomp_datx).")
      message("  -> Falling back to pure-R reference path. Rebuild hiercpd_core.cpp for Rcpp acceleration.")
      ok <- FALSE
    }
  }

  needed <- c(
    "set_precomp_datx",
    "eval_shape_at_times_cpp",
    "compute_f_all_cpp",
    "compute_mu_all_cpp",
    "compute_mu_and_variance_cpp",
    "log_lik_series_cpp",
    "log_lik_interval_cpp",
    "lower_la_pair_cpp",
    "upper_la_pair_cpp",
    "split_scan_upper_cpp",
    "split_seq_scan_upper_cpp",
    "profile_atom_scores_cpp",
    "merge_eval_cpp",
    "compute_atom_loglik_cpp",
    "compute_H_at_times_cpp",
    "rtrunc_invgamma_lower1_cpp"
  )

  if (!ok || !all(vapply(needed, .hiercpd_has_compiled_symbol, logical(1)))) {
    message("15_rcpp_bridge.R: compiled symbols incomplete; using R reference path.")
    return(FALSE)
  }

  TRUE
}

if (.hiercpd_try_compile_cpp()) {
  .hiercpd_rcpp_enabled <- TRUE

  # -------------------------------------------------------------------------
  # [PERF] D-cube cache management.
  #   set_precomp_datx() copies the cube into the compiled DLL's G_DATX exactly
  #   once per distinct cube. We track a cheap O(1) fingerprint (dims + corner
  #   values). x_global is deterministic in T for this model, so equal dims =>
  #   equal D; the corner values additionally guard against pathological reuse
  #   across data sets that happen to share dimensions.
  # -------------------------------------------------------------------------
  .hiercpd_datx_state <- new.env(parent = emptyenv())
  .hiercpd_datx_state$set <- FALSE
  .hiercpd_datx_state$fp  <- NULL

  .hiercpd_ensure_datx <- function(precomp) {
    D <- precomp$D_at_x
    if (is.null(D)) return(invisible(FALSE))
    n <- length(D)
    fp <- c(dim(D), D[1L], D[n])         # O(1) fingerprint
    if (!.hiercpd_datx_state$set || !identical(.hiercpd_datx_state$fp, fp)) {
      set_precomp_datx(D)
      .hiercpd_datx_state$set <- TRUE
      .hiercpd_datx_state$fp  <- fp
    }
    invisible(TRUE)
  }

  eval_shape_at_times <- function(theta_j, atom, t_indices, precomp) {
    if (length(t_indices) == 0) return(numeric(0))
    if (is.null(atom)) return(rep(0, length(t_indices)))
    .hiercpd_ensure_datx(precomp)
    as.numeric(eval_shape_at_times_cpp(
      as.numeric(theta_j), atom, as.integer(t_indices),
      as.numeric(precomp$x_global)
    ))
  }

  compute_f_all_timepoints <- function(j, state, params, precomp, model) {
    c_j <- state$cluster[j]
    K_c <- state$K[c_j]
    .hiercpd_ensure_datx(precomp)
    as.numeric(compute_f_all_cpp(
      as.integer(state$S_lower[j, ]),
      K_c,
      as.numeric(params$theta[j, ]),
      params$atoms[[c_j]],
      as.numeric(precomp$x_global)
    ))
  }

  compute_mu_all <- function(j, state, params, precomp, model) {
    c_j <- state$cluster[j]
    K_c <- state$K[c_j]
    .hiercpd_ensure_datx(precomp)
    as.numeric(compute_mu_all_cpp(
      params$alpha,
      params$beta[c_j],
      as.integer(state$S_lower[j, ]),
      K_c,
      as.numeric(params$theta[j, ]),
      as.numeric(params$gamma[[j]]),
      params$atoms[[c_j]],
      as.numeric(precomp$x_global)
    ))
  }

  log_lik_series_continuous <- function(j, Y, state, params, precomp, model) {
    c_j <- state$cluster[j]
    K_c <- state$K[c_j]
    .hiercpd_ensure_datx(precomp)
    ll <- log_lik_series_cpp(
      as.integer(state$S_lower[j, ]),
      K_c,
      as.numeric(params$theta[j, ]),
      as.numeric(params$gamma[[j]]),
      params$atoms[[c_j]],
      as.numeric(Y[j, ]),
      params$alpha,
      params$beta[c_j],
      as.numeric(params$v[[j]]),  # [v_j scalar] length-1 vector
      as.numeric(params$phi[j, ]),
      as.integer(params$xi[j, ]),
      as.numeric(precomp$x_global)
    )
    if (!is.finite(ll)) ll <- -1e300
    ll
  }

  log_lik_interval_continuous <- function(j, t_indices, state_seq, Y, state,
                                          params, precomp, model,
                                          var_state_override = NULL) {
    if (length(t_indices) == 0) return(0)
    c_j <- state$cluster[j]
    .hiercpd_ensure_datx(precomp)
    ll <- log_lik_interval_cpp(
      as.integer(t_indices),
      as.integer(state_seq),
      as.numeric(Y[j, ]),
      as.numeric(params$theta[j, ]),
      as.numeric(params$gamma[[j]]),
      params$atoms[[c_j]],
      params$alpha,
      params$beta[c_j],
      as.numeric(params$v[[j]]),  # [v_j scalar] length-1 vector
      as.numeric(params$phi[j, ]),
      as.integer(params$xi[j, ]),
      as.numeric(precomp$x_global),
      if (is.null(var_state_override)) 0L else as.integer(var_state_override)
    )
    if (!is.finite(ll)) ll <- -1e300
    ll
  }

  compute_atom_loglik <- function(cc, k, atom, j_in_c, t_in_k_upper,
                                  state, params, Y, precomp, model) {
    if (length(j_in_c) == 0 || length(t_in_k_upper) == 0) return(0)
    .hiercpd_ensure_datx(precomp)

    theta_mat <- params$theta[j_in_c, , drop = FALSE]
    gamma_vec <- vapply(j_in_c, function(jj) get_state_intercept(jj, k, params), numeric(1))
    v_vec <- vapply(j_in_c, function(jj) get_state_variance(jj, k, params), numeric(1))

    ll <- compute_atom_loglik_cpp(
      atom,
      theta_mat,
      gamma_vec,
      as.integer(t_in_k_upper),
      state$S_lower[j_in_c, , drop = FALSE],
      as.integer(k),
      Y[j_in_c, , drop = FALSE],
      params$alpha,
      params$beta[cc],
      v_vec,
      params$phi[j_in_c, , drop = FALSE],
      params$xi[j_in_c, , drop = FALSE],
      as.numeric(precomp$x_global)
    )
    if (!is.finite(ll)) ll <- -1e300
    ll
  }

  # [PERF] Normalized shape profile H_theta(x_t) = theta^T D(x_t) theta / max(SZ, eps).
  #   Overrides the pure-R reference in 02_shape_function.R; reads the cached
  #   G_DATX cube. Used by .collect_state_series_data() in 08 to build the
  #   (shape_beta, shape_gamma) design columns without an R per-time loop.
  compute_H_at_times <- function(theta_j, t_indices, precomp) {
    if (length(t_indices) == 0) return(numeric(0))
    .hiercpd_ensure_datx(precomp)
    as.numeric(compute_H_at_times_cpp(as.numeric(theta_j), as.integer(t_indices)))
  }

  # [PERF] Batch exact-inverse-CDF truncated InvGamma on [1, Inf).
  #   Overrides the pure-R reference in 00_utils.R. Used by update_phi_continuous
  #   in 08 to draw all outlier-point phi at once (no per-(j,t) rejection loop).
  rtrunc_invgamma_lower1 <- function(shape, scale) {
    as.numeric(rtrunc_invgamma_lower1_cpp(as.numeric(shape), as.numeric(scale)))
  }

  # [sigma^2_c] The C++ independent-prior intercept sampler is no longer valid
  # (current model uses a per-cluster sigma^2_c conjugate draw in 08).
  # Returning NULL forces the R conjugate path.
  sample_state_intercepts_continuous_fast <- function(j, state, params, Y, precomp, model, sigma2_gamma) {
    return(NULL)
  }

  # Lower / Upper Local Adjusting (06_local_adjusting.R) now call
  # lower_la_pair_cpp / upper_la_pair_cpp DIRECTLY (guarded internally by
  # .hiercpd_rcpp_enabled + exists()), using the Normal likelihood with v_j in
  # the variance component. They are NOT overridden here; 06 invokes
  # .hiercpd_ensure_datx(precomp) before each direct kernel call.
  #
  # update_atoms() B2 (08_param_update_continuous.R) likewise calls
  # compute_atoms_loglik_batch_cpp DIRECTLY when present (one call evaluates all
  # N pool atoms for a state vs N per-atom calls). It is intentionally NOT in the
  # `needed` gate above: an older C++ without it must still keep all other
  # acceleration, so 08 guards it with exists() and falls back to the per-atom
  # compute_atom_loglik path.

  message("15_rcpp_bridge.R: Rcpp acceleration ACTIVE (decoupled kernels, cached D-cube, DIS upper/lower LA).")
} else {
  message("15_rcpp_bridge.R: running in pure-R reference mode.")
}
