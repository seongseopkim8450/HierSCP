# Interval Adjusting: split pass and merge pass on the cluster-level regimes

# One Interval Adjusting round: split pass, then merge pass, per cluster
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

# Categorical draw that stops on invalid weights
sample_categorical_log_strict <- function(log_weights) {
  if (!length(log_weights) || anyNA(log_weights) || any(log_weights == Inf)) {
    stop("IA categorical draw received invalid log-weights.", call. = FALSE)
  }
  if (!any(is.finite(log_weights))) {
    stop("IA categorical draw has zero finite mass.", call. = FALSE)
  }
  m <- max(log_weights[is.finite(log_weights)])
  w <- exp(log_weights - m)
  w[!is.finite(w)] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    stop("IA categorical normalization failed.", call. = FALSE)
  }
  sample.int(length(w), size = 1L, prob = w / sw)
}

# Prefix sums of the sufficient statistics over 1:T, computed once per merge pass
.merge_integrated_cache <- function(cc, state, params, Y, precomp, model) {
  jj <- which(state$cluster == cc)
  if (!length(jj)) stop("Integrated merge requires a nonempty cluster.", call. = FALSE)
  list(cc = cc, T_len = model$T_len, j_in_c = jj,
       prior = .split_mag_prior(params), signs = .split_sign_pairs(),
       prefix = .split_prefix_stats(cc, 1L, model$T_len + 1L, jj,
                                    state, params, Y, precomp))
}

.merge_empty_stats <- function() {
  list(Sgg = numeric(4L), Sgb = numeric(4L), Sbb = numeric(4L),
       Sgr = numeric(4L), Sbr = numeric(4L), Srr = 0, SlogW = 0, nobs = 0L)
}

.merge_add_stats <- function(a, b) {
  Map(`+`, a, b)
}

# Sufficient statistics of the block first:last
.merge_interval_stats <- function(cache, first, last) {
  if (first > last) return(.merge_empty_stats())
  if (first < 1L || last > cache$T_len) {
    stop("Integrated merge interval is outside the data range.", call. = FALSE)
  }
  out <- .merge_empty_stats()
  for (nm in c("Sgg", "Sgb", "Sbb", "Sgr", "Sbr")) {
    out[[nm]] <- cache$prefix[[nm]][, last] -
      if (first > 1L) cache$prefix[[nm]][, first - 1L] else numeric(4L)
  }
  for (nm in c("Srr", "SlogW")) {
    out[[nm]] <- cache$prefix[[nm]][last] -
      if (first > 1L) cache$prefix[[nm]][first - 1L] else 0
  }
  out$nobs <- (last - first + 1L) * length(cache$j_in_c)
  out
}

# Log evidence of a block with its atom integrated out under F0; theta = (g, b)
.merge_integrated_evidence <- function(stats, prior, signs = .split_sign_pairs()) {
  m0 <- c(prior$m_g, prior$m_b)
  P0 <- diag(c(prior$p_g, prior$p_b))
  mu <- matrix(0, 4L, 2L)
  Sigma <- vector("list", 4L)
  lev <- numeric(4L)
  if (stats$nobs == 0L) {
    mu[,] <- rep(m0, each = 4L)
    Sigma <- rep(list(diag(1 / c(prior$p_g, prior$p_b))), 4L)
    # empty block: evidence 1
    return(list(log_m = 0, log_sign_evidence = lev, mu = mu,
                Sigma = Sigma, signs = signs, nobs = 0L))
  }
  if (any(!is.finite(unlist(stats, use.names = FALSE)))) {
    stop("Nonfinite sufficient statistics in integrated merge.", call. = FALSE)
  }
  logC <- -0.5 * stats$nobs * log(2 * pi) + 0.5 * stats$SlogW
  q0 <- sum(m0 * (P0 %*% m0))
  logdet0 <- log(prior$p_g) + log(prior$p_b)
  for (di in seq_len(4L)) {
    P <- P0 + matrix(c(stats$Sgg[di], stats$Sgb[di],
                       stats$Sgb[di], stats$Sbb[di]), 2L, 2L)
    R <- chol(P)
    S <- chol2inv(R)
    rhs <- as.numeric(P0 %*% m0) + c(stats$Sgr[di], stats$Sbr[di])
    md <- as.numeric(backsolve(R, forwardsolve(t(R), rhs)))
    mu[di, ] <- md
    Sigma[[di]] <- S
    lev[di] <- logC + 0.5 * (logdet0 - 2 * sum(log(diag(R)))) -
      0.5 * (stats$Srr + q0 - sum(rhs * md)) +
      .split_biv_pos_logmass(md, S) - prior$log_Z0
  }
  lm <- .ia_logsumexp(lev) - log(4)
  if (!is.finite(lm)) stop("Integrated merge evidence is not finite.", call. = FALSE)
  list(log_m = lm, log_sign_evidence = lev, mu = mu, Sigma = Sigma,
       signs = signs, nobs = stats$nobs)
}

# Merge weights of regime k: structural weights times m(B + I_k) / m(B), where B is
# the data already allocated to the candidate pool entry
.merge_integrated_evaluate <- function(cc, k, state, params, Y, precomp, model,
                                       cache = NULL) {
  K <- state$K[cc]
  if (length(k) != 1L || is.na(k) || k != as.integer(k) || k < 1L || k > K) {
    stop("Invalid current regime in integrated merge.", call. = FALSE)
  }
  z <- params$z_state[[cc]]
  pool <- params$atom_pool[[cc]]
  if (length(z) != K || !length(pool) || anyNA(z) ||
      any(z != as.integer(z)) || any(z < 1L | z > length(pool))) {
    stop("Atom-pool allocation is inconsistent before integrated merge.", call. = FALSE)
  }
  tau <- state$tau_upper[[cc]]
  if (length(tau) != K || anyNA(tau) || tau[1L] != 1L ||
      any(diff(c(tau, model$T_len + 1L)) <= 0)) {
    stop("Upper boundaries are inconsistent before integrated merge.", call. = FALSE)
  }
  if (is.null(cache)) cache <- .merge_integrated_cache(cc, state, params, Y, precomp, model)
  if (cache$cc != cc || cache$T_len != model$T_len ||
      !identical(cache$j_in_c, which(state$cluster == cc))) {
    stop("Integrated merge cache belongs to a different cluster/data range.", call. = FALSE)
  }
  ends <- c(tau[-1L] - 1L, model$T_len)
  block_stats <- lapply(seq_len(K), function(kk) {
    .merge_interval_stats(cache, tau[kk], ends[kk])
  })
  candidate_pool <- c(if (k > 1L) z[k - 1L] else NA_integer_, z[k],
                      if (k < K) z[k + 1L] else NA_integer_)
  touched <- unique(candidate_pool[!is.na(candidate_pool)])
  base <- augmented <- setNames(vector("list", length(touched)), as.character(touched))
  for (n in touched) {
    st <- .merge_empty_stats()
    # every other regime with pool entry n
    for (kk in which(z == n & seq_len(K) != k)) {
      st <- .merge_add_stats(st, block_stats[[kk]])
    }
    key <- as.character(n)
    base[[key]] <- .merge_integrated_evidence(st, cache$prior, cache$signs)
    augmented[[key]] <- .merge_integrated_evidence(
      .merge_add_stats(st, block_stats[[k]]), cache$prior, cache$signs)
  }
  pred <- rep(-Inf, 3L)
  for (d in which(!is.na(candidate_pool))) {
    key <- as.character(candidate_pool[d])
    pred[d] <- augmented[[key]]$log_m - base[[key]]$log_m
  }
  ps <- params$pi_star[[cc]]
  if (length(ps) != K || any(!is.finite(ps))) {
    stop("Upper persistence vector is inconsistent before integrated merge.", call. = FALSE)
  }
  pk <- min(max(ps[k], 1e-10), 1 - 1e-10)
  pm <- if (k > 1L) min(max(ps[k - 1L], 1e-10), 1 - 1e-10) else 0.5
  lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0
  if (!is.finite(lambda0) || lambda0 <= 0) stop("Invalid lambda0_K.", call. = FALSE)
  nk <- ends[k] - tau[k] + 1L
  bonus <- log(K) - log(lambda0)
  # structural weights relative to p_K(K): merge left, keep, merge right
  lp <- c(nk * log(pm) + log1p(-pm) + bonus,
          max(nk - model$m_min, 0) * log(pk) + log1p(-pk),
          log1p(-pk) + bonus)
  Kmin <- if (!is.null(model$K_min)) model$K_min else 1L
  allowed <- c(k > 1L && K > Kmin, TRUE, k < K && K > Kmin)
  # k - 1 and k + 1 would become adjacent with the same pool entry
  if (k > 1L && k < K && z[k - 1L] == z[k + 1L]) allowed[c(1L, 3L)] <- FALSE
  lp[!allowed] <- -Inf
  lw <- lp + pred
  if (anyNA(lw) || any(lw == Inf) || !is.finite(lw[2L])) {
    stop("Invalid integrated merge action weights.", call. = FALSE)
  }
  list(log_predictive = pred, log_prior = lp, log_weights = lw,
       allowed = allowed, touched_pool = touched, candidate_pool = candidate_pool,
       base = base, augmented = augmented)
}

# Sign pair from the four evidences, then (g, b) from the truncated bivariate normal
.draw_atom_from_integrated_merge <- function(evidence) {
  di <- sample_categorical_log_strict(evidence$log_sign_evidence)
  gb <- rbivnorm_pos_exact(evidence$mu[di, ], evidence$Sigma[[di]])
  canonicalize_atom(list(gamma1 = evidence$signs[di, 1L],
                         gamma2 = evidence$signs[di, 2L],
                         shape_beta = gb[2L], shape_gamma = gb[1L]))
}

# Log predictive ratios of the three merge actions
.merge_eval <- function(cc, k, I_k, j_in_c, state, params, Y, precomp, model) {
  .merge_integrated_evaluate(cc, k, state, params, Y, precomp, model)$log_predictive
}

# Merge pass over k = 2, ..., K_c; after a merge the relabeled regime k is visited again
merge_pass_once <- function(cc, state, params, Y, precomp, model) {
  if (state$K[cc] < 2L || !any(state$cluster == cc)) {
    return(list(state = state, params = params))
  }
  cache <- .merge_integrated_cache(cc, state, params, Y, precomp, model)
  k <- 2L
  max_steps <- 2L * state$K[cc]
  for (step in seq_len(max_steps)) {
    if (k > state$K[cc]) break
    terminal <- k == state$K[cc]
    res <- merge_step(k, cc, state, params, Y, precomp, model, cache = cache)
    state <- res$state; params <- res$params
    if (terminal) break
    if (!res$merged) k <- k + 1L
  }
  list(state = state, params = params)
}

# Split pass over the regimes present at its start
split_pass_once <- function(cc, state, params, Y, precomp, model) {
  tau_start <- state$tau_upper[[cc]]
  if (length(tau_start) < 1) {
    return(list(state = state, params = params))
  }

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

# Sign pairs (delta1, delta2) in a fixed order
.split_sign_pairs <- function() {
  matrix(c( 1L,  1L,
             1L, -1L,
            -1L,  1L,
            -1L, -1L), ncol = 2L, byrow = TRUE)
}

# Magnitude prior as a Gaussian kernel for theta = (g, b): means, precisions, log Z0
.split_mag_prior <- function(params) {
  h <- params$hyper
  a_b <- if (!is.null(h$a_shape_beta))  h$a_shape_beta  else 2.0
  r_b <- if (!is.null(h$b_shape_beta))  h$b_shape_beta  else 1.5
  a_g <- if (!is.null(h$a_shape_gamma)) h$a_shape_gamma else 1.5
  r_g <- if (!is.null(h$b_shape_gamma)) h$b_shape_gamma else 2.0
  if (!is.finite(a_b) || a_b <= 1 || !is.finite(a_g) || a_g <= 1) {
    stop("Integrated split requires a_shape_beta>1 and a_shape_gamma>1.", call. = FALSE)
  }
  if (!is.finite(r_b) || r_b <= 0 || !is.finite(r_g) || r_g <= 0) {
    stop("Integrated split requires positive magnitude rates.", call. = FALSE)
  }
  m_g <- (a_g - 1) / r_g
  m_b <- (a_b - 1) / r_b
  p_g <- r_g^2 / (a_g - 1)
  p_b <- r_b^2 / (a_b - 1)
  list(a_g = a_g, r_g = r_g, a_b = a_b, r_b = r_b,
       m_g = m_g, m_b = m_b, p_g = p_g, p_b = p_b,
       log_Z0 = pnorm(sqrt(a_g - 1), log.p = TRUE) +
                pnorm(sqrt(a_b - 1), log.p = TRUE))
}

.split_biv_pos_logmass <- function(mu, Sig) {
  if (!exists("biv_pos_logmass", mode = "function", inherits = TRUE)) {
    stop("biv_pos_logmass helper is not loaded; source 04_sampling.R first.",
         call. = FALSE)
  }
  biv_pos_logmass(mu, Sig)
}

.split_biv_pos_mass <- function(mu, Sig) {
  lp <- .split_biv_pos_logmass(mu, Sig)
  if (!is.finite(lp)) return(0)
  if (lp < log(.Machine$double.xmin)) return(0)
  exp(lp)
}

.ia_logsumexp <- function(x) {
  if (anyNA(x) || any(x == Inf)) stop("IA: invalid log-weight", call. = FALSE)
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Prefix sums of the statistics A, s and S for each sign pair, from t = tau_k
.split_prefix_stats_R <- function(cc, tau_k, g_right, j_in_c,
                                  state, params, Y, precomp) {
  tt <- seq.int(tau_k, g_right - 1L)
  L <- length(tt)
  if (L < 1L) stop("Integrated split has empty prospective-left range.", call. = FALSE)
  signs <- .split_sign_pairs()
  Sgg <- Sgb <- Sbb <- Sgr <- Sbr <- matrix(0, 4L, L)
  Srr <- SlogW <- numeric(L)

  for (j in j_in_c) {
    gamma_j <- params$gamma[[j]]
    gamma_j <- if (is.null(gamma_j) || !length(gamma_j) || !is.finite(gamma_j[1L])) 0 else gamma_j[1L]
    r <- as.numeric(Y[j, tt] - params$alpha - params$beta[cc] - gamma_j)
    sig2 <- as.numeric(get_effective_variance_all(j, state, params)[tt])
    if (any(!is.finite(sig2)) || any(sig2 <= 0) || any(!is.finite(r))) {
      stop("Invalid residual/variance in integrated split.", call. = FALSE)
    }
    w <- 1 / sig2
    x <- as.numeric(precomp$x_global[tt])
    H <- as.numeric(compute_H_at_times(params$theta[j, ], tt, precomp))
    if (length(H) != L || any(!is.finite(H))) stop("Invalid H_j values in integrated split.", call. = FALSE)

    Srr <- Srr + w * r^2
    SlogW <- SlogW + log(w)
    for (di in seq_len(4L)) {
      d1 <- signs[di, 1L]; d2 <- signs[di, 2L]
      zg <- d1 * x
      zb <- (d1 + d2) / 2 * x - d2 * H
      Sgg[di, ] <- Sgg[di, ] + w * zg^2
      Sgb[di, ] <- Sgb[di, ] + w * zg * zb
      Sbb[di, ] <- Sbb[di, ] + w * zb^2
      Sgr[di, ] <- Sgr[di, ] + w * zg * r
      Sbr[di, ] <- Sbr[di, ] + w * zb * r
    }
  }

  for (di in seq_len(4L)) {
    Sgg[di, ] <- cumsum(Sgg[di, ]); Sgb[di, ] <- cumsum(Sgb[di, ])
    Sbb[di, ] <- cumsum(Sbb[di, ]); Sgr[di, ] <- cumsum(Sgr[di, ])
    Sbr[di, ] <- cumsum(Sbr[di, ])
  }
  list(Sgg = Sgg, Sgb = Sgb, Sbb = Sbb, Sgr = Sgr, Sbr = Sbr,
       Srr = cumsum(Srr), SlogW = cumsum(SlogW))
}

# Same, with the compiled kernel when available
.split_prefix_stats <- function(cc, tau_k, g_right, j_in_c,
                                state, params, Y, precomp) {
  use_cpp <- isTRUE(get0(".hiercpd_rcpp_enabled", ifnotfound = FALSE)) &&
             exists("split_integrated_prefix_stats_cpp", mode = "function")
  if (use_cpp) {
    if (exists(".hiercpd_ensure_datx", mode = "function")) .hiercpd_ensure_datx(precomp)
    gamma_vec <- vapply(j_in_c, function(jj) {
      z <- params$gamma[[jj]]
      if (is.null(z) || !length(z) || !is.finite(z[1L])) 0 else as.numeric(z[1L])
    }, numeric(1))
    v_vec <- vapply(j_in_c, function(jj) {
      z <- params$v[[jj]]
      if (is.null(z) || !length(z) || !is.finite(z[1L]) || z[1L] <= 0) 1 else as.numeric(z[1L])
    }, numeric(1))
    xi_sub <- params$xi[j_in_c, , drop = FALSE]
    xi_sub[is.na(xi_sub)] <- 0L
    storage.mode(xi_sub) <- "integer"
    ans <- tryCatch(split_integrated_prefix_stats_cpp(
      params$theta[j_in_c, , drop = FALSE], gamma_vec,
      Y[j_in_c, , drop = FALSE], params$alpha, params$beta[cc], v_vec,
      params$phi[j_in_c, , drop = FALSE], xi_sub,
      as.integer(tau_k), as.integer(g_right), as.numeric(precomp$x_global)
    ), error = function(e) NULL)
    if (!is.null(ans)) return(ans)
  }
  .split_prefix_stats_R(cc, tau_k, g_right, j_in_c, state, params, Y, precomp)
}

# Evidence m_delta(h) and posterior (mu, Sigma) for every candidate h and sign pair
.build_integrated_split_cache <- function(cc, tau_k, G_split, j_in_c,
                                          state, params, Y, precomp) {
  if (!length(G_split)) stop("Empty split candidate set.", call. = FALSE)
  prior <- .split_mag_prior(params)
  stats <- .split_prefix_stats(cc, tau_k, max(G_split), j_in_c,
                               state, params, Y, precomp)
  signs <- .split_sign_pairs()
  nc <- length(G_split)
  log_ev <- matrix(-Inf, nc, 4L)
  mu <- array(NA_real_, dim = c(nc, 4L, 2L))
  s11 <- s12 <- s22 <- matrix(NA_real_, nc, 4L)
  log_new <- rep(-Inf, nc)
  quad0 <- prior$p_g * prior$m_g^2 + prior$p_b * prior$m_b^2

  for (ci in seq_len(nc)) {
    sl <- as.integer(G_split[ci] - tau_k)
    if (sl < 1L || sl > length(stats$Srr)) next
    nobs <- length(j_in_c) * sl
    log_like_const <- -0.5 * nobs * log(2 * pi) + 0.5 * stats$SlogW[sl]

    for (di in seq_len(4L)) {
      P11 <- prior$p_g + stats$Sgg[di, sl]
      P12 <- stats$Sgb[di, sl]
      P22 <- prior$p_b + stats$Sbb[di, sl]
      detP <- P11 * P22 - P12^2
      if (!is.finite(detP) || detP <= 0) next

      S11 <- P22 / detP; S12 <- -P12 / detP; S22 <- P11 / detP
      rhs_g <- prior$p_g * prior$m_g + stats$Sgr[di, sl]
      rhs_b <- prior$p_b * prior$m_b + stats$Sbr[di, sl]
      mu_g <- S11 * rhs_g + S12 * rhs_b
      mu_b <- S12 * rhs_g + S22 * rhs_b
      Sig <- matrix(c(S11, S12, S12, S22), 2L, 2L)
      log_ppos <- .split_biv_pos_logmass(c(mu_g, mu_b), Sig)
      if (is.na(log_ppos) || log_ppos == Inf) {
        stop("Integrated-split orthant log-probability is invalid.",
             call. = FALSE)
      }

      quad_post <- rhs_g * mu_g + rhs_b * mu_b
      log_ev[ci, di] <- log_like_const +
        0.5 * (log(prior$p_g) + log(prior$p_b) - log(detP)) +
        0.5 * (quad_post - quad0) - 0.5 * stats$Srr[sl] -
        prior$log_Z0 + log_ppos
      mu[ci, di, ] <- c(mu_g, mu_b)
      s11[ci, di] <- S11; s12[ci, di] <- S12; s22[ci, di] <- S22
    }
    log_new[ci] <- .ia_logsumexp(log_ev[ci, ]) - log(4)
  }

  list(G_split = as.integer(G_split), signs = signs,
       log_sign_evidence = log_ev, log_new = log_new,
       mu = mu, s11 = s11, s12 = s12, s22 = s22)
}

# Atom for a confirmed split at g
.draw_atom_from_integrated_split <- function(cache, g) {
  ci <- match(as.integer(g), cache$G_split)
  if (is.na(ci)) stop("Confirmed split location not found in integrated cache.", call. = FALSE)
  lev <- cache$log_sign_evidence[ci, ]
  di <- sample_categorical_log_strict(lev)
  m <- as.numeric(cache$mu[ci, di, ])
  S <- matrix(c(cache$s11[ci, di], cache$s12[ci, di],
                cache$s12[ci, di], cache$s22[ci, di]), 2L, 2L)
  if (!exists("rbivnorm_pos_exact", mode = "function", inherits = TRUE)) {
    stop("rbivnorm_pos_exact helper is not loaded; source 04_sampling.R first.",
         call. = FALSE)
  }
  gb <- rbivnorm_pos_exact(m, S)
  canonicalize_atom(list(
    gamma1 = cache$signs[di, 1L], gamma2 = cache$signs[di, 2L],
    shape_beta = gb[2L], shape_gamma = gb[1L]
  ))
}

# Alias of merge_pass_once()
same_shape_merge_sweep <- function(cc, state, params, Y, precomp, model) {
  merge_pass_once(cc, state, params, Y, precomp, model)
}

# Merge step for regime k: the candidate atoms are integrated out, an action is drawn,
# and every touched pool entry is redrawn given its data
merge_step <- function(k, cc, state, params, Y, precomp, model, cache = NULL) {
  ev <- .merge_integrated_evaluate(cc, k, state, params, Y, precomp, model, cache)
  decision <- sample_categorical_log_strict(ev$log_weights)
  if (!ev$allowed[decision]) stop("An inadmissible merge action was selected.", call. = FALSE)
  action <- c("left", "keep", "right")[decision]
  destination <- ev$candidate_pool[decision]
  merged <- decision != 2L
  if (merged) {
    res <- execute_merge(cc, k, action, state, params, model)
    state <- res$state; params <- res$params
  }
  for (n in ev$touched_pool) {
    key <- as.character(n)
    posterior <- if (n == destination) ev$augmented[[key]] else ev$base[[key]]
    params$atom_pool[[cc]][[n]] <- .draw_atom_from_integrated_merge(posterior)
  }
  # regime-level atoms from the pool
  z <- params$z_state[[cc]]
  for (kk in seq_len(state$K[cc])) {
    params$atoms[[cc]][[kk]] <- params$atom_pool[[cc]][[z[kk]]]
  }
  list(state = state, params = params, merged = merged,
       decision = decision, action = action, evaluation = ev)
}

# Split scan of regime k over G_split = {tau_k + m_min, ..., tau_{k+1} - m_min}
split_step <- function(k, cc, state, params, Y, precomp, model) {
  K_c <- state$K[cc]; T_len <- model$T_len; m_min <- model$m_min
  no_op <- list(state = state, params = params, did_split = FALSE)
  tau_U <- state$tau_upper[[cc]]
  if (k > length(tau_U)) return(no_op)
  tau_k <- tau_U[k]
  tau_kp1 <- if (k + 1L <= length(tau_U)) tau_U[k + 1L] else T_len + 1L
  if (is.na(tau_k) || is.na(tau_kp1)) return(no_op)
  seg_len <- tau_kp1 - tau_k
  if (seg_len < 2L * m_min) return(no_op)

  g_left <- tau_k + m_min
  g_right <- tau_kp1 - m_min
  if (g_left > g_right) return(no_op)
  G_split <- seq.int(g_left, g_right)
  j_in_c <- which(state$cluster == cc)
  if (!length(j_in_c)) return(no_op)

  ps_cur <- params$pi_star[[cc]]
  ps_k <- if (k <= length(ps_cur)) min(max(ps_cur[k], 1e-10), 1 - 1e-10) else 0.5
  has_left <- k > 1L
  ps_km1 <- if (has_left && (k - 1L) <= length(ps_cur)) {
    min(max(ps_cur[k - 1L], 1e-10), 1 - 1e-10)
  } else 0.5

  atom_km1 <- if (has_left) params$atoms[[cc]][[k - 1L]] else NULL
  lambda0 <- if (!is.null(params$hyper$lambda0_K)) params$hyper$lambda0_K else 2.0
  log_poisson_split_penalty <- log(lambda0) - log(K_c + 1L)

  # the prospective left atom is integrated out; nothing is drawn yet
  cache <- .build_integrated_split_cache(
    cc, tau_k, G_split, j_in_c, state, params, Y, precomp)

  scan <- .split_seq_scan(
    cc, k, tau_k, tau_kp1, G_split, j_in_c,
    cache, atom_km1,
    ps_k, ps_km1, has_left, m_min, log_poisson_split_penalty,
    state, params, Y, precomp, model)

  if (scan$action == "nosplit") return(no_op)
  if (scan$action == "absorb") {
    res <- .absorb_left_into_prev(cc, k, scan$g, state, params, model)
    return(list(state = res$state, params = res$params, did_split = FALSE))
  }

  # atom drawn only after the split is confirmed
  new_left_atom <- .draw_atom_from_integrated_split(cache, scan$g)
  res <- execute_split(cc, k, scan$g, state, params, model,
                       new_left_atom = new_left_atom)
  list(state = res$state, params = res$params, did_split = TRUE)
}

# Absorb [tau_k, g - 1] into regime k - 1
.absorb_left_into_prev <- function(cc, k, g, state, params, model) {
  if (k <= 1L) return(list(state = state, params = params))
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

# Sequential scan: absorb, confirm or defer at each candidate; defer carries the suffix mass
.split_seq_scan <- function(cc, k, tau_k, tau_kp1, G_split, j_in_c,
                            integrated_cache, atom_km1,
                            ps_k, ps_km1, has_left, m_min,
                            log_poisson_split_penalty,
                            state, params, Y, precomp, model) {
  n_cand <- length(G_split)
  if (!n_cand) return(list(action = "nosplit", g = NA_integer_))
  if (!identical(as.integer(G_split), integrated_cache$G_split)) {
    stop("Integrated split cache does not match candidate set.", call. = FALSE)
  }

  # likelihood under the atom of regime k - 1, for absorb
  km1_params <- params
  if (has_left && !is.null(atom_km1)) {
    km1_params$atoms[[cc]] <- params$atoms[[cc]]
    km1_params$atoms[[cc]][[k]] <- atom_km1
  }
  seg_start <- tau_k
  n_seg <- tau_kp1 - tau_k
  ll_km1 <- numeric(n_seg)
  ll_right <- numeric(n_seg)
  for (loc in seq_len(n_seg)) {
    t1 <- seg_start + loc - 1L
    sr <- 0; sm <- 0
    for (j in j_in_c) {
      sr <- sr + log_lik_interval(j, t1, k, Y, state, params, precomp, model)
      if (has_left) sm <- sm + log_lik_interval(j, t1, k, Y, state, km1_params, precomp, model)
    }
    ll_right[loc] <- sr
    ll_km1[loc] <- sm
  }
  cum_right <- cumsum(ll_right)
  cum_km1 <- cumsum(ll_km1)
  total_right <- cum_right[n_seg]

  pk <- min(max(ps_k, 1e-10), 1 - 1e-10)
  pkm1 <- min(max(ps_km1, 1e-10), 1 - 1e-10)
  log_pi_k <- log(pk); log_1mpi_k <- log1p(-pk)
  log_pi_m <- log(pkm1); log_1mpi_m <- log1p(-pkm1)

  log_confirm <- rep(-Inf, n_cand)
  log_absorb <- rep(-Inf, n_cand)
  for (ci in seq_len(n_cand)) {
    g <- G_split[ci]
    sl <- g - tau_k
    if (sl < 1L || sl > n_seg - 1L) next
    LL_right <- total_right - cum_right[sl]
    log_confirm[ci] <- max(sl - m_min, 0L) * log_pi_k + log_1mpi_k +
      log_poisson_split_penalty + integrated_cache$log_new[ci] + LL_right
    if (has_left) {
      log_absorb[ci] <- sl * log_pi_m + log_1mpi_m +
        cum_km1[sl] + LL_right
    }
  }

  log_nosplit <- max(n_seg - m_min, 0L) * log_pi_k + total_right
  # suffix sums: total weight of the outcomes to the right
  log_suffix <- rep(-Inf, n_cand + 1L)
  log_suffix[n_cand + 1L] <- log_nosplit
  for (ci in rev(seq_len(n_cand))) {
    log_suffix[ci] <- .ia_logsumexp(c(
      log_absorb[ci], log_confirm[ci], log_suffix[ci + 1L]))
  }
  if (!is.finite(log_suffix[1L])) stop("IA split scan has no finite terminal weight.", call. = FALSE)

  for (ci in seq_len(n_cand)) {
    cand <- c(log_absorb[ci], log_confirm[ci], log_suffix[ci + 1L])
    sel <- sample_categorical_log_strict(cand)  # 1 absorb, 2 confirm, 3 defer
    if (sel == 1L) return(list(action = "absorb", g = G_split[ci]))
    if (sel == 2L) return(list(action = "split",  g = G_split[ci]))
    if (ci == n_cand) return(list(action = "nosplit", g = NA_integer_))
  }
  stop("IA split scan ended without a terminal outcome.", call. = FALSE)
}

# Merge regime k into its left or right neighbor
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

# New regime k on [tau_k, g - 1]; the old regime k becomes k + 1
execute_split <- function(cc, k, g, state, params, model, new_left_atom) {
  if (is.null(new_left_atom)) stop("execute_split requires the post-confirmation atom draw.", call. = FALSE)
  K_c <- state$K[cc]

  su <- state$S_upper[cc, ]
  su[su >= k] <- su[su >= k] + 1L
  if (g > 1L) {
    idxL <- seq_len(g - 1L)
    sel <- idxL[su[idxL] == (k + 1L)]
    if (length(sel)) su[sel] <- k
  }
  state$S_upper[cc, ] <- su
  state$K[cc] <- K_c + 1L
  state$tau_upper[[cc]] <- extract_changepoints(su)

  j_in_c <- which(state$cluster == cc)
  for (j in j_in_c) {
    state$S_lower[j, ] <- state$S_upper[cc, ]
    state$tau_lower[[j]] <- extract_changepoints(state$S_lower[j, ])
  }

  old_atoms <- params$atoms[[cc]]
  new_K <- K_c + 1L
  new_atoms <- vector("list", new_K)
  if (k > 1L) for (kk in seq_len(k - 1L)) new_atoms[[kk]] <- old_atoms[[kk]]
  new_atoms[[k]] <- canonicalize_atom(new_left_atom)
  for (kk in k:K_c) new_atoms[[kk + 1L]] <- old_atoms[[kk]]
  params$atoms[[cc]] <- new_atoms

  # new atom in an unoccupied pool entry, chosen with probability proportional to its weight
  if (!is.null(params$z_state) && !is.null(params$atom_pool)) {
    N <- if (!is.null(model$N_atom)) model$N_atom else 20L
    oz <- params$z_state[[cc]]
    if (is.null(oz) || length(oz) != K_c) {
      stop("Atom-pool allocation is inconsistent before split.", call. = FALSE)
    }
    pool <- params$atom_pool[[cc]]
    if (is.null(pool) || length(pool) != N) {
      stop("Atom pool is missing or has wrong size at confirmed split.", call. = FALSE)
    }
    free <- setdiff(seq_len(N), unique(oz))
    if (!length(free)) stop("No unoccupied atom-pool slot available for confirmed split.", call. = FALSE)
    if (!is.null(params$atom_pi[[cc]]) && length(params$atom_pi[[cc]]) == N) {
      ww <- pmax(params$atom_pi[[cc]][free], 1e-300)
      slot <- free[sample.int(length(free), 1L, prob = ww / sum(ww))]
    } else {
      slot <- free[sample.int(length(free), 1L)]
    }
    pool[[slot]] <- canonicalize_atom(new_left_atom)
    nz <- integer(new_K)
    if (k > 1L) nz[seq_len(k - 1L)] <- oz[seq_len(k - 1L)]
    nz[k] <- slot
    nz[(k + 1L):new_K] <- oz[k:K_c]
    params$atom_pool[[cc]] <- pool
    params$z_state[[cc]] <- nz
  }

  # stick mass of regime k shared equally by k and k + 1
  op <- params$pi_weights[[cc]]
  np <- numeric(new_K)
  if (k > 1L) for (kk in seq_len(k - 1L)) np[kk] <- op[kk]
  np[k] <- op[k] / 2
  np[k + 1L] <- np[k]
  if (k < K_c) for (kk in (k + 1L):K_c) np[kk + 1L] <- op[kk]
  np <- np / sum(np)
  params$pi_weights[[cc]] <- np
  params$pi_star[[cc]] <- compute_pi_star(np)

  # lower stay probabilities: copy of regime k, or the prior mean when k was terminal
  if (!is.null(params$pi_star_lower[[cc]])) {
    old_lower <- params$pi_star_lower[[cc]]
    if (length(old_lower) != K_c) stop("Lower persistence length mismatch before split.", call. = FALSE)
    if (k < K_c) {
      new_lower <- old_lower[k]
    } else {
      aL <- if (!is.null(params$hyper$a_pi_lower)) params$hyper$a_pi_lower else 50
      bL <- if (!is.null(params$hyper$b_pi_lower)) params$hyper$b_pi_lower else 1
      new_lower <- aL / (aL + bL)
    }
    lower <- append(old_lower, new_lower, after = k - 1L)
    lower[new_K] <- 1
    params$pi_star_lower[[cc]] <- lower
  }

  list(state = state, params = params)
}

# Relabel the regimes to 1:K after a merge; atoms, pool labels and stay probabilities follow
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
    sl <- cummax(mapped)  # nondecreasing
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

  op <- params$pi_weights[[cc]]
  np <- numeric(K_new)
  for (i in seq_along(old_labels)) {
    np[i] <- if (old_labels[i] <= length(op)) op[old_labels[i]] else 1 / K_new
  }
  np <- np / sum(np)
  params$pi_weights[[cc]] <- np
  params$pi_star[[cc]]    <- compute_pi_star(np)

  # lower stay probabilities of the surviving regimes
  if (!is.null(params$pi_star_lower[[cc]])) {
    old_lower <- params$pi_star_lower[[cc]]
    if (length(old_lower) < max_old) stop("Lower persistence length mismatch before relabel.", call. = FALSE)
    lower <- old_lower[old_labels]
    lower[K_new] <- 1
    params$pi_star_lower[[cc]] <- lower
  }

  list(state = state, params = params)
}

cat("07_interval_adjusting.R loaded (v11: integrated split and shared-pool merge; all candidate atoms redrawn).\n")
