
# --- Log-sum-exp (numerically stable) ---
log_sum_exp <- function(log_vals) {
  max_val <- max(log_vals[is.finite(log_vals)])
  if (!is.finite(max_val)) return(-Inf)
  max_val + log(sum(exp(log_vals - max_val)))
}

# --- Softmax (numerically stable) ---
softmax <- function(log_vals) {
  log_vals[!is.finite(log_vals)] <- -1e300
  max_val <- max(log_vals)
  if (!is.finite(max_val)) {
    n <- length(log_vals)
    return(rep(1/n, n))
  }
  w <- exp(log_vals - max_val)
  w[!is.finite(w)] <- 0
  s <- sum(w)
  if (s <= 0 || !is.finite(s)) return(rep(1/length(w), length(w)))
  w / s
}

# --- Categorical sampling from log-weights ---
sample_categorical_log <- function(log_weights) {
  probs <- softmax(log_weights)
  probs[is.na(probs)] <- 0
  s <- sum(probs)
  if (s <= 0 || !is.finite(s)) probs <- rep(1/length(probs), length(probs))
  sample.int(length(probs), size = 1, prob = probs)
}

# --- InvGamma sampling (shape-scale parameterization) ---
# p(x) ∝ x^{-(a+1)} exp(-b/x)
rinvgamma <- function(n, shape, scale) {
  1 / rgamma(n, shape = shape, rate = scale)
}


rtrunc_invgamma_lower1 <- function(shape, scale) {
  n <- length(shape)
  out <- rep(1.0, n)
  ok <- is.finite(shape) & is.finite(scale) & shape > 0 & scale > 0
  if (any(ok)) {
    a <- shape[ok]; b <- scale[ok]; gscale <- 1 / b
    F1 <- pgamma(1, shape = a, scale = gscale)        # P(Y <= 1), vectorized
    res <- rep(1.0, length(a))
    valid <- is.finite(F1) & F1 > 1e-300
    if (any(valid)) {
      u <- runif(sum(valid), 0, F1[valid])
      y <- qgamma(u, shape = a[valid], scale = gscale[valid])
      x <- 1 / y
      x[!is.finite(x) | x < 1] <- 1.0
      res[valid] <- x
    }
    out[ok] <- res
  }
  out
}

dinvgamma_log <- function(x, shape, scale) {
  -(shape + 1) * log(x) - scale / x + shape * log(scale) - lgamma(shape)
}

# --- Normal log-density ---
dnorm_log <- function(x, mean, sd) {
  -0.5 * log(2 * pi) - log(sd) - 0.5 * ((x - mean) / sd)^2
}

# --- Normal log-density (variance parameterization) ---
dnorm_log_var <- function(x, mean, var) {
  -0.5 * log(2 * pi) - 0.5 * log(var) - 0.5 * (x - mean)^2 / var
}

# --- Poisson log-density ---
dpois_log <- function(y, lambda) {
  y * log(lambda) - lambda - lfactorial(y)
}

# --- Half-normal sampling (truncated to positive) ---
rhalfnorm <- function(n, sd) {
  abs(rnorm(n, 0, sd))
}


# --- Canonical DP atom helpers ---
normalize_sign_carrier <- function(x, default = 1) {
  s <- sign(x)
  if (!is.finite(s) || s == 0) s <- default
  as.integer(ifelse(s >= 0, 1L, -1L))
}

canonicalize_atom <- function(atom) {
  if (is.null(atom)) return(NULL)
  sb <- atom$shape_beta
  sg <- atom$shape_gamma
  if (is.null(sb) || !is.finite(sb) || sb <= 0) sb <- 1.0
  if (is.null(sg) || !is.finite(sg) || sg < 0) sg <- 0.0
  list(
    gamma1 = normalize_sign_carrier(atom$gamma1, default = 1),
    gamma2 = normalize_sign_carrier(atom$gamma2, default = 1),
    shape_beta = as.numeric(max(sb, 1e-6)),
    shape_gamma = as.numeric(sg)
  )
}




sample_atom_from_base <- function(params) {
  hyper <- params$hyper

  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 2.0
  b_beta  <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 0.5
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 0.5
  b_gamma <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 1.0

  sb <- rgamma(1, shape = a_beta, rate = b_beta)
  if (!is.finite(sb) || sb < 0.01) sb <- 0.01
  sg <- rgamma(1, shape = a_gamma, rate = b_gamma)
  if (!is.finite(sg) || sg < 0) sg <- 0

  canonicalize_atom(list(
    gamma1 = sample(c(-1, 1), size = 1),
    gamma2 = sample(c(-1, 1), size = 1),
    shape_beta = sb,
    shape_gamma = sg
  ))
}

atom_log_prior <- function(atom, params) {
  if (is.null(atom)) return(-Inf)
  atom <- canonicalize_atom(atom)

  # Uniform prior on 4 sign combinations
  lp_sign <- -log(4)

  hyper <- params$hyper
  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 2.0
  b_beta  <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 0.5
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 0.5
  b_gamma <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 1.0

  sb <- atom$shape_beta
  sg <- atom$shape_gamma
  if (!is.finite(sb) || sb <= 0) return(-Inf)
  if (!is.finite(sg) || sg < 0) return(-Inf)

  # Gamma log-density: (a-1)*log(x) - b*x + a*log(b) - lgamma(a)
  lp_beta <- (a_beta - 1) * log(sb) - b_beta * sb + a_beta * log(b_beta) - lgamma(a_beta)

  # For shape_gamma: handle sg=0 gracefully when a_gamma < 1
  if (sg < 1e-300) {
    # Gamma(a<1) density → ∞ at 0; use a small positive floor
    sg_safe <- 1e-300
    lp_gamma <- (a_gamma - 1) * log(sg_safe) - b_gamma * sg_safe + a_gamma * log(b_gamma) - lgamma(a_gamma)
  } else {
    lp_gamma <- (a_gamma - 1) * log(sg) - b_gamma * sg + a_gamma * log(b_gamma) - lgamma(a_gamma)
  }

  if (!is.finite(lp_beta)) lp_beta <- -1e300
  if (!is.finite(lp_gamma)) lp_gamma <- -1e300

  lp_sign + lp_beta + lp_gamma
}


# --- Resize state-intercept vector gamma_{j,1:K} ---
resize_state_intercepts <- function(gamma_vec, K_new, fill_value = NULL) {
  if (is.null(gamma_vec) || length(gamma_vec) == 0) gamma_vec <- 0
  gamma_vec <- as.numeric(gamma_vec)
  K_old <- length(gamma_vec)
  if (K_new <= 0) return(numeric(0))
  if (K_old == K_new) return(gamma_vec)
  if (K_old > K_new) return(gamma_vec[seq_len(K_new)])
  if (is.null(fill_value)) fill_value <- gamma_vec[K_old]
  c(gamma_vec, rep(fill_value, K_new - K_old))
}

# --- Truncated normal sampling (positive only) ---
rtnorm_pos <- function(n, mean, sd) {
  samples <- numeric(n)
  for (i in 1:n) {
    repeat {
      x <- rnorm(1, mean, sd)
      if (x > 0) { samples[i] <- x; break }
    }
  }
  samples
}

safe_log <- function(x) {
  ifelse(x > 0, log(x), -Inf)
}

compute_pi_star <- function(pi_weights) {
  K <- length(pi_weights)
  pi_star <- numeric(K)
  cumsum_pi <- cumsum(pi_weights)
  for (k in 1:K) {
    remaining <- 1 - ifelse(k > 1, cumsum_pi[k - 1], 0)
    if (remaining > 1e-300) {
      pi_star[k] <- pi_weights[k] / remaining
    } else {
      pi_star[k] <- 0.5  # fallback
    }
    pi_star[k] <- min(max(pi_star[k], 1e-10), 1 - 1e-10)
  }
  pi_star
}

# --- Sample stick-breaking weights ---
sample_stick_breaking <- function(K, varsigma) {
  v <- rbeta(K, 1, varsigma)
  v[K] <- 1  # last stick gets everything remaining
  pi_weights <- numeric(K)
  pi_weights[1] <- v[1]
  for (k in 2:K) {
    pi_weights[k] <- v[k] * prod(1 - v[1:(k - 1)])
  }
  pi_weights
}

# --- Dirichlet sampling ---
rdirichlet <- function(alpha_vec) {
  x <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
  x / sum(x)
}

cat("00_utils.R loaded.\n")
