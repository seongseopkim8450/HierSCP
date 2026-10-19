# Numerical helpers, samplers and the base measure F0 of the shape atoms

log_sum_exp <- function(log_vals) {
  max_val <- max(log_vals[is.finite(log_vals)])
  if (!is.finite(max_val)) return(-Inf)
  max_val + log(sum(exp(log_vals - max_val)))
}

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

# Categorical draw from unnormalized log-weights
sample_categorical_log <- function(log_weights) {
  probs <- softmax(log_weights)
  probs[is.na(probs)] <- 0
  s <- sum(probs)
  if (s <= 0 || !is.finite(s)) probs <- rep(1/length(probs), length(probs))
  sample.int(length(probs), size = 1, prob = probs)
}

# Inverse-gamma in the shape-scale parameterization
rinvgamma <- function(n, shape, scale) {
  1 / rgamma(n, shape = shape, rate = scale)
}

# InvGamma(shape, scale) truncated to [1, Inf), by inversion for z = 1/phi on (0, 1]
rtrunc_invgamma_lower1 <- function(shape, scale) {
  n <- length(shape)
  if (n == 0L) return(numeric(0))
  if (length(scale) != n || any(!is.finite(shape)) ||
      any(!is.finite(scale)) || any(shape <= 0) || any(scale <= 0)) {
    stop("rtrunc_invgamma_lower1: invalid parameters.")
  }

  log_Z <- pgamma(1, shape = shape, rate = scale,
                  lower.tail = TRUE, log.p = TRUE)
  if (any(!is.finite(log_Z))) {
    stop("rtrunc_invgamma_lower1: invalid truncation probability.")
  }

  log_u <- log(runif(n)) + log_Z
  z <- qgamma(log_u, shape = shape, rate = scale,
              lower.tail = TRUE, log.p = TRUE)
  if (any(!is.finite(z)) || any(z <= 0) || any(z > 1)) {
    stop("rtrunc_invgamma_lower1: invalid precision draw.")
  }

  phi <- 1 / z
  if (any(!is.finite(phi)) || any(phi < 1)) {
    stop("rtrunc_invgamma_lower1: invalid scale draw.")
  }
  phi
}

dinvgamma_log <- function(x, shape, scale) {
  -(shape + 1) * log(x) - scale / x + shape * log(scale) - lgamma(shape)
}

dnorm_log <- function(x, mean, sd) {
  -0.5 * log(2 * pi) - log(sd) - 0.5 * ((x - mean) / sd)^2
}

dnorm_log_var <- function(x, mean, var) {
  -0.5 * log(2 * pi) - 0.5 * log(var) - 0.5 * (x - mean)^2 / var
}

dpois_log <- function(y, lambda) {
  y * log(lambda) - lambda - lfactorial(y)
}

rhalfnorm <- function(n, sd) {
  abs(rnorm(n, 0, sd))
}

# Sign coded as +1 or -1
normalize_sign_carrier <- function(x, default = 1) {
  s <- sign(x)
  if (!is.finite(s) || s == 0) s <- default
  as.integer(ifelse(s >= 0, 1L, -1L))
}

# Atom (delta1, delta2, b, g) stored as (gamma1, gamma2, shape_beta, shape_gamma)
canonicalize_atom <- function(atom) {
  if (is.null(atom)) return(NULL)
  sb <- atom$shape_beta
  sg <- atom$shape_gamma
  if (is.null(sb) || !is.finite(sb) || sb <= 0) sb <- 1.0
  if (is.null(sg) || !is.finite(sg) || sg <= 0) sg <- 1.0
  list(
    gamma1 = normalize_sign_carrier(atom$gamma1, default = 1),
    gamma2 = normalize_sign_carrier(atom$gamma2, default = 1),
    shape_beta = as.numeric(sb),
    shape_gamma = as.numeric(sg)
  )
}

# Magnitude prior N+((a - 1)/r, (a - 1)/r^2); its normalizer Phi(sqrt(a - 1)) is free of r
shape_magnitude_prior_params <- function(a, r) {
  if (length(a) != 1L || !is.finite(a) || a <= 1) {
    stop("Shape-magnitude prior requires a > 1.", call. = FALSE)
  }
  if (length(r) != 1L || !is.finite(r) || r <= 0) {
    stop("Shape-magnitude prior requires rate r > 0.", call. = FALSE)
  }
  mu <- (a - 1) / r
  var <- (a - 1) / r^2
  sd <- sqrt(var)
  prec <- 1 / var
  log_Z <- pnorm(sqrt(a - 1), log.p = TRUE)
  list(mu = mu, var = var, sd = sd, prec = prec, log_Z = log_Z)
}

rshape_magnitude_prior <- function(n, a, r) {
  pr <- shape_magnitude_prior_params(a, r)
  n <- as.integer(n)
  if (n < 1L) return(numeric(0))
  F0 <- pnorm((0 - pr$mu) / pr$sd)
  u <- F0 + runif(n) * (1 - F0)
  x <- qnorm(u, mean = pr$mu, sd = pr$sd)
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop("Positive-truncated Gaussian draw failed.", call. = FALSE)
  }
  x
}

dshape_magnitude_prior_log <- function(x, a, r) {
  pr <- shape_magnitude_prior_params(a, r)
  out <- rep(-Inf, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- dnorm(x[ok], mean = pr$mu, sd = pr$sd, log = TRUE) - pr$log_Z
  out
}

# Draw from F0: uniform sign pair and N+ magnitudes
sample_atom_from_base <- function(params) {
  hyper <- params$hyper
  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 2.0
  r_beta  <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 1.5
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 1.5
  r_gamma <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 2.0

  sb <- rshape_magnitude_prior(1L, a_beta,  r_beta)
  sg <- rshape_magnitude_prior(1L, a_gamma, r_gamma)

  canonicalize_atom(list(
    gamma1 = sample(c(-1L, 1L), size = 1L),
    gamma2 = sample(c(-1L, 1L), size = 1L),
    shape_beta = sb,
    shape_gamma = sg
  ))
}

# log F0 density of an atom
atom_log_prior <- function(atom, params) {
  if (is.null(atom)) return(-Inf)
  hyper <- params$hyper

  a_beta  <- if (!is.null(hyper$a_shape_beta))  hyper$a_shape_beta  else 2.0
  r_beta  <- if (!is.null(hyper$b_shape_beta))  hyper$b_shape_beta  else 1.5
  a_gamma <- if (!is.null(hyper$a_shape_gamma)) hyper$a_shape_gamma else 1.5
  r_gamma <- if (!is.null(hyper$b_shape_gamma)) hyper$b_shape_gamma else 2.0

  sb <- atom$shape_beta
  sg <- atom$shape_gamma
  if (is.null(sb) || is.null(sg) || !is.finite(sb) || sb <= 0 || !is.finite(sg) || sg <= 0) return(-Inf)

  -log(4) +
    dshape_magnitude_prior_log(sb, a_beta,  r_beta) +
    dshape_magnitude_prior_log(sg, a_gamma, r_gamma)
}

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

# Positive-truncated normal by rejection
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

# Stay probabilities pi*_k = pi_k / (mass of regimes k, ..., K); regime K is terminal
compute_pi_star <- function(pi_weights) {
  if (!is.numeric(pi_weights) || !length(pi_weights) ||
      anyNA(pi_weights) || any(!is.finite(pi_weights)) ||
      any(pi_weights < 0)) {
    stop("Stick masses must be a nonempty finite nonnegative numeric vector.",
         call. = FALSE)
  }
  K <- length(pi_weights)
  # tail sums taken directly; 1 - cumsum() loses small tails
  remaining <- rev(cumsum(rev(pi_weights)))
  if (any(!is.finite(remaining)) || any(remaining <= 0)) {
    stop("Stick masses must have a finite positive remaining tail at every regime.",
         call. = FALSE)
  }
  pi_star <- pi_weights / remaining
  pi_star <- pmin(pmax(pi_star, 1e-10), 1 - 1e-10)
  pi_star[K] <- 1 - 1e-10
  pi_star
}

# Stick-breaking masses with v_k ~ Beta(1, varsigma) and v_K = 1
sample_stick_breaking <- function(K, varsigma) {
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) ||
      K < 1 || K != as.integer(K)) {
    stop("Stick initialization requires a positive integer regime count.", call. = FALSE)
  }
  if (!is.numeric(varsigma) || length(varsigma) != 1L ||
      !is.finite(varsigma) || varsigma <= 0) {
    stop("Stick initialization requires a finite positive concentration.", call. = FALSE)
  }
  K <- as.integer(K)
  v <- rbeta(K, 1, varsigma)
  if (any(!is.finite(v)) || any(v < 0 | v > 1)) {
    stop("Invalid Beta draw in stick initialization.", call. = FALSE)
  }
  # rbeta() can return exactly 0 or 1
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  v[K] <- 1
  pi_weights <- numeric(K)
  pi_weights[1] <- v[1]
  if (K >= 2L) for (k in 2:K) {
    pi_weights[k] <- v[k] * prod(1 - v[1:(k - 1)])
  }
  pi_weights
}

rdirichlet <- function(alpha_vec) {
  x <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
  x / sum(x)
}

cat("00_utils.R loaded.\n")

# Strict version: stops on invalid weights instead of using uniform weights
sample_categorical_log_strict <- function(log_weights) {
  if (!is.numeric(log_weights) || !length(log_weights) ||
      anyNA(log_weights) || any(log_weights == Inf)) {
    stop("Invalid categorical log-weights.")
  }
  finite <- is.finite(log_weights)
  if (!any(finite)) stop("Categorical distribution has zero finite mass.")
  m <- max(log_weights[finite])
  w <- numeric(length(log_weights))
  w[finite] <- exp(log_weights[finite] - m)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) stop("Invalid categorical normalizer.")
  sample.int(length(w), 1L, prob = w / s)
}
