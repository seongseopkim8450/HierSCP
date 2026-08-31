###############################################################################
# 04_sampling.R
# Stepping-out Slice Sampling (Neal 2003) and 
# Elliptical Slice Sampling (Murray et al. 2010)
###############################################################################

#' Stepping-out Slice Sampling (Neal 2003)
#'
#' @param x0 Current value
#' @param log_target Function: log target density (unnormalized)
#' @param w Initial bracket width
#' @param m Maximum number of stepping-out steps
#' @param lower Lower bound for x (default -Inf)
#' @param upper Upper bound for x (default Inf)
#' @return New sample
slice_sample <- function(x0, log_target, w = 1, m = 10, lower = -Inf, upper = Inf) {
  safe_eval <- function(x) {
    val <- tryCatch(log_target(x), error = function(e) NA_real_)
    if (length(val) != 1 || is.na(val) || !is.finite(val)) return(-Inf)
    val
  }

  if (length(x0) != 1 || is.na(x0) || !is.finite(x0)) {
    x0 <- if (is.finite(lower) && is.finite(upper)) (lower + upper) / 2 else 0
  }
  if (is.finite(lower)) x0 <- max(x0, lower)
  if (is.finite(upper)) x0 <- min(x0, upper)
  if (!is.finite(w) || w <= 0) w <- 1
  if (!is.finite(m) || m < 1) m <- 10

  log_y0 <- safe_eval(x0)
  if (!is.finite(log_y0)) return(x0)

  log_z <- log_y0 - rexp(1)

  u <- runif(1)
  L <- x0 - w * u
  R <- L + w
  if (!is.finite(L)) L <- x0 - w
  if (!is.finite(R)) R <- x0 + w

  J <- floor(m * runif(1))
  K <- (m - 1) - J

  while (J > 0 && L >= lower) {
    ll <- safe_eval(L)
    if (!(ll > log_z)) break
    L <- L - w
    if (!is.finite(L)) break
    J <- J - 1
  }
  while (K > 0 && R <= upper) {
    rr <- safe_eval(R)
    if (!(rr > log_z)) break
    R <- R + w
    if (!is.finite(R)) break
    K <- K - 1
  }

  L <- max(L, lower)
  R <- min(R, upper)
  if (!is.finite(L) || !is.finite(R) || L > R) return(x0)

  max_iter <- 200
  for (iter in 1:max_iter) {
    x1 <- runif(1, L, R)
    if (is.na(x1) || !is.finite(x1)) return(x0)
    log_y1 <- safe_eval(x1)

    if (log_y1 >= log_z) {
      return(x1)
    }

    if (x1 < x0) {
      L <- x1
    } else {
      R <- x1
    }
    if (!is.finite(L) || !is.finite(R) || L > R) return(x0)
  }

  x0
}

#' Slice sampling for log-transformed positive parameters
#' Samples on log scale, returns on original scale
#'
#' @param x0 Current value (positive)
#' @param log_target Log target density on original scale
#' @param w Bracket width on log scale
#' @param m Max stepping-out steps
#' @return New sample (positive)
slice_sample_positive <- function(x0, log_target, w = 1, m = 10) {
  if (length(x0) != 1 || is.na(x0) || !is.finite(x0) || x0 <= 0) x0 <- 1

  log_x0 <- log(x0)

  log_target_log <- function(log_x) {
    x <- exp(log_x)
    if (!is.finite(x) || x <= 0) return(-Inf)
    val <- tryCatch(log_target(x), error = function(e) NA_real_)
    if (length(val) != 1 || is.na(val) || !is.finite(val)) return(-Inf)
    val + log_x  # Jacobian: d/d(log x) = x
  }

  log_x1 <- slice_sample(log_x0, log_target_log, w = w, m = m,
                         lower = -20, upper = 20)
  x1 <- exp(log_x1)
  if (!is.finite(x1) || x1 <= 0) return(x0)
  x1
}

#' Elliptical Slice Sampling (Murray et al. 2010)
#' 
#' Samples from posterior ∝ N(mu, Sigma) * L(theta)
#' by operating on centered variable g = theta - mu ~ N(0, Sigma)
#'
#' @param theta_current Current parameter vector
#' @param prior_mean Prior mean vector (mu)
#' @param prior_cov_diag Diagonal of prior covariance (vector)
#' @param log_likelihood Function: log-likelihood given theta
#' @return New theta vector
elliptical_slice_sample <- function(theta_current, prior_mean, prior_cov_diag,
                                     log_likelihood) {
  
  p <- length(theta_current)
  
  # Center: g = theta - mu
  g_current <- theta_current - prior_mean
  
  # Sample auxiliary from prior: nu ~ N(0, Sigma)
  nu <- rnorm(p, mean = 0, sd = sqrt(prior_cov_diag))
  
  # Current log-likelihood
  log_L_current <- log_likelihood(theta_current)
  
  # Define slice level
  log_u <- log_L_current + log(runif(1))
  
  # Draw initial angle
  angle <- runif(1, 0, 2 * pi)
  angle_min <- angle - 2 * pi
  angle_max <- angle
  
  max_iter <- 200
  for (iter in 1:max_iter) {
    # Propose on ellipse
    g_proposed <- g_current * cos(angle) + nu * sin(angle)
    theta_proposed <- g_proposed + prior_mean
    
    log_L_proposed <- log_likelihood(theta_proposed)
    
    if (log_L_proposed > log_u) {
      return(theta_proposed)
    }
    
    # Shrink bracket
    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    angle <- runif(1, angle_min, angle_max)
  }
  
  # Fallback: return current
  theta_current
}

cat("04_sampling.R loaded.\n")
