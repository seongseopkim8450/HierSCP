# Slice samplers, truncated normal samplers and bivariate normal orthant probabilities

# Univariate slice sampler with stepping out and shrinkage (Neal, 2003)
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

# Slice sampling of a positive parameter on the log scale
slice_sample_positive <- function(x0, log_target, w = 1, m = 10) {
  if (length(x0) != 1 || is.na(x0) || !is.finite(x0) || x0 <= 0) x0 <- 1

  log_x0 <- log(x0)

  log_target_log <- function(log_x) {
    x <- exp(log_x)
    if (!is.finite(x) || x <= 0) return(-Inf)
    val <- tryCatch(log_target(x), error = function(e) NA_real_)
    if (length(val) != 1 || is.na(val) || !is.finite(val)) return(-Inf)
    val + log_x  # Jacobian of x = exp(u)
  }

  log_x1 <- slice_sample(log_x0, log_target_log, w = w, m = m,
                         lower = -20, upper = 20)
  x1 <- exp(log_x1)
  if (!is.finite(x1) || x1 <= 0) return(x0)
  x1
}

# Elliptical slice sampling (Murray et al., 2010) with a diagonal Gaussian prior
elliptical_slice_sample <- function(theta_current, prior_mean, prior_cov_diag,
                                     log_likelihood) {
  p <- length(theta_current)

  g_current <- theta_current - prior_mean

  nu <- rnorm(p, mean = 0, sd = sqrt(prior_cov_diag))

  log_L_current <- log_likelihood(theta_current)

  log_u <- log_L_current + log(runif(1))

  angle <- runif(1, 0, 2 * pi)
  angle_min <- angle - 2 * pi
  angle_max <- angle

  max_iter <- 200
  for (iter in 1:max_iter) {
    g_proposed <- g_current * cos(angle) + nu * sin(angle)
    theta_proposed <- g_proposed + prior_mean

    log_L_proposed <- log_likelihood(theta_proposed)

    if (log_L_proposed > log_u) {
      return(theta_proposed)
    }

    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    angle <- runif(1, angle_min, angle_max)
  }

  # no point accepted within max_iter
  theta_current
}

# 64-point Gauss-Legendre rule on [0, 1]
.biv_gl64 <- local({
  n <- 64L
  k <- seq_len(n - 1L)
  off <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- off
  J[cbind(k + 1L, k)] <- off
  ee <- eigen(J, symmetric = TRUE)
  ord <- order(ee$values)
  nodes <- ee$values[ord]
  weights <- 2 * ee$vectors[1L, ord]^2
  list(y = (nodes + 1) / 2, w = weights / 2)
})

# N(mean, sd^2) truncated to (lower, Inf); exponential proposal in the tail (Robert, 1995)
rtnorm_lower_exact <- function(mean = 0, sd = 1, lower = 0) {
  if (length(mean) != 1L || length(sd) != 1L || length(lower) != 1L ||
      !is.finite(mean) || !is.finite(sd) || sd <= 0 || is.na(lower)) {
    stop("rtnorm_lower_exact: invalid arguments.", call. = FALSE)
  }
  if (lower == -Inf) return(rnorm(1L, mean, sd))
  if (!is.finite(lower)) stop("rtnorm_lower_exact: finite lower bound required.", call. = FALSE)

  a <- (lower - mean) / sd
  if (a <= 0) {
    repeat {
      z <- rnorm(1L)
      if (is.finite(z) && z > a) return(mean + sd * z)
    }
  }

  lambda <- 0.5 * (a + sqrt(a * a + 4))
  repeat {
    z <- a + rexp(1L, rate = lambda)
    if (log(runif(1L)) <= -0.5 * (z - lambda)^2) {
      return(mean + sd * z)
    }
  }
}

# log Phi_2(a, b; rho) by adaptive quadrature, for probabilities that underflow
.biv_logcdf_adaptive <- function(a, b, rho) {
  if (length(a) != 1L || length(b) != 1L || length(rho) != 1L ||
      !is.finite(a) || !is.finite(b) || !is.finite(rho) || abs(rho) >= 1) {
    stop(".biv_logcdf_adaptive: invalid standardized parameters.", call. = FALSE)
  }

  if (b < a) {
    tmp <- a; a <- b; b <- tmp
  }

  cond_sd <- sqrt(1 - rho^2)

  # Phi_2(a, b; rho) = int_{z < a} phi(z) Phi((b - rho z) / sqrt(1 - rho^2)) dz;
  # with z = a - exp(u) the integrand is bounded after removing its log-mode
  lower_u <- min(-40, -log1p(abs(a)) - 10)
  upper_u <- log(max(50, abs(a) + 50))

  log_integrand <- function(u) {
    t <- exp(u)
    z <- a - t
    arg <- (b - rho * z) / cond_sd
    dnorm(z, log = TRUE) + pnorm(arg, log.p = TRUE) + u
  }

  objective <- function(u) {
    lv <- log_integrand(u)
    if (length(lv) != 1L || !is.finite(lv)) return(.Machine$double.xmax / 100)
    -lv
  }

  opt <- tryCatch(
    optimize(objective, interval = c(lower_u, upper_u), tol = 1e-10),
    error = function(e) NULL
  )

  # grid check in case optimize() misses the mode
  ug <- seq(lower_u, upper_u, length.out = 65L)
  lg <- log_integrand(ug)
  lg[is.na(lg)] <- -Inf
  m_grid <- max(lg)
  m_opt <- if (!is.null(opt) && is.finite(opt$objective)) -opt$objective else -Inf
  m <- max(m_grid, m_opt)
  if (!is.finite(m)) {
    stop(".biv_logcdf_adaptive: could not locate a finite log-integrand mode.",
         call. = FALSE)
  }

  scaled_integrand <- function(u) {
    lv <- log_integrand(u) - m
    ans <- exp(lv)
    ans[!is.finite(ans)] <- 0
    ans
  }

  integ <- tryCatch(
    integrate(
      scaled_integrand,
      lower = lower_u, upper = upper_u,
      subdivisions = 400L,
      rel.tol = 1e-10, abs.tol = 1e-12,
      stop.on.error = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(integ) || length(integ$value) != 1L ||
      !is.finite(integ$value) || integ$value <= 0) {
    stop(".biv_logcdf_adaptive: adaptive quadrature failed.", call. = FALSE)
  }

  lp <- m + log(integ$value)
  if (!is.finite(lp)) {
    stop(".biv_logcdf_adaptive: non-finite log probability.", call. = FALSE)
  }

  # rounding can give lp slightly above 0
  min(lp, 0)
}

# log P(X > 0) for X ~ N2(mu, Sig)
biv_pos_logmass <- function(mu, Sig) {
  if (length(mu) != 2L || any(!is.finite(mu)) ||
      !is.matrix(Sig) || !identical(dim(Sig), c(2L, 2L)) ||
      any(!is.finite(Sig))) {
    stop("biv_pos_logmass: invalid bivariate Gaussian parameters.", call. = FALSE)
  }

  Sig <- 0.5 * (Sig + t(Sig))
  tryCatch(chol(Sig), error = function(e) {
    stop("biv_pos_logmass: covariance is not positive definite.", call. = FALSE)
  })

  s <- sqrt(diag(Sig))
  if (any(!is.finite(s)) || any(s <= 0)) {
    stop("biv_pos_logmass: invalid marginal standard deviations.", call. = FALSE)
  }

  rho <- Sig[1L, 2L] / (s[1L] * s[2L])
  if (!is.finite(rho)) stop("biv_pos_logmass: invalid correlation.", call. = FALSE)
  if (abs(rho) >= 1) {
    if (abs(rho) > 1 + 1e-10) {
      stop("biv_pos_logmass: correlation outside [-1,1].", call. = FALSE)
    }
    rho <- sign(rho) * (1 - 1e-12)
  }

  # P(X > 0) = Phi_2(mu_1/s_1, mu_2/s_2; rho)
  a <- mu[1L] / s[1L]
  b <- mu[2L] / s[2L]

  if (abs(rho) < 1e-14) {
    return(pnorm(a, log.p = TRUE) + pnorm(b, log.p = TRUE))
  }

  # TVPACK first; log-scale quadrature if its result underflows or fails
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    p <- tryCatch(
      as.numeric(mvtnorm::pmvnorm(
        lower = c(0, 0), upper = c(Inf, Inf),
        mean = as.numeric(mu), sigma = Sig,
        algorithm = mvtnorm::TVPACK()
      )),
      error = function(e) NA_real_
    )

    if (length(p) == 1L && is.finite(p)) {
      if (p >= 1 && p <= 1 + 1e-10) return(0)
      if (p > 1e-250 && p < 1) return(log(p))
    }
  }

  .biv_logcdf_adaptive(a, b, rho)
}

# Exact draw from N2(mu, Sig) truncated to (0, Inf)^2
rbivnorm_pos_exact <- function(mu, Sig, max_tries = 1000000L) {
  if (length(mu) != 2L || any(!is.finite(mu)) ||
      !is.matrix(Sig) || !identical(dim(Sig), c(2L, 2L)) ||
      any(!is.finite(Sig))) {
    stop("rbivnorm_pos_exact: invalid bivariate Gaussian parameters.", call. = FALSE)
  }

  Sig <- 0.5 * (Sig + t(Sig))
  tryCatch(chol(Sig), error = function(e) {
    stop("rbivnorm_pos_exact: covariance is not positive definite.", call. = FALSE)
  })

  s <- sqrt(diag(Sig))
  rho <- Sig[1L, 2L] / (s[1L] * s[2L])
  if (!is.finite(rho)) stop("rbivnorm_pos_exact: invalid correlation.", call. = FALSE)
  if (abs(rho) >= 1) {
    if (abs(rho) > 1 + 1e-10) {
      stop("rbivnorm_pos_exact: correlation outside [-1,1].", call. = FALSE)
    }
    rho <- sign(rho) * (1 - 1e-12)
  }
  cond_sd <- sqrt(1 - rho^2)
  ell <- -mu / s

  # mode at the corner: exponential proposals, acceptance exp(-y'Qy/2)
  Q <- matrix(c(1, -rho, -rho, 1), 2L, 2L) / (1 - rho^2)
  lambda_corner <- as.numeric(Q %*% ell)
  if (all(is.finite(lambda_corner)) && all(lambda_corner > 0.5)) {
    for (attempt in seq_len(as.integer(max_tries))) {
      y <- c(rexp(1L, rate = lambda_corner[1L]),
             rexp(1L, rate = lambda_corner[2L]))
      logacc <- -0.5 * as.numeric(crossprod(y, Q %*% y))
      if (log(runif(1L)) <= logacc) {
        u <- ell + y
        out <- mu + s * u
        if (all(is.finite(out)) && all(out > 0)) return(as.numeric(out))
        stop("rbivnorm_pos_exact: corner-envelope draw violated positivity numerically.",
             call. = FALSE)
      }
    }
    stop("rbivnorm_pos_exact: corner-envelope sampler exceeded max_tries.",
         call. = FALSE)
  }

  # otherwise one coordinate from its marginal by rejection, then the other given it
  envelope_logM <- function(li, lj) {
    if (rho > 1e-14) {
      0
    } else if (rho < -1e-14) {
      pnorm((rho * li - lj) / cond_sd, log.p = TRUE)
    } else {
      pnorm(-lj, log.p = TRUE)
    }
  }

  logden <- numeric(2L)
  for (i in 1:2) {
    j <- 3L - i
    logden[i] <- pnorm(-ell[i], log.p = TRUE) + envelope_logM(ell[i], ell[j])
  }
  first <- if (logden[1L] <= logden[2L]) 1L else 2L
  second <- 3L - first

  li <- ell[first]
  lj <- ell[second]
  logM <- envelope_logM(li, lj)

  for (attempt in seq_len(as.integer(max_tries))) {
    ui <- rtnorm_lower_exact(0, 1, li)
    logq <- pnorm((rho * ui - lj) / cond_sd, log.p = TRUE)
    logacc <- min(0, logq - logM)

    if (log(runif(1L)) <= logacc) {
      uj <- rtnorm_lower_exact(rho * ui, cond_sd, lj)
      u <- numeric(2L)
      u[first] <- ui
      u[second] <- uj
      out <- mu + s * u
      if (all(is.finite(out)) && all(out > 0)) return(as.numeric(out))
      stop("rbivnorm_pos_exact: accepted draw violated positivity numerically.",
           call. = FALSE)
    }
  }

  stop(
    "rbivnorm_pos_exact: exact accept-reject sampler exceeded max_tries; ",
    "the selected positive-orthant conditional is exceptionally extreme.",
    call. = FALSE
  )
}

cat("04_sampling.R loaded.\n")
