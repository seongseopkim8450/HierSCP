# Cosine basis and D_{mm'}(x) = int_0^x int_s^1 phi_m(t) phi_m'(t) dt ds,
# so that H_j(x) = theta' D(x) theta / theta'theta

# phi_m and D(x) on a regular grid of [0, 1]
precompute_basis <- function(M, n_grid = 500) {
  x_grid <- seq(0, 1, length.out = n_grid)

  phi_mat <- matrix(0, nrow = M + 1, ncol = n_grid)
  phi_mat[1, ] <- 1
  for (mm in 1:M) {
    phi_mat[mm + 1, ] <- sqrt(2) * cos(mm * pi * x_grid)
  }

  D_array <- array(0, dim = c(M + 1, M + 1, n_grid))

  for (mi in 0:M) {
    for (mj in mi:M) {
      D_vals <- compute_D_analytical(mi, mj, x_grid)
      D_array[mi + 1, mj + 1, ] <- D_vals
      if (mi != mj) {
        D_array[mj + 1, mi + 1, ] <- D_vals
      }
    }
  }

  list(
    M = M,
    n_grid = n_grid,
    x_grid = x_grid,
    phi_mat = phi_mat,
    D_array = D_array
  )
}

# D_{mm'}(x) in closed form
compute_D_analytical <- function(m, mp, x) {
  n <- length(x)
  result <- numeric(n)

  if (m == 0 && mp == 0) {
    result <- x - x^2 / 2
  } else if (m == 0 && mp > 0) {
    mp_pi <- mp * pi
    result <- sqrt(2) / mp_pi^2 * (cos(mp_pi * x) - 1)
  } else if (m > 0 && mp == 0) {
    m_pi <- m * pi
    result <- sqrt(2) / m_pi^2 * (cos(m_pi * x) - 1)
  } else if (m == mp) {
    two_m_pi <- 2 * m * pi
    result <- x - x^2 / 2 + (cos(two_m_pi * x) - 1) / two_m_pi^2
  } else {
    # phi_m phi_m' = cos((m - m') pi t) + cos((m + m') pi t)
    d_minus <- (m - mp) * pi
    d_plus  <- (m + mp) * pi

    result <- (cos(d_minus * x) - 1) / d_minus^2 +
              (cos(d_plus * x) - 1) / d_plus^2
  }

  result
}

eval_basis <- function(m, x) {
  if (m == 0) return(rep(1, length(x)))
  sqrt(2) * cos(m * pi * x)
}

# Z(x) = sum_m theta_m phi_m(x)
eval_Z <- function(theta_j, x_local, precomp) {
  M <- precomp$M
  n <- length(x_local)
  phi_local <- matrix(0, nrow = M + 1, ncol = n)
  phi_local[1, ] <- 1
  for (mm in 1:M) {
    phi_local[mm + 1, ] <- sqrt(2) * cos(mm * pi * x_local)
  }
  as.numeric(crossprod(theta_j, phi_local))
}

# theta' D(x) theta at arbitrary x
eval_double_integral <- function(theta_j, x_local, precomp) {
  M <- precomp$M
  n <- length(x_local)
  result <- numeric(n)

  for (i in seq_along(x_local)) {
    xi <- x_local[i]
    D_xi <- matrix(0, M + 1, M + 1)
    for (mi in 0:M) {
      for (mj in mi:M) {
        val <- compute_D_analytical(mi, mj, xi)
        D_xi[mi + 1, mj + 1] <- val
        if (mi != mj) D_xi[mj + 1, mi + 1] <- val
      }
    }
    result[i] <- as.numeric(t(theta_j) %*% D_xi %*% theta_j)
  }
  result
}

# theta' D(x) theta with D(x) interpolated on the grid
eval_double_integral_fast <- function(theta_j, x_local, precomp) {
  M <- precomp$M

  theta_outer <- tcrossprod(theta_j)

  n <- length(x_local)
  result <- numeric(n)

  for (i in seq_along(x_local)) {
    xi <- x_local[i]
    idx <- findInterval(xi, precomp$x_grid, all.inside = TRUE)
    dx <- precomp$x_grid[2] - precomp$x_grid[1]
    w <- (xi - precomp$x_grid[idx]) / dx

    D_lo <- precomp$D_array[, , idx]
    D_hi <- precomp$D_array[, , min(idx + 1, precomp$n_grid)]
    D_xi <- (1 - w) * D_lo + w * D_hi

    result[i] <- sum(theta_outer * D_xi)
  }
  result
}

cat("01_precompute.R loaded.\n")
