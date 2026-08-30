###############################################################################
# 01_precompute.R
# Karhunen-Loève basis precomputation for shape-restricted functions
###############################################################################

#' Precompute all basis-related quantities
#' 
#' @param M Number of K-L basis terms (0 to M)
#' @param n_grid Number of grid points for precomputation
#' @return List with basis evaluation functions and precomputed arrays
precompute_basis <- function(M, n_grid = 500) {
  
  x_grid <- seq(0, 1, length.out = n_grid)
  
  # --- Basis function values: phi[m+1, i] = phi_m(x_grid[i]) ---
  # phi_0(x) = 1, phi_m(x) = sqrt(2)*cos(m*pi*x) for m >= 1
  phi_mat <- matrix(0, nrow = M + 1, ncol = n_grid)
  phi_mat[1, ] <- 1
  for (mm in 1:M) {
    phi_mat[mm + 1, ] <- sqrt(2) * cos(mm * pi * x_grid)
  }
  
  # --- D_{mm'}(x) = int_0^x int_s^1 phi_m(t)*phi_{m'}(t) dt ds ---
  # Precompute analytically for all (m, m') pairs at grid points
  D_array <- array(0, dim = c(M + 1, M + 1, n_grid))
  
  for (mi in 0:M) {
    for (mj in mi:M) {  # symmetric
      D_vals <- compute_D_analytical(mi, mj, x_grid)
      D_array[mi + 1, mj + 1, ] <- D_vals
      if (mi != mj) {
        D_array[mj + 1, mi + 1, ] <- D_vals  # symmetry
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

#' Analytical computation of D_{mm'}(x) = int_0^x int_s^1 phi_m(t)*phi_{m'}(t) dt ds
#' 
#' @param m First basis index (0-based)
#' @param mp Second basis index (0-based)
#' @param x Vector of evaluation points in [0,1]
#' @return Vector of D values
compute_D_analytical <- function(m, mp, x) {
  n <- length(x)
  result <- numeric(n)
  
  if (m == 0 && mp == 0) {
    # D_{00}(x) = x - x^2/2
    result <- x - x^2 / 2
    
  } else if (m == 0 && mp > 0) {
    # D_{0,m'}(x) = sqrt(2)/(m'*pi)^2 * [cos(m'*pi*x) - 1]
    mp_pi <- mp * pi
    result <- sqrt(2) / mp_pi^2 * (cos(mp_pi * x) - 1)
    
  } else if (m > 0 && mp == 0) {
    # Same as D_{0,m} by symmetry of phi product
    m_pi <- m * pi
    result <- sqrt(2) / m_pi^2 * (cos(m_pi * x) - 1)
    
  } else if (m == mp) {
    # D_{mm}(x) = x - x^2/2 + [cos(2m*pi*x) - 1]/(2m*pi)^2
    two_m_pi <- 2 * m * pi
    result <- x - x^2 / 2 + (cos(two_m_pi * x) - 1) / two_m_pi^2
    
  } else {
    # m > 0, m' > 0, m != m'
    # phi_m * phi_{m'} = 2*cos(m*pi*t)*cos(m'*pi*t) 
    #                   = cos((m-m')*pi*t) + cos((m+m')*pi*t)
    d_minus <- (m - mp) * pi  # could be negative, but cos is even
    d_plus  <- (m + mp) * pi
    
    result <- (cos(d_minus * x) - 1) / d_minus^2 + 
              (cos(d_plus * x) - 1) / d_plus^2
  }
  
  result
}

#' Evaluate basis functions at arbitrary points
#' 
#' @param m Basis index (0-based)
#' @param x Evaluation points
#' @return Vector of phi_m(x) values
eval_basis <- function(m, x) {
  if (m == 0) return(rep(1, length(x)))
  sqrt(2) * cos(m * pi * x)
}

#' Evaluate Z_j(x) = sum_{m=0}^{M} theta_{j,m} * phi_m(x) at given points
#' 
#' @param theta_j Vector of K-L coefficients (length M+1)
#' @param x_local Vector of local time points in [0,1]
#' @param precomp Precomputed basis object
#' @return Vector of Z values at x_local
eval_Z <- function(theta_j, x_local, precomp) {
  M <- precomp$M
  # Build basis matrix at x_local
  n <- length(x_local)
  phi_local <- matrix(0, nrow = M + 1, ncol = n)
  phi_local[1, ] <- 1
  for (mm in 1:M) {
    phi_local[mm + 1, ] <- sqrt(2) * cos(mm * pi * x_local)
  }
  # Z(x) = theta^T phi(x)
  as.numeric(crossprod(theta_j, phi_local))
}

#' Compute double integral at arbitrary x using interpolation
#' 
#' @param theta_j Vector of K-L coefficients (length M+1)
#' @param x_local Vector of local time points in [0,1]
#' @param precomp Precomputed basis object
#' @return Vector of double integral values: theta^T D(x) theta
eval_double_integral <- function(theta_j, x_local, precomp) {
  M <- precomp$M
  n <- length(x_local)
  result <- numeric(n)
  
  for (i in seq_along(x_local)) {
    xi <- x_local[i]
    # Compute D_{mm'}(xi) for all pairs
    D_xi <- matrix(0, M + 1, M + 1)
    for (mi in 0:M) {
      for (mj in mi:M) {
        val <- compute_D_analytical(mi, mj, xi)
        D_xi[mi + 1, mj + 1] <- val
        if (mi != mj) D_xi[mj + 1, mi + 1] <- val
      }
    }
    # theta^T D(x) theta
    result[i] <- as.numeric(t(theta_j) %*% D_xi %*% theta_j)
  }
  result
}

#' Vectorized: compute double integral using precomputed grid + interpolation
#' 
#' @param theta_j Vector of K-L coefficients
#' @param x_local Vector of local time points
#' @param precomp Precomputed basis object
#' @return Vector of theta^T D(x) theta values
eval_double_integral_fast <- function(theta_j, x_local, precomp) {
  M <- precomp$M
  
  # Compute theta_m * theta_{m'} outer product (only need to do once per theta update)
  theta_outer <- tcrossprod(theta_j)  # (M+1) x (M+1)
  
  n <- length(x_local)
  result <- numeric(n)
  
  for (i in seq_along(x_local)) {
    xi <- x_local[i]
    # Find closest grid index for interpolation
    idx <- findInterval(xi, precomp$x_grid, all.inside = TRUE)
    # Linear interpolation weight
    dx <- precomp$x_grid[2] - precomp$x_grid[1]
    w <- (xi - precomp$x_grid[idx]) / dx
    
    # Interpolate D array
    D_lo <- precomp$D_array[, , idx]
    D_hi <- precomp$D_array[, , min(idx + 1, precomp$n_grid)]
    D_xi <- (1 - w) * D_lo + w * D_hi
    
    # theta^T D(xi) theta
    result[i] <- sum(theta_outer * D_xi)
  }
  result
}

cat("01_precompute.R loaded.\n")
