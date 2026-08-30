
if (!exists("compute_atom_loglik")) {
  compute_atom_loglik <- function(cc, k, atom, j_in_c, t_in_k_upper,
                                  state, params, Y, precomp, model) {
    if (length(j_in_c) == 0 || length(t_in_k_upper) == 0) return(0)
    ll <- 0
    for (j in j_in_c) {
      t_jk <- t_in_k_upper[ state$S_lower[j, t_in_k_upper] == k ]
      if (length(t_jk) == 0) next
      theta_j <- params$theta[j, ]
      gamma_j <- get_state_intercept(j, k, params)
      f_jk <- eval_shape_at_times(theta_j, atom, t_jk, precomp)
      mu <- params$alpha + params$beta[cc] + gamma_j + f_jk
      v_j <- params$v[[j]]; if (is.null(v_j)||!is.finite(v_j[1])||v_j[1]<=0) v_j <- 1 else v_j <- v_j[1]
      for (i in seq_along(t_jk)) {
        t <- t_jk[i]
        phi_t <- params$phi[j,t]; if (!is.finite(phi_t)||phi_t<1) phi_t <- 1
        xi_t <- params$xi[j,t]; if (is.na(xi_t)) xi_t <- 0
        sig2 <- v_j * phi_t^xi_t; if (!is.finite(sig2)||sig2<=0) sig2 <- 1
        ll <- ll + dnorm_log_var(Y[j,t], mu[i], sig2)
      }
    }
    if (!is.finite(ll)) ll <- -1e300
    ll
  }
}
