

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace arma;


static arma::cube G_DATX;
static bool G_DATX_SET = false;

// [[Rcpp::export]]
Rcpp::IntegerVector set_precomp_datx(arma::cube D_at_x) {
  G_DATX = D_at_x;        // single copy, once per run (or on dimension change)
  G_DATX_SET = true;
  return Rcpp::IntegerVector::create((int)D_at_x.n_rows,
                                     (int)D_at_x.n_cols,
                                     (int)D_at_x.n_slices);
}

inline double dnorm_log_var_s(double x, double mu, double sig2) {
  const double LOG2PI = 1.8378770664093453;
  return -0.5 * LOG2PI - 0.5 * std::log(sig2) - 0.5 * (x - mu) * (x - mu) / sig2;
}

inline double safe_sign(double x) {
  return (x >= 0.0) ? 1.0 : -1.0; // exact match to R code: sign(0) -> 1
}

inline double clamp_prob(double x) {
  if (!std::isfinite(x)) return 0.5;
  if (x < 1e-10) return 1e-10;
  if (x > 1.0 - 1e-10) return 1.0 - 1e-10;
  return x;
}

inline double last_finite_or_zero(const NumericVector& x) {
  double out = 0.0;
  bool found = false;
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) {
      out = x[i];
      found = true;
    }
  }
  return found ? out : 0.0;
}


inline double logsumexp_pair(double a, double b) {
  if (!std::isfinite(a)) return b;
  if (!std::isfinite(b)) return a;
  double m = (a > b) ? a : b;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// [DIS] Direct categorical sampling from log probabilities.
// Replaces sequential hazard scan. No path dependency.
inline int sample_categorical_log_cpp(const NumericVector& log_probs) {
  const int n = log_probs.size();
  if (n <= 0) return 1;
  if (n == 1) return 1;

  // Find max for numerical stability
  double max_lp = R_NegInf;
  for (int i = 0; i < n; ++i) {
    if (R_finite(log_probs[i]) && log_probs[i] > max_lp) max_lp = log_probs[i];
  }
  if (!R_finite(max_lp)) {
    return static_cast<int>(std::floor(R::runif(0.0, 1.0) * n)) + 1;
  }

  // Compute normalized probabilities via softmax
  double sum_p = 0.0;
  NumericVector probs(n);
  for (int i = 0; i < n; ++i) {
    if (R_finite(log_probs[i])) {
      probs[i] = std::exp(log_probs[i] - max_lp);
    } else {
      probs[i] = 0.0;
    }
    sum_p += probs[i];
  }
  if (sum_p <= 0.0) {
    return static_cast<int>(std::floor(R::runif(0.0, 1.0) * n)) + 1;
  }

  // Draw from categorical distribution
  double u = R::runif(0.0, 1.0) * sum_p;
  double cum = 0.0;
  for (int i = 0; i < n; ++i) {
    cum += probs[i];
    if (u <= cum) return i + 1; // 1-based
  }
  return n; // fallback
}

inline double get_gamma_val(const NumericVector& gamma_j, int k_1based) {
  if (gamma_j.size() == 0) return 0.0;
  if (k_1based >= 1 && k_1based <= gamma_j.size() && R_finite(gamma_j[k_1based - 1])) {
    return gamma_j[k_1based - 1];
  }
  return last_finite_or_zero(gamma_j);
}


inline double get_v_val(const NumericVector& v_j, int k_1based) {
  if (v_j.size() == 0) return 1.0;
  // Always return the single scalar value
  if (R_finite(v_j[0]) && v_j[0] > 0.0) return v_j[0];
  return 1.0;
}

struct AtomPars {
  bool present;
  double gamma1;
  double gamma2;
  double shape_beta;
  double shape_gamma;
  AtomPars() : present(false), gamma1(0.1), gamma2(0.1), shape_beta(1.0), shape_gamma(0.0) {}
};

inline AtomPars parse_atom(SEXP atom_sexp) {
  AtomPars out;
  if (Rf_isNull(atom_sexp) || !Rf_isNewList(atom_sexp)) return out;

  List atom(atom_sexp);
  out.present = true;
  if (atom.containsElementNamed("gamma1")) out.gamma1 = as<double>(atom["gamma1"]);
  
  if (atom.containsElementNamed("gamma2")) out.gamma2 = as<double>(atom["gamma2"]);
  
  if (atom.containsElementNamed("shape_beta")) out.shape_beta = as<double>(atom["shape_beta"]);
  if (atom.containsElementNamed("shape_gamma")) out.shape_gamma = as<double>(atom["shape_gamma"]);
  
  // 방어 코드: 부호 제약 보장
  if (!std::isfinite(out.shape_beta) || out.shape_beta <= 0.0) out.shape_beta = 1.0;
  if (!std::isfinite(out.shape_gamma) || out.shape_gamma < 0.0) out.shape_gamma = 0.0;
  
  return out;
  
}

inline double theta_sq_sum(const NumericVector& theta_j) {
  double out = 0.0;
  for (int m = 1; m < theta_j.size(); ++m) out += theta_j[m] * theta_j[m];
  return out;
}

inline double compute_DI_single(const NumericVector& theta_j,
                                int t0) {
  const int Mp1 = theta_j.size();
  double DI = 0.0;
  for (int mi = 1; mi < Mp1; ++mi) {
    const double th_mi = theta_j[mi];
    DI += th_mi * th_mi * G_DATX(mi, mi, t0);
    for (int mj = mi + 1; mj < Mp1; ++mj) {
      DI += 2.0 * th_mi * theta_j[mj] * G_DATX(mi, mj, t0);
    }
  }
  return DI;
}


inline double shape_value_single(const NumericVector& theta_j,
                                 const AtomPars& atom,
                                 int t0,
                                 const NumericVector& x_global) {
  if (!atom.present) return 0.0;
  if (theta_j.size() <= 1) return 0.0;

  const double delta1 = safe_sign(atom.gamma1);
  const double delta2 = safe_sign(atom.gamma2);
  const double SZ = theta_sq_sum(theta_j);
  const double DI = compute_DI_single(theta_j, t0);
  
  
  const double SZ_safe = (SZ < 1e-12) ? 1e-12 : SZ;
  const double H_val = DI / SZ_safe;

  const double linear_coef = delta1 * atom.shape_gamma + atom.shape_beta * (delta1 + delta2) / 2.0;
  
  return linear_coef * x_global[t0] - delta2 * atom.shape_beta * H_val;
}

// [[Rcpp::export]]
NumericVector dnorm_log_var_cpp(NumericVector x, NumericVector mean, NumericVector var) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = dnorm_log_var_s(x[i], mean[i], var[i]);
  return out;
}

// [[Rcpp::export]]
NumericVector eval_shape_at_times_cpp(NumericVector theta_j,
                                      SEXP atom_sexp,
                                      IntegerVector t_indices,
                                      NumericVector x_global) {
  AtomPars atom = parse_atom(atom_sexp);
  NumericVector out(t_indices.size(), 0.0);
  if (!atom.present || theta_j.size() <= 1) return out;

  for (int i = 0; i < t_indices.size(); ++i) {
    int t0 = t_indices[i] - 1;
    if (t0 < 0 || t0 >= x_global.size()) continue;
    out[i] = shape_value_single(theta_j, atom, t0, x_global);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector compute_f_all_cpp(IntegerVector S_lower_j,
                                int K_c,
                                NumericVector theta_j,
                                List atoms_c,
                                NumericVector x_global) {
  int T_len = S_lower_j.size();
  NumericVector f_vals(T_len, 0.0);
  if (theta_j.size() <= 1) return f_vals;

  std::vector<AtomPars> atoms(std::max(K_c, 0));
  for (int k = 1; k <= K_c; ++k) {
    if (k <= atoms_c.size()) atoms[k - 1] = parse_atom(atoms_c[k - 1]);
  }

  for (int t = 0; t < T_len; ++t) {
    int k = S_lower_j[t];
    if (k < 1 || k > K_c) continue;
    f_vals[t] = shape_value_single(theta_j, atoms[k - 1], t, x_global);
  }
  return f_vals;
}

// [[Rcpp::export]]
NumericVector compute_mu_all_cpp(double alpha,
                                 double beta_c,
                                 IntegerVector S_lower_j,
                                 int K_c,
                                 NumericVector theta_j,
                                 NumericVector gamma_j,
                                 List atoms_c,
                                 NumericVector x_global) {
  int T_len = S_lower_j.size();
  NumericVector mu(T_len);
  NumericVector f_vals = compute_f_all_cpp(S_lower_j, K_c, theta_j, atoms_c, x_global);
  const double ab = alpha + beta_c;

  for (int t = 0; t < T_len; ++t) {
    int k = S_lower_j[t];
    double gamma_val = get_gamma_val(gamma_j, k);
    mu[t] = ab + gamma_val + f_vals[t];
  }
  return mu;
}

// [[Rcpp::export]]
List compute_mu_and_variance_cpp(IntegerVector S_lower_j,
                                 int K_c,
                                 NumericVector theta_j,
                                 NumericVector gamma_j,
                                 List atoms_c,
                                 double alpha,
                                 double beta_c,
                                 NumericVector v_j,
                                 NumericVector phi_j,
                                 IntegerVector xi_j,
                                 NumericVector x_global) {
  int T_len = S_lower_j.size();
  NumericVector mu = compute_mu_all_cpp(alpha, beta_c, S_lower_j, K_c, theta_j, gamma_j,
                                        atoms_c, x_global);
  NumericVector sigma2(T_len, 1.0);

  for (int t = 0; t < T_len; ++t) {
    int k = S_lower_j[t];
    double v_val = get_v_val(v_j, k);
    double phi_val = (t < phi_j.size() && R_finite(phi_j[t]) && phi_j[t] >= 1.0) ? phi_j[t] : 1.0;
    int xi_val = (t < xi_j.size()) ? xi_j[t] : 0;
    double s2 = v_val * std::pow(phi_val, static_cast<double>(xi_val));
    if (!std::isfinite(s2) || s2 <= 0.0) s2 = 1.0;
    sigma2[t] = s2;
  }

  return List::create(Named("mu") = mu, Named("sigma2") = sigma2);
}

// [[Rcpp::export]]
double log_lik_series_cpp(IntegerVector S_lower_j,
                          int K_c,
                          NumericVector theta_j,
                          NumericVector gamma_j,
                          List atoms_c,
                          NumericVector Y_j,
                          double alpha,
                          double beta_c,
                          NumericVector v_j,
                          NumericVector phi_j,
                          IntegerVector xi_j,
                          NumericVector x_global) {
  List mu_sig = compute_mu_and_variance_cpp(S_lower_j, K_c, theta_j, gamma_j, atoms_c,
                                            alpha, beta_c, v_j, phi_j, xi_j,
                                            x_global);
  NumericVector mu = mu_sig["mu"];
  NumericVector sigma2 = mu_sig["sigma2"];

  double ll = 0.0;
  for (int t = 0; t < Y_j.size(); ++t) ll += dnorm_log_var_s(Y_j[t], mu[t], sigma2[t]);
  if (!std::isfinite(ll)) ll = -1e300;
  return ll;
}

// [[Rcpp::export]]
double log_lik_interval_cpp(IntegerVector t_indices,
                            IntegerVector state_seq,
                            NumericVector Y_j,
                            NumericVector theta_j,
                            NumericVector gamma_j,
                            List atoms_c,
                            double alpha,
                            double beta_c,
                            NumericVector v_j,
                            NumericVector phi_j,
                            IntegerVector xi_j,
                            NumericVector x_global,
                            int var_state_override = 0) {
  int n = t_indices.size();
  if (n == 0) return 0.0;

  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    int t0 = t_indices[i] - 1;
    if (t0 < 0 || t0 >= Y_j.size()) continue;

    int k = state_seq[i];
    double gamma_val = get_gamma_val(gamma_j, k);
    AtomPars atom;
    if (k >= 1 && k <= atoms_c.size()) atom = parse_atom(atoms_c[k - 1]);
    double f_val = shape_value_single(theta_j, atom, t0, x_global);
    double mu_val = alpha + beta_c + gamma_val + f_val;

    int k_for_var = (var_state_override > 0) ? var_state_override : k;
    double v_val = get_v_val(v_j, k_for_var);
    double phi_val = (t0 < phi_j.size() && R_finite(phi_j[t0]) && phi_j[t0] >= 1.0) ? phi_j[t0] : 1.0;
    int xi_val = (t0 < xi_j.size()) ? xi_j[t0] : 0;
    double sig2 = v_val * std::pow(phi_val, static_cast<double>(xi_val));
    if (!std::isfinite(sig2) || sig2 <= 0.0) sig2 = 1.0;

    ll += dnorm_log_var_s(Y_j[t0], mu_val, sig2);
  }

  if (!std::isfinite(ll)) ll = -1e300;
  return ll;
}

// [DIS] Direct Independence Sampling for lower-level changepoint.
// [[Rcpp::export]]
int lower_la_pair_cpp(NumericVector theta_j,
                      double gamma_k,
                      SEXP atom_k_sexp,
                      double gamma_kp1,
                      SEXP atom_kp1_sexp,
                      IntegerVector G_L,
                      NumericVector Y_j,
                      double alpha,
                      double beta_c,
                      double v_k,
                      double v_kp1,
                      NumericVector phi_j,
                      IntegerVector xi_j,
                      double pi_star_k,
                      int anchor_tau,                 
                      NumericVector x_global) {
  int n_free = G_L.size();
  if (n_free == 0) return 1;

  AtomPars atom_k = parse_atom(atom_k_sexp);
  AtomPars atom_kp1 = parse_atom(atom_kp1_sexp);
  const double log_pi = std::log(clamp_prob(pi_star_k));
  const double log_1mpi = std::log(1.0 - clamp_prob(pi_star_k));
  const double ab = alpha + beta_c;

  NumericVector ll_k(n_free), ll_kp1(n_free);
  for (int idx = 0; idx < n_free; ++idx) {
    int t0 = G_L[idx] - 1;
    double mu_k = ab + gamma_k + shape_value_single(theta_j, atom_k, t0, x_global);
    double mu_kp1 = ab + gamma_kp1 + shape_value_single(theta_j, atom_kp1, t0, x_global);

    double phi_val = (t0 < phi_j.size() && R_finite(phi_j[t0]) && phi_j[t0] >= 1.0) ? phi_j[t0] : 1.0;
    int xi_val = (t0 < xi_j.size()) ? xi_j[t0] : 0;
    double sig2_k = v_k * std::pow(phi_val, static_cast<double>(xi_val));
    double sig2_kp1 = v_kp1 * std::pow(phi_val, static_cast<double>(xi_val));
    if (!std::isfinite(sig2_k) || sig2_k <= 0.0) sig2_k = 1.0;
    if (!std::isfinite(sig2_kp1) || sig2_kp1 <= 0.0) sig2_kp1 = 1.0;

    ll_k[idx] = dnorm_log_var_s(Y_j[t0], mu_k, sig2_k);
    ll_kp1[idx] = dnorm_log_var_s(Y_j[t0], mu_kp1, sig2_kp1);
  }

  NumericVector cum_k(n_free), cum_kp1(n_free), log_probs(n_free);
  cum_k[0] = ll_k[0];
  cum_kp1[0] = ll_kp1[0];
  for (int i = 1; i < n_free; ++i) {
    cum_k[i] = cum_k[i - 1] + ll_k[i];
    cum_kp1[i] = cum_kp1[i - 1] + ll_kp1[i];
  }
  double total_kp1 = cum_kp1[n_free - 1];


  const int g_plus = G_L[n_free - 1];
  for (int i = 0; i < n_free; ++i) {
    int tau = G_L[i];
    double ll_before = (i > 0) ? cum_k[i - 1] : 0.0;
    double ll_after = total_kp1 - ((i > 0) ? cum_kp1[i - 1] : 0.0);
    double log_prior = static_cast<double>(anchor_tau - tau) * log_pi;
    if (tau != (g_plus + 1)) log_prior += log_1mpi;
    log_probs[i] = log_prior + ll_before + ll_after;
  }

  return sample_categorical_log_cpp(log_probs);
}


// [DIS] Direct Independence Sampling for upper-level changepoint.
// [[Rcpp::export]]
int upper_la_pair_cpp(NumericMatrix theta_mat,
                      NumericVector gamma_k_vec,
                      SEXP atom_k_sexp,
                      NumericVector gamma_kp1_vec,
                      SEXP atom_kp1_sexp,
                      IntegerVector G_U,
                      NumericMatrix Y_sub,
                      double alpha,
                      double beta_c,
                      NumericVector v_k_vec,
                      NumericVector v_kp1_vec,
                      NumericMatrix phi_sub,
                      IntegerMatrix xi_sub,
                      double pi_star_k,
                      int anchor_tau,                 // [신규] anchor τ^U
                      NumericVector x_global) {
  int n_free = G_U.size();
  if (n_free == 0) return 1;
  int n_series = theta_mat.nrow();

  AtomPars atom_k = parse_atom(atom_k_sexp);
  AtomPars atom_kp1 = parse_atom(atom_kp1_sexp);
  const double log_pi = std::log(clamp_prob(pi_star_k));
  const double log_1mpi = std::log(1.0 - clamp_prob(pi_star_k));
  const double ab = alpha + beta_c;

  NumericVector agg_ll_k(n_free, 0.0), agg_ll_kp1(n_free, 0.0);
  for (int s = 0; s < n_series; ++s) {
    NumericVector theta_j = theta_mat(s, _);
    double gamma_k = (s < gamma_k_vec.size()) ? gamma_k_vec[s] : 0.0;
    double gamma_kp1 = (s < gamma_kp1_vec.size()) ? gamma_kp1_vec[s] : 0.0;
    double v_k = (s < v_k_vec.size() && std::isfinite(v_k_vec[s]) && v_k_vec[s] > 0.0) ? v_k_vec[s] : 1.0;
    double v_kp1 = (s < v_kp1_vec.size() && std::isfinite(v_kp1_vec[s]) && v_kp1_vec[s] > 0.0) ? v_kp1_vec[s] : 1.0;

    for (int idx = 0; idx < n_free; ++idx) {
      int t0 = G_U[idx] - 1;
      double mu_k = ab + gamma_k + shape_value_single(theta_j, atom_k, t0, x_global);
      double mu_kp1 = ab + gamma_kp1 + shape_value_single(theta_j, atom_kp1, t0, x_global);

      double phi_val = (s < phi_sub.nrow() && t0 < phi_sub.ncol() && R_finite(phi_sub(s, t0)) && phi_sub(s, t0) >= 1.0)
        ? phi_sub(s, t0) : 1.0;
      int xi_val = (s < xi_sub.nrow() && t0 < xi_sub.ncol()) ? xi_sub(s, t0) : 0;
      double sig2_k = v_k * std::pow(phi_val, static_cast<double>(xi_val));
      double sig2_kp1 = v_kp1 * std::pow(phi_val, static_cast<double>(xi_val));
      if (!std::isfinite(sig2_k) || sig2_k <= 0.0) sig2_k = 1.0;
      if (!std::isfinite(sig2_kp1) || sig2_kp1 <= 0.0) sig2_kp1 = 1.0;

      agg_ll_k[idx] += dnorm_log_var_s(Y_sub(s, t0), mu_k, sig2_k);
      agg_ll_kp1[idx] += dnorm_log_var_s(Y_sub(s, t0), mu_kp1, sig2_kp1);
    }
  }

  NumericVector cum_k(n_free), cum_kp1(n_free), log_probs(n_free);
  cum_k[0] = agg_ll_k[0];
  cum_kp1[0] = agg_ll_kp1[0];
  for (int i = 1; i < n_free; ++i) {
    cum_k[i] = cum_k[i - 1] + agg_ll_k[i];
    cum_kp1[i] = cum_kp1[i - 1] + agg_ll_kp1[i];
  }
  double total_kp1 = cum_kp1[n_free - 1];

  
  const int g_plus = G_U[n_free - 1];
  for (int i = 0; i < n_free; ++i) {
    int tau = G_U[i];
    double ll_before = (i > 0) ? cum_k[i - 1] : 0.0;
    double ll_after = total_kp1 - ((i > 0) ? cum_kp1[i - 1] : 0.0);
    double log_prior = static_cast<double>(anchor_tau - tau) * log_pi;
    if (tau != (g_plus + 1)) log_prior += log_1mpi;
    log_probs[i] = log_prior + ll_before + ll_after;
  }

  return sample_categorical_log_cpp(log_probs);
}



// [[Rcpp::export]]
int split_scan_upper_cpp(NumericMatrix theta_mat,
                         NumericVector gamma_left_vec,    // gamma^*_{j} for left block (scalar per series)
                         SEXP atom_left_sexp,             // NEW LEFT atom (proposed) -> split (p)
                         NumericVector gamma_right_vec,   // gamma^*_{j} for right block
                         SEXP atom_right_sexp,            // OLD atom (current state k)
                         SEXP atom_km1_sexp,              // PREVIOUS atom (state k-1) -> absorb (p-)
                         NumericVector gamma_km1_vec,     // gamma^*_{j} for k-1 (absorbed left block)
                         IntegerVector G_split,           // candidate g (first index of right block)
                         int tau_k,                       // segment start (1-based)
                         int tau_kp1,                     // segment end + 1 (1-based)
                         NumericMatrix Y_sub,
                         double alpha,
                         double beta_c,
                         NumericVector v_left_vec,        // v for left block under NEW atom (= series v)
                         NumericVector v_right_vec,       // v for right block (= series v)
                         NumericVector v_km1_vec,         // v for left block under k-1 atom (= series v)
                         NumericMatrix phi_sub,
                         IntegerMatrix xi_sub,
                         double pi_star_k,                // persistence of state k   (split base)
                         double pi_star_km1,              // persistence of state k-1 (absorb base)
                         int m_min,
                         double log_poisson_split_penalty,
                         int has_left_neighbor,           // 1 if k>1 (absorb admissible), else 0
                         NumericVector x_global) {
  const int n_cand = G_split.size();
  const int n_series = theta_mat.nrow();
  
  const int idx_nosplit = 2 * n_cand + 1;   // 1-based index of the nosplit option
  if (n_cand == 0) return idx_nosplit;

  AtomPars atom_left  = parse_atom(atom_left_sexp);    // p  : new left atom
  AtomPars atom_right = parse_atom(atom_right_sexp);   // right block (old k)
  AtomPars atom_km1   = parse_atom(atom_km1_sexp);     // p- : absorbed into k-1

  const double pk   = clamp_prob(pi_star_k);
  const double pkm1 = clamp_prob(pi_star_km1);
  const double log_pi_k    = std::log(pk);
  const double log_1mpi_k  = std::log(1.0 - pk);
  const double log_pi_km1  = std::log(pkm1);
  const double log_1mpi_km1= std::log(1.0 - pkm1);
  const double ab = alpha + beta_c;

  const int seg_start = tau_k;
  const int seg_end   = tau_kp1 - 1;
  const int n_seg = seg_end - seg_start + 1;
  if (n_seg <= 0) return idx_nosplit;

  // Per-time aggregated log-likelihood:
  //   ll_new[loc]   : left block under NEW atom        (for split p)
  //   ll_km1[loc]   : left block under k-1 atom         (for absorb p-)
  //   ll_right[loc] : block under OLD atom k            (right block in both; whole seg in nosplit)
  std::vector<double> ll_new(n_seg, 0.0), ll_km1(n_seg, 0.0), ll_right(n_seg, 0.0);

  for (int s = 0; s < n_series; ++s) {
    NumericVector theta_j = theta_mat(s, _);
    double g_new   = (s < gamma_left_vec.size())  ? gamma_left_vec[s]  : 0.0;
    double g_right = (s < gamma_right_vec.size()) ? gamma_right_vec[s] : 0.0;
    double g_km1   = (s < gamma_km1_vec.size())   ? gamma_km1_vec[s]   : 0.0;
    double v_n = (s < v_left_vec.size()  && std::isfinite(v_left_vec[s])  && v_left_vec[s]  > 0.0) ? v_left_vec[s]  : 1.0;
    double v_r = (s < v_right_vec.size() && std::isfinite(v_right_vec[s]) && v_right_vec[s] > 0.0) ? v_right_vec[s] : 1.0;
    double v_m = (s < v_km1_vec.size()   && std::isfinite(v_km1_vec[s])   && v_km1_vec[s]   > 0.0) ? v_km1_vec[s]   : 1.0;

    for (int loc = 0; loc < n_seg; ++loc) {
      int t0 = (seg_start + loc) - 1;   // 0-based
      if (t0 < 0 || t0 >= Y_sub.ncol()) continue;

      double f_new = shape_value_single(theta_j, atom_left,  t0, x_global);
      double f_r   = shape_value_single(theta_j, atom_right, t0, x_global);
      double f_m   = has_left_neighbor ? shape_value_single(theta_j, atom_km1, t0, x_global) : 0.0;

      double mu_new = ab + g_new   + f_new;
      double mu_r   = ab + g_right + f_r;
      double mu_m   = ab + g_km1   + f_m;

      double phi_val = (s < phi_sub.nrow() && t0 < phi_sub.ncol() && R_finite(phi_sub(s, t0)) && phi_sub(s, t0) >= 1.0)
        ? phi_sub(s, t0) : 1.0;
      int xi_val = (s < xi_sub.nrow() && t0 < xi_sub.ncol()) ? xi_sub(s, t0) : 0;
      double base_pow = std::pow(phi_val, static_cast<double>(xi_val));
      double sig2_n = v_n * base_pow; if (!std::isfinite(sig2_n) || sig2_n <= 0.0) sig2_n = 1.0;
      double sig2_r = v_r * base_pow; if (!std::isfinite(sig2_r) || sig2_r <= 0.0) sig2_r = 1.0;
      double sig2_m = v_m * base_pow; if (!std::isfinite(sig2_m) || sig2_m <= 0.0) sig2_m = 1.0;

      double y = Y_sub(s, t0);
      ll_new[loc]   += dnorm_log_var_s(y, mu_new, sig2_n);
      ll_right[loc] += dnorm_log_var_s(y, mu_r,   sig2_r);
      if (has_left_neighbor) ll_km1[loc] += dnorm_log_var_s(y, mu_m, sig2_m);
    }
  }

  // Prefix cumulative sums.
  std::vector<double> cum_new(n_seg), cum_km1(n_seg), cum_right(n_seg);
  cum_new[0]   = ll_new[0];
  cum_km1[0]   = ll_km1[0];
  cum_right[0] = ll_right[0];
  for (int i = 1; i < n_seg; ++i) {
    cum_new[i]   = cum_new[i - 1]   + ll_new[i];
    cum_km1[i]   = cum_km1[i - 1]   + ll_km1[i];
    cum_right[i] = cum_right[i - 1] + ll_right[i];
  }
  const double total_right = cum_right[n_seg - 1];   

  NumericVector log_probs(idx_nosplit);   

  for (int c = 0; c < n_cand; ++c) {
    int g = G_split[c];                 // first index of right block (1-based)
    int split_loc = g - seg_start;      // # of left-block points = local prefix length
    if (split_loc < 1 || split_loc > n_seg - 1) {
      log_probs[c]            = R_NegInf;   // split
      log_probs[n_cand + c]   = R_NegInf;   // absorb
      continue;
    }
    // left block = [tau_k, g-1] = local prefix 0..split_loc-1  -> cum[...][split_loc-1]
    double LL_right = total_right - cum_right[split_loc - 1];
    int n_left = g - tau_k;             // length of the left block

    // ---- SPLIT (p): left = NEW atom, base = pi*_k, m_min-discounted ----
    {
      double LL_left = cum_new[split_loc - 1];
      int eff = n_left - m_min; if (eff < 0) eff = 0;
      double lp = static_cast<double>(eff) * log_pi_k + log_1mpi_k
                  + log_poisson_split_penalty + LL_left + LL_right;
      log_probs[c] = std::isfinite(lp) ? lp : R_NegInf;
    }

    // ---- ABSORB (p-): left = k-1 atom, base = pi*_{k-1}, NO m_min discount ----
    if (has_left_neighbor) {
      double LL_left = cum_km1[split_loc - 1];
      int eff = n_left; if (eff < 0) eff = 0;     // image: exponent = g - tau_k (no m_min)
      double lp = static_cast<double>(eff) * log_pi_km1 + log_1mpi_km1
                  + LL_left + LL_right;
      log_probs[n_cand + c] = std::isfinite(lp) ? lp : R_NegInf;
    } else {
      log_probs[n_cand + c] = R_NegInf;           // no k-1 neighbor -> absorb inadmissible
    }
  }

  // ---- NOSPLIT: whole segment under OLD atom k, base = pi*_k ----
  {
    int eff = n_seg - m_min; if (eff < 0) eff = 0;
    double lp = static_cast<double>(eff) * log_pi_k + total_right;
    log_probs[idx_nosplit - 1] = std::isfinite(lp) ? lp : R_NegInf;
  }

  return sample_categorical_log_cpp(log_probs);
}



int split_seq_scan_upper_cpp(NumericMatrix theta_mat,
                             NumericVector gamma_left_vec,
                             SEXP atom_left_sexp,
                             NumericVector gamma_right_vec,
                             SEXP atom_right_sexp,
                             SEXP atom_km1_sexp,
                             NumericVector gamma_km1_vec,
                             IntegerVector G_split,
                             int tau_k,
                             int tau_kp1,
                             NumericMatrix Y_sub,
                             double alpha,
                             double beta_c,
                             NumericVector v_left_vec,
                             NumericVector v_right_vec,
                             NumericVector v_km1_vec,
                             NumericMatrix phi_sub,
                             IntegerMatrix xi_sub,
                             double pi_star_k,
                             double pi_star_km1,
                             int m_min,
                             double log_poisson_split_penalty,
                             int has_left_neighbor,
                             double log_q_atom,        
                             NumericVector x_global) {
  const int n_cand = G_split.size();
  const int n_series = theta_mat.nrow();
  const int idx_nosplit = 2 * n_cand + 1;
  if (n_cand == 0) return idx_nosplit;

  AtomPars atom_left  = parse_atom(atom_left_sexp);
  AtomPars atom_right = parse_atom(atom_right_sexp);
  AtomPars atom_km1   = parse_atom(atom_km1_sexp);

  const double pk   = clamp_prob(pi_star_k);
  const double pkm1 = clamp_prob(pi_star_km1);
  const double log_pi_k     = std::log(pk);
  const double log_1mpi_k   = std::log(1.0 - pk);
  const double log_pi_km1   = std::log(pkm1);
  const double log_1mpi_km1 = std::log(1.0 - pkm1);
  const double ab = alpha + beta_c;

  const int seg_start = tau_k;
  const int seg_end   = tau_kp1 - 1;
  const int n_seg = seg_end - seg_start + 1;
  if (n_seg <= 0) return idx_nosplit;

  
  std::vector<double> ll_new(n_seg, 0.0), ll_km1(n_seg, 0.0), ll_right(n_seg, 0.0);
  for (int s = 0; s < n_series; ++s) {
    NumericVector theta_j = theta_mat(s, _);
    double g_new   = (s < gamma_left_vec.size())  ? gamma_left_vec[s]  : 0.0;
    double g_right = (s < gamma_right_vec.size()) ? gamma_right_vec[s] : 0.0;
    double g_km1   = (s < gamma_km1_vec.size())   ? gamma_km1_vec[s]   : 0.0;
    double v_n = (s < v_left_vec.size()  && std::isfinite(v_left_vec[s])  && v_left_vec[s]  > 0.0) ? v_left_vec[s]  : 1.0;
    double v_r = (s < v_right_vec.size() && std::isfinite(v_right_vec[s]) && v_right_vec[s] > 0.0) ? v_right_vec[s] : 1.0;
    double v_m = (s < v_km1_vec.size()   && std::isfinite(v_km1_vec[s])   && v_km1_vec[s]   > 0.0) ? v_km1_vec[s]   : 1.0;
    for (int loc = 0; loc < n_seg; ++loc) {
      int t0 = (seg_start + loc) - 1;
      if (t0 < 0 || t0 >= Y_sub.ncol()) continue;
      double f_new = shape_value_single(theta_j, atom_left,  t0, x_global);
      double f_r   = shape_value_single(theta_j, atom_right, t0, x_global);
      double f_m   = has_left_neighbor ? shape_value_single(theta_j, atom_km1, t0, x_global) : 0.0;
      double mu_new = ab + g_new   + f_new;
      double mu_r   = ab + g_right + f_r;
      double mu_m   = ab + g_km1   + f_m;
      double phi_val = (s < phi_sub.nrow() && t0 < phi_sub.ncol() && R_finite(phi_sub(s, t0)) && phi_sub(s, t0) >= 1.0)
        ? phi_sub(s, t0) : 1.0;
      int xi_val = (s < xi_sub.nrow() && t0 < xi_sub.ncol()) ? xi_sub(s, t0) : 0;
      double base_pow = std::pow(phi_val, static_cast<double>(xi_val));
      double sig2_n = v_n * base_pow; if (!std::isfinite(sig2_n) || sig2_n <= 0.0) sig2_n = 1.0;
      double sig2_r = v_r * base_pow; if (!std::isfinite(sig2_r) || sig2_r <= 0.0) sig2_r = 1.0;
      double sig2_m = v_m * base_pow; if (!std::isfinite(sig2_m) || sig2_m <= 0.0) sig2_m = 1.0;
      double y = Y_sub(s, t0);
      ll_new[loc]   += dnorm_log_var_s(y, mu_new, sig2_n);
      ll_right[loc] += dnorm_log_var_s(y, mu_r,   sig2_r);
      if (has_left_neighbor) ll_km1[loc] += dnorm_log_var_s(y, mu_m, sig2_m);
    }
  }
  std::vector<double> cum_new(n_seg), cum_km1(n_seg), cum_right(n_seg);
  cum_new[0] = ll_new[0]; cum_km1[0] = ll_km1[0]; cum_right[0] = ll_right[0];
  for (int i = 1; i < n_seg; ++i) {
    cum_new[i]   = cum_new[i - 1]   + ll_new[i];
    cum_km1[i]   = cum_km1[i - 1]   + ll_km1[i];
    cum_right[i] = cum_right[i - 1] + ll_right[i];
  }
  const double total_right = cum_right[n_seg - 1];
  const int eff_ns_raw = n_seg - m_min;
  const double eff_ns = (eff_ns_raw < 0) ? 0.0 : static_cast<double>(eff_ns_raw);
  const double logp_nosplit = eff_ns * log_pi_k + total_right;

  
  auto logp0 = [&](int g) -> double {
    int sl = g - seg_start;
    if (sl < 1 || sl > n_seg - 1) return R_NegInf;
    double LL_left  = cum_new[sl - 1];
    double LL_right = total_right - cum_right[sl - 1];
    int eff_raw = (g - tau_k) - m_min; double eff = (eff_raw < 0) ? 0.0 : (double)eff_raw;
    double v = eff * log_pi_k + log_1mpi_k + log_poisson_split_penalty
               + log_q_atom + LL_left + LL_right;
    return std::isfinite(v) ? v : R_NegInf;
  };
  auto logpm = [&](int g) -> double {
    if (!has_left_neighbor) return R_NegInf;
    int sl = g - seg_start;
    if (sl < 1 || sl > n_seg - 1) return R_NegInf;
    double LL_left  = cum_km1[sl - 1];
    double LL_right = total_right - cum_right[sl - 1];
    double eff = (double)(g - tau_k);   // 지수 g-tau_k (m_min 할인 없음)
    double v = eff * log_pi_km1 + log_1mpi_km1 + LL_left + LL_right;
    return std::isfinite(v) ? v : R_NegInf;
  };
  
  auto logp_plus = [&](int ci_g) -> double {
    std::vector<double> terms;
    terms.reserve(2 * n_cand + 1);
    for (int d = ci_g + 1; d < n_cand; ++d) {
      int h = G_split[d];
      int sl = h - seg_start;
      if (sl < 1 || sl > n_seg - 1) continue;
      double LL_right_h = total_right - cum_right[sl - 1];
      
      double t_split = eff_ns * log_pi_k + log_1mpi_k + cum_new[sl - 1] + LL_right_h;
      if (std::isfinite(t_split)) terms.push_back(t_split);
      if (has_left_neighbor) {
        double eff_a = (double)(h - tau_k);
        double t_abs = eff_a * log_pi_km1 + log_1mpi_km1 + cum_km1[sl - 1] + LL_right_h;
        if (std::isfinite(t_abs)) terms.push_back(t_abs);
      }
    }
    if (std::isfinite(logp_nosplit)) terms.push_back(logp_nosplit);   
    if (terms.empty()) return R_NegInf;
    double M = terms[0];
    for (double t : terms) if (t > M) M = t;
    double acc = 0.0;
    for (double t : terms) acc += std::exp(t - M);
    return M + std::log(acc);
  };

  
  for (int ci = 0; ci < n_cand; ++ci) {
    int g = G_split[ci];
    bool terminal = (ci == n_cand - 1);
    double lp0 = logp0(g);
    double lpm = logpm(g);
    
    NumericVector cand(3);
    cand[0] = lpm;
    cand[1] = lp0;
    cand[2] = terminal ? logp_nosplit : logp_plus(ci);
    
    bool any_ok = std::isfinite(cand[0]) || std::isfinite(cand[1]) || std::isfinite(cand[2]);
    if (!any_ok) {
      if (terminal) return idx_nosplit; else continue;
    }
    int sel = sample_categorical_log_cpp(cand);   // 0=absorb,1=split,2=defer|nosplit
    if (sel == 1) return (ci + 1);                 // SPLIT at G_split[ci]
    if (sel == 0) return (n_cand + ci + 1);        // ABSORB at G_split[ci]
    
    if (terminal) return idx_nosplit;              // nosplit 종착
    // else defer → 다음 후보 g 로
  }
  return idx_nosplit;   // 끝까지 이연 → nosplit
}


// [[Rcpp::export]]
NumericVector merge_eval_cpp(NumericMatrix theta_mat,
                             IntegerVector I_k,             // 1-based time indices of the interval
                             NumericVector gamma_km1_vec, SEXP atom_km1_sexp, NumericVector v_km1_vec,
                             NumericVector gamma_k_vec,   SEXP atom_k_sexp,   NumericVector v_k_vec,
                             NumericVector gamma_kp1_vec, SEXP atom_kp1_sexp, NumericVector v_kp1_vec,
                             NumericMatrix Y_sub,
                             double alpha, double beta_c,
                             NumericMatrix phi_sub,
                             IntegerMatrix xi_sub,
                             int has_km1, int has_kp1,      // 1 if the neighbor label exists
                             NumericVector x_global) {
  NumericVector out(3);
  const int n = I_k.size();
  const int n_series = theta_mat.nrow();
  if (n == 0) { out[0]=0.0; out[1]=0.0; out[2]=0.0; return out; }

  AtomPars a_km1 = parse_atom(atom_km1_sexp);
  AtomPars a_k   = parse_atom(atom_k_sexp);
  AtomPars a_kp1 = parse_atom(atom_kp1_sexp);
  const double ab = alpha + beta_c;

  double ll_km1 = 0.0, ll_k = 0.0, ll_kp1 = 0.0;

  for (int s = 0; s < n_series; ++s) {
    NumericVector theta_j = theta_mat(s, _);
    double g_m  = (s < gamma_km1_vec.size()) ? gamma_km1_vec[s] : 0.0;
    double g_k  = (s < gamma_k_vec.size())   ? gamma_k_vec[s]   : 0.0;
    double g_p  = (s < gamma_kp1_vec.size()) ? gamma_kp1_vec[s] : 0.0;
    double v_m  = (s < v_km1_vec.size() && std::isfinite(v_km1_vec[s]) && v_km1_vec[s] > 0.0) ? v_km1_vec[s] : 1.0;
    double v_k  = (s < v_k_vec.size()   && std::isfinite(v_k_vec[s])   && v_k_vec[s]   > 0.0) ? v_k_vec[s]   : 1.0;
    double v_p  = (s < v_kp1_vec.size() && std::isfinite(v_kp1_vec[s]) && v_kp1_vec[s] > 0.0) ? v_kp1_vec[s] : 1.0;

    for (int i = 0; i < n; ++i) {
      int t0 = I_k[i] - 1;
      if (t0 < 0 || t0 >= Y_sub.ncol()) continue;
      double y = Y_sub(s, t0);

      double phi_val = (s < phi_sub.nrow() && t0 < phi_sub.ncol() && R_finite(phi_sub(s, t0)) && phi_sub(s, t0) >= 1.0)
        ? phi_sub(s, t0) : 1.0;
      int xi_val = (s < xi_sub.nrow() && t0 < xi_sub.ncol()) ? xi_sub(s, t0) : 0;
      double base_pow = std::pow(phi_val, static_cast<double>(xi_val));

      double f_k = shape_value_single(theta_j, a_k, t0, x_global);
      double mu_k = ab + g_k + f_k;
      double s2k = v_k * base_pow; if (!std::isfinite(s2k) || s2k <= 0.0) s2k = 1.0;
      ll_k += dnorm_log_var_s(y, mu_k, s2k);

      if (has_km1) {
        double f_m = shape_value_single(theta_j, a_km1, t0, x_global);
        double mu_m = ab + g_m + f_m;
        double s2m = v_m * base_pow; if (!std::isfinite(s2m) || s2m <= 0.0) s2m = 1.0;
        ll_km1 += dnorm_log_var_s(y, mu_m, s2m);
      }
      if (has_kp1) {
        double f_p = shape_value_single(theta_j, a_kp1, t0, x_global);
        double mu_p = ab + g_p + f_p;
        double s2p = v_p * base_pow; if (!std::isfinite(s2p) || s2p <= 0.0) s2p = 1.0;
        ll_kp1 += dnorm_log_var_s(y, mu_p, s2p);
      }
    }
  }

  out[0] = has_km1 ? ll_km1 : R_NegInf;
  out[1] = ll_k;
  out[2] = has_kp1 ? ll_kp1 : R_NegInf;
  return out;
}


// [[Rcpp::export]]
double compute_atom_loglik_cpp(SEXP atom_sexp,
                               NumericMatrix theta_mat,
                               NumericVector gamma_vec,
                               IntegerVector t_in_k_upper,
                               IntegerMatrix S_lower_sub,
                               int k,
                               NumericMatrix Y_sub,
                               double alpha,
                               double beta_c,
                               NumericVector v_vec,
                               NumericMatrix phi_sub,
                               IntegerMatrix xi_sub,
                               NumericVector x_global) {
  int n_series = theta_mat.nrow();
  AtomPars atom = parse_atom(atom_sexp);
  double log_lik = 0.0;

  for (int s = 0; s < n_series; ++s) {
    NumericVector theta_j = theta_mat(s, _);
    double gamma_jk = (s < gamma_vec.size()) ? gamma_vec[s] : 0.0;
    double v_jk = (s < v_vec.size() && std::isfinite(v_vec[s]) && v_vec[s] > 0.0) ? v_vec[s] : 1.0;

    for (int i = 0; i < t_in_k_upper.size(); ++i) {
      int t0 = t_in_k_upper[i] - 1;
      if (t0 < 0 || t0 >= S_lower_sub.ncol()) continue;
      if (S_lower_sub(s, t0) != k) continue;

      double mu_val = alpha + beta_c + gamma_jk + shape_value_single(theta_j, atom, t0, x_global);
      double phi_val = (s < phi_sub.nrow() && t0 < phi_sub.ncol() && R_finite(phi_sub(s, t0)) && phi_sub(s, t0) >= 1.0)
        ? phi_sub(s, t0) : 1.0;
      int xi_val = (s < xi_sub.nrow() && t0 < xi_sub.ncol()) ? xi_sub(s, t0) : 0;
      double sig2 = v_jk * std::pow(phi_val, static_cast<double>(xi_val));
      if (!std::isfinite(sig2) || sig2 <= 0.0) sig2 = 1.0;

      log_lik += dnorm_log_var_s(Y_sub(s, t0), mu_val, sig2);
    }
  }

  if (!std::isfinite(log_lik)) log_lik = -1e300;
  return log_lik;
}



// [[Rcpp::export]]
NumericVector compute_atoms_loglik_batch_cpp(NumericMatrix atom_mat,
                                             NumericMatrix theta_mat,
                                             NumericVector gamma_vec,
                                             IntegerVector t_in_k_upper,
                                             IntegerMatrix S_lower_sub,
                                             int k,
                                             NumericMatrix Y_sub,
                                             double alpha,
                                             double beta_c,
                                             NumericVector v_vec,
                                             NumericMatrix phi_sub,
                                             IntegerMatrix xi_sub,
                                             NumericVector x_global) {
  const int N = atom_mat.nrow();
  const int n_series = theta_mat.nrow();
  NumericVector out(N, 0.0);
  if (N == 0) return out;

  // Atom-only precomputation (matches safe_sign + linear_coef in
  // shape_value_single; independent of series/time, so done once).
  std::vector<double> lin(N), c2(N);   // lin = linear_coef, c2 = delta2*shape_beta
  for (int n = 0; n < N; ++n) {
    double delta1 = safe_sign(atom_mat(n, 0));
    double delta2 = safe_sign(atom_mat(n, 1));
    double b      = atom_mat(n, 2);
    double g      = atom_mat(n, 3);
    
    if (!std::isfinite(b) || b <= 0.0) b = 1.0;
    if (!std::isfinite(g) || g <  0.0) g = 0.0;
    lin[n] = delta1 * g + b * (delta1 + delta2) / 2.0;
    c2[n]  = delta2 * b;
  }

  for (int s = 0; s < n_series; ++s) {
    NumericVector theta_j = theta_mat(s, _);
    const bool theta_ok = (theta_j.size() > 1);
    double gamma_jk = (s < gamma_vec.size()) ? gamma_vec[s] : 0.0;
    double v_jk = (s < v_vec.size() && std::isfinite(v_vec[s]) && v_vec[s] > 0.0) ? v_vec[s] : 1.0;

    // SZ depends only on theta_j (atom- and time-independent within a series).
    const double SZ = theta_sq_sum(theta_j);
    const double SZ_safe = (SZ < 1e-12) ? 1e-12 : SZ;
    const double base_mu_s = alpha + beta_c + gamma_jk;

    for (int i = 0; i < t_in_k_upper.size(); ++i) {
      int t0 = t_in_k_upper[i] - 1;
      if (t0 < 0 || t0 >= S_lower_sub.ncol()) continue;
      if (S_lower_sub(s, t0) != k) continue;

      // H = DI/SZ_safe : computed ONCE, reused for every atom (== the value
      // shape_value_single would compute for each atom at this (s,t0)).
      const double H_val = theta_ok ? (compute_DI_single(theta_j, t0) / SZ_safe) : 0.0;
      const double x_t = x_global[t0];

      double phi_val = (s < phi_sub.nrow() && t0 < phi_sub.ncol() && R_finite(phi_sub(s, t0)) && phi_sub(s, t0) >= 1.0)
        ? phi_sub(s, t0) : 1.0;
      int xi_val = (s < xi_sub.nrow() && t0 < xi_sub.ncol()) ? xi_sub(s, t0) : 0;
      double sig2 = v_jk * std::pow(phi_val, static_cast<double>(xi_val));
      if (!std::isfinite(sig2) || sig2 <= 0.0) sig2 = 1.0;

      const double yv = Y_sub(s, t0);
      for (int n = 0; n < N; ++n) {
        // shape_value_single(theta_j, atom_n, t0, x_global) with reused H_val:
        double shape_n = theta_ok ? (lin[n] * x_t - c2[n] * H_val) : 0.0;
        double mu_n = base_mu_s + shape_n;
        out[n] += dnorm_log_var_s(yv, mu_n, sig2);
      }
    }
  }

  for (int n = 0; n < N; ++n) {
    if (!std::isfinite(out[n])) out[n] = -1e300;
  }
  return out;
}


// [[Rcpp::export]]
NumericVector sample_state_intercepts_independent_cpp(
    int j_1based,
    IntegerVector S_lower_j,
    int K_c,
    NumericVector Y_j,
    double alpha,
    double beta_c,
    NumericVector f_j,
    NumericVector sigma2_j,
    double sigma2_gamma)
{
    NumericVector gamma_new(K_c, 0.0);
    NumericVector obs_prec(K_c, 0.0);
    NumericVector obs_mean(K_c, 0.0);

    // 1. 각 구간(k)별로 우도의 정밀도와 평균을 누적
    for (int t = 0; t < S_lower_j.size(); ++t) {
        int k = S_lower_j[t] - 1; // 0-based index
        if (k >= 0 && k < K_c) {
            double inv_sig2 = 1.0 / sigma2_j[t];
            obs_prec[k] += inv_sig2;
            
            double resid = Y_j[t] - alpha - beta_c - f_j[t];
            obs_mean[k] += resid * inv_sig2;
        }
    }

    // 2. 완벽한 독립 사전분포 기반 사후분포 추출 (행렬 역산 O(K^3) -> O(K)로 최적화)
    for (int k = 0; k < K_c; ++k) {
        double V_post = 1.0 / (1.0 / sigma2_gamma + obs_prec[k]);
        double m_post = V_post * obs_mean[k];
        gamma_new[k] = R::rnorm(m_post, std::sqrt(V_post));
    }
    return gamma_new;
}



NumericVector compute_H_at_times_cpp(NumericVector theta_j,
                                     IntegerVector t_indices) {
  NumericVector out(t_indices.size(), 0.0);
  if (theta_j.size() <= 1) return out;
  const double SZ = theta_sq_sum(theta_j);
  const double SZ_safe = (SZ < 1e-12) ? 1e-12 : SZ;
  const int Tslices = (int)G_DATX.n_slices;
  for (int i = 0; i < t_indices.size(); ++i) {
    int t0 = t_indices[i] - 1;
    if (t0 < 0 || t0 >= Tslices) continue;
    out[i] = compute_DI_single(theta_j, t0) / SZ_safe;
  }
  return out;
}


NumericVector rtrunc_invgamma_lower1_cpp(NumericVector shape, NumericVector scale) {
  int n = shape.size();
  NumericVector out(n, 1.0);
  for (int i = 0; i < n; ++i) {
    double a = shape[i], b = scale[i];
    if (!std::isfinite(a) || !std::isfinite(b) || a <= 0.0 || b <= 0.0) { out[i] = 1.0; continue; }
    double gscale = 1.0 / b;                          // gamma scale = 1/rate
    double F1 = R::pgamma(1.0, a, gscale, 1, 0);      // P(Y <= 1)
    if (!std::isfinite(F1) || F1 <= 1e-300) { out[i] = 1.0; continue; }  // ~no truncated mass
    double u = R::runif(0.0, F1);
    double y = R::qgamma(u, a, gscale, 1, 0);
    if (!std::isfinite(y) || y <= 0.0) { out[i] = 1.0; continue; }
    double x = 1.0 / y;
    out[i] = (std::isfinite(x) && x >= 1.0) ? x : 1.0;
  }
  return out;
}



// [[Rcpp::export]]
NumericVector profile_atom_scores_cpp(NumericVector x_all,
                                      NumericVector H_all,
                                      NumericVector rc_all,
                                      NumericVector w_all,
                                      IntegerVector seg_start,  // 0-based 시작 인덱스
                                      IntegerVector seg_len,
                                      NumericVector W_j,
                                      double delta1,
                                      double delta2,
                                      NumericVector cand_b,
                                      NumericVector cand_g) {
  const int n_series = seg_start.size();
  const int n_cand = cand_b.size();
  const double half_sum = (delta1 + delta2) / 2.0;

  
  double T_rcrc = 0.0, T_Arc = 0.0, T_Brc = 0.0, T_AA = 0.0, T_AB = 0.0, T_BB = 0.0;

  for (int s = 0; s < n_series; ++s) {
    const int st = seg_start[s];
    const int L  = seg_len[s];
    if (L <= 0) continue;
    const double Wj = (W_j[s] > 0 && std::isfinite(W_j[s])) ? W_j[s] : 1.0;

    double sw_zg = 0.0, sw_zb = 0.0;
    for (int t = 0; t < L; ++t) {
      const double w = w_all[st + t];
      const double xv = x_all[st + t];
      const double Hv = H_all[st + t];
      const double zg = delta1 * xv;
      const double zb = half_sum * xv - delta2 * Hv;
      sw_zg += w * zg;
      sw_zb += w * zb;
    }
    const double cg = sw_zg / Wj;   
    const double cb = sw_zb / Wj;   

    for (int t = 0; t < L; ++t) {
      const double w  = w_all[st + t];
      const double xv = x_all[st + t];
      const double Hv = H_all[st + t];
      const double rc = rc_all[st + t];
      const double zg = delta1 * xv;
      const double zb = half_sum * xv - delta2 * Hv;
      const double Ac = zg - cg;
      const double Bc = zb - cb;
      T_rcrc += w * rc * rc;
      T_Arc  += w * Ac * rc;
      T_Brc  += w * Bc * rc;
      T_AA   += w * Ac * Ac;
      T_AB   += w * Ac * Bc;
      T_BB   += w * Bc * Bc;
    }
  }

  NumericVector log_scores(n_cand);
  for (int i = 0; i < n_cand; ++i) {
    const double g = cand_g[i];
    const double b = cand_b[i];
    const double quad = T_rcrc - 2.0*g*T_Arc - 2.0*b*T_Brc
                        + g*g*T_AA + 2.0*g*b*T_AB + b*b*T_BB;
    log_scores[i] = -0.5 * quad;
  }
  return log_scores;
}
