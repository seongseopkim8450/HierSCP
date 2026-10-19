# OECD CLI application: data, model fit, figures and convergence diagnostics

library(dplyr)
library(lubridate)
library(RcppArmadillo)
library(ggplot2)
library(Rcpp)

# YOUR PATH to the code: MUST CHANGE IT TO YOUR OWN
setwd("C://Users/limit/Desktop/HierSCP_code_2/")  
###########

source("corefunctions_code/00_utils.R")
source("corefunctions_code/01_precompute.R")
source("corefunctions_code/02_shape_function.R")
source("corefunctions_code/03_state_management.R")
source("corefunctions_code/04_sampling.R")
source("corefunctions_code/05_likelihood.R")
source("corefunctions_code/06_local_adjusting.R")
source("corefunctions_code/07_interval_adjusting.R")
source("corefunctions_code/08_param_update_continuous.R")
source("corefunctions_code/10_initialization.R")
source("corefunctions_code/11_mcmc_main.R")
sourceCpp("corefunctions_code/hiercpd_core.cpp")
source("corefunctions_code/15_rcpp_bridge.R")
source("corefunctions_code/compute_atom_loglik_fallback.R")

MODEL_DIR <- "realdata_analysis"  # folder with OECD_CLI.csv
csv_path  <- "OECD_CLI.csv"

setwd(MODEL_DIR)
cli_data<-read.csv("OECD_CLI.csv")

unique(cli_data$Reference.area)

# Sample window
start_ym  <- c(2000, 1)
end_ym    <- c(2024, 12)

# Economies and clusters
areas         <- c("USA","CAN","GBR","DEU","FRA","ITA","JPN","KOR","CHN")
group_members <- list(c(1,2), c(3,4,5,6), c(7,8,9))
group_names   <- c("North America","Major European Countries","East Asian Countries")

# Sampler settings
SEED <- 2026
K_init <- 5
K_min <- 2
K_max <- 12
M      <- 8
m_min  <- 9 # 9 months(3 quarters)
lambda0_K <- 1
n_iter <- 100000; n_burnin <- 50000; n_thin <- 2
N_launch <- 1; n_warmup <- 100; ia_every <- 20; verbose_mcmc <- 150

OVERLAY_RAW <- TRUE  # observed CLI behind the fitted curves

# Output files
DATA_RDS  <- "cli_data.rds"
FIT_RDS   <- "cli_fit.rds"
PDF_UPPER <- "HierSCP_CLI_upper.pdf"
PDF_LOWER <- "HierSCP_CLI_lower_G%d.pdf"

# Steps to run
DO_BUILD <- TRUE
DO_FIT   <- TRUE
DO_PLOT  <- TRUE


## month index t (t = 1 is 2000-01)
t_of <- function(year, month) (year - 2000) * 12 + month

# State of saved draw b: paths, regime counts and clusters
reconstruct_state_from_sample <- function(b, samples, model) {
  J <- as.integer(model$J)
  C <- as.integer(model$C)
  T_len <- as.integer(model$T_len)
  
  if (length(b) != 1L || !is.finite(b) || b != as.integer(b)) {
    stop("reconstruct_state_from_sample: b must be a single integer index.", call. = FALSE)
  }
  b <- as.integer(b)
  B <- if (!is.null(samples$n_saved)) as.integer(samples$n_saved) else dim(samples$S_lower)[1L]
  if (b < 1L || b > B) {
    stop("reconstruct_state_from_sample: sample index b is out of range.", call. = FALSE)
  }
  
  required <- c("S_lower", "S_upper", "K", "cluster")
  miss <- required[vapply(required, function(nm) is.null(samples[[nm]]), logical(1))]
  if (length(miss)) {
    stop(sprintf("reconstruct_state_from_sample: missing stored field(s): %s",
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  
  if (!identical(dim(samples$S_lower)[2:3], c(J, T_len))) {
    stop("reconstruct_state_from_sample: S_lower dimensions do not match model.", call. = FALSE)
  }
  if (!identical(dim(samples$S_upper)[2:3], c(C, T_len))) {
    stop("reconstruct_state_from_sample: S_upper dimensions do not match model.", call. = FALSE)
  }
  
  # explicit dims keep the shape when J or C is 1
  S_lower <- array(samples$S_lower[b, , , drop = FALSE], dim = c(J, T_len))
  S_upper <- array(samples$S_upper[b, , , drop = FALSE], dim = c(C, T_len))
  storage.mode(S_lower) <- "integer"
  storage.mode(S_upper) <- "integer"
  
  K <- as.integer(samples$K[b, seq_len(C), drop = TRUE])
  cluster <- as.integer(samples$cluster[b, seq_len(J), drop = TRUE])
  
  if (length(K) != C || anyNA(K) || any(K < 1L)) {
    stop("reconstruct_state_from_sample: invalid stored K vector.", call. = FALSE)
  }
  if (length(cluster) != J || anyNA(cluster) || any(cluster < 1L) || any(cluster > C)) {
    stop("reconstruct_state_from_sample: invalid stored cluster vector.", call. = FALSE)
  }
  if (anyNA(S_lower) || anyNA(S_upper)) {
    stop("reconstruct_state_from_sample: stored state arrays contain NA.", call. = FALSE)
  }
  
  state <- list(
    S_lower = S_lower,
    S_upper = S_upper,
    K = K,
    cluster = cluster
  )
  state$tau_upper <- lapply(seq_len(C), function(cc) extract_changepoints(state$S_upper[cc, ]))
  state$tau_lower <- lapply(seq_len(J), function(j)  extract_changepoints(state$S_lower[j, ]))
  state
}

# Parameters of saved draw b (for summaries; not a restart point)
reconstruct_params_from_sample <- function(b, samples, model) {
  J <- as.integer(model$J)
  C <- as.integer(model$C)
  T_len <- as.integer(model$T_len)
  
  if (length(b) != 1L || !is.finite(b) || b != as.integer(b)) {
    stop("reconstruct_params_from_sample: b must be a single integer index.", call. = FALSE)
  }
  b <- as.integer(b)
  B <- if (!is.null(samples$n_saved)) as.integer(samples$n_saved) else length(samples$alpha)
  if (b < 1L || b > B) {
    stop("reconstruct_params_from_sample: sample index b is out of range.", call. = FALSE)
  }
  
  required <- c("alpha", "beta", "theta", "gamma_store", "atoms_store")
  miss <- required[vapply(required, function(nm) is.null(samples[[nm]]), logical(1))]
  if (length(miss)) {
    stop(sprintf("reconstruct_params_from_sample: missing stored field(s): %s",
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  
  alpha <- as.numeric(samples$alpha[b])
  beta <- as.numeric(samples$beta[b, seq_len(C), drop = TRUE])
  if (length(alpha) != 1L || !is.finite(alpha)) {
    stop("reconstruct_params_from_sample: invalid stored alpha.", call. = FALSE)
  }
  if (length(beta) != C || any(!is.finite(beta))) {
    stop("reconstruct_params_from_sample: invalid stored beta vector.", call. = FALSE)
  }
  
  theta_dim <- dim(samples$theta)
  if (length(theta_dim) != 3L || theta_dim[2L] != J) {
    stop("reconstruct_params_from_sample: theta dimensions do not match model.", call. = FALSE)
  }
  theta <- array(samples$theta[b, , , drop = FALSE], dim = c(J, theta_dim[3L]))
  storage.mode(theta) <- "double"
  if (any(!is.finite(theta))) {
    stop("reconstruct_params_from_sample: stored theta contains non-finite values.", call. = FALSE)
  }
  
  gamma <- samples$gamma_store[[b]]
  if (is.null(gamma) || length(gamma) != J) {
    stop("reconstruct_params_from_sample: stored gamma list has wrong length.", call. = FALSE)
  }
  gamma <- lapply(seq_len(J), function(j) {
    gj <- as.numeric(gamma[[j]])
    if (!length(gj) || any(!is.finite(gj))) {
      stop(sprintf("reconstruct_params_from_sample: invalid gamma for series %d.", j),
           call. = FALSE)
    }
    gj
  })
  
  atoms <- samples$atoms_store[[b]]
  if (is.null(atoms) || length(atoms) != C) {
    stop("reconstruct_params_from_sample: stored atom list has wrong cluster dimension.",
         call. = FALSE)
  }
  
  params <- list(
    alpha = alpha,
    beta = beta,
    theta = theta,
    gamma = gamma,
    atoms = atoms
  )
  
  # optional fields, when stored
  if (!is.null(samples$v)) {
    vv <- as.numeric(samples$v[b, seq_len(J), drop = TRUE])
    if (length(vv) == J && all(is.finite(vv)) && all(vv > 0)) {
      params$v <- lapply(vv, function(x) as.numeric(x))
    }
  }
  if (!is.null(samples$phi)) {
    ph <- array(samples$phi[b, , , drop = FALSE], dim = c(J, T_len))
    storage.mode(ph) <- "double"
    params$phi <- ph
  }
  if (!is.null(samples$xi_store)) {
    xi <- array(samples$xi_store[b, , , drop = FALSE], dim = c(J, T_len))
    storage.mode(xi) <- "integer"
    params$xi <- xi
  }
  if (!is.null(samples$theta_out)) params$theta_out <- as.numeric(samples$theta_out[b])
  if (!is.null(samples$nu)) params$nu <- as.numeric(samples$nu[b])
  if (!is.null(samples$sigma2_gamma)) params$sigma2_gamma <- as.numeric(samples$sigma2_gamma[b])
  
  params
}


# Plotting helpers (190 mm wide, 9 pt, vector PDF)
.FIG_WIDTH  <- 7.48
.FIG_PT     <- 9

# Okabe-Ito colors, one per shape class
.SP_COL  <- c("++" = "#E69F00",
              "+-" = "#CC79A7",
              "-+" = "#009E73",
              "--" = "#0072B2")
.SP_NAME <- c("++" = "(+,+) increasing, convex", "+-" = "(+,-) increasing, concave",
              "-+" = "(-,+) decreasing, convex", "--" = "(-,-) decreasing, concave")
.RAW_COL <- "grey60"
.CI_COL  <- "grey35"  # changepoint intervals

# Diverging scale for the jump at a changepoint: red up, indigo down, square-root scaled
.JUMP_RAMP <- colorRampPalette(c("#110847", "#24187A", "#3F33A3", "#7068C4", "#A9A6DA",
                                 "#CFCFCF",
                                 "#F5B3A3", "#EB6A4E", "#CB1F1F", "#8E0A12", "#550008"),
                               space = "Lab")
.JUMP_PAL  <- .JUMP_RAMP(101)
.jump_u    <- function(d, lim) sign(d) * sqrt(pmin(1, abs(d) / lim))
.jump_col  <- function(d, lim) .JUMP_PAL[round((.jump_u(d, lim) + 1) * 50) + 1]

COUNTRY_NAMES <- c(USA = "United States", CAN = "Canada", GBR = "United Kingdom",
                   DEU = "Germany", FRA = "France", ITA = "Italy",
                   JPN = "Japan", KOR = "Korea", CHN = "China")

.sp_key <- function(d1, d2) paste0(if (d1 > 0) "+" else "-", if (d2 > 0) "+" else "-")
.sp_lab <- function(d1, d2) sprintf("(%s,%s)", if (d1 > 0) "+" else "-", if (d2 > 0) "+" else "-")
.alpha_col <- function(hex, a) adjustcolor(hex, alpha.f = a)
.cp_pick   <- function(primary, fallback) if (!is.null(primary) && length(primary) > 0) primary else fallback
.spread <- function(y, gap) {  # keep direct labels apart
  o <- order(y); ys <- y[o]
  if (length(ys) > 1) for (i in 2:length(ys)) if (ys[i] - ys[i - 1] < gap) ys[i] <- ys[i - 1] + gap
  out <- y; out[o] <- ys; out
}

# Regimes of one panel: sign pairs, boundaries and 95% intervals
.seg_struct <- function(post, g, cp_list, T_len, min_seg_len = 2) {
  K  <- post$K_hat[g]
  d1 <- d2 <- integer(K)
  for (k in seq_len(K)) {
    sc <- post$shape_class[[g]][[k]]
    d1[k] <- if (!is.null(sc$d1_sign)) as.integer(sc$d1_sign) else 1L
    d2[k] <- if (!is.null(sc$d2_sign)) as.integer(sc$d2_sign) else 1L
  }
  starts <- rep(NA_real_, K); starts[1] <- 1
  cp_loc <- cp_lo <- cp_hi <- rep(NA_real_, K)
  if (!is.null(cp_list) && length(cp_list) > 0) {
    for (info in cp_list) {
      if (is.null(info) || is.null(info$k)) next
      k <- as.integer(info$k); if (k < 2L || k > K) next
      loc <- if (!is.null(info$point_est)) info$point_est else info$median
      starts[k] <- loc; cp_loc[k] <- loc
      if (!is.null(info$lo)) cp_lo[k] <- as.numeric(info$lo)
      if (!is.null(info$hi)) cp_hi[k] <- as.numeric(info$hi)
    }
  }
  for (k in seq_len(K)[-1]) if (is.na(starts[k])) starts[k] <- starts[k - 1]
  starts <- pmin(pmax(round(starts), 1), T_len)
  for (k in seq_len(K)[-1]) if (starts[k] < starts[k - 1]) starts[k] <- starts[k - 1]
  ends <- c(starts[-1] - 1, T_len); ends <- pmin(pmax(ends, starts), T_len)

  segs <- lapply(seq_len(K), function(k)
    list(d1 = d1[k], d2 = d2[k], start = starts[k], end = ends[k],
         cp = cp_loc[k], lo = cp_lo[k], hi = cp_hi[k]))
  repeat {  # merge slivers created by rounding
    if (length(segs) <= 1) break
    w   <- vapply(segs, function(s) s$end - s$start + 1, numeric(1))
    bad <- which(w < min_seg_len)
    if (length(bad) == 0) break
    k <- bad[1]
    if (k == 1) {
      segs[[2]]$start <- segs[[1]]$start
      segs[[2]]$cp <- segs[[2]]$lo <- segs[[2]]$hi <- NA_real_
      segs[[1]] <- NULL
    } else {
      segs[[k - 1]]$end <- segs[[k]]$end
      segs[[k]] <- NULL
    }
  }
  f <- function(nm) vapply(segs, function(s) as.numeric(s[[nm]]), numeric(1))
  list(K = length(segs), starts = f("start"), ends = f("end"), d1 = f("d1"), d2 = f("d2"),
       cp_loc = f("cp"), cp_lo = f("lo"), cp_hi = f("hi"))
}

.seg_of_x <- function(x, starts, K) {
  se <- c(starts, Inf); idx <- rep(K, length(x))
  for (k in seq_len(K)) idx[x >= se[k] & x < se[k + 1]] <- k
  idx[x < starts[1]] <- 1L
  idx
}

# Jump of the fitted curve at each changepoint: both sides extrapolated to the
# boundary by local line fits that skip the h months next to it
.panel_jumps <- function(x, mu, S, h = 2, w = 4) {
  t_all <- seq_len(max(S$ends))
  m_t <- approx(x, mu, xout = t_all, rule = 2)$y
  lev_fit <- function(tt, x0) {
    tt <- tt[tt >= 1 & tt <= length(m_t)]
    if (length(tt) == 0) return(NA_real_)
    if (length(tt) == 1) return(m_t[tt])
    b <- stats::coef(stats::lm(m_t[tt] ~ tt)); unname(b[1] + b[2] * x0)
  }
  K <- S$K; before <- after <- jump <- rep(NA_real_, K)
  for (k in seq_len(K)[-1]) {
    tau <- if (is.finite(S$cp_loc[k])) round(S$cp_loc[k]) else S$starts[k]
    x0  <- tau - 0.5
    Lk1 <- tau - S$starts[k - 1]; Lk <- S$ends[k] - tau + 1
    hl  <- min(h, max(0, floor((Lk1 - 2) / 2))); hr <- min(h, max(0, floor((Lk - 2) / 2)))
    a1 <- max(S$starts[k - 1] + hl, tau - hl - w); b1 <- tau - hl - 1
    a2 <- tau + hr;                               b2 <- min(S$ends[k] - hr, tau + hr + w - 1)
    left  <- if (a1 <= b1) a1:b1 else integer(0)
    right <- if (a2 <= b2) a2:b2 else integer(0)
    before[k] <- lev_fit(left, x0); after[k] <- lev_fit(right, x0)
    jump[k] <- after[k] - before[k]
  }
  list(before = before, after = after, jump = jump)
}

# One panel: curve and band by shape class, changepoints colored by the jump
.draw_panel <- function(x, mu, lo, hi, S, meta, raw_list = NULL, raw_labels = NULL,
                        panel_label = "", ylab = "CLI", xlab_on = TRUE) {
  T_len <- meta$T_len
  J  <- .panel_jumps(x, mu, S)
  jl <- if (!is.null(meta$jump_lim)) meta$jump_lim else max(1, abs(J$jump), na.rm = TRUE)
  yr <- range(c(mu, lo, hi, unlist(raw_list), meta$ref_line), na.rm = TRUE)
  pad <- 0.04 * diff(yr); strip_h <- 0.11 * diff(yr); tri_h <- 0.07 * diff(yr)
  ylim <- c(yr[1] - pad, yr[2] + pad + tri_h + strip_h)
  plot(NA, xlim = c(0.5, T_len + 0.5), ylim = ylim, xlab = "", ylab = "",
       axes = FALSE, xaxs = "i", yaxs = "i")
  usr <- par("usr"); strip_lo <- usr[4] - strip_h; tri_y <- strip_lo - tri_h / 2

  # CLI = 100 and the observed series
  if (!is.null(meta$ref_line)) abline(h = meta$ref_line, lty = 3, col = "grey30", lwd = 0.8)
  if (!is.null(raw_list)) for (z in raw_list) lines(seq_len(T_len), z, col = .RAW_COL, lwd = 0.6)

  # changepoints: solid for class, dashed for intensity transitions
  for (k in seq_len(S$K)[-1]) {
    if (!is.finite(S$cp_loc[k])) next
    class_change <- (S$d1[k] != S$d1[k - 1]) || (S$d2[k] != S$d2[k - 1])
    segments(S$cp_loc[k], usr[3], S$cp_loc[k], strip_lo, col = .alpha_col(.jump_col(J$jump[k], jl), 0.75),
             lty = if (class_change) 1 else 2, lwd = if (class_change) 1.0 else 0.9)
    if (is.finite(S$cp_lo[k]) && is.finite(S$cp_hi[k])) {
      yb <- usr[3] + 0.035 * diff(usr[3:4])
      segments(S$cp_lo[k], yb, S$cp_hi[k], yb, col = .alpha_col(.CI_COL, 0.7), lwd = 2.0, lend = 1)
    }
  }

  # curve and band regime by regime, not joined across changepoints
  seg_x <- .seg_of_x(x, S$starts, S$K)
  key_k <- vapply(seq_len(S$K), function(k) .sp_key(S$d1[k], S$d2[k]), character(1))
  for (k in seq_len(S$K)) {
    in_k <- (seg_x == k)
    sel  <- in_k & is.finite(lo) & is.finite(hi)
    if (sum(sel) >= 2) {
      xx <- x[sel]
      polygon(c(xx, rev(xx)), c(lo[sel], rev(hi[sel])),
              col = .alpha_col(.SP_COL[key_k[k]], 0.30), border = NA)
    }
    sel <- in_k & is.finite(mu)
    if (sum(sel) >= 2) lines(x[sel], mu[sel], col = .SP_COL[key_k[k]], lwd = 1.7)
  }

  # jump arrows and markers
  for (k in seq_len(S$K)[-1]) {
    if (!is.finite(S$cp_loc[k]) || !is.finite(J$jump[k])) next
    cx <- S$cp_loc[k] - 0.5; col_k <- .jump_col(J$jump[k], jl)
    if (abs(J$jump[k]) > 0.025 * diff(usr[3:4])) {
      arrows(cx, J$before[k], cx, J$after[k], length = 0.05, angle = 25,
             col = col_k, lwd = 2.4)
    }
    points(S$cp_loc[k], tri_y, pch = if (J$jump[k] >= 0) 24 else 25,
           bg = col_k, col = "grey20", lwd = 0.4,
           cex = 0.55 + 0.55 * min(1, abs(J$jump[k]) / jl))
  }

  # class strip, one label per run of equal classes
  for (k in seq_len(S$K)) {
    xr <- if (k < S$K) S$starts[k + 1] - 0.5 else T_len + 0.5
    rect(S$starts[k] - 0.5, strip_lo, xr, usr[4], col = .alpha_col(.SP_COL[key_k[k]], 0.75),
         border = "white", lwd = 0.6)
  }
  run_start <- 1
  for (k in seq_len(S$K)) {
    last_in_run <- (k == S$K) || key_k[k + 1] != key_k[k]
    if (last_in_run) {
      a <- S$starts[run_start]; b <- S$ends[k]
      if (b - a + 1 >= 16) text((a + b) / 2, (strip_lo + usr[4]) / 2, .sp_lab(S$d1[k], S$d2[k]),
                                cex = 0.8, font = 2, col = "black")
      run_start <- k + 1
    }
  }
  abline(h = strip_lo, col = "grey40", lwd = 0.5)

  # member labels
  if (!is.null(raw_list) && !is.null(raw_labels)) {
    yl <- .spread(vapply(raw_list, function(z) z[T_len], numeric(1)), 0.07 * diff(usr[3:4]))
    text(T_len + 2, yl, raw_labels, adj = 0, cex = 0.68, col = "grey30", xpd = NA)
  }

  yrs_all <- unique(meta$months %/% 12)
  minor_at <- vapply(yrs_all, function(Y0) which(meta$months %/% 12 == Y0)[1], integer(1))
  axis(1, at = minor_at, labels = FALSE, tcl = -0.15, lwd = 0, lwd.ticks = 0.6)
  axis(1, at = meta$tick_at, labels = if (xlab_on) meta$tick_years else FALSE,
       tcl = -0.35, lwd = 0, lwd.ticks = 0.8, mgp = c(2, 0.35, 0), cex.axis = 0.9)
  yt <- pretty(yr); yt <- yt[yt >= yr[1] - pad & yt <= yr[2] + pad]
  axis(2, at = yt, las = 1, tcl = -0.3, mgp = c(2, 0.5, 0), cex.axis = 0.9, lwd = 0, lwd.ticks = 0.8)
  mtext(ylab, side = 2, line = 2.3, cex = 0.85)
  mtext(panel_label, side = 3, line = 0.15, adj = 0, cex = 0.9, font = 2)
  box(lwd = 0.6)
}

# Legend and color bar for the jump
.draw_legend <- function(jump_lim = 1) {
  op <- par(mar = c(0, 0.5, 0, 0.5)); on.exit(par(op)); plot.new()
  lab <- c(unname(.SP_NAME), "95% pointwise credible band", "observed CLI",
           "class transition", "intensity transition", "95% interval of changepoint",
           "CLI = 100 (long-term trend)")
  fill <- c(.alpha_col(.SP_COL, 0.75), .alpha_col("grey40", 0.30), NA, NA, NA, NA, NA)
  bord <- c(rep("grey30", 4), NA, NA, NA, NA, NA, NA)
  lty  <- c(rep(NA, 5), 1, 1, 2, 1, 3)
  lwd  <- c(rep(NA, 5), 0.8, 1.3, 1.1, 2.0, 0.8)
  col  <- c(rep(NA, 5), .RAW_COL, "grey25", "grey25", .alpha_col(.CI_COL, 0.7), "grey30")
  legend(x = 0, y = 0.5, xjust = 0, yjust = 0.5, legend = lab, fill = fill, border = bord,
         lty = lty, lwd = lwd, col = col, ncol = 3, bty = "n", cex = 0.82,
         x.intersp = 0.6, seg.len = 1.6)

  x0 <- 0.775; x1 <- 0.975; y0 <- 0.42; y1 <- 0.62
  nb <- length(.JUMP_PAL)
  rasterImage(as.raster(matrix(.JUMP_PAL, nrow = 1)), x0, y0, x1, y1, interpolate = FALSE)
  rect(x0, y0, x1, y1, border = "grey30", lwd = 0.5)
  tp <- c(1, 3)[c(1, 3) < 0.8 * jump_lim]
  tk <- c(-jump_lim, -rev(tp), 0, tp, jump_lim)
  tx <- x0 + (.jump_u(tk, jump_lim) + 1) / 2 * (x1 - x0)
  segments(tx, y0, tx, y0 - 0.04, lwd = 0.5)
  lab_tk <- formatC(tk, format = "fg", digits = 3, flag = "+"); lab_tk[tk == 0] <- "0"
  text(tx, y0 - 0.13, lab_tk, cex = 0.72)
  text((x0 + x1) / 2, y1 + 0.13, "jump at changepoint (CLI points)", cex = 0.78)
  points(x0 - 0.012, (y0 + y1) / 2, pch = 25, bg = .JUMP_PAL[1], col = "grey20", lwd = 0.4, cex = 0.9)
  points(x1 + 0.012, (y0 + y1) / 2, pch = 24, bg = .JUMP_PAL[nb], col = "grey20", lwd = 0.4, cex = 0.9)
}

.open_pdf <- function(file, height) {
  if (!is.null(file)) pdf(file, width = .FIG_WIDTH, height = height, pointsize = .FIG_PT,
                          family = "Helvetica", useDingbats = FALSE)
}

# Cluster-level figure, one panel per cluster
plot_hierscp_upper <- function(post, meta, save_pdf = NULL, min_seg_len = 2) {
  C <- length(post$K_hat); T_len <- meta$T_len; Y <- meta$Y_for_overlay
  .open_pdf(save_pdf, height = 2.05 * C + 0.95)
  layout(matrix(seq_len(C + 1), ncol = 1), heights = c(rep(1, C), 0.95 / 2.05))
  for (g in seq_len(C)) {
    par(mar = c(if (g == C) 1.6 else 0.6, 3.4, 1.4, 2.4))
    S  <- .seg_struct(post, g, .cp_pick(post$pma_upper[[g]], post$cp_summary[[g]]),
                      T_len, min_seg_len = min_seg_len)
    uc <- post$upper_curves[[g]]; members <- meta$group_members[[g]]
    raw_list <- if (!is.null(Y)) lapply(members, function(j) Y[j, ]) else NULL
    lbl <- sprintf("(%s) %s", letters[g], meta$group_names[g])
    .draw_panel(seq_len(T_len), uc$mean, uc$lo, uc$hi, S, meta,
                raw_list = raw_list, raw_labels = meta$series_names[members],
                panel_label = lbl, ylab = "CLI", xlab_on = (g == C))
  }
  .draw_legend(meta$jump_lim)
  if (!is.null(save_pdf)) dev.off()
}

# Fitted curve of series j
.lower_curve <- function(post, j, T_len) {
  sc <- if (!is.null(post$smooth_curves) && j <= length(post$smooth_curves)) post$smooth_curves[[j]] else NULL
  if (!is.null(sc) && !is.null(sc$t_fine)) list(x = sc$t_fine, mu = sc$mean, lo = sc$lo, hi = sc$hi)
  else list(x = seq_len(T_len), mu = post$mu_mean[j, ], lo = post$mu_lo[j, ], hi = post$mu_hi[j, ])
}

# Series-level figure for one cluster
plot_hierscp_lower <- function(post, meta, group, save_pdf = NULL, min_seg_len = 2) {
  T_len <- meta$T_len; Y <- meta$Y_for_overlay
  series <- meta$group_members[[group]]; n <- length(series)
  .open_pdf(save_pdf, height = 1.75 * n + 0.95)
  layout(matrix(seq_len(n + 1), ncol = 1), heights = c(rep(1, n), 0.95 / 1.75))
  for (i in seq_along(series)) {
    j <- series[i]
    par(mar = c(if (i == n) 1.6 else 0.6, 3.4, 1.4, 2.4))
    S <- .seg_struct(post, group, .cp_pick(post$pma_lower[[j]], post$lower_cp_summary[[j]]),
                     T_len, min_seg_len = min_seg_len)
    cv <- .lower_curve(post, j, T_len)
    raw_list <- if (!is.null(Y)) list(Y[j, ]) else NULL
    nm <- meta$series_names[j]
    full <- if (!is.null(meta$country_names) && !is.na(meta$country_names[nm])) meta$country_names[nm] else nm
    .draw_panel(cv$x, cv$mu, cv$lo, cv$hi, S, meta, raw_list = raw_list,
                panel_label = sprintf("(%s) %s", letters[i], full), ylab = "CLI",
                xlab_on = (i == n))
  }
  .draw_legend(meta$jump_lim)
  if (!is.null(save_pdf)) dev.off()
}

# Common color limit: largest |jump| over all panels, rounded up
.global_jump_lim <- function(post, meta, min_seg_len = 2) {
  T_len <- meta$T_len; js <- numeric(0)
  for (g in seq_along(post$K_hat)) {
    S <- .seg_struct(post, g, .cp_pick(post$pma_upper[[g]], post$cp_summary[[g]]), T_len, min_seg_len)
    js <- c(js, .panel_jumps(seq_len(T_len), post$upper_curves[[g]]$mean, S)$jump)
    for (j in meta$group_members[[g]]) {
      S <- .seg_struct(post, g, .cp_pick(post$pma_lower[[j]], post$lower_cp_summary[[j]]), T_len, min_seg_len)
      cv <- .lower_curve(post, j, T_len)
      js <- c(js, .panel_jumps(cv$x, cv$mu, S)$jump)
    }
  }
  max(1, ceiling(max(abs(js), na.rm = TRUE)))
}

# Step 1: monthly panel
if (DO_BUILD) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  df <- df[, c("REF_AREA","TIME_PERIOD","OBS_VALUE")]
  df <- df[!is.na(df$OBS_VALUE) & nzchar(df$TIME_PERIOD), ]
  yr <- as.integer(substr(df$TIME_PERIOD, 1, 4)); mo <- as.integer(substr(df$TIME_PERIOD, 6, 7))
  df$midx <- yr * 12L + (mo - 1L)
  m0 <- start_ym[1]*12L + (start_ym[2]-1L); m1 <- end_ym[1]*12L + (end_ym[2]-1L)
  months <- m0:m1; T_len <- length(months)
  Y <- matrix(NA_real_, length(areas), T_len, dimnames = list(areas, NULL))
  for (i in seq_along(areas)) {
    sub <- df[df$REF_AREA == areas[i] & df$midx >= m0 & df$midx <= m1, ]
    Y[i, match(sub$midx, months)] <- sub$OBS_VALUE
  }
  if (anyNA(Y)) {
    bad <- areas[apply(Y, 1, anyNA)]
    stop("panel not balanced in this window; series with gaps: ", paste(bad, collapse = ", "))
  }
  yrs <- months %/% 12; mos <- months %% 12 + 1
  date_frac  <- yrs + (mos - 1)/12
  tick_years <- seq(ceiling(min(yrs)/5)*5, floor(max(yrs)), by = 5)
  tick_at    <- vapply(tick_years, function(Y0) which.min(abs(date_frac - Y0)), integer(1))
  meta <- list(areas = areas, group_members = group_members, group_names = group_names,
               series_names = areas, months = months, date_frac = date_frac,
               tick_at = tick_at, tick_years = tick_years, T_len = T_len, window = c(m0, m1))
  saveRDS(list(Y = Y, meta = meta), DATA_RDS)
  cat(sprintf("[STEP1] %s : %d series x %d months (%d-%02d..%d-%02d)\n",
              DATA_RDS, nrow(Y), T_len, m0%/%12, m0%%12+1, m1%/%12, m1%%12+1))
}

# Step 2: model fit
if (DO_FIT) {
  dat <- readRDS(DATA_RDS)

  Y_raw <- dat$Y
  meta <- dat$meta

  # y = (CLI - 100) / mean_t |CLI_t - 100|
  cli_center <- 100
  cli_scale  <- apply(Y_raw, 1, function(x) mean(abs(x - cli_center)))
  Y <- (Y_raw - cli_center) / cli_scale
  set.seed(SEED)

  model <- create_model_spec(
    Y, C = length(meta$group_members), group_members = meta$group_members,
    K_init = K_init, K_min = K_min, K_max = K_max, M = M, m_min = m_min,
    type = "continuous", fixed_clusters = TRUE, fixed_cluster_sizes = TRUE,
    lambda0_K = lambda0_K)
  model$prestandardized <- TRUE  # data already on the model scale
  res <- run_mcmc(Y, model, n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin,
                  N_launch = N_launch, n_warmup = n_warmup, ia_every = ia_every,
                  verbose = verbose_mcmc)
  post <- compute_conditional_posterior_summary(res, Y, verbose = TRUE)
  saveRDS(list(model = model, res = res, post = post, meta = meta, Y = Y, Y_raw = Y_raw,
               cli_center = cli_center, cli_scale = cli_scale), FIT_RDS)
  cat(sprintf("[STEP2] %s  K_hat=[%s]\n", FIT_RDS, paste(post$K_hat, collapse = ",")))
}

# Step 3: figures
if (DO_PLOT) {
  fit  <- readRDS(FIT_RDS)
  post <- fit$post; meta <- fit$meta
  meta$Y_for_overlay <- if (OVERLAY_RAW) fit$Y_raw else NULL  # observed CLI
  meta$ref_line      <- fit$cli_center  # long-term trend
  meta$country_names <- COUNTRY_NAMES
  means <- rep(fit$cli_center, nrow(fit$Y_raw)); sds <- fit$cli_scale

  J <- nrow(post$mu_mean)  # curves back to the CLI scale
  for (j in seq_len(J)) {
    post$mu_mean[j, ] <- post$mu_mean[j, ]*sds[j] + means[j]
    post$mu_lo[j, ]   <- post$mu_lo[j, ]  *sds[j] + means[j]
    post$mu_hi[j, ]   <- post$mu_hi[j, ]  *sds[j] + means[j]
    if (length(post$smooth_curves) >= j && !is.null(post$smooth_curves[[j]])) {
      sc <- post$smooth_curves[[j]]
      post$smooth_curves[[j]]$mean <- sc$mean*sds[j] + means[j]
      post$smooth_curves[[j]]$lo   <- sc$lo  *sds[j] + means[j]
      post$smooth_curves[[j]]$hi   <- sc$hi  *sds[j] + means[j]
    }
  }
  for (g in seq_along(post$upper_curves)) {  # cluster curves
    mem <- meta$group_members[[g]]
    post$upper_curves[[g]]$mean <- colMeans(post$mu_mean[mem, , drop = FALSE])
    mg <- mean(means[mem]); sg <- mean(sds[mem]); uc <- post$upper_curves[[g]]
    post$upper_curves[[g]]$lo <- uc$lo*sg + mg
    post$upper_curves[[g]]$hi <- uc$hi*sg + mg
  }
  meta$jump_lim <- .global_jump_lim(post, meta)  # one color scale for all figures

  plot_hierscp_upper(post, meta, save_pdf = PDF_UPPER)
  for (g in seq_along(meta$group_members))
    plot_hierscp_lower(post, meta, group = g, save_pdf = sprintf(PDF_LOWER, g))
  cat(sprintf("[STEP3] wrote %s and %s (jump colour limit = %g CLI points)\n", PDF_UPPER,
              paste(sprintf(PDF_LOWER, seq_along(meta$group_members)), collapse = ", "), meta$jump_lim))

  # changepoint table
  .md <- function(t){t<-max(1,min(meta$T_len,round(t)));mm<-meta$months[t];sprintf("%d-%02d",mm%/%12,mm%%12+1)}
  cat(sprintf("\nWindow %s..%s, T=%d, n_saved=%d\n", .md(1), .md(meta$T_len),
              meta$T_len, fit$res$samples$n_saved))
  for (g in seq_along(post$K_hat)) {
    S  <- .seg_struct(post, g, .cp_pick(post$pma_upper[[g]], post$cp_summary[[g]]), meta$T_len)
    JP <- .panel_jumps(seq_len(meta$T_len), post$upper_curves[[g]]$mean, S)
    cat(sprintf("\n[G%d %s] K=%d\n", g, meta$group_names[g], S$K))
    for (k in seq_len(S$K)[-1]) if (is.finite(S$cp_loc[k]))
      cat(sprintf("  CP %s [95%%: %s, %s]   %s -> %s   jump %+.2f (%.2f -> %.2f)\n", .md(S$cp_loc[k]),
                  if (is.finite(S$cp_lo[k])) .md(S$cp_lo[k]) else "NA",
                  if (is.finite(S$cp_hi[k])) .md(S$cp_hi[k]) else "NA",
                  .sp_lab(S$d1[k-1], S$d2[k-1]), .sp_lab(S$d1[k], S$d2[k]),
                  JP$jump[k], JP$before[k], JP$after[k]))
  }
}

# Convergence diagnostics

# install.packages("posterior")
library(posterior)

compute_mu_draws <- function(res, series = NULL, times = NULL, draws = NULL,
                             cli_center = NULL, cli_scale = NULL,
                             series_names = NULL, verbose = TRUE) {
  smp <- res$samples; model <- res$model; precomp <- res$precomp
  J <- model$J; T_len <- model$T_len; B <- smp$n_saved
  if (is.null(series)) series <- seq_len(J)        # which series (1..J)
  if (is.null(times))  times  <- seq_len(T_len)    # which months (1..T)
  if (is.null(draws))  draws  <- seq_len(B)        # which saved draws (1..B)
  
  nm  <- if (!is.null(series_names)) series_names[series] else as.character(series)
  out <- array(NA_real_, dim = c(length(draws), length(series), length(times)),
               dimnames = list(draw = draws, series = nm, t = times))
  
  for (i in seq_along(draws)) {
    b  <- draws[i]
    pb <- reconstruct_params_from_sample(b, smp, model)
    sb <- reconstruct_state_from_sample(b, smp, model)
    for (s in seq_along(series)) {
      mu_bs <- tryCatch(compute_mu_all(series[s], sb, pb, precomp, model),
                        error = function(e) rep(NA_real_, T_len))
      out[i, s, ] <- mu_bs[times]
    }
    if (verbose && i %% 1000 == 0) cat(sprintf("  %d / %d draws\n", i, length(draws)))
  }
  
  # back to the CLI scale: mu * (mean |CLI - 100|) + 100
  if (!is.null(cli_scale)) {
    for (s in seq_along(series)) out[, s, ] <- out[, s, ] * cli_scale[series[s]] + cli_center
  }
  out
}





##post <- fit$post
smp  <- fit$res$samples
d3   <- function(x) c(bulk = ess_bulk(x), tail = ess_tail(x), rhat = rhat(x))
tab  <- rbind(alpha = d3(smp$alpha), nu = d3(smp$nu), sigma2_gamma = d3(smp$sigma2_gamma))
round(tab, 3)

# ESS and R-hat of mu_{j,t} over all saved draws
# The same over the draws B_j of the fitted curves
mu_d <- compute_mu_draws(fit$res)
E3 <- array(NA_real_, c(3, 9, 300), dimnames = list(c("bulk", "tail", "rhat"), NULL, NULL))
for (j in 1:9) {
  idx <- post$B_j_list[[j]]
  E3[, j, ] <- apply(mu_d[idx, j, ], 2, d3)
}
c(bulk_ge400 = mean(E3["bulk", , ] >= 400), tail_ge400 = mean(E3["tail", , ] >= 400),
  rhat_lt101 = mean(E3["rhat", , ] < 1.01))

c(bulk_med=median(E3["bulk", , ]),tail_med=median(E3["tail", , ]),Rhat_med=median(E3["rhat", , ]))
