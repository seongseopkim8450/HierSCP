
library(dplyr)
library(lubridate)
library(RcppArmadillo)
library(ggplot2)

setwd("C://Users/User/Desktop/HierSCP_code/") # folder PATH

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
source("corefunctions_code/15_rcpp_bridge.R")
source("corefunctions_code/compute_atom_loglik_fallback.R")

MODEL_DIR <- "C://Users/User/Desktop/HierSCP_code/realdata_data_code" # YOUR DIR                       
csv_path  <- "OECD_CLI.csv"            

setwd(MODEL_DIR)
cli_data<-read.csv("OECD_CLI.csv")

unique(cli_data$Reference.area)


start_ym  <- c(2000, 1)                
end_ym    <- c(2026, 2)                

areas         <- c("USA","CAN","GBR","DEU","FRA","ITA","JPN","KOR","CHN")
group_members <- list(c(1,2), c(3,4,5,6), c(7,8,9))   # series-index -> cluster
group_names   <- c("North America","Major European Countries","East Asian Countries")

## 
SEED <- 2026                        
K_init <- 2
K_min <- 2
K_max <- 10   
M      <- 12                          
m_min  <- 12                           
lambda0_K <- 1.0                       
n_iter <- 50000; n_burnin <- 25000; n_thin <- 2    
N_launch <- 2; n_warmup <- 5000; ia_every <- 20; verbose_mcmc <- 15

## plotting -------------------------------------------------------------------
OVERLAY_RAW <- TRUE                    

## outputs --------------------------------------------------------------------
DATA_RDS  <- "cli_data.rds"
FIT_RDS   <- "cli_fit.rds"
PDF_UPPER <- "HierSCP_CLI_upper.pdf"
PDF_LOWER <- "HierSCP_CLI_lower_G%d.pdf"

## run switches ---------------------------------------------------------------
DO_BUILD <- TRUE
DO_FIT   <- TRUE                      
DO_PLOT  <- TRUE
## ============================================================================


.sp_key  <- function(d1, d2) paste0(if (d1 > 0) "+" else "-", if (d2 > 0) "+" else "-")
.sp_lab  <- function(d1, d2) sprintf("(%s,%s)", if (d1 > 0) "+" else "-",
                                     if (d2 > 0) "+" else "-")
.SP_COL  <- c("++" = "#E69F00",   # inc-convex   (orange)
              "+-" = "#CC79A7",   # inc-concave  (purple)
              "-+" = "#009E73",   # dec-convex   (green)
              "--" = "#0072B2")   # dec-concave  (blue)
.SP_NAME <- c("++" = "(+,+) inc-convex", "+-" = "(+,-) inc-concave",
              "-+" = "(-,+) dec-convex", "--" = "(-,-) dec-concave")
.CP_COL  <- "#E41A1C"             # changepoints: strong red, distinct from curve
.RAW_COL <- "grey55"

.alpha_col <- function(hex, a) adjustcolor(hex, alpha.f = a)
.cp_pick   <- function(primary, fallback) if (!is.null(primary) && length(primary) > 0) primary else fallback

.seg_struct <- function(post, g, cp_list, T_len, min_seg_len = 2) {
  K  <- post$K_hat[g]
  d1 <- d2 <- integer(K); b <- gmag <- numeric(K)
  for (k in seq_len(K)) {
    sc <- post$shape_class[[g]][[k]]; dm <- post$delta_mean[[g]][[k]]
    d1[k] <- if (!is.null(sc$d1_sign)) as.integer(sc$d1_sign) else 1L
    d2[k] <- if (!is.null(sc$d2_sign)) as.integer(sc$d2_sign) else 1L
    b[k]    <- if (!is.null(dm$shape_beta_mean)  && is.finite(dm$shape_beta_mean))  abs(dm$shape_beta_mean)  else NA_real_
    gmag[k] <- if (!is.null(dm$shape_gamma_mean) && is.finite(dm$shape_gamma_mean)) abs(dm$shape_gamma_mean) else NA_real_
  }
  
  starts <- rep(NA_real_, K); starts[1] <- 1
  cp_loc <- rep(NA_real_, K)
  if (!is.null(cp_list) && length(cp_list) > 0) {
    for (info in cp_list) {
      if (is.null(info) || is.null(info$k)) next
      k <- as.integer(info$k); if (k < 2L || k > K) next
      loc <- if (!is.null(info$point_est)) info$point_est else info$median
      starts[k] <- loc; cp_loc[k] <- loc
    }
  }
  for (k in seq_len(K)[-1]) if (is.na(starts[k])) starts[k] <- starts[k - 1]
  starts <- pmin(pmax(round(starts), 1), T_len)
  for (k in seq_len(K)[-1]) if (starts[k] < starts[k - 1]) starts[k] <- starts[k - 1]  
  ends <- c(starts[-1] - 1, T_len); ends <- pmin(pmax(ends, starts), T_len)
  
  segs <- lapply(seq_len(K), function(k)
    list(d1 = d1[k], d2 = d2[k], b = b[k], g = gmag[k],
         start = starts[k], end = ends[k], cp = cp_loc[k]))
  repeat {
    if (length(segs) <= 1) break
    w   <- vapply(segs, function(s) s$end - s$start + 1, numeric(1))
    bad <- which(w < min_seg_len)
    if (length(bad) == 0) break
    k <- bad[1]
    if (k == 1) {                       
      segs[[2]]$start <- segs[[1]]$start
      segs[[2]]$cp    <- NA_real_
      segs[[1]] <- NULL
    } else {                            
      segs[[k - 1]]$end <- segs[[k]]$end
      segs[[k]] <- NULL
    }
  }
  
  Kf <- length(segs)
  list(K      = Kf,
       starts = vapply(segs, `[[`, numeric(1), "start"),
       ends   = vapply(segs, `[[`, numeric(1), "end"),
       d1     = vapply(segs, `[[`, numeric(1), "d1"),
       d2     = vapply(segs, `[[`, numeric(1), "d2"),
       b      = vapply(segs, `[[`, numeric(1), "b"),
       g      = vapply(segs, `[[`, numeric(1), "g"),
       cp_loc = vapply(segs, `[[`, numeric(1), "cp"))
}

.seg_of_x <- function(x, starts, K) {
  se <- c(starts, Inf); idx <- rep(K, length(x))
  for (k in seq_len(K)) idx[x >= se[k] & x < se[k + 1]] <- k
  idx[x < starts[1]] <- 1L
  idx
}


.draw_panel <- function(x, mu, lo, hi, S, meta,
                        raw_list = NULL, raw_col = .RAW_COL,
                        title = "", ylab = "CLI",
                        cex_axis = 1.25, cex_lab = 1.35,
                        cex_main = 1.45, cex_sp = 1.20) {
  T_len <- meta$T_len
  yvals <- c(mu, lo, hi, unlist(raw_list))
  yr <- range(yvals, na.rm = TRUE); pad <- 0.06 * diff(yr)
  ylim <- c(yr[1] - pad, yr[2] + 1.7 * pad)
  
  plot(NA, xlim = c(1, T_len), ylim = ylim, xlab = "", ylab = "",
       main = "", axes = FALSE, xaxs = "i")
  
  seg_x  <- .seg_of_x(x, S$starts, S$K)
  key_k  <- vapply(seq_len(S$K), function(k) .sp_key(S$d1[k], S$d2[k]), character(1))
  
  
  for (k in seq_len(S$K)) {
    xr <- if (k < S$K) S$starts[k + 1] else T_len + 0.5
    rect(S$starts[k] - 0.5, ylim[1], xr - 0.5, ylim[2],
         col = .alpha_col(.SP_COL[key_k[k]], 0.11), border = NA)
  }
  
  
  ok <- is.finite(lo) & is.finite(hi)
  for (k in seq_len(S$K)) {
    sel <- ok & (seg_x == k)
    if (sum(sel) >= 2) {
      xx <- x[sel]; polygon(c(xx, rev(xx)), c(lo[sel], rev(hi[sel])),
                            col = .alpha_col(.SP_COL[key_k[k]], 0.28), border = NA)
    }
  }
  
  
  if (!is.null(raw_list)) for (z in raw_list)
    lines(seq_len(T_len), z, col = .alpha_col(raw_col, 0.45), lwd = 0.9)
  
  n <- length(x)
  col_x <- .SP_COL[ key_k[seg_x] ]
  segments(x[-n], mu[-n], x[-1], mu[-1], col = col_x[-1], lwd = 3.0)
  
  for (k in seq_len(S$K)[-1]) if (is.finite(S$cp_loc[k]))
    abline(v = S$cp_loc[k], col = .CP_COL, lwd = 2.4)
  
  reg_start <- 1
  for (k in seq_len(S$K)) {
    last_in_run <- (k == S$K) || !(S$d1[k + 1] == S$d1[k] && S$d2[k + 1] == S$d2[k])
    if (last_in_run) {
      mid <- (S$starts[reg_start] + S$ends[k]) / 2
      text(mid, ylim[2] - 0.55 * pad, .sp_lab(S$d1[k], S$d2[k]),
           col = .SP_COL[key_k[k]], font = 2, cex = cex_sp, xpd = NA)
      reg_start <- k + 1
    }
  }
  
  title(main = title, cex.main = cex_main)
  title(ylab = ylab, cex.lab = cex_lab, line = 2.6)
  axis(1, at = meta$tick_at, labels = meta$tick_years, cex.axis = cex_axis)
  axis(2, las = 1, cex.axis = cex_axis); box()
}

# ---- legend -----------------------------------------------------------------
.draw_legend <- function(cex = 1.35) {
  op <- par(mar = c(0.4, 1, 0.4, 1)); on.exit(par(op)); plot.new()
  fills <- c(.alpha_col(.SP_COL["++"], .55), .alpha_col(.SP_COL["+-"], .55),
             .alpha_col(.SP_COL["-+"], .55), .alpha_col(.SP_COL["--"], .55),
             .alpha_col("grey40", 0.35), NA, NA)
  legend("center", ncol = 4, bty = "n", cex = cex,
         legend = c(.SP_NAME[["++"]], .SP_NAME[["+-"]], .SP_NAME[["-+"]], .SP_NAME[["--"]],
                    "95% pointwise CI", "changepoint", "member series"),
         fill   = fills,
         border = c("grey30", "grey30", "grey30", "grey30", NA, NA, NA),
         lty    = c(NA, NA, NA, NA, NA, 1, 1),
         lwd    = c(NA, NA, NA, NA, NA, 2.6, 1.2),
         col    = c(NA, NA, NA, NA, NA, .CP_COL, .RAW_COL))
}

# ---- upper level (group curves) ---------------------------------------------
plot_hierscp_upper <- function(post, meta, save_pdf = NULL,
                               width = 11, height = NULL, min_seg_len = 2) {
  C <- length(post$K_hat); T_len <- meta$T_len; Y <- meta$Y_for_overlay
  if (is.null(height)) height <- 3.0 * C + 1.3
  if (!is.null(save_pdf)) pdf(save_pdf, width = width, height = height)
  layout(matrix(seq_len(C + 1), ncol = 1), heights = c(rep(3, C), 0.9))
  par(mar = c(3.4, 4.6, 2.8, 1.2), mgp = c(2.6, 0.8, 0))
  for (g in seq_len(C)) {
    S  <- .seg_struct(post, g, .cp_pick(post$pma_upper[[g]], post$cp_summary[[g]]),
                      T_len, min_seg_len = min_seg_len)
    uc <- post$upper_curves[[g]]; members <- meta$group_members[[g]]
    raw_list <- if (!is.null(Y)) lapply(members, function(j) Y[j, ]) else NULL
    ttl <- sprintf("Upper level - Group %d (%s):  %s", g, meta$group_names[g],
                   paste(meta$series_names[members], collapse = ", "))
    .draw_panel(seq_len(T_len), uc$mean, uc$lo, uc$hi, S, meta,
                raw_list = raw_list, title = ttl, ylab = "CLI (group)")
  }
  .draw_legend(); if (!is.null(save_pdf)) dev.off()
}

# ---- lower level (series curves) --------------------------------------------
plot_hierscp_lower <- function(post, meta, group = NULL, save_pdf = NULL,
                               width = 11, height = NULL, min_seg_len = 2) {
  C <- length(post$K_hat); T_len <- meta$T_len; Y <- meta$Y_for_overlay
  groups <- if (is.null(group)) seq_len(C) else group
  series <- unlist(meta$group_members[groups]); n <- length(series)
  if (is.null(height)) height <- 2.9 * n + 1.3
  if (!is.null(save_pdf)) pdf(save_pdf, width = width, height = height)
  layout(matrix(seq_len(n + 1), ncol = 1), heights = c(rep(3, n), 0.9))
  par(mar = c(3.4, 4.6, 2.8, 1.2), mgp = c(2.6, 0.8, 0))
  s2g <- integer(0); for (g in seq_len(C)) for (j in meta$group_members[[g]]) s2g[j] <- g
  for (j in series) {
    g <- s2g[j]
    S <- .seg_struct(post, g, .cp_pick(post$pma_lower[[j]], post$lower_cp_summary[[j]]),
                     T_len, min_seg_len = min_seg_len)
    sc <- if (!is.null(post$smooth_curves) && j <= length(post$smooth_curves)) post$smooth_curves[[j]] else NULL
    if (!is.null(sc) && !is.null(sc$t_fine)) { x <- sc$t_fine; mu <- sc$mean; lo <- sc$lo; hi <- sc$hi }
    else { x <- seq_len(T_len); mu <- post$mu_mean[j, ]; lo <- post$mu_lo[j, ]; hi <- post$mu_hi[j, ] }
    raw_list <- if (!is.null(Y)) list(Y[j, ]) else NULL
    ttl <- sprintf("Lower level - Series %d: %s  (Group %d, %s)",
                   j, meta$series_names[j], g, meta$group_names[g])
    .draw_panel(x, mu, lo, hi, S, meta, raw_list = raw_list, raw_col = "grey45",
                title = ttl, ylab = "CLI")
  }
  .draw_legend(); if (!is.null(save_pdf)) dev.off()
}

## ============================ STEP 1: BUILD PANEL ===========================
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


## ============================ STEP 2: FIT + SUMMARY =========================
if (DO_FIT) {
  #src <- c("00_utils.R","01_precompute.R","02_shape_function.R","03_state_management.R",
  #         "04_sampling.R","05_likelihood.R","06_local_adjusting.R","07_interval_adjusting.R",
  #         "08_param_update_continuous.R","10_initialization.R","11_mcmc_main.R",
  #         "12_post_summary.R","compute_atom_loglik_fallback.R","15_rcpp_bridge.R")
  #for (f in src) source(file.path(MODEL_DIR, f))     # 15 auto-compiles hiercpd_core.cpp
  
  dat <- readRDS(DATA_RDS)
  Y <- dat$Y
  meta <- dat$meta
  
  Y <- apply(Y, 1, function(x){(x-mean(x))/sd(x)})
  Y <- t(Y)
  set.seed(SEED)
  model <- create_model_spec(
    Y, C = length(meta$group_members), group_members = meta$group_members,
    K_init = K_init, K_min = K_min, K_max = K_max, M = M, m_min = m_min,
    type = "continuous", fixed_clusters = TRUE, fixed_cluster_sizes = TRUE,
    lambda0_K = lambda0_K)
  res <- run_mcmc(Y, model, n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin,
                  N_launch = N_launch, n_warmup = n_warmup, ia_every = ia_every,
                  verbose = verbose_mcmc)
  post <- compute_conditional_posterior_summary(res, Y, cp_window = 3, verbose = TRUE)
  saveRDS(list(model = model, res = res, post = post, meta = meta, Y = Y), FIT_RDS)
  cat(sprintf("[STEP2] %s  K_hat=[%s]\n", FIT_RDS, paste(post$K_hat, collapse = ",")))
}


if (DO_PLOT) {
  fit  <- readRDS(FIT_RDS)
  post <- fit$post; meta <- fit$meta; Y <- fit$Y
  meta$Y_for_overlay <- if (OVERLAY_RAW) Y else NULL
  si <- fit$res$precomp$std_info; means <- si$means; sds <- si$sds
  
  J <- nrow(post$mu_mean)                     # un-standardise curves -> original CLI scale
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
  for (g in seq_along(post$upper_curves)) {   # group curve
    mem <- meta$group_members[[g]]
    post$upper_curves[[g]]$mean <- colMeans(post$mu_mean[mem, , drop = FALSE])
    mg <- mean(means[mem]); sg <- mean(sds[mem]); uc <- post$upper_curves[[g]]
    post$upper_curves[[g]]$lo <- uc$lo*sg + mg
    post$upper_curves[[g]]$hi <- uc$hi*sg + mg
  }
  
  plot_hierscp_upper(post, meta, save_pdf = PDF_UPPER)
  for (g in seq_along(meta$group_members))
    plot_hierscp_lower(post, meta, group = g, save_pdf = sprintf(PDF_LOWER, g), height = 9.5)
  cat(sprintf("[STEP3] wrote %s and %s\n", PDF_UPPER, sprintf(PDF_LOWER, seq_along(meta$group_members))))
  
  .md <- function(t){t<-max(1,min(meta$T_len,round(t)));mm<-meta$months[t];sprintf("%d-%02d",mm%/%12,mm%%12+1)}
  .sp <- function(d1,d2) sprintf("(%s,%s)", if(d1>0)"+"else"-", if(d2>0)"+"else"-")
  cat(sprintf("\nWindow %s..%s, T=%d, n_saved=%d\n", .md(1), .md(meta$T_len),
              meta$T_len, fit$res$samples$n_saved))
  for (g in seq_along(post$K_hat)) {
    S <- .seg_struct(post, g, .cp_pick(post$pma_upper[[g]], post$cp_summary[[g]]), meta$T_len)
    cat(sprintf("\n[G%d %s] K=%d (after sliver merge)\n", g, meta$group_names[g], S$K))
    for (k in seq_len(S$K)[-1]) if (is.finite(S$cp_loc[k]))
      cat(sprintf("  CP %s   %s -> %s\n", .md(S$cp_loc[k]),
                  .sp(S$d1[k-1], S$d2[k-1]), .sp(S$d1[k], S$d2[k])))
  }
}






