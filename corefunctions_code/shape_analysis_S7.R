

## YOUR PATH: MUST CHANGE IT TO YOUR OWN
setwd("C://Users/limit/Desktop/HierSCP_code_2/corefunctions_code")
####
RESULT_SHOW = FALSE #### whether to show, tablulate, plot the results 

.spt_or <- function(a, b) if (!is.null(a)) a else b


# Scenario label "S1", "S2", "S3" from the stored scenario name.
.spt_scenario_label <- function(raw) {
  s <- .spt_or(raw$scenario, .spt_or(raw$sim$truth$scenario, raw$sim_m$scenario$name))
  if (is.null(s) || !length(s) || is.na(s[1])) return(NA_character_)
  s <- as.character(s[1])
  m <- regmatches(s, regexpr("^S[0-9]+", s))
  if (length(m)) m else s
}


# One changepoint value from a summary item, exactly as extract_point_estimates():
# point_est, then mean, then median, rounded to an integer.
.spt_cp_value <- function(info) {
  if (is.null(info)) return(NA_integer_)
  for (nm in c("point_est", "mean", "median")) {
    v <- info[[nm]]
    if (!is.null(v) && length(v) == 1L && is.finite(v)) return(as.integer(round(v)))
  }
  NA_integer_
}


.spt_est_regime_path <- function(cp_items, T_len) {
  path <- rep(1L, T_len)
  out_empty <- list(path = path, n_cp = 0L, mismatch = 0L)
  if (is.null(cp_items) || !length(cp_items)) return(out_empty)

  taus <- integer(0)
  ks   <- integer(0)
  for (i in seq_along(cp_items)) {
    info <- cp_items[[i]]
    if (is.null(info)) next
    tau <- .spt_cp_value(info)
    if (is.na(tau) || tau < 2L || tau > T_len) next
    k <- if (!is.null(info$k) && length(info$k) == 1L && is.finite(info$k)) {
      as.integer(info$k)
    } else NA_integer_
    taus <- c(taus, tau)
    ks   <- c(ks, k)
  }
  if (!length(taus)) return(out_empty)

  o <- order(taus, ks)
  taus <- taus[o]
  ks   <- ks[o]

  # Order-based labels: identical to cps_to_segvec(extract_point_estimates(.)).
  ut <- unique(taus)
  path_order <- rep(1L, T_len)
  for (i in seq_along(ut)) path_order[ut[i]:T_len] <- as.integer(i + 1L)

  if (anyNA(ks)) {
    path <- path_order
  } else {
    for (i in seq_along(taus)) path[taus[i]:T_len] <- ks[i]
  }

  list(path = path, n_cp = length(ut), mismatch = sum(path != path_order))
}


.spt_est_pair <- function(sc) {
  if (is.null(sc)) return(c(NA_integer_, NA_integer_))
  d1 <- NA_integer_
  d2 <- NA_integer_
  if (!is.null(sc$d1_sign)) d1 <- as.integer(sc$d1_sign)
  if (!is.null(sc$d2_sign)) d2 <- as.integer(sc$d2_sign)
  if (is.na(d1) && !is.null(sc$delta1)) d1 <- as.integer(sc$delta1)
  if (is.na(d2) && !is.null(sc$delta2)) d2 <- as.integer(sc$delta2)
  if ((is.na(d1) || is.na(d2)) && !is.null(sc$shape_key)) {
    map <- list(inc_convex  = c(+1L, +1L), inc_concave = c(+1L, -1L),
                dec_concave = c(-1L, -1L), dec_convex  = c(-1L, +1L))
    if (!is.null(map[[sc$shape_key]])) {
      d1 <- map[[sc$shape_key]][1]
      d2 <- map[[sc$shape_key]][2]
    }
  }
  c(as.integer(sign(d1)), as.integer(sign(d2)))
}


.spt_est_prob <- function(dm, sc) {
  p <- c(NA_real_, NA_real_)
  fb <- 0L
  if (!is.null(dm)) {
    if (!is.null(dm$d1_prob_pos) && is.finite(dm$d1_prob_pos)) p[1] <- as.numeric(dm$d1_prob_pos)
    if (!is.null(dm$d2_prob_pos) && is.finite(dm$d2_prob_pos)) p[2] <- as.numeric(dm$d2_prob_pos)
  }
  for (i in 1:2) {
    if (is.na(p[i]) && !is.null(sc)) {
      conf <- sc[[c("d1_confidence", "d2_confidence")[i]]]
      sgn  <- sc[[c("d1_sign", "d2_sign")[i]]]
      if (!is.null(conf) && !is.null(sgn) && is.finite(conf) && is.finite(sgn)) {
        p[i] <- if (sgn > 0) conf else 1 - conf
        fb <- fb + 1L
      }
    }
  }
  list(p = p, fallback = fb)
}


# True sign pairs of the regimes of cluster g (rows = regimes, cols = delta1, delta2).
.spt_true_pairs <- function(sim_m, g) {
  st <- sim_m$shape_truth[[g]]
  if (is.null(st)) stop("sim_m$shape_truth is missing for cluster ", g, call. = FALSE)
  if (is.data.frame(st) || is.matrix(st)) {
    out <- cbind(as.integer(st[, "delta1"]), as.integer(st[, "delta2"]))
  } else {
    out <- do.call(rbind, lapply(st, function(x) {
      c(as.integer(.spt_or(x$d1, x$delta1)), as.integer(.spt_or(x$d2, x$delta2)))
    }))
  }
  storage.mode(out) <- "integer"
  out
}


spt_metrics_one <- function(raw, cp_source = c("lower_cp_summary", "pma_lower")) {
  cp_source <- match.arg(cp_source)
  post  <- raw$post
  sim_m <- raw$sim_m
  if (is.null(sim_m)) {
    if (exists("adapt_sim_for_metrics", mode = "function")) {
      sim_m <- adapt_sim_for_metrics(raw$sim)
    } else {
      stop("raw$sim_m is missing and adapt_sim_for_metrics() is not loaded.", call. = FALSE)
    }
  }
  if (is.null(post)) stop("raw$post is missing.", call. = FALSE)

  J     <- as.integer(sim_m$J)
  T_len <- as.integer(sim_m$T_len)
  C     <- as.integer(sim_m$C)

  S_true <- sim_m$S_lower_true
  if (is.null(S_true) || !identical(dim(S_true), c(J, T_len))) {
    stop("sim_m$S_lower_true is missing or has the wrong dimension.", call. = FALSE)
  }

  # Cluster indices: estimate side (post) and truth side (sim_m) must agree.
  cl_true <- as.integer(sim_m$cluster_assign)
  gm_est  <- .spt_or(post$group_members, sim_m$group_members)
  cl_est  <- integer(J)
  for (g in seq_along(gm_est)) cl_est[as.integer(gm_est[[g]])] <- g
  if (any(cl_est != cl_true)) {
    stop("Cluster indexing of post$group_members and the truth differ.", call. = FALSE)
  }

  K_hat  <- as.integer(post$K_hat)
  K_true <- as.integer(sim_m$K_true)
  cp_lists <- post[[cp_source]]
  if (is.null(cp_lists)) cp_lists <- vector("list", J)

  # Regime-level estimated pairs and probabilities, per cluster.
  est_pairs <- vector("list", C)
  est_probs <- vector("list", C)
  n_prob_fallback <- 0L
  for (g in seq_len(C)) {
    Kg <- K_hat[g]
    ep <- matrix(NA_integer_, Kg, 2)
    pp <- matrix(NA_real_, Kg, 2)
    for (k in seq_len(Kg)) {
      sc <- if (length(post$shape_class[[g]]) >= k) post$shape_class[[g]][[k]] else NULL
      dm <- if (length(post$delta_mean[[g]])  >= k) post$delta_mean[[g]][[k]]  else NULL
      ep[k, ] <- .spt_est_pair(sc)
      pr <- .spt_est_prob(dm, sc)
      pp[k, ] <- pr$p
      n_prob_fallback <- n_prob_fallback + pr$fallback
    }
    est_pairs[[g]] <- ep
    est_probs[[g]] <- pp
  }
  true_pairs <- lapply(seq_len(C), function(g) .spt_true_pairs(sim_m, g))

  # Time-point accumulation.
  hit1 <- hit2 <- hitb <- numeric(J)
  sq1  <- sq2  <- numeric(J)
  n_mismatch <- 0L
  n_clamp    <- 0L
  n_na       <- 0L
  for (j in seq_len(J)) {
    g  <- cl_true[j]
    rp <- .spt_est_regime_path(cp_lists[[j]], T_len)
    kh <- rp$path
    n_mismatch <- n_mismatch + rp$mismatch
    if (any(kh > K_hat[g])) {               # cannot occur in normal output
      n_clamp <- n_clamp + sum(kh > K_hat[g])
      kh <- pmin(kh, K_hat[g])
    }
    kt <- as.integer(S_true[j, ])

    dh <- est_pairs[[g]][kh, , drop = FALSE]
    ph <- est_probs[[g]][kh, , drop = FALSE]
    d0 <- true_pairs[[g]][kt, , drop = FALSE]

    n_na <- n_na + sum(is.na(dh)) + sum(is.na(ph))
    m1 <- !is.na(dh[, 1]) & dh[, 1] == d0[, 1]
    m2 <- !is.na(dh[, 2]) & dh[, 2] == d0[, 2]
    hit1[j] <- sum(m1)
    hit2[j] <- sum(m2)
    hitb[j] <- sum(m1 & m2)
    sq1[j]  <- sum((ph[, 1] - as.numeric(d0[, 1] == 1L))^2)
    sq2[j]  <- sum((ph[, 2] - as.numeric(d0[, 2] == 1L))^2)
  }
  JT <- J * T_len

  # RISE, as in compute_all_metrics().
  rise <- rep(NA_real_, J)
  if (!is.null(post$mu_mean) && !is.null(sim_m$mu_true)) {
    for (j in seq_len(J)) {
      rise[j] <- sqrt(mean((post$mu_mean[j, ] - sim_m$mu_true[j, ])^2, na.rm = TRUE))
    }
  }

  # Regime-level legacy quantities (clusters with K_hat = K only).
  n_km <- 0L
  N_reg <- 0L
  c1 <- c2 <- cb <- 0L
  s1 <- s2 <- 0
  for (g in seq_len(C)) {
    if (is.na(K_hat[g]) || K_hat[g] != K_true[g] || is.null(post$shape_class[[g]])) next
    n_km <- n_km + 1L
    for (k in seq_len(K_true[g])) {
      e  <- est_pairs[[g]][k, ]
      p  <- est_probs[[g]][k, ]
      tr <- true_pairs[[g]][k, ]
      m1 <- !is.na(e[1]) && e[1] == tr[1]
      m2 <- !is.na(e[2]) && e[2] == tr[2]
      c1 <- c1 + m1
      c2 <- c2 + m2
      cb <- cb + (m1 && m2)
      s1 <- s1 + (p[1] - as.numeric(tr[1] == 1L))^2
      s2 <- s2 + (p[2] - as.numeric(tr[2] == 1L))^2
      N_reg <- N_reg + 1L
    }
  }

  row <- data.frame(
    scenario    = .spt_scenario_label(raw),
    seed        = as.integer(.spt_or(raw$seed, NA_integer_)),
    SCApt_joint = sum(hitb) / JT,
    SCApt_d1    = sum(hit1) / JT,
    SCApt_d2    = sum(hit2) / JT,
    Bpt_d1      = sum(sq1) / JT,
    Bpt_d2      = sum(sq2) / JT,
    RISE        = mean(rise, na.rm = TRUE),
    n_Kmatch    = n_km,
    C           = C,
    N_reg       = N_reg,
    SA_reg_d1   = if (N_reg > 0) c1 / N_reg else NA_real_,
    SA_reg_d2   = if (N_reg > 0) c2 / N_reg else NA_real_,
    SA_reg_both = if (N_reg > 0) cb / N_reg else NA_real_,
    B_reg_d1    = if (N_reg > 0) s1 / N_reg else NA_real_,
    B_reg_d2    = if (N_reg > 0) s2 / N_reg else NA_real_,
    c_reg_d1 = c1, c_reg_d2 = c2, c_reg_both = cb,
    sum_reg_d1 = s1, sum_reg_d2 = s2,
    n_label_mismatch = n_mismatch,
    n_prob_fallback  = n_prob_fallback,
    n_clamped        = n_clamp,
    n_na             = n_na,
    cp_source        = cp_source,
    stringsAsFactors = FALSE
  )

  by_cluster <- do.call(rbind, lapply(seq_len(C), function(g) {
    idx <- which(cl_true == g)
    nT  <- length(idx) * T_len
    data.frame(
      scenario    = row$scenario,
      seed        = row$seed,
      cluster     = g,
      K_true      = K_true[g],
      K_hat       = K_hat[g],
      SCApt_joint = sum(hitb[idx]) / nT,
      SCApt_d1    = sum(hit1[idx]) / nT,
      SCApt_d2    = sum(hit2[idx]) / nT,
      Bpt_d1      = sum(sq1[idx]) / nT,
      Bpt_d2      = sum(sq2[idx]) / nT,
      RISE        = mean(rise[idx], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  list(row = row, by_cluster = by_cluster)
}


# Table S7.2 layout: means over replications by scenario.
spt_table <- function(per_rep, digits = 3) {
  scs <- sort(unique(per_rep$scenario))
  out <- do.call(rbind, lapply(scs, function(s) {
    d <- per_rep[per_rep$scenario == s, , drop = FALSE]
    data.frame(
      scenario    = s,
      R           = nrow(d),
      RISE        = mean(d$RISE, na.rm = TRUE),
      SCApt_joint = mean(d$SCApt_joint),
      SCApt_d1    = mean(d$SCApt_d1),
      SCApt_d2    = mean(d$SCApt_d2),
      Bpt_d1      = mean(d$Bpt_d1),
      Bpt_d2      = mean(d$Bpt_d2),
      stringsAsFactors = FALSE
    )
  }))
  num <- vapply(out, is.numeric, logical(1)) & names(out) != "R"
  out[num] <- lapply(out[num], round, digits = digits)
  out
}


spt_check_table <- function(per_rep, digits = 3) {
  scs <- sort(unique(per_rep$scenario))
  out <- do.call(rbind, lapply(scs, function(s) {
    d <- per_rep[per_rep$scenario == s, , drop = FALSE]
    N <- sum(d$N_reg)
    data.frame(
      scenario          = s,
      R                 = nrow(d),
      Kmatch_rate       = sum(d$n_Kmatch) / sum(d$C),
      N_reg_total       = N,
      SA_reg_d1_pooled  = if (N > 0) sum(d$c_reg_d1) / N else NA_real_,
      SA_reg_d2_pooled  = if (N > 0) sum(d$c_reg_d2) / N else NA_real_,
      B_reg_d1_pooled   = if (N > 0) sum(d$sum_reg_d1) / N else NA_real_,
      B_reg_d2_pooled   = if (N > 0) sum(d$sum_reg_d2) / N else NA_real_,
      B_reg_d1_repmean  = mean(d$B_reg_d1, na.rm = TRUE),
      B_reg_d2_repmean  = mean(d$B_reg_d2, na.rm = TRUE),
      n_label_mismatch  = sum(d$n_label_mismatch),
      n_prob_fallback   = sum(d$n_prob_fallback),
      stringsAsFactors = FALSE
    )
  }))
  num <- vapply(out, is.numeric, logical(1)) &
    !names(out) %in% c("R", "N_reg_total", "n_label_mismatch", "n_prob_fallback")
  out[num] <- lapply(out[num], round, digits = digits)
  out
}


# Per-cluster means over replications, by scenario.
spt_cluster_table <- function(by_cluster, digits = 3) {
  key <- unique(by_cluster[, c("scenario", "cluster")])
  key <- key[order(key$scenario, key$cluster), , drop = FALSE]
  out <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
    d <- by_cluster[by_cluster$scenario == key$scenario[i] &
                      by_cluster$cluster == key$cluster[i], , drop = FALSE]
    data.frame(
      scenario    = key$scenario[i],
      cluster     = key$cluster[i],
      R           = nrow(d),
      Kmatch_rate = mean(d$K_hat == d$K_true),
      RISE        = mean(d$RISE, na.rm = TRUE),
      SCApt_joint = mean(d$SCApt_joint),
      SCApt_d1    = mean(d$SCApt_d1),
      SCApt_d2    = mean(d$SCApt_d2),
      Bpt_d1      = mean(d$Bpt_d1),
      Bpt_d2      = mean(d$Bpt_d2),
      stringsAsFactors = FALSE
    )
  }))
  num <- vapply(out, is.numeric, logical(1)) & !names(out) %in% c("cluster", "R")
  out[num] <- lapply(out[num], round, digits = digits)
  out
}

spt_latex_rows <- function(table, digits = 3,
                           labels = c(S1 = "S1 (benign)", S2 = "S2 (contaminated)",
                                      S3 = "S3 (misspecified)")) {
  f <- function(v) formatC(v, format = "f", digits = digits)
  rows <- vapply(seq_len(nrow(table)), function(i) {
    s  <- table$scenario[i]
    lb <- if (!is.na(labels[s])) labels[[s]] else s
    paste0(lb, " & ", f(table$RISE[i]), " & ", f(table$SCApt_joint[i]), " & ",
           f(table$SCApt_d1[i]), " & ", f(table$SCApt_d2[i]), " & ",
           f(table$Bpt_d1[i]), " & ", f(table$Bpt_d2[i]), " \\\\")
  }, character(1))
  cat(rows, sep = "\n")
  invisible(rows)
}

spt_metrics_folder <- function(paths, pattern = "_raw\\.rds$", recursive = FALSE,
                               cp_source = c("lower_cp_summary", "pma_lower"),
                               digits = 3, verbose = TRUE) {
  cp_source <- match.arg(cp_source)
  files <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, pattern = pattern, full.names = TRUE, recursive = recursive)
    else if (file.exists(p)) p
    else { warning("Not found: ", p); character(0) }
  }))
  if (!length(files)) stop("No files matching '", pattern, "' were found.", call. = FALSE)

  res <- vector("list", length(files))
  for (i in seq_along(files)) {
    raw <- readRDS(files[i])
    res[[i]] <- tryCatch({
      o <- spt_metrics_one(raw, cp_source = cp_source)
      o$row$file <- basename(files[i])
      o
    }, error = function(e) {
      warning(sprintf("%s: %s", basename(files[i]), conditionMessage(e)))
      NULL
    })
    rm(raw)
    if (verbose && (i %% 10L == 0L || i == length(files))) {
      message(sprintf("  processed %d / %d files", i, length(files)))
    }
  }
  res <- Filter(Negate(is.null), res)
  per_rep    <- do.call(rbind, lapply(res, `[[`, "row"))
  by_cluster <- do.call(rbind, lapply(res, `[[`, "by_cluster"))
  per_rep    <- per_rep[order(per_rep$scenario, per_rep$seed), , drop = FALSE]
  rownames(per_rep) <- NULL

  out <- list(
    table      = spt_table(per_rep, digits = digits),
    check      = spt_check_table(per_rep, digits = digits),
    by_cluster = spt_cluster_table(by_cluster, digits = digits),
    per_rep    = per_rep,
    per_rep_cluster = by_cluster
  )

  if (verbose) {
    if (sum(per_rep$n_label_mismatch) > 0L) {
      message("Note: ", sum(per_rep$n_label_mismatch),
              " time points used index-based regime labels that differ from the ",
              "order-based labels (a changepoint summary item was missing).")
    }
    if (sum(per_rep$n_prob_fallback) > 0L) {
      message("Note: ", sum(per_rep$n_prob_fallback),
              " regime sign probabilities were recovered from shape_class confidences.")
    }
    if (sum(per_rep$n_na) > 0L) {
      message("Note: ", sum(per_rep$n_na), " missing sign or probability entries were found.")
    }
  }
  out
}


if(RESULT_SHOW){
folders <- c(S1 = "S1_replicate_50rep",
             S2 = "S2_replicate_50rep",
             S3 = "S3_replicate_50rep")

CP_SOURCE <- "lower_cp_summary"

res <- spt_metrics_folder(folders, cp_source = CP_SOURCE)

cat("\n== Table S7.2: RISE, SCApt (joint, d1, d2), time-point Brier (B1, B2); means over replications\n")
print(res$table)

cat("\n== LaTeX rows for Table S7.2\n")
spt_latex_rows(res$table)

cat("\n== Regime-level values on clusters with K_hat = K (current Table S7.2 definition), with coverage\n")
print(res$check)

cat("\n== By cluster\n")
print(res$by_cluster)

write.csv(res$table,      "table_S7_2_shape.csv",          row.names = FALSE)
write.csv(res$check,      "table_S7_2_regime_level.csv",   row.names = FALSE)
write.csv(res$by_cluster, "table_S7_2_by_cluster.csv",     row.names = FALSE)
write.csv(res$per_rep,    "table_S7_2_per_replication.csv", row.names = FALSE)
saveRDS(res, "shape_metrics_S7.rds")
}
