# One-at-a-time sensitivity analysis on S2 (seeds 2026-2045): m_min = 8, 12 (default 10),
# lambda0 = 0.5, 2 (default 1), m0 = 0.01, 0.10 (default 0.03); driver run_sensitivity_S2.R

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Source only the function definitions and constants of a script
.sens_source_functions <- function(path, envir = globalenv()) {
  is_const <- function(x) {
    if (is.atomic(x) || is.null(x)) return(TRUE)
    if (is.call(x) && is.name(x[[1]]) &&
        as.character(x[[1]]) %in% c("c", "list", ":", "-", "+", "(")) {
      return(all(vapply(as.list(x)[-1], is_const, logical(1))))
    }
    FALSE
  }
  exprs <- parse(path, keep.source = FALSE)
  n <- 0L
  for (e in exprs) {
    if (is.call(e) && is.name(e[[1]]) && as.character(e[[1]]) %in% c("<-", "=") &&
        is.name(e[[2]])) {
      rhs <- e[[3]]
      is_fun <- is.call(rhs) && identical(rhs[[1]], as.name("function"))
      if (is_fun || is_const(rhs)) {
        eval(e, envir)
        n <- n + 1L
      }
    }
  }
  invisible(n)
}

SENS_SETTINGS_S2 <- list(
  default     = list(),
  mmin_8      = list(model_m_min = 8L),
  mmin_12     = list(model_m_min = 12L),
  lambda0_0.5 = list(lambda0_K = 0.5),
  lambda0_2   = list(lambda0_K = 2.0),
  m0_0.01     = list(hyper = list(m0 = 0.01)),
  m0_0.10     = list(hyper = list(m0 = 0.10))
)

# labels
SENS_LABELS_S2 <- c(
  default     = "Default",
  mmin_8      = "$m_{\\min}=8$",
  mmin_12     = "$m_{\\min}=12$",
  lambda0_0.5 = "$\\lambda_0=0.5$",
  lambda0_2   = "$\\lambda_0=2$",
  m0_0.01     = "$m_0=0.01$",
  m0_0.10     = "$m_0=0.10$"
)

# model$hyper_override replaces entries of params$hyper after initialization (m0)
if (!isTRUE(attr(initialize_model, "sens_wrapper"))) {
  .initialize_model_base <- initialize_model
  initialize_model <- function(Y, model) {
    init <- .initialize_model_base(Y, model)
    ov <- model$hyper_override
    for (nm in names(ov)) init$params$hyper[[nm]] <- ov[[nm]]
    init
  }
  attr(initialize_model, "sens_wrapper") <- TRUE
}

# One replication under one setting, as in run_one_hierscp_rep()
run_one_sens_rep <- function(seed, setting_name, setting = SENS_SETTINGS_S2[[setting_name]],
                             scenario = "S2",
                             sim_J = 9L, sim_T_len = 120L, sim_C = 3L,
                             sim_cp_upper = list(c(50L), c(42L), c(40L, 85L)),
                             sim_M_basis = 8L, sim_m_min = 10L,
                             model_M = 8L, model_m_min = 10L,
                             model_K_init = 2, model_K_min = 2, model_K_max = 6,
                             lambda0_K = 1.0,
                             n_iter = 20000L, n_burnin = 10000L, n_thin = 1L,
                             n_warmup = 100L, ia_every = 20L, N_launch = 2L,
                             verbose_mcmc = 0L, f1_margin = 10, raw_path = NULL) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  
  sim <- simulate_shape_driven(
    scenario = scenario, J = sim_J, T_len = sim_T_len, C = sim_C,
    cp_upper = sim_cp_upper, M_basis = sim_M_basis, m_min = sim_m_min,
    seed = seed, return_scale = "standardized", attach_model_objects = FALSE
  )
  sim_m <- adapt_sim_for_metrics(sim)
  
  m_min_use   <- as.integer(setting$model_m_min %||% model_m_min)
  lambda0_use <- as.numeric(setting$lambda0_K %||% lambda0_K)
  
  model <- create_model_spec(
    sim$Y, C = sim_C,
    K_init = model_K_init, K_max = model_K_max, K_min = model_K_min,
    M = model_M, m_min = m_min_use, type = "continuous",
    fixed_clusters = TRUE, fixed_cluster_sizes = TRUE,
    lambda0_K = lambda0_use, group_members = sim_m$group_members
  )
  model$hyper_override <- setting$hyper
  
  result <- run_mcmc(
    sim$Y, model,
    n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin,
    n_warmup = n_warmup, ia_every = ia_every,
    N_launch = N_launch, verbose = verbose_mcmc
  )
  
  applied <- list(
    m_min   = result$model$m_min,
    lambda0 = result$final_params$hyper$lambda0_K,
    m0      = result$final_params$hyper$m0
  )
  
  post <- compute_conditional_posterior_summary(result, sim$Y, verbose = FALSE)
  rm(result)
  
  if (!is.null(raw_path)) {
    saveRDS(list(sim = sim, sim_m = sim_m, post = post, model = model,
                 scenario = scenario, seed = seed,
                 setting = setting_name, setting_spec = setting, applied = applied),
            raw_path)
  }
  
  m <- compute_all_metrics(post, sim, f1_margin = f1_margin, verbose = FALSE)
  row <- extract_hierscp_summary(m, sim_m, scenario = scenario, seed = seed)
  row$setting <- setting_name
  row$applied_m_min   <- applied$m_min
  row$applied_lambda0 <- applied$lambda0
  row$applied_m0      <- applied$m0
  row
}

# Parallel runs over (setting, seed); finished jobs are skipped
run_sensitivity_S2 <- function(seeds, settings, workers, out_folder = "S2_sensitivity", sim_source,
                               core_dir = getwd(), sens_script = file.path(core_dir, "sensitivity_S2.R"), ...) {
  core_dir    <- normalizePath(core_dir)
  sim_source  <- normalizePath(sim_source)
  sens_script <- normalizePath(sens_script)
  for (s in names(settings)) dir.create(file.path(out_folder, s), showWarnings = FALSE, recursive = TRUE)
  out_dir <- normalizePath(out_folder)
  
  jobs <- expand.grid(setting = names(settings), seed = seeds, stringsAsFactors = FALSE)
  jobs <- jobs[order(jobs$seed, match(jobs$setting, names(settings))), ]
  jobs$out_path <- file.path(out_dir, jobs$setting, sprintf("S2_%s_seed%d.rds", jobs$setting, jobs$seed))
  jobs <- jobs[!file.exists(jobs$out_path), ]
  job_list <- lapply(seq_len(nrow(jobs)), function(i) as.list(jobs[i, ]))
  message(sprintf("[sensitivity] %d jobs to run | %d workers", length(job_list), workers))
  if (!length(job_list)) return(invisible(collect_sensitivity_rows(out_dir)))
  
  dots    <- list(...)
  src_fun <- .sens_source_functions
  
  run_job <- function(jb) {
    log_path <- sub("\\.rds$", ".log", jb$out_path)
    con <- file(log_path, open = "wt")
    stamp <- function(txt) {
      cat(sprintf("[%s] %s seed %d | %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  jb$setting, jb$seed, txt), file = con)
      flush(con)
    }
    t0 <- Sys.time()
    stamp("START")
    sink(con, type = "output")
    on.exit({ sink(); close(con) }, add = TRUE)
    
    out <- tryCatch({
      setwd(core_dir)
      suppressMessages({
        library(Rcpp); library(RcppArmadillo); library(mvtnorm)
        for (f in sort(list.files(core_dir, pattern = "^[0-9].*\\.R$"))) source(f)
        src_fun(sim_source)
        source(sens_script)
      })
      do.call(run_one_sens_rep,
              c(list(seed = jb$seed, setting_name = jb$setting,
                     setting = settings[[jb$setting]],
                     raw_path = sub("\\.rds$", "_raw.rds", jb$out_path)), dots))
    }, error = function(e) { stamp(paste("FAILED:", conditionMessage(e))); NULL })
    
    if (!is.null(out)) {
      saveRDS(out, jb$out_path)
      stamp(sprintf("DONE %.1f min | m_min=%s lambda0=%s m0=%s",
                    as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    out$applied_m_min[1], out$applied_lambda0[1], out$applied_m0[1]))
    }
    invisible(NULL)
  }
  
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  future.apply::future_lapply(job_list, run_job, future.seed = TRUE,
                              future.globals = c("core_dir", "sim_source", "sens_script",
                                                 "settings", "dots", "src_fun"),
                              future.packages = character(0), future.chunk.size = 1L)
  invisible(collect_sensitivity_rows(out_dir))
}

collect_sensitivity_rows <- function(out_dir) {
  fs <- list.files(out_dir, pattern = "^S2_.*_seed[0-9]+\\.rds$", recursive = TRUE, full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, readRDS))
}

# Check that the data sets equal those of the main S2 run
sens_check_default_data <- function(default_folder, seeds = 2026:2045, scenario = "S2") {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  out <- do.call(rbind, lapply(seeds, function(s) {
    fp <- list.files(default_folder, pattern = sprintf("_seed%d_raw\\.rds$", s), full.names = TRUE)
    if (!length(fp)) return(data.frame(seed = s, file = NA_character_, same_data = NA))
    Y0 <- readRDS(fp[1])$sim$Y
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")  # as in run_one_sens_rep()
    Y1 <- simulate_shape_driven(scenario = scenario, seed = s, return_scale = "standardized",
                                attach_model_objects = FALSE)$Y
    data.frame(seed = s, file = basename(fp[1]),
               same_data = isTRUE(all.equal(Y0, Y1, tolerance = 0)))
  }))
  message(if (all(out$same_data %in% TRUE)) {
    "All data sets are identical to the main S2 run: the main-run default can be reused."
  } else {
    "Some data sets differ from the main S2 run: run the default setting as well."
  })
  out
}

# Metrics of the saved fits
.sens_rows_from_raw <- function(files, setting_name) {
  do.call(rbind, lapply(files, function(fp) {
    r  <- readRDS(fp)
    m  <- compute_all_metrics(r$post, r$sim, f1_margin = 10, verbose = FALSE)
    sr <- extract_hierscp_summary(m, r$sim_m, scenario = "S2", seed = r$seed)
    sa <- sr[sr$level == "series_all", , drop = FALSE]
    sp <- spt_metrics_one(r)$row
    data.frame(setting = setting_name, seed = r$seed,
               ARI = sa$ARI, Cover = sa$Cover, NVI = sa$NVI, F1 = sa$F1, dK = sa$dK,
               Kmatch = mean(r$post$K_hat == r$sim_m$K_true),
               RISE = sp$RISE, SCApt_joint = sp$SCApt_joint, SCApt_d1 = sp$SCApt_d1,
               SCApt_d2 = sp$SCApt_d2, Bpt_d1 = sp$Bpt_d1, Bpt_d2 = sp$Bpt_d2,
               m_min = r$applied$m_min %||% r$model$m_min,
               lambda0 = r$applied$lambda0 %||% r$model$lambda0_K,
               m0 = r$applied$m0 %||% 0.03,
               stringsAsFactors = FALSE)
  }))
}

# Summary by setting and paired differences from the default
summarize_sensitivity_S2 <- function(out_folder = "S2_sensitivity", default_folder = NULL,
                                     seeds = 2026:2045,
                                     setting_order = names(SENS_SETTINGS_S2), digits = 3) {
  seed_pat <- function(s) sprintf("_seed%d_raw\\.rds$", s)
  
  per_rep <- list()
  for (s in setting_order) {
    d <- file.path(out_folder, s)
    fs <- list.files(d, pattern = "_raw\\.rds$", full.names = TRUE)
    if (s == "default" && !length(fs) && !is.null(default_folder)) {
      fs <- unlist(lapply(seeds, function(x) list.files(default_folder, pattern = seed_pat(x), full.names = TRUE)))
    }
    if (length(fs)) per_rep[[s]] <- .sens_rows_from_raw(fs, s)
  }
  per_rep <- do.call(rbind, per_rep)
  per_rep <- per_rep[per_rep$seed %in% seeds, , drop = FALSE]
  rownames(per_rep) <- NULL
  
  sets <- intersect(setting_order, unique(per_rep$setting))
  table <- do.call(rbind, lapply(sets, function(s) {
    d <- per_rep[per_rep$setting == s, , drop = FALSE]
    data.frame(setting = s, R = nrow(d),
               ARI = median(d$ARI), Cover = median(d$Cover), NVI = median(d$NVI),
               F1 = median(d$F1), dK = median(d$dK), Kmatch_rate = mean(d$Kmatch),
               RISE = mean(d$RISE), SCApt_joint = mean(d$SCApt_joint),
               SCApt_d1 = mean(d$SCApt_d1), SCApt_d2 = mean(d$SCApt_d2),
               Bpt_d1 = mean(d$Bpt_d1), Bpt_d2 = mean(d$Bpt_d2),
               m_min = paste(unique(d$m_min), collapse = "/"),
               lambda0 = paste(unique(d$lambda0), collapse = "/"),
               m0 = paste(unique(d$m0), collapse = "/"),
               stringsAsFactors = FALSE)
  }))
  num <- vapply(table, is.numeric, logical(1)) & names(table) != "R"
  table[num] <- lapply(table[num], round, digits = digits)
  
  paired <- NULL
  if ("default" %in% sets) {
    base <- per_rep[per_rep$setting == "default", , drop = FALSE]
    mets <- c("ARI", "Cover", "NVI", "F1", "dK", "RISE", "SCApt_joint", "SCApt_d1",
              "SCApt_d2", "Bpt_d1", "Bpt_d2")
    paired <- do.call(rbind, lapply(setdiff(sets, "default"), function(s) {
      d  <- per_rep[per_rep$setting == s, , drop = FALSE]
      mm <- merge(d, base, by = "seed", suffixes = c("", ".def"))
      v  <- vapply(mets, function(x) mean(mm[[x]] - mm[[paste0(x, ".def")]]), numeric(1))
      data.frame(setting = s, n_pairs = nrow(mm), t(round(v, digits)), stringsAsFactors = FALSE)
    }))
  }
  
  list(table = table, paired = paired, per_rep = per_rep)
}

# LaTeX table rows
sens_latex_rows <- function(table, digits = 3, labels = SENS_LABELS_S2) {
  f <- function(v) formatC(v, format = "f", digits = digits)
  rows <- vapply(seq_len(nrow(table)), function(i) {
    s  <- table$setting[i]
    lb <- if (!is.na(labels[s])) labels[[s]] else s
    paste(lb, f(table$ARI[i]), f(table$Cover[i]), f(table$NVI[i]), f(table$F1[i]),
          f(table$dK[i]), f(table$RISE[i]), f(table$SCApt_joint[i]), f(table$SCApt_d1[i]),
          f(table$SCApt_d2[i]), f(table$Bpt_d1[i]), paste0(f(table$Bpt_d2[i]), " \\\\"),
          sep = " & ")
  }, character(1))
  cat(rows, sep = "\n")
  invisible(rows)
}





