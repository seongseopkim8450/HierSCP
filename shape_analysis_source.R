library(stringr)

setwd("C://Users/User/Desktop/HierSCP_code/shape_analysis")

.est_lower_states <- function(po, J, T_len, cp_field) {
  S <- matrix(NA_integer_, J, T_len)
  for (j in seq_len(J)) {
    cps <- po$lower_cp_summary[[j]]
    tau <- if (is.null(cps) || !length(cps)) integer(0)
    else vapply(cps, function(z) {
      v <- z[[cp_field]]
      if (is.null(v) || !length(v)) NA_integer_ else as.integer(round(v[1]))
    }, integer(1))
    tau <- sort(tau[is.finite(tau) & tau >= 2L & tau <= T_len])
    st <- rep(1L, T_len)
    for (i in seq_along(tau)) st[tau[i]:T_len] <- i + 1L
    S[j, ] <- st
  }
  S
}


##
compute_shape_metrics <- function(rds_path, cp_field = "median") {
  
  x  <- readRDS(rds_path)
  tr <- x$sim$truth
  po <- x$post
  
  J <- tr$J; T_len <- tr$T_len; C <- tr$C
  clus <- tr$cluster_assign
  
  ## RISE
  mu0 <- tr$mu_std; muhat <- po$mu_mean
  RISE <- mean(sqrt(rowMeans((muhat - mu0)^2)))
  
  ## est signs
  est <- lapply(seq_len(C), function(cc) {
    sc <- po$shape_class[[cc]]
    if (is.null(sc) || !length(sc)) return(NULL)
    data.frame(
      d1 = vapply(sc, function(s) as.integer(s$d1_sign),       integer(1)),
      d2 = vapply(sc, function(s) as.integer(s$d2_sign),       integer(1)),
      p1 = vapply(sc, function(s) as.numeric(s$d1_confidence), numeric(1)),
      p2 = vapply(sc, function(s) as.numeric(s$d2_confidence), numeric(1))
    )
  })
  
  ## True signs
  tru <- lapply(seq_len(C), function(cc) {
    a <- tr$atoms[[cc]]
    data.frame(d1 = as.integer(a$delta1), d2 = as.integer(a$delta2))
  })
  
  K_true <- vapply(tru, nrow, integer(1))
  K_hat  <- vapply(est, function(e) if (is.null(e)) 0L else nrow(e), integer(1))
  Ceq    <- which(K_hat == K_true & K_hat > 0L)
  
  ## SCApt 
  S_true <- tr$S_lower
  S_hat  <- .est_lower_states(po, J, T_len, cp_field)
  pt_both <- pt_slope <- pt_curv <- 0L; ptden <- 0L
  for (j in seq_len(J)) {
    e <- est[[clus[j]]]; t0 <- tru[[clus[j]]]
    if (is.null(e)) next
    kt <- S_true[j, ]; kh <- S_hat[j, ]
    ok <- is.finite(kt) & is.finite(kh) &
      kt >= 1 & kt <= nrow(t0) & kh >= 1 & kh <= nrow(e)
    if (!any(ok)) next
    ok_slope <- e$d1[kh[ok]] == t0$d1[kt[ok]]
    ok_curv  <- e$d2[kh[ok]] == t0$d2[kt[ok]]
    pt_slope <- pt_slope + sum(ok_slope)
    pt_curv  <- pt_curv  + sum(ok_curv)
    pt_both  <- pt_both  + sum(ok_slope & ok_curv)
    ptden    <- ptden    + sum(ok)
  }
  SCApt       <- if (ptden) pt_both /ptden else NA_real_
  SCApt_slope <- if (ptden) pt_slope/ptden else NA_real_
  SCApt_curv  <- if (ptden) pt_curv /ptden else NA_real_
  
  ## Brier 
  d1h <- d2h <- d1t <- d2t <- integer(0); p1v <- p2v <- numeric(0)
  for (cc in Ceq) {
    d1h <- c(d1h, est[[cc]]$d1); d1t <- c(d1t, tru[[cc]]$d1); p1v <- c(p1v, est[[cc]]$p1)
    d2h <- c(d2h, est[[cc]]$d2); d2t <- c(d2t, tru[[cc]]$d2); p2v <- c(p2v, est[[cc]]$p2)
  }
  if (length(d1h)) {
    p1p <- ifelse(d1h == 1L, p1v, 1 - p1v)      # P(delta1 = +1 | y)
    p2p <- ifelse(d2h == 1L, p2v, 1 - p2v)      # P(delta2 = +1 | y)
    y1 <- as.integer(d1t == 1L); y2 <- as.integer(d2t == 1L)
    sq_slope <- sum((p1p - y1)^2)
    sq_curv  <- sum((p2p - y2)^2)
    nB <- length(y1)
  } else { sq_slope <- sq_curv <- 0; nB <- 0L }
  
  data.frame(
    scenario = as.character(x$scenario), seed = as.integer(x$seed),
    RISE = RISE,
    SCApt = SCApt, SCApt_slope = SCApt_slope, SCApt_curv = SCApt_curv,
    sq_slope = sq_slope, sq_curv = sq_curv, nB = nB,
    K_hat = paste(K_hat, collapse=","), K_true = paste(K_true, collapse=","),
    stringsAsFactors = FALSE
  )
}


summarize_shape_metrics <- function(files, cp_field = "median") {
  
  per_rep <- do.call(rbind, lapply(files, compute_shape_metrics, cp_field = cp_field))
  
  num <- c("RISE", "SCApt", "SCApt_slope", "SCApt_curv")
  
  summary <- do.call(rbind, lapply(split(per_rep, per_rep$scenario), function(d)
    cbind(scenario = d$scenario[1], n_rep = nrow(d),
          as.data.frame(lapply(d[num], mean, na.rm = TRUE)),
          Brier_slope = sum(d$sq_slope)/sum(d$nB),
          Brier_curv  = sum(d$sq_curv) /sum(d$nB),
          stringsAsFactors = FALSE)))
  rownames(summary) <- NULL
  
  list(per_rep = per_rep, summary = summary)
}



S1_file<-list.files("S1_replicate_50rep/",pattern = "_raw\\.rds$", full.names = TRUE)
S2_file<-list.files("S2_replicate_50rep/",pattern = "_raw\\.rds$", full.names = TRUE)
S3_file<-list.files("S3_replicate_50rep/",pattern = "_raw\\.rds$", full.names = TRUE)

file_raw_extract<-function(list_f){
 
  raw_fname<-list_f[str_ends(list_f,"_raw\\.rds")] 
  
  return(raw_fname)
}



s1_raw_file <- file_raw_extract(S1_file)
s2_raw_file <- file_raw_extract(S2_file)
s3_raw_file <- file_raw_extract(S3_file)


summarize_shape_s1<-summarize_shape_metrics(files = s1_raw_file,cp_field = "median")
summarize_shape_s2<-summarize_shape_metrics(files = s2_raw_file,cp_field = "median")
summarize_shape_s3<-summarize_shape_metrics(files = s3_raw_file,cp_field = "median")


print(summarize_shape_s1$summary)
print(summarize_shape_s2$summary)
print(summarize_shape_s3$summary)








