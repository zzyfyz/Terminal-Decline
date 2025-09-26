library(MASS)
library(survival)
library(tidyr)
library(dplyr)

dirg <- "C:/Yizhou/Term/Terminal-Decline/Unbalanced"
setwd(dirg)

num_simulations <- 300


## ---------- REST EXACT CLUSTER SIZES (set once, used every sim) ----------
sizes_control  <- c(7, 8, 8, 11, 14, 15, 37, 38, 48)
sizes_treat    <- c(1, 3, 6, 17, 64, 102)
cluster_info <- dplyr::tibble(
  cluster  = 1:(length(sizes_control) + length(sizes_treat)),
  size     = c(sizes_control, sizes_treat),
  treatment= c(rep(0L, length(sizes_control)), rep(1L, length(sizes_treat)))
)
stopifnot(sum(cluster_info$size) == 379L, nrow(cluster_info) == 15L)

# Subject -> cluster mapping (fixed across reps)
subject_cluster <- rep(cluster_info$cluster, times = cluster_info$size)
treatment_vec  <- cluster_info$treatment[subject_cluster]
n <- length(subject_cluster)

# Simulation loop
for (sim in 1:num_simulations) {
  set.seed(123 + sim)
  
  ## ---------- USE REST MAPPING ----------
  J_clusters <- nrow(cluster_info)
  # subject_cluster, treatment_vec, n are already defined above and fixed
  
  ## Study design
  study_duration <- 9          # months
  measurement_interval <- 3     # every 3 months: 0,3,6,9,12
  
  ## Fixed-effect parameters (yours)
  alpha00 <- 5
  alpha01 <- 0.7
  alpha02 <- -0.03
  alpha03 <- 2
  alpha04 <- -0.6
  alpha05 <- 0.01
  alpha06 <- 0.8
  alpha07 <- -0.5
  alpha08 <- -0.5
  
  alpha10 <- -5.4
  alpha11 <- -0.02
  alpha12 <- 0.01
  
  b <- 0.3
  c <- -0.2
  
  gamma   <- 1.4
  sigma_u <- 0.5
  sigma_b <- 1
  sigma_e <- 1
  
  ## Covariates
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n, 70, 10)
  
  ## Cluster & subject random effects
  ui <- rnorm(J_clusters, mean = 0, sd = sigma_u)  # cluster RE
  bi <- rnorm(n,          mean = 0, sd = sigma_b)  # subject RE
  
  ## Weibull survival with frailty & subject RE
  linear  <- alpha10 + alpha11*x1 + alpha12*x2 + b*bi + c*ui[subject_cluster]
  lambda  <- exp(linear)
  U       <- runif(n)
  survival_times <- (-log(U) / lambda)^(1 / gamma)
  
  ## Censoring
  censoring_rate  <- 1 / (4 * study_duration)
  censoring_times <- rexp(n, rate = censoring_rate)
  
  ## Observed time & status
  observed_times <- pmin(survival_times, censoring_times, study_duration)
  status <- as.integer(survival_times <= pmin(censoring_times, study_duration))
  
  ## tuning
  dev3   <- 1.0    # +/- jitter around 3
  dev6   <- 1.0    # +/- jitter around 6
  dev9   <- 0.25   # +/- jitter around 9
  eps    <- 1e-6
  gap_min <- 0.0   # require at least this gap from previous kept time (e.g., 0.5)
  
  measurement_times <- t(vapply(seq_len(n), function(i) {
    Tobs  <- observed_times[i]                      # min(survival, censoring, 9)
    died  <- (status[i] == 1L)
    surv9 <- (!died) && (abs(Tobs - study_duration) < 1e-8)  # alive through month 9
    
    # scheduled with jitter (0 fixed)
    t1 <- 0
    t2r <- 3 + runif(1, -dev3, dev3)
    t3r <- 6 + runif(1, -dev6, dev6)
    t4r <- 9 + runif(1, -dev9, dev9)
    
    # clamp to [0, 9] for drawing; we'll enforce observed-time after
    t2 <- pmax(0, pmin(t2r, study_duration))
    t3 <- pmax(0, pmin(t3r, study_duration))
    t4 <- pmax(0, pmin(t4r, study_duration))
    
    # apply ">= observed_time -> NA" for time_2 and time_3
    if (t2 >= Tobs - eps) t2 <- NA_real_
    if (t3 >= Tobs - eps) t3 <- NA_real_
    
    # previous available time (to enforce ordering for time_4)
    prev <- max(c(t1, t2, t3), na.rm = TRUE)
    
    # LAST MEASUREMENT RULES
    if (surv9 && (t4r > study_duration)) {
      # survivor to 9 and simulated > 9  -> keep exactly 9
      t4 <- study_duration
      # ensure ordering; if no room, drop
      if (t4 <= prev + gap_min + eps) t4 <- NA_real_
    } else {
      # everyone else: do NOT truncate to observed time; just obey the NA rule
      if (t4 >= Tobs - eps) t4 <- NA_real_
      # enforce ordering if kept
      if (!is.na(t4) && t4 <= prev + gap_min + eps) t4 <- NA_real_
    }
    
    c(t1, t2, t3, t4)
  }, FUN.VALUE = numeric(4)))
  
  measurement_times_df <- as.data.frame(measurement_times)
  colnames(measurement_times_df) <- paste0("time_", 1:4)
  
  ## QoL at measurement times (backward-time structure)
  qol_values <- t(sapply(seq_len(n), function(i) {
    out <- rep(NA_real_, 4)
    for (j in 1:4) {
      t_ij <- measurement_times[i, j]
      if (!is.na(t_ij)) {
        bt <- observed_times[i] - t_ij
        out[j] <- alpha00 +
          (alpha03 / (1 + exp(alpha04 * bt))) +
          treatment_vec[i] * (alpha05 + alpha06 * exp(alpha07 * bt + alpha08)) +
          alpha01 * x1[i] + alpha02 * x2[i] +
          bi[i] + ui[subject_cluster[i]] + rnorm(1, 0, sigma_e)
      }
    }
    out
  }))
  qol_values_df <- as.data.frame(qol_values)
  colnames(qol_values_df) <- paste0("qol_", 1:4)
  
  ## Mask (1 observed, 0 missing)
  mask_df <- as.data.frame(!is.na(measurement_times_df)) * 1L
  colnames(mask_df) <- paste0("mask_", 1:4)
  
  ## Final dataset
  final_data <- data.frame(subject = seq_len(n)) %>%
    cbind(measurement_times_df, qol_values_df) %>%
    dplyr::left_join(
      data.frame(
        subject = seq_len(n),
        x1 = x1, x2 = x2,
        cluster = subject_cluster,
        treatment = treatment_vec,
        observed_time  = observed_times,
        survival_time  = survival_times,
        censoring_time = censoring_times,
        status = status
      ),
      by = "subject"
    )
  
  ## Save
  write.csv(final_data, file = sprintf("sim.data.%d.csv", sim - 1), row.names = FALSE)
  write.csv(mask_df,   file = sprintf("mask.%d.csv",      sim - 1), row.names = FALSE)
}