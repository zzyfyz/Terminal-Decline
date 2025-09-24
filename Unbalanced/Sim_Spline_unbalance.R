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
  study_duration <- 12          # months
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
  
  alpha10 <- -4.8
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
  censoring_rate  <- 1 / (20 * study_duration)
  censoring_times <- rexp(n, rate = censoring_rate)
  
  ## Observed time & status
  observed_times <- pmin(survival_times, censoring_times, study_duration)
  status <- as.integer(survival_times <= pmin(censoring_times, study_duration))
  
  ## Measurement times (≤ 5 per subject; jitter around schedule)
  measurement_times <- t(sapply(seq_len(n), function(i) {
    sched <- seq(0, study_duration, by = measurement_interval)
    actual <- ifelse(sched > 0, sched + runif(length(sched), -1.5, 1.5), sched)
    actual <- actual[actual < observed_times[i]]
    actual[actual > study_duration] <- study_duration
    if (length(actual) < 5) actual <- c(actual, rep(NA_real_, 5 - length(actual)))
    actual[1:5]
  }))
  measurement_times_df <- as.data.frame(measurement_times)
  colnames(measurement_times_df) <- paste0("time_", 1:5)
  
  ## QoL at measurement times (backward-time structure)
  qol_values <- t(sapply(seq_len(n), function(i) {
    out <- rep(NA_real_, 5)
    for (j in 1:5) {
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
  colnames(qol_values_df) <- paste0("qol_", 1:5)
  
  ## Mask (1 observed, 0 missing)
  mask_df <- as.data.frame(!is.na(measurement_times_df)) * 1L
  colnames(mask_df) <- paste0("mask_", 1:5)
  
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
