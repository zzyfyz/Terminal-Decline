library(MASS)
library(survival)
library(tidyr)
library(dplyr)

dirg <- "C:/Yizhou/Term/Terminal-Decline/Piecewise"
setwd(dirg)

num_simulations <- 300

# ---- helper: simulate one survival time under piecewise baseline hazard ----
rpexp_one <- function(U, eta, breaks, h0) {
  # U: uniform(0,1), eta: linear predictor (log HR), breaks: cutpoints (increasing), h0: baseline hazards per interval
  # Returns a single survival time (no administrative censoring applied here)
  # We solve H0(t) * exp(eta) = -log(U)
  need <- -log(U) / exp(eta)  # required baseline cumulative hazard
  K <- length(h0)
  # Traverse intervals to consume cumulative baseline hazard
  for (k in seq_len(K)) {
    left <- breaks[k]
    right <- breaks[k + 1]
    width <- right - left
    incr <- h0[k] * width  # baseline cumhaz over this full interval
    if (need <= incr || !is.finite(incr)) {
      # Finish inside this interval
      t <- left + need / h0[k]
      return(t)
    } else {
      need <- need - incr
    }
  }
  # Fallback (should not reach if last break is Inf); return large time
  return(breaks[length(breaks)] + need / tail(h0, 1))
}

# ---- main simulation loop ----
for (sim in 1:num_simulations) {
  set.seed(123 + sim)  
  
  # Parameters
  cluster <- 10
  cluster_subj <- 25
  n <- cluster * cluster_subj
  study_duration <- 12          # months
  measurement_interval <- 3     # months
  
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
  
  b <- 0.3     # scaling for subject RE in hazard
  c <- -0.2    # scaling for cluster RE in hazard
  
  sigma_u <- 0.5
  sigma_b <- 1
  sigma_e <- 1
  
  # Fixed effects covariates
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n, 70, 10)
  subject_cluster <- rep(1:cluster, each = cluster_subj)
  treatment_clusters <- sample(1:cluster, size = cluster/2, replace = FALSE)
  treatment <- ifelse(subject_cluster %in% treatment_clusters, 0, 1)
  
  # Random effects
  ui <- rnorm(cluster, mean = 0, sd = sigma_u)
  bi <- rnorm(n, mean = 0, sd = sigma_b)
  
  # --------- PIECEWISE-CONSTANT BASELINE HAZARD SETUP ----------
  # Choose cutpoints (months) and hazards (per month) for baseline
  # Example shape: higher early hazard, dip mid, uptick late
  breaks <- c(0, 4, 8, Inf)     # 5 intervals up to 12m, last is open
  h0     <- c(2, 4, 6)  # baseline hazards per interval
  
  # Linear predictor for PH model (log HR multiplier)
  eta <- alpha10 + alpha11 * x1 + alpha12 * x2 + b * bi + c * ui[subject_cluster]
  
  # Simulate survival times via inverse transform with piecewise baseline
  U <- runif(n)
  survival_times <- vapply(
    seq_len(n),
    function(i) rpexp_one(U[i], eta[i], breaks = breaks, h0 = h0),
    numeric(1)
  )
  
  # Random censoring (independent)
  censoring_rate <- 1 / (20 * study_duration)
  censoring_times <- rexp(n, rate = censoring_rate)
  
  # Observed times and status, apply administrative study cap
  observed_times <- pmin(survival_times, censoring_times, study_duration)
  status <- as.numeric(survival_times <= pmin(censoring_times, study_duration))
  
  # Generate randomized measurement times for each subject
  measurement_times <- t(sapply(1:n, function(i) {
    scheduled_times <- seq(0, study_duration, by = measurement_interval)
    actual_times <- ifelse(
      scheduled_times > 0,
      scheduled_times + runif(length(scheduled_times), -1.5, 1.5),
      scheduled_times
    )
    actual_times <- actual_times[actual_times < observed_times[i]]  # only before observed time
    actual_times[actual_times > study_duration] <- study_duration
    if (length(actual_times) < 5) {
      actual_times <- c(actual_times, rep(NA, 5 - length(actual_times)))
    }
    actual_times[1:5]
  }))
  
  measurement_times_df <- as.data.frame(measurement_times)
  colnames(measurement_times_df) <- c("time_1", "time_2", "time_3", "time_4", "time_5")
  
  # QoL at measurement times (uses true survival_time for backward time)
  qol_values <- t(sapply(1:n, function(i) {
    out <- numeric(5)
    for (j in 1:5) {
      t <- measurement_times[i, j]
      if (!is.na(t)) {
        backward_time <- survival_times[i] - t
        out[j] <- alpha00 +
          (alpha03 / (1 + exp(alpha04 * backward_time))) +
          treatment[i] * (alpha05 + alpha06 * exp(alpha07 * backward_time + alpha08)) +
          alpha01 * x1[i] + alpha02 * x2[i] +
          bi[i] + ui[subject_cluster[i]] + rnorm(1, 0, sigma_e)
      } else {
        out[j] <- NA
      }
    }
    out
  }))
  
  qol_values_df <- as.data.frame(qol_values)
  colnames(qol_values_df) <- c("qol_1", "qol_2", "qol_3", "qol_4", "qol_5")
  
  # Mask (1 observed, 0 missing)
  mask_df <- as.data.frame(!is.na(measurement_times_df) * 1)
  colnames(mask_df) <- c("mask_1", "mask_2", "mask_3", "mask_4", "mask_5")
  
  # Final dataset
  final_data <- data.frame(subject = 1:n) %>%
    cbind(measurement_times_df, qol_values_df) %>%
    left_join(
      data.frame(
        subject = 1:n, x1, x2, cluster = subject_cluster, treatment = treatment,
        observed_time = observed_times, survival_time = survival_times,
        censoring_time = censoring_times, status = status
      ),
      by = "subject"
    )
  
  # Save
  write.csv(final_data, paste0("sim.data.", sim - 1, ".csv"), row.names = FALSE)
  write.csv(mask_df, paste0("mask.", sim - 1, ".csv"), row.names = FALSE)
}