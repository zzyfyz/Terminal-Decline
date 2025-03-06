library(MASS)
library(survival)
library(tidyr)
library(dplyr)

dirg <- "C:/Yizhou/Term/Terminal-Decline/Spline"
setwd(dirg)

num_simulations <- 100

# Simulation loop
for (sim in 1:num_simulations) {
  set.seed(123 + 99)  
  
  # Parameters
  cluster <- 10
  cluster_subj <- 25
  n <- cluster * cluster_subj
  study_duration <- 12  # Study duration in months
  measurement_interval <- 3  # Intended interval (every 6 months)
  
  alpha00 <- 6
  alpha01 <- 0.05
  alpha02 <- -0.003
  alpha03 <- -1
  alpha04 <- 0.6 
  alpha05 <- 0.02
  alpha06 <- 0.3 
  alpha07 <- -0.5
  alpha08 <- -0.5
  alpha11 <- -0.02
  alpha12 <- 0.001
  
  b <- 0.3
  c <- 0.2
  lambda0 <- 0.02
  gamma <- 1.3
  sigma_u <- 0.5
  sigma_b <- 1
  sigma_e <- 0.5
  
  # Fixed effects covariates
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n, 70, 10)
  subject_cluster <- rep(1:cluster, each = cluster_subj)
  treatment_clusters <- sample(1:cluster, size = cluster/2, replace = FALSE)
  treatment <- ifelse(subject_cluster %in% treatment_clusters, 0, 1)
  
  # Random effects
  ui <- rnorm(cluster, mean = 0, sd = sigma_u)
  bi <- rnorm(n, mean = 0, sd = sigma_b)
  
  # Cox frailty model for survival data assuming Weibull distribution
  linear <- alpha11 * x1 + alpha12 * x2 + b * bi + c * ui[subject_cluster]
  lambda <- lambda0 * exp(linear)
  U <- runif(n)
  
  # Simulate survival times
  survival_times <- (-log(U) / lambda)^(1 / gamma)
  censoring_times <- study_duration
  observed_times <- pmin(survival_times, censoring_times)
  status <- as.numeric(survival_times <= censoring_times)
  
  # Generate randomized measurement times for each subject
  measurement_times <- t(sapply(1:n, function(i) {
    scheduled_times <- seq(0, study_duration, by = measurement_interval)
    actual_times <- ifelse(scheduled_times > 0, scheduled_times + runif(length(scheduled_times), -1.5, 1.5), scheduled_times)
    actual_times <- actual_times[actual_times < survival_times[i]]  # Only keep valid times
    actual_times[actual_times > study_duration] <- study_duration

    
    if (length(actual_times) < 5) {
      actual_times <- c(actual_times, rep(NA, 5 - length(actual_times)))
    }
    return(actual_times[1:5]) 
  }))
  
  # Convert to data frame with fixed column names
  measurement_times_df <- as.data.frame(measurement_times)
  colnames(measurement_times_df) <- c("time_1", "time_2", "time_3", "time_4", "time_5")
  
  # Generate QoL values at measurement times
  qol_values <- t(sapply(1:n, function(i) {
    qol_measurements <- numeric(5)  # 
    for (j in 1:5) {
      t <- measurement_times[i, j]
      if (!is.na(t)) {
        backward_time <- survival_times[i] - t
        qol_measurements[j] <- alpha00 + 
          (alpha03 / (1 + exp(alpha04 * backward_time))) + 
          treatment[i] * (alpha05 + alpha06 * exp(alpha07 * backward_time + alpha08)) + 
          alpha01 * x1[i] + alpha02 * x2[i] + 
          bi[i] + ui[subject_cluster[i]] + rnorm(1, mean = 0, sd = sigma_e)
      } else {
        qol_measurements[j] <- NA
      }
    }
    return(qol_measurements)
  }))
  
  # Convert QoL values to data frame with fixed column names
  qol_values_df <- as.data.frame(qol_values)
  colnames(qol_values_df) <- c("qol_1", "qol_2", "qol_3", "qol_4", "qol_5")
  
  # Create mask dataset (1 if observed, 0 if missing)
  mask <- !is.na(measurement_times_df) * 1  # Convert TRUE/FALSE to 1/0
  mask_df <- as.data.frame(mask)
  colnames(mask_df) <- c("mask_1", "mask_2", "mask_3", "mask_4", "mask_5")
  
  # Create final dataset
  final_data <- data.frame(subject = 1:n) %>%
    cbind(measurement_times_df, qol_values_df) %>%
    left_join(data.frame(subject = 1:n, x1, x2, cluster = subject_cluster, treatment = treatment, 
                         observed_time = observed_times, survival_time = survival_times, status = status),
              by = "subject")
  
  #table(final_data$status)
  summary(qol_values_df)
  
  # Save simulation data
  write.csv(final_data, paste0("sim.data.", sim-1, ".csv"), row.names = FALSE)
  write.csv(mask_df, paste0("mask.", sim-1, ".csv"), row.names = FALSE)
}

