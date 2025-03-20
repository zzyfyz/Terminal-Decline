library(MASS)
library(survival)
library(tidyr)
library(dplyr)

dirg <- "C:/Yizhou/Term/Terminal-Decline/Spline"
setwd(dirg)

num_simulations <- 300

# Simulation loop
for (sim in 1:num_simulations) {
  set.seed(123 + sim)  
  
  # Parameters
  cluster <- 10
  cluster_subj <- 25
  n <- cluster * cluster_subj
  study_duration <- 12  # Study duration in months
  measurement_interval <- 3  # Intended interval (every 6 months)
  
  alpha00 <- 5
  alpha01 <- 0.7
  alpha02 <- -0.03
  alpha03 <- 2
  alpha04 <- -0.6 
  alpha05 <- 0.01
  alpha06 <- 0.8 
  alpha07 <- -0.5
  alpha08 <- -0.5
  
  alpha10 <- -6
  alpha11 <- -0.02
  alpha12 <- 0.01
  
  b <- 0.3
  c <- -0.2
  
  gamma <- 1.4
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
  
  # Cox frailty model for survival data assuming Weibull distribution
  linear <- alpha10 + alpha11 * x1 + alpha12 * x2 + b * bi + c * ui[subject_cluster]
  lambda <- exp(linear)
  U <- runif(n)
  
  # Simulate survival times
  survival_times <- (-log(U) / lambda)^(1 / gamma)
  
  # Generate random censoring times from an exponential distribution
  censoring_rate <- 1 / (20 * study_duration)
  censoring_times <- rexp(n, rate = censoring_rate)
  
  # Determine observed times and status
  observed_times <- pmin(survival_times, censoring_times, study_duration)
  status <- as.numeric(survival_times <= pmin(censoring_times, study_duration))
  
  # Generate randomized measurement times for each subject
  measurement_times <- t(sapply(1:n, function(i) {
    scheduled_times <- seq(0, study_duration, by = measurement_interval)
    actual_times <- ifelse(scheduled_times > 0, scheduled_times + runif(length(scheduled_times), -1.5, 1.5), scheduled_times)
    actual_times <- actual_times[actual_times < observed_times[i]]  # Only keep valid times
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
        backward_time <- observed_times[i] - t
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
                         observed_time = observed_times, survival_time = survival_times, censoring_time = censoring_times, status = status),
              by = "subject")
  
  # Save simulation data
  write.csv(final_data, paste0("sim.data.", sim-1, ".csv"), row.names = FALSE)
  write.csv(mask_df, paste0("mask.", sim-1, ".csv"), row.names = FALSE)
}


#alpha00 <- 7
#alpha01 <- 0.7
#alpha02 <- -0.03
#alpha03 <- 2
#alpha04 <- -0.6 
#alpha05 <- 0.01
#alpha06 <- 0.8 
#alpha07 <- -0.5
#alpha08 <- -0.5
#alpha11 <- -0.02
#alpha12 <- 0.01

# Define a sequence of backward times from 0 to study duration
#time_points <- seq(0, study_duration, length.out = 100)  # 100 points for smooth curve

# Compute the true expected QoL for both treatment groups
#true_qol <- data.frame(
  #backward_time = time_points,
  #control = alpha00 + (alpha03 / (1 + exp(alpha04 * time_points))) + alpha01*0.5 + alpha02 * 70,
  #treatment = alpha00 + (alpha03 / (1 + exp(alpha04 * time_points))) + 
    #(alpha05 + alpha06 * exp(alpha07 * time_points + alpha08)) + alpha01*0.5 + alpha02 * 70
#)

# Convert to long format for ggplot
#true_qol_long <- true_qol %>%
  #tidyr::pivot_longer(cols = c("control", "treatment"), 
                      #names_to = "group", values_to = "expected_qol")

# Plot the true mean QoL trajectory
#ggplot(true_qol_long, aes(x = backward_time, y = expected_qol, color = group)) +
  #geom_line(size = 1.2) +
  #labs(title = "True Mean QoL Trajectory",
       #x = "Backward Time from Death",
       #y = "Expected QoL",
       #color = "Group") +
  #scale_color_manual(values = c("control" = "blue", "treatment" = "red"),
                     #labels = c("Control", "Treatment")) +
  #theme_minimal()
